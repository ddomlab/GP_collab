# -----------------------------
# Nested 5-Fold CV with Stan GP Model (updated)
# -----------------------------

# 0) Load libraries
library(rstan)
library(dplyr)
library(jsonlite)
library(pheatmap)
library(Matrix)

# 1) Read data
data <- read.csv("training_data.csv", stringsAsFactors = FALSE)
data <- data[ , -1]  # drop first column if it's an index

# 2) Extract response and predictors
response  <- data$log.Rg..nm.
predictors <- data %>% select(-Rg1..nm., -log.Rg..nm., -substructure.cluster)

# 3.a) Compute generalized Tanimoto distance for fingerprint columns
fp_cols <- grep("^Monomer_ECFP6_count_bit", names(predictors), value = TRUE)
if (length(fp_cols) == 0) stop("No fingerprint predictors found.")
fp_data <- as.matrix(predictors[, fp_cols])
n       <- nrow(fp_data)

tanimoto_similarity_counts <- function(x, y) {
  num <- sum(pmin(x, y))
  den <- sum(pmax(x, y))
  if (den == 0) 0 else num / den
}

tanimoto_dist <- matrix(0, n, n)
for (i in seq_len(n)) {
  for (j in i:n) {
    sim          <- tanimoto_similarity_counts(fp_data[i,], fp_data[j,])
    d_val        <- 1 - sim
    tanimoto_dist[i, j] <- d_val
    tanimoto_dist[j, i] <- d_val
  }
}

# 3.b) Heatmap of Tanimoto distance
rownames(data) <- paste0("mol", seq_len(nrow(data)))
dimnames(tanimoto_dist) <- list(rownames(data), rownames(data))

desired_lvls <- c("PPV","Fluorene","Thiophene")
cluster_fact <- factor(data$substructure.cluster, levels = desired_lvls)
ord_idx    <- order(cluster_fact)
counts     <- table(cluster_fact)
boundaries <- cumsum(counts)

mat_ord <- tanimoto_dist[ord_idx, ord_idx]
ann     <- data.frame(cluster = cluster_fact, row.names = rownames(data))
ann_ord <- ann[ord_idx, , drop = FALSE]

my_cols <- list(
  cluster = c(
    PPV       = "#FF66CC",
    Fluorene  = "#FF9999",
    Thiophene = "#66CCFF"
  )
)

pheatmap(
  mat_ord,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  annotation_row = ann_ord,
  annotation_col = ann_ord,
  gaps_row       = boundaries,
  gaps_col       = boundaries,
  border_color   = "white",
  show_rownames  = FALSE,
  show_colnames  = FALSE
)

# 4) Continuous predictors
other_cols <- setdiff(names(predictors), fp_cols)
if (length(other_cols) == 0) stop("No continuous predictors found.")
X_cont <- as.matrix(predictors[, other_cols])

# 5) Stan model code (unchanged)
stan_model_code <- '
data {
  int<lower=1> N;
  int<lower=1> P;
  vector[N] y;
  matrix[N, P] X;
  matrix[N, N] D_fp;
}
parameters {
  real intercept;
  real<lower=0> sigma_cont;
  vector<lower=0>[P] l_cont;
  real<lower=0> sigma_fp;
  real<lower=0> l_fp;
  real<lower=0> sigma_noise;
  vector[N] f;
}
transformed parameters {
  matrix[N,N] K_cont;
  matrix[N,N] K_fp;
  matrix[N,N] K;
  K_cont = rep_matrix(0.0, N, N);
  for (j in 1:P) {
    for (i in 1:N) {
      for (k in i:N) {
        real sq = square(X[i,j] - X[k,j]);
        real v  = sigma_cont^2 * exp(-0.5 * sq / square(l_cont[j]));
        K_cont[i,k] = K_cont[i,k] + v;
        if (i != k) K_cont[k,i] = K_cont[i,k];
      }
    }
  }
  K_fp = rep_matrix(0.0, N, N);
  for (i in 1:N) {
    for (k in i:N) {
      real v = sigma_fp^2 * exp(- D_fp[i,k] / (2 * square(l_fp)));
      K_fp[i,k] = v;
      if (i != k) K_fp[k,i] = v;
    }
  }
  K = K_cont + K_fp;
}
model {
  intercept ~ normal(0, 5);
  sigma_cont ~ normal(0, 1);
  l_cont ~ inv_gamma(5, 5);
  sigma_fp ~ normal(0, 1);
  l_fp ~ inv_gamma(5, 5);
  sigma_noise ~ normal(0, 1);
  f ~ multi_normal(rep_vector(intercept, N), K + diag_matrix(rep_vector(square(sigma_noise), N)));
  y ~ normal(f, sigma_noise);
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    y_rep[n] = normal_rng(f[n], sigma_noise);
    log_lik[n] = normal_lpdf(y[n] | f[n], sigma_noise);
  }
}
'

# 6) Compile Stan model
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stan_mod <- stan_model(model_code = stan_model_code)

# 7) Read 5-fold CV splits from JSON
cv_all <- fromJSON("five_fold_cv.json")
seeds  <- names(cv_all)

# 8) Outer loop: CV over seeds and folds
all_results <- list()
all_coefs   <- list()

for (seed in seeds) {
  cv        <- cv_all[[seed]]
  train_idx <- lapply(cv$train, function(x) x + 1)
  test_idx  <- lapply(cv$test,  function(x) x + 1)
  
  seed_res <- list()
  coef_res <- list()
  
  for (fold in seq_along(train_idx)) {
    tr <- train_idx[[fold]]
    te <- test_idx[[fold]]
    
    # Subset
    X_tv <- X_cont[tr, , drop = FALSE]
    y_tv <- response[tr]
    D_fp_tv <- tanimoto_dist[tr, tr]
    X_te <- X_cont[te, , drop = FALSE]
    y_te <- response[te]
    D_fp_te <- tanimoto_dist[te, tr]
    
    # Scale train + val
    ctr <- colMeans(X_tv, na.rm = TRUE)
    scl <- apply(X_tv, 2, sd, na.rm = TRUE)
    scl[scl == 0] <- 1
    X_tv_s <- sweep(sweep(X_tv, 2, ctr, "-"), 2, scl, "/")
    
    y_ctr <- mean(y_tv, na.rm = TRUE)
    y_scl <- sd(y_tv,   na.rm = TRUE)
    if (y_scl == 0) y_scl <- 1
    y_tv_s <- (y_tv - y_ctr) / y_scl
    
    # Stan data
    stan_data <- list(
      N    = length(y_tv_s),
      P    = ncol(X_tv_s),
      X    = X_tv_s,
      y    = as.numeric(y_tv_s),
      D_fp = D_fp_tv
    )
    
    # Fit model reproducibly
    fit <- sampling(
      stan_mod,
      data   = stan_data,
      chains = 4,
      iter   = 2000,
      warmup = 1000,
      seed   = as.integer(seed)
    )
    
    post <- extract(fit)
    
    
    # Posterior means
    alpha <- mean(post$intercept)
    sc    <- mean(post$sigma_cont)
    lc    <- apply(post$l_cont, 2, mean)
    s_fp  <- mean(post$sigma_fp)
    l_fp  <- mean(post$l_fp)
    s_n   <- mean(post$sigma_noise)
    f_bar <- apply(post$f, 2, mean)
    
    # Reconstruct train covariance
    N_tr <- length(y_tv_s)
    K_tr <- matrix(0, N_tr, N_tr)
    for (j in seq_len(ncol(X_tv_s))) {
      for (i in seq_len(N_tr)) {
        for (k in seq_len(N_tr)) {
          sq <- (X_tv_s[i,j] - X_tv_s[k,j])^2
          K_tr[i,k] <- K_tr[i,k] + sc^2 * exp(-0.5 * sq / lc[j]^2)
        }
      }
    }
    K_tr <- K_tr + s_fp^2 * exp(- D_fp_tv / (2 * l_fp^2)) +
      diag(s_n^2 + 1e-8, N_tr)
    K_tr_pd <- as.matrix(nearPD(K_tr)$mat)
    L       <- chol(K_tr_pd)
    
    # Cross-covariance test
    N_te <- length(te)
    Kc_te <- matrix(0, N_te, N_tr)
    for (j in seq_len(ncol(X_tv_s))) {
      for (i in seq_len(N_te)) {
        for (k in seq_len(N_tr)) {
          x_scaled <- (X_te[i,j] - ctr[j]) / scl[j]
          Kc_te[i,k] <- Kc_te[i,k] +
            sc^2 * exp(-0.5 * (x_scaled - X_tv_s[k,j])^2 / lc[j]^2)
        }
      }
    }
    Kfp_te <- s_fp^2 * exp(- D_fp_te / (2 * l_fp^2))
    Kc_te  <- Kc_te + Kfp_te
    
    # Predictive mean
    w       <- backsolve(L, f_bar - alpha, transpose = TRUE)
    w       <- backsolve(L, w,            transpose = FALSE)
    mu_te_s <- alpha + Kc_te %*% w
    mu_te   <- mu_te_s * y_scl + y_ctr
    
    # Metrics
    rmse <- sqrt(mean((mu_te - y_te)^2))
    mae  <- mean(abs(mu_te - y_te))
    seed_res[[fold]] <- data.frame(seed = seed, fold = fold,
                                   RMSE = rmse, MAE = mae)
    
    # Coef summaries
    sum_stats <- summary(fit,
                         pars = c('intercept','sigma_cont','l_cont','sigma_fp','l_fp','sigma_noise')
    )$summary
    coef_df <- data.frame(
      parameter = rownames(sum_stats),
      mean      = sum_stats[, 'mean'],
      sd        = sum_stats[, 'sd'],
      lower     = sum_stats[, '2.5%'],
      upper     = sum_stats[, '97.5%'],
      seed      = seed,
      fold      = fold,
      row.names = NULL
    )
    coef_res[[fold]] <- coef_df
  }
  
  all_results[[seed]] <- bind_rows(seed_res)
  all_coefs[[seed]]   <- bind_rows(coef_res)
}

# 9) Combine and save
final_results <- bind_rows(all_results)
final_coefs   <- bind_rows(all_coefs)

print(final_results)
print(aggregate(cbind(RMSE, MAE) ~ seed, data = final_results, FUN = mean))
print(final_coefs)

write.csv(final_results, "gp_cv_results.csv", row.names = FALSE)
write.csv(final_coefs,   "gp_cv_coefs.csv",   row.names = FALSE)

save.image("gp_model.RData")
