# -----------------------------
# Fully updated Nested 5-Fold CV with Stan GP Model
# -----------------------------
rm(list = ls())
gc
# 0) Libraries
library(rstan)
library(dplyr)
library(stringr)
library(jsonlite)
library(pheatmap)
library(Matrix)
library(ggplot2)

# 1) Read & clean data
data <- read.csv("training_data.csv", stringsAsFactors = FALSE)
if (names(data)[1] %in% c("X", "index", "ID")) {
  data <- data[ , -1]
} else {
  message("First column '", names(data)[1], "' left in place.")
}

# 2) Response & predictors
response   <- data$log.Rg..nm.
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

# 4) Continuous predictors matrix
other_cols <- setdiff(names(predictors), fp_cols)
if (length(other_cols) == 0) stop("No continuous predictors.")
X_cont <- as.matrix(predictors[, other_cols])


# 5) Stan model code (common sigma)
stan_model_code <- '
data {
  int<lower=1> N;
  int<lower=1> P;
  vector[N] y;
  matrix[N,N] D_fp;
  array[P] matrix[N,N] D2_cont;
}
parameters {
  real intercept;
  real<lower=0> sigma;         // common kernel amplitude
  vector<lower=0>[P] l_cont;   // length-scales for each continuous feature
  real<lower=0> l_fp;          // length-scale for fingerprint kernel
  real<lower=0> sigma_noise;
  vector[N] f;
}
transformed parameters {
  matrix[N,N] K_cont = rep_matrix(0.0, N, N);
  for (j in 1:P) {
    K_cont += square(sigma)
           * exp(-0.5 * D2_cont[j] / square(l_cont[j]));
  }
  matrix[N,N] K_fp = square(sigma)
                   * exp(- D_fp / (2 * square(l_fp)));
  matrix[N,N] K    = K_cont + K_fp;
}
model {
  intercept   ~ normal(0,2);
  sigma       ~ normal(0,1);
  l_cont      ~ inv_gamma(5,5);
  l_fp        ~ inv_gamma(5,5);
  sigma_noise ~ normal(0,1);

  f ~ multi_normal(rep_vector(intercept,N),
                   K + diag_matrix(rep_vector(square(sigma_noise), N)));
  y ~ normal(f, sigma_noise);
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    y_rep[n]   = normal_rng(f[n], sigma_noise);
    log_lik[n] = normal_lpdf(y[n] | f[n], sigma_noise);
  }
}
'
start <- Sys.time()

# 6) Compile Stan model
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
stan_mod <- stan_model(model_code = stan_model_code)

# 7) Load CV splits
cv_all <- fromJSON("five_fold_cv.json")
seeds  <- names(cv_all)

# 8) CV loops
all_results <- list()
all_coefs   <- list()
all_ppi_tr  <- list()
all_ppi_te  <- list()

for (seed in names(cv_all)) {
  cv        <- cv_all[[seed]]
  train_idx <- lapply(cv$train, function(x) x + 1)
  test_idx  <- lapply(cv$test,  function(x) x + 1)
  
  seed_res <- vector("list", length(train_idx))
  coef_res <- vector("list", length(train_idx))
  ppi_tr   <- vector("list", length(train_idx))
  ppi_te   <- vector("list", length(train_idx))
  
  for (fold in seq_along(train_idx)) {
    tr <- train_idx[[fold]]
    te <- test_idx[[fold]]
    
    # --- subset & scale ---
    X_tv   <- X_cont[tr, , drop = FALSE];  y_tv <- response[tr]
    X_te   <- X_cont[te, , drop = FALSE];  y_te <- response[te]
    D_fp_tv <- tanimoto_dist[tr, tr]
    D_fp_te <- tanimoto_dist[te, tr]
    
    ctr <- colMeans(X_tv, na.rm = TRUE)
    scl <- apply(X_tv, 2, sd, na.rm = TRUE); scl[scl==0] <- 1
    X_tv_s <- sweep(sweep(X_tv, 2, ctr, "-"), 2, scl, "/")
    X_te_s <- sweep(sweep(X_te, 2, ctr, "-"), 2, scl, "/")
    
    y_ctr <- mean(y_tv, na.rm = TRUE)
    y_scl <- sd(  y_tv, na.rm = TRUE); if(y_scl==0) y_scl <- 1
    y_tv_s <- (y_tv - y_ctr) / y_scl
    
    N_tr <- length(y_tv_s)
    P    <- ncol(X_tv_s)
    
    # --- precompute D2_cont ---
    D2_cont <- array(0, dim = c(P, N_tr, N_tr))
    for (j in 1:P) D2_cont[j,,] <- as.matrix(dist(X_tv_s[,j]))^2
    
    # --- Stan fit: save y_rep too (for train‐set PPI) ---
    stan_data <- list(N=N_tr, P=P, y=as.numeric(y_tv_s),
                      D_fp=D_fp_tv, D2_cont=D2_cont)
    fit <- sampling(
      stan_mod,
      data    = stan_data,
      pars    = c("intercept","sigma","l_cont","l_fp",
                  "sigma_noise","f","y_rep"),
      include = TRUE,
      chains  = 4, iter = 2000, warmup = 1000,
      control = list(adapt_delta=0.8, max_treedepth=10),
      seed    = as.integer(seed)
    )
    
    # --- posterior summaries for parameters + f ---
    sum_stats <- summary(
      fit,
      pars  = c("intercept","sigma","l_cont","l_fp","sigma_noise","f"),
      probs = c(0.025, 0.5, 0.975)
    )$summary
    
    alpha <-  sum_stats["intercept","50%"]
    sigma <-  sum_stats["sigma","50%"]
    lc    <-  sum_stats[grep("^l_cont", rownames(sum_stats)),"50%"]
    l_fp  <-  sum_stats["l_fp","50%"]
    s_n   <-  sum_stats["sigma_noise","50%"]
    f_bar <-  sum_stats[grep("^f\\[", rownames(sum_stats)),"50%"]
    
    # --- TRAIN‐SET posterior predictive intervals ---
    y_rep_draws <- extract(fit, pars="y_rep")$y_rep  # draws × N_tr
    ppi_tr_mat  <- apply(y_rep_draws, 2,
                         quantile, probs=c(0.025,0.5,0.975))
    ppi_tr[[fold]] <- data.frame(
      obs = tr,
      y_train = y_tv,
      y_ctr  = y_ctr,
      y_scl  = y_scl,
      lower  = ppi_tr_mat[1,],
      median = ppi_tr_mat[2,],
      upper  = ppi_tr_mat[3,],
      seed, fold, row.names = NULL
    )
    
    # --- vectorized predictions for test ‐ no change ---
    X_tv_k <- sweep(X_tv_s, 2, lc, "/")
    X_te_k <- sweep(X_te_s, 2, lc, "/")
    A2     <- rowSums(X_te_k^2);    B2 <- rowSums(X_tv_k^2)
    D2_te  <- outer(A2, rep(1,N_tr)) +
      outer(rep(1,length(te)), B2) -
      2 * (X_te_k %*% t(X_tv_k))
    Kc_cont_te <- sigma^2 * exp(-0.5 * D2_te)
    K_fp_te    <- sigma^2 * exp(- D_fp_te / (2*l_fp^2))
    Kc_te      <- Kc_cont_te + K_fp_te
    
    # --- training kernel + Cholesky in Matrix form ---
    K_cont_tr <- sigma^2 * Reduce(`+`,
                                  lapply(1:P, function(j)
                                    exp(-0.5 * D2_cont[j,,] / lc[j]^2)))
    K_fp_tr <- sigma^2 * exp(- D_fp_tv/(2*l_fp^2))
    K_tr    <- K_cont_tr + K_fp_tr + diag(s_n^2+1e-8, N_tr)
    K_tr_pd <- nearPD(K_tr, corr=FALSE)$mat
    L       <- chol(K_tr_pd)
    
    w       <- backsolve(L, f_bar-alpha, transpose=TRUE)
    w       <- backsolve(L, w,          transpose=FALSE)
    mu_te_s <- alpha + Kc_te %*% w
    mu_te   <- mu_te_s * y_scl + y_ctr
    
    # --- TEST‐SET posterior predictive intervals ---
    post <- extract(fit,
                    pars=c("intercept","sigma","l_cont","l_fp",
                           "sigma_noise"),
                    permuted=TRUE)
    n_draws <- length(post$sigma_noise)
    n_test  <- length(te)
    y_pred_draws <- matrix(NA, n_draws, n_test)
    for (s in seq_len(n_draws)) {
      α_s   <- post$intercept[s]
      σ_s   <- post$sigma[s]
      lc_s  <- post$l_cont[s,]
      lfp_s <- post$l_fp[s]
      σn_s  <- post$sigma_noise[s]
      # rebuild Kc_te for draw s
      X_tv_k_s <- sweep(X_tv_s, 2, lc_s, "/")
      X_te_k_s <- sweep(X_te_s, 2, lc_s, "/")
      A2_s   <- rowSums(X_te_k_s^2); B2_s <- rowSums(X_tv_k_s^2)
      D2_te_s<- outer(A2_s, rep(1,N_tr))+outer(rep(1,n_test), B2_s)-
        2*(X_te_k_s %*% t(X_tv_k_s))
      Kc_cont_te_s <- σ_s^2 * exp(-0.5 * D2_te_s)
      K_fp_te_s    <- σ_s^2 * exp(- D_fp_te/(2*lfp_s^2))
      Kc_te_s      <- Kc_cont_te_s + K_fp_te_s
      w_s <- backsolve(L, f_bar-α_s, transpose=TRUE)
      w_s <- backsolve(L, w_s,       transpose=FALSE)
      mu_te_s <- α_s + Kc_te_s %*% w_s
      y_pred_draws[s,] <- rnorm(n_test, mean=mu_te_s, sd=σn_s)
    }
    ppi_te_mat <- apply(y_pred_draws, 2,
                        quantile, probs=c(0.025,0.5,0.975))
    ppi_te[[fold]] <- data.frame(
      obs    = te,
      y_test = y_te,
      y_ctr  = y_ctr,
      y_scl  = y_scl,
      lower  = ppi_te_mat[1,],
      median = ppi_te_mat[2,],
      upper  = ppi_te_mat[3,],
      seed,fold, row.names = NULL
    )
    
    # --- metrics & summaries ---
    rmse <- sqrt(mean((mu_te - y_te)^2))
    mae  <- mean(abs(mu_te - y_te))
    seed_res[[fold]] <- data.frame(seed, fold, RMSE=rmse, MAE=mae)
    
    pars_sum <- sum_stats[ c("intercept","sigma",
                             grep("^l_cont", rownames(sum_stats), value=TRUE),
                             "l_fp","sigma_noise"), ]
    coef_res[[fold]] <- data.frame(
      parameter = rownames(pars_sum),
      mean      = pars_sum[,"50%"],
      sd        = pars_sum[,"sd"],
      lower     = pars_sum[,"2.5%"],
      upper     = pars_sum[,"97.5%"],
      seed, fold, row.names = NULL
    )
    
    # --- cleanup + break ---
    rm(fit, sum_stats, f_bar,
       D2_cont, K_cont_tr, K_fp_tr, K_tr, K_tr_pd, L,
       Kc_cont_te, K_fp_te, Kc_te,
       y_rep_draws, post, y_pred_draws)
    gc()
  }
  
  all_results[[seed]] <- bind_rows(seed_res)
  all_coefs[[seed]]   <- bind_rows(coef_res)
  all_ppi_tr[[seed]]  <- bind_rows(ppi_tr)
  all_ppi_te[[seed]]  <- bind_rows(ppi_te)
  gc()
}

end <- Sys.time()

load("GP_regression_full_results.RData")

print("Total time taken:\n")
print(end - start)

# 9) Aggregate & save
final_results <- bind_rows(all_results)
cont_names <- names(predictors[, other_cols])
final_coefs <- bind_rows(all_coefs) %>%
  mutate(parameter = as.character(parameter)) %>%
  mutate(parameter = case_when(
    parameter == "intercept"   ~ "intercept",
    parameter == "sigma"       ~ "sigma",
    parameter == "l_fp"        ~ "l_fp",
    parameter == "sigma_noise" ~ "sigma_noise",
    str_detect(parameter, "^l_cont\\[\\d+\\]$") ~ {
      idx <- str_extract(parameter, "(?<=\\[)\\d+(?=\\])") %>% as.integer()
      paste0("l_", cont_names[idx])
    },
    TRUE ~ parameter
  ))

final_ppi_tr  <- bind_rows(all_ppi_tr) %>%
  mutate(
    lower_orig  = lower  * y_scl + y_ctr,
    median_orig = median * y_scl + y_ctr,
    upper_orig  = upper  * y_scl + y_ctr,
    true        = y_train
  )
final_ppi_te <- bind_rows(all_ppi_te) %>%
  mutate(
    lower_orig  = lower  * y_scl + y_ctr,
    median_orig = median * y_scl + y_ctr,
    upper_orig  = upper  * y_scl + y_ctr,
    true        = y_test
  )

# 10) Results summary and plots

# Overall CV metrics (RMSE & MAE)
overall_metrics <- final_results %>%
  summarise(
    RMSE_mean = mean(RMSE),
    RMSE_sd   = sd(RMSE),
    MAE_mean  = mean(MAE),
    MAE_sd    = sd(MAE)
  )

# Metrics by seed (to check variability across random splits)
seed_metrics <- final_results %>%
  group_by(seed) %>%
  summarise(
    RMSE_mean = mean(RMSE),
    MAE_mean  = mean(MAE)
  )

# Posterior‐coefficient summary
coef_summary <- final_coefs %>%
  group_by(parameter) %>%
  summarise(
    mean_estimate = mean(mean),     # average of the fold‐wise medians
    sd_estimate   = mean(sd),       # average of the fold‐wise posterior sds
    lower_avg     = mean(lower),    # average lower 95% bound
    upper_avg     = mean(upper)     # average upper 95% bound
  ) %>%
  arrange(parameter)

# Test‐set PPI performance
ppi_te_summary <- final_ppi_te %>%
  summarise(
    coverage    = mean(true >= lower_orig & true <= upper_orig),
    avg_width   = mean(upper_orig - lower_orig),
    width_sd    = sd(upper_orig - lower_orig)
  )

# Train‐set PPI performance
ppi_tr_summary <- final_ppi_tr %>%
  summarise(
    coverage    = mean(true >= lower_orig & true <= upper_orig),
    avg_width   = mean(upper_orig - lower_orig),
    width_sd    = sd(upper_orig - lower_orig)
  )

#  Test‐set PPI performance by seed
ppi_te_summary_by_seed <- final_ppi_te %>%
  group_by(seed) %>%
  summarise(
    coverage    = mean(true >= lower_orig & true <= upper_orig),
    avg_width   = mean(upper_orig - lower_orig),
    width_sd    = sd(upper_orig - lower_orig)
  )

# Train‐set PPI performance by seed
ppi_tr_summary_by_seed <- final_ppi_tr %>%
  group_by(seed) %>%
  summarise(
    coverage    = mean(true >= lower_orig & true <= upper_orig),
    avg_width   = mean(upper_orig - lower_orig),
    width_sd    = sd(upper_orig - lower_orig)
  )


# Print everything
print(overall_metrics)
print(seed_metrics)
print(coef_summary)
print(ppi_te_summary)
print(ppi_tr_summary)
print(ppi_te_summary_by_seed)
print(ppi_tr_summary_by_seed)

# Plots

# Histogram of Test set PPI widths
final_ppi_te %>%
  mutate(width = upper_orig - lower_orig) %>%
  ggplot(aes(x = width)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "white") +
  labs(x = "PPI width", y = "Count",
       title = "Distribution of Test set PPI Widths") +
  theme_minimal()

# Scatter plot of predicted vs true values with 95% PPIs
final_ppi_te %>%
  ggplot(aes(x = true, y = median_orig)) +
  geom_errorbar(aes(ymin = lower_orig, ymax = upper_orig), width = 0) +
  geom_point(alpha = 0.3, size = 2, color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "True y", y = "Predicted median",
       title = "Predicted vs True with 95% PPIs") +
  theme_minimal()

# Average PPI for each observation in test set
avg_ppi <- final_ppi_te %>%
  group_by(obs, true = y_test) %>%
  summarise(
    avg_median = mean(median_orig),
    avg_lower  = mean(lower_orig),
    avg_upper  = mean(upper_orig),
    .groups = "drop"
  )

# Scatter plot of averaged predictions vs true
ggplot(avg_ppi, aes(x = true, y = avg_median)) +
  geom_errorbar(aes(ymin = avg_lower, ymax = avg_upper), width = 0) +
  geom_point(size = 2, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x     = "True y",
    y     = "Average predicted median",
    title = "True vs. Averaged Predicted Median\nwith Averaged 95% PPIs"
  ) +
  theme_minimal()

#save.image("GP_regression_full_results.RData")
