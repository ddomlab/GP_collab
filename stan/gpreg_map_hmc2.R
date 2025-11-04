# ================================================================
# MAP Cross-Validation for Gaussian Processes (Stan version)
# + Optional HMC refit per seed-fold with MAP warm start
# Flexible kernels for continuous covariates + fingerprints
# CONT_KERNEL: "none" | "rbf" | "matern32" | "matern52"
# FP_KERNEL  : "none" | "rbf" | "matern32" | "matern52" | "tanimoto" (raw D_fp as kernel)
# COMBINE    : "additive" | "multiplicative_all" | "modified_multiplicative"
# ================================================================
rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(Matrix)
  library(stringr)
  library(ggplot2)
})

# --- rstan opts (don't require rstan_options in older envs) ---
if ("rstan_options" %in% ls("package:rstan")) rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# -----------------------------
# User knobs (can be overridden via Sys.setenv from a driver)
# -----------------------------
DATA_CSV   <- Sys.getenv("DATA_CSV",   unset = "training_data_imputed.csv")
TANI_PATH  <- Sys.getenv("TANI_PATH",  unset = "tanimoto_dist_imputed.rds")

CONT_KERNEL <- Sys.getenv("CONT_KERNEL", unset = "matern32")           # none|rbf|matern32|matern52
FP_KERNEL   <- Sys.getenv("FP_KERNEL",   unset = "rbf")               # none|rbf|matern32|matern52|tanimoto
COMBINE     <- Sys.getenv("COMBINE",     unset = "multiplicative_all") # additive|multiplicative_all|modified_multiplicative
# --- custom seeds and K folds ---
SEEDS <- c(6, 13, 42, 69, 420, 1234567890, 473129)
K_FOLDS <- 5

VERBOSE     <- as.logical(Sys.getenv("VERBOSE", unset = "TRUE"))

# ---- HMC knobs ----
DO_HMC         <- as.logical(Sys.getenv("DO_HMC", unset = "TRUE"))
HMC_MAX_SPLITS <- as.integer(Sys.getenv("HMC_MAX_SPLITS", unset = "10"))  # cap total seed-fold HMC refits
CHAINS         <- as.integer(Sys.getenv("HMC_CHAINS", unset = "4"))
ITER           <- as.integer(Sys.getenv("HMC_ITER", unset = "2000"))
WARMUP         <- as.integer(Sys.getenv("HMC_WARMUP", unset = "1000"))
ADAPT_DELTA    <- as.numeric(Sys.getenv("HMC_ADAPT_DELTA", unset = "0.8"))
MAX_TD         <- as.integer(Sys.getenv("HMC_MAX_TD", unset = "10"))
HMC_PRED_NSAMP <- as.integer(Sys.getenv("HMC_PRED_NSAMP", unset = "500")) # posterior draws to use for prediction per split

# ---- Diagnostics & aggregation knobs ----
PLOT_NSPLITS  <- as.integer(Sys.getenv("PLOT_NSPLITS", unset = "3"))     # how many random splits to plot
PLOT_SHOW_ALL <- as.logical(Sys.getenv("PLOT_SHOW_ALL", unset = "FALSE"))# TRUE -> plot all seed-folds
PLOT_PARAMS   <- strsplit(Sys.getenv("PLOT_PARAMS", unset = "sigma,sigma_noise,l_fp"), ",")[[1]]

# ---- Tiny diagnostic draws (for later plotting in Rmd) ----
SAVE_DIAG_DRAWS      <- as.logical(Sys.getenv("SAVE_DIAG_DRAWS", unset = "FALSE"))
SAVE_DIAG_DRAWS_ALL  <- as.logical(Sys.getenv("SAVE_DIAG_DRAWS_ALL", unset = "FALSE")) # if TRUE, save for ALL HMC splits; else only the plotted subset
DIAG_DRAWS_MAX       <- as.integer(Sys.getenv("DIAG_DRAWS_MAX", unset = "2000"))       # cap total draws saved per split across all chains
# default DIAG_PARAMS = PLOT_PARAMS to keep things consistent, but can be overridden
DIAG_PARAMS <- strsplit(
  Sys.getenv("DIAG_PARAMS", unset = paste(PLOT_PARAMS, collapse = ",")), ","
)[[1]]

# ---- Optional seed limiter (quality-of-life) ----
MAX_SEEDS     <- as.integer(Sys.getenv("MAX_SEEDS", unset = "9999"))   # limit number of seeds processed (CV+HMC)

# --- suffix to label outputs ---
make_suffix <- function(cont, fp, comb) {
  s <- paste0("cont-", cont, "__fp-", fp, "__combine-", comb)
  gsub("[^A-Za-z0-9_\\-\\.]+", "-", s)
}
SUFFIX    <- make_suffix(CONT_KERNEL, FP_KERNEL, COMBINE)
OUT_RDATA <- file.path(paste0("/share/statistics/rchakra6/gp_reg/hmc/GP_CV_MAP_Stan_imputed_split-10_", SUFFIX, ".RData"))
HMC_RDATA <- file.path(paste0("/share/statistics/rchakra6/gp_reg/hmc/GP_HMC_Stan_imputed_split-10_",    SUFFIX, ".RData"))
PLOTS_PDF <- file.path(paste0("/share/statistics/rchakra6/gp_reg/hmc/GP_HMC_Stan_imputed_split-10_",    SUFFIX, "_plots.pdf"))

cat(sprintf(">>> Using kernels: CONT=%s  FP=%s  COMBINE=%s\n", CONT_KERNEL, FP_KERNEL, COMBINE))
cat(sprintf(">>> Output files:  CV=%s  HMC=%s  PDF=%s\n", OUT_RDATA, HMC_RDATA, PLOTS_PDF))

# -----------------------------
# Load data
# -----------------------------
dat <- read.csv(DATA_CSV, stringsAsFactors = FALSE)
if (names(dat)[1] %in% c("X","index","ID")) dat <- dat[, -1]

y <- dat$log.Rg..nm.
X <- dat %>% select(-Rg1..nm., -log.Rg..nm., -canonical_name)

fp_cols    <- grep("^Monomer_ECFP6_count_bit", names(X), value = TRUE)
other_cols <- setdiff(names(X), fp_cols)

if (CONT_KERNEL != "none" && length(other_cols) == 0) {
  stop("CONT_KERNEL != 'none' but no continuous predictors found.")
}
if (FP_KERNEL != "none" && length(fp_cols) == 0) {
  stop("FP_KERNEL != 'none' but no fingerprint predictors found.")
}

X_cont     <- if (length(other_cols) > 0) as.matrix(X[, other_cols, drop = FALSE]) else NULL
cont_names <- if (!is.null(X_cont)) colnames(X_cont) else character(0)

# -----------------------------
# Fingerprint distance (Tanimoto)
# -----------------------------
if (length(fp_cols) > 0) {
  fp_mat <- as.matrix(X[, fp_cols, drop = FALSE])
  n <- nrow(fp_mat)
  if (file.exists(TANI_PATH)) {
    D_fp_full <- readRDS(TANI_PATH)
    if (!is.matrix(D_fp_full) || any(dim(D_fp_full) != c(n,n))) {
      stop("Cached TANI_PATH dims mismatch.")
    }
  } else {
    tanimoto_similarity_counts <- function(a, b) {
      num <- sum(pmin(a, b)); den <- sum(pmax(a, b)); if (den == 0) 0 else num/den
    }
    D_fp_full <- matrix(0, n, n)
    for (i in 1:n) for (j in i:n) {
      sim <- tanimoto_similarity_counts(fp_mat[i,], fp_mat[j,])
      d   <- 1 - sim
      D_fp_full[i,j] <- d; D_fp_full[j,i] <- d
    }
    saveRDS(D_fp_full, TANI_PATH)
  }
} else {
  D_fp_full <- NULL
}

# -----------------------------
# Helpers
# -----------------------------
chol_retry <- function(M, base = 1e-8, tries = 8) {
  for (t in 0:(tries-1)) {
    jitter <- base * (10^t)
    MM <- M + diag(jitter, nrow(M))
    out <- tryCatch(chol(MM), error = function(e) NULL)
    if (!is.null(out)) return(list(L = out, jitter = jitter))
  }
  stop("Cholesky failed after jitter escalation.")
}

build_D2_cont <- function(Xs) {
  P <- ncol(Xs); N <- nrow(Xs)
  out <- array(0, dim = c(P, N, N))
  for (j in 1:P) out[j,,] <- as.matrix(dist(Xs[, j]))^2
  out
}
build_D2_cross <- function(xnew, xtr) {
  A2 <- rowSums(xnew^2); B2 <- rowSums(xtr^2)
  outer(A2, rep(1, nrow(xtr))) + outer(rep(1, nrow(xnew)), B2) - 2 * (xnew %*% t(xtr))
}

make_cv_splits <- function(n, k, seed) {
  set.seed(seed)
  idx_all <- sample.int(n)  # random permutation of rows 1..n

  folds <- split(idx_all, rep(1:k, length.out = n))

  train_list <- vector("list", k)
  test_list  <- vector("list", k)

  for (f in seq_len(k)) {
    test_idx  <- folds[[f]]
    train_idx <- setdiff(idx_all, test_idx)

    train_list[[f]] <- train_idx
    test_list[[f]]  <- test_idx
  }

  list(train = train_list, test = test_list)
}


# shape guards
to_sq <- function(M, n) {
  if (is.null(dim(M))) return(matrix(as.numeric(M), n, n))
  d <- dim(M)
  if (length(d) == 2 && d[1] == n && d[2] == n) return(as.matrix(M))
  if (prod(d) == n*n) return(matrix(as.numeric(M), n, n))
  stop(sprintf("to_sq: cannot coerce to %dx%d (got %s)", n, n, paste(d, collapse="x")))
}
to_rect <- function(M, nrow, ncol) {
  if (is.null(dim(M))) return(matrix(as.numeric(M), nrow, ncol))
  d <- dim(M)
  if (length(d) == 2 && d[1] == nrow && d[2] == ncol) return(as.matrix(M))
  if (prod(d) == nrow*ncol) return(matrix(as.numeric(M), nrow, ncol))
  stop(sprintf("to_rect: cannot coerce to %dx%d (got %s)", nrow, ncol, paste(d, collapse="x")))
}

# Continuous kernels (R helpers, no amplitude)
k_rbf_from_D2 <- function(D2, l) exp(-0.5 * D2 / (l^2))
k_m32_from_D  <- function(D, l) { R <- pmax(0, D)/l; (1 + sqrt(3)*R) * exp(-sqrt(3)*R) }
k_m52_from_D  <- function(D, l) { R <- pmax(0, D)/l; (1 + sqrt(5)*R + 5*R^2/3) * exp(-sqrt(5)*R) }

build_Kc_perdim_R <- function(kind, D2_cont, lc) {
  P <- dim(D2_cont)[1]
  Ks <- vector("list", P)
  if (kind == "rbf") {
    for (j in 1:P) Ks[[j]] <- k_rbf_from_D2(D2_cont[j,,], lc[j])
  } else if (kind == "matern32") {
    for (j in 1:P) Ks[[j]] <- k_m32_from_D(sqrt(D2_cont[j,,]), lc[j])
  } else if (kind == "matern52") {
    for (j in 1:P) Ks[[j]] <- k_m52_from_D(sqrt(D2_cont[j,,]), lc[j])
  } else stop("Unknown CONT_KERNEL in R helper")
  Ks
}

# FP kernel (R helper; raw D_fp if "tanimoto")
k_fp_from <- function(kind, D_fp, l_fp) {
  if (kind == "rbf")        exp( -(D_fp^2) / (2*l_fp^2) )
  else if (kind == "matern32") { R <- D_fp / l_fp; (1 + sqrt(3)*R) * exp(-sqrt(3)*R) }
  else if (kind == "matern52") { R <- D_fp / l_fp; (1 + sqrt(5)*R + 5*R^2/3) * exp(-sqrt(5)*R) }
  else if (kind == "tanimoto") D_fp  # raw distance as "kernel" (non-PSD!)
  else stop("Unknown FP_KERNEL in R helper")
}

# -----------------------------
# Stan model (handles "none" via codes)
# -----------------------------
stan_code <- '
data {
  int<lower=1> N;
  int<lower=1> P;                      // If CONT_CODE==0 (none), pass P=1 with dummy matrices
  vector[N] y;
  matrix[N,N] D_fp;                    // If FP_CODE==0 (none), pass zeros
  array[P] matrix[N,N] D2_cont;        // If CONT_CODE==0, provide one zero matrix
  int<lower=0,upper=3> CONT_CODE;      // 0 none, 1 rbf, 2 m32, 3 m52
  int<lower=0,upper=4> FP_CODE;        // 0 none, 1 rbf, 2 m32, 3 m52, 4 tanimoto raw
  int<lower=1,upper=3> COMBINE_CODE;   // 1 additive, 2 mult_all, 3 mod_mult
}
parameters {
  real<lower=0> sigma;
  vector<lower=0>[P] l_cont;           // ignored if CONT_CODE==0
  real<lower=0> l_fp;                  // ignored if FP_CODE==0
  real<lower=0> sigma_noise;
}
transformed parameters {
  matrix[N,N] Kc;
  matrix[N,N] Kf;
  matrix[N,N] K_core;
  matrix[N,N] K;

  // ----- Kc (continuous, no amplitude) -----
  if (CONT_CODE == 0) {
    Kc = rep_matrix(0.0, N, N);
  } else if (COMBINE_CODE == 2) {
    // multiplicative_all: product over dims
    Kc = rep_matrix(1.0, N, N);
    for (j in 1:P) {
      matrix[N,N] Kj;
      if (CONT_CODE == 1) {
        Kj = exp(-0.5 * D2_cont[j] ./ square(l_cont[j]));
      } else if (CONT_CODE == 2) {
        matrix[N,N] R = sqrt(D2_cont[j]) ./ l_cont[j];
        Kj = (1 + sqrt(3) * R) .* exp(-sqrt(3) * R);
      } else { // m52
        matrix[N,N] R = sqrt(D2_cont[j]) ./ l_cont[j];
        Kj = (1 + sqrt(5) * R + (5.0/3.0) * square(R)) .* exp(-sqrt(5) * R);
      }
      Kc = Kc .* Kj;
    }
  } else {
    // additive or modified_multiplicative (average later)
    Kc = rep_matrix(0.0, N, N);
    for (j in 1:P) {
      matrix[N,N] Kj;
      if (CONT_CODE == 1) {
        Kj = exp(-0.5 * D2_cont[j] ./ square(l_cont[j]));
      } else if (CONT_CODE == 2) {
        matrix[N,N] R = sqrt(D2_cont[j]) ./ l_cont[j];
        Kj = (1 + sqrt(3) * R) .* exp(-sqrt(3) * R);
      } else { // m52
        matrix[N,N] R = sqrt(D2_cont[j]) ./ l_cont[j];
        Kj = (1 + sqrt(5) * R + (5.0/3.0) * square(R)) .* exp(-sqrt(5) * R);
      }
      Kc += Kj;
    }
    if (COMBINE_CODE == 3) Kc /= P; // modified_multiplicative: average within cont
  }

  // ----- Kf (fingerprints, no amplitude) -----
  if (FP_CODE == 0) {
    Kf = rep_matrix(0.0, N, N);
  } else if (FP_CODE == 1) { // rbf
    Kf = exp( -square(D_fp) / (2 * square(l_fp)) );
  } else if (FP_CODE == 2) { // m32
    matrix[N,N] R = D_fp ./ l_fp;
    Kf = (1 + sqrt(3) * R) .* exp(-sqrt(3) * R);
  } else if (FP_CODE == 3) { // m52
    matrix[N,N] R = D_fp ./ l_fp;
    Kf = (1 + sqrt(5) * R + (5.0/3.0) * square(R)) .* exp(-sqrt(5) * R);
  } else { // FP_CODE == 4, tanimoto raw distance as "kernel"
    Kf = D_fp;
  }

  // ----- Combine core (no amplitude) -----
  if (COMBINE_CODE == 1) {
    K_core = Kc + Kf;                        // additive
  } else {
    if (CONT_CODE != 0 && FP_CODE != 0)      K_core = Kc .* Kf;  // multiplicative flavors
    else if (CONT_CODE != 0)                 K_core = Kc;        // only cont present
    else                                      K_core = Kf;        // only fp present
  }

  // Scale and jitter
  K = square(sigma) * K_core;
  K += diag_matrix(rep_vector(1e-6, N));
}
model {
  // Priors (l_* sampled even if unused; harmless)
  l_cont ~ inv_gamma(5,5);
  l_fp   ~ inv_gamma(5,5);
  sigma  ~ normal(0,1);
  sigma_noise ~ normal(0,1);

  {
    matrix[N,N] Ky = K + diag_matrix(rep_vector(square(sigma_noise), N));
    y ~ multi_normal(rep_vector(0, N), Ky);
  }
}
'

stan_mod <- stan_model(model_code = stan_code)

# -----------------------------
# Codes for kernels/combine (for Stan)
# -----------------------------
CONT_CODE <- switch(CONT_KERNEL,
                    "none"=0L, "rbf"=1L, "matern32"=2L, "matern52"=3L, stop("Bad CONT_KERNEL"))
FP_CODE   <- switch(FP_KERNEL,
                    "none"=0L, "rbf"=1L, "matern32"=2L, "matern52"=3L, "tanimoto"=4L, stop("Bad FP_KERNEL"))
COMBINE_CODE <- switch(COMBINE,
                       "additive"=1L, "multiplicative_all"=2L, "modified_multiplicative"=3L, stop("Bad COMBINE"))

if (CONT_CODE == 0L && FP_CODE == 0L) stop("Both CONT_KERNEL and FP_KERNEL are 'none' (degenerate).")

# -----------------------------
# CV Loop (MAP via optimizing) + Predictions (RMSE/MAE)
# -----------------------------
results   <- list()
fold_info <- list()
if (length(SEEDS) > MAX_SEEDS) {
  SEEDS <- SEEDS[seq_len(MAX_SEEDS)]
}


for (seed in SEEDS) {

  splits <- make_cv_splits(n = nrow(dat), k = K_FOLDS, seed = seed)
  train_idx <- splits$train  # already 1-based
  test_idx  <- splits$test   # already 1-based

  for (f in seq_len(K_FOLDS)) {
    tr <- train_idx[[f]]
    te <- test_idx[[f]]
    Ntr <- length(tr); Nte <- length(te)
    
    # Standardize (train-only)
    X_tr <- if (!is.null(X_cont)) X_cont[tr,,drop=FALSE] else NULL
    X_te <- if (!is.null(X_cont)) X_cont[te,,drop=FALSE] else NULL
    y_tr <- y[tr]; y_te <- y[te]
    
    if (!is.null(X_tr)) {
      ctr <- colMeans(X_tr)
      scl <- apply(X_tr, 2, sd); scl[scl==0] <- 1
      X_tr_s <- sweep(sweep(X_tr, 2, ctr, "-"), 2, scl, "/")
      X_te_s <- sweep(sweep(X_te, 2, ctr, "-"), 2, scl, "/")
    } else {
      X_tr_s <- NULL; X_te_s <- NULL
    }
    
    y_ctr <- mean(y_tr); y_scl <- sd(y_tr); if (y_scl == 0) y_scl <- 1
    y_tr_s <- as.numeric((y_tr - y_ctr)/y_scl)
    
    # Distances and dummies
    if (!is.null(X_tr_s) && CONT_KERNEL != "none") {
      D2_cont_tr <- build_D2_cont(X_tr_s)
      P_tr <- dim(D2_cont_tr)[1]
    } else {
      P_tr <- 1L
      D2_cont_tr <- array(0, dim = c(1L, Ntr, Ntr))
    }
    
    if (!is.null(D_fp_full) && FP_KERNEL != "none") {
      D_fp_tr <- D_fp_full[tr, tr, drop=FALSE]
      D_fp_te <- D_fp_full[te, tr, drop=FALSE]
      D_fp_tt <- D_fp_full[te, te, drop=FALSE]
    } else {
      D_fp_tr <- matrix(0, Ntr, Ntr)
      D_fp_te <- matrix(0, Nte, Ntr)
      D_fp_tt <- matrix(0, Nte, Nte)
    }
    
    # Stan data
    stan_data <- list(
      N = Ntr,
      P = P_tr,
      y = y_tr_s,
      D_fp = D_fp_tr,
      D2_cont = D2_cont_tr,
      CONT_CODE = CONT_CODE,
      FP_CODE = FP_CODE,
      COMBINE_CODE = COMBINE_CODE
    )
    
    # MAP via optimizing
    init_fun <- function() {
      list(sigma = 0.5,
           l_cont = rep(1.0, P_tr),
           l_fp = 0.5,
           sigma_noise = 0.2)
    }
    fit <- rstan::optimizing(
      stan_mod, data = stan_data, as_vector = FALSE,
      init = init_fun, seed = as.integer(seed)
    )
    
    par <- fit$par
    sigma_hat <- as.numeric(par$sigma)
    noise_hat <- as.numeric(par$sigma_noise)
    lc_hat    <- as.numeric(par$l_cont)
    lfp_hat   <- as.numeric(par$l_fp)
    
    # ------------ Predict on test (in R) ------------
    # K_tr
    if (!is.null(X_tr_s) && CONT_KERNEL != "none") {
      Ks_tr <- build_Kc_perdim_R(CONT_KERNEL, D2_cont_tr, lc_hat[seq_len(ifelse(CONT_KERNEL=="none",1L,P_tr))])
      if (COMBINE == "multiplicative_all") {
        Kc_tr <- Reduce(`*`, Ks_tr)
      } else if (COMBINE == "modified_multiplicative") {
        Kc_tr <- Reduce(`+`, Ks_tr) / length(Ks_tr)
      } else {
        Kc_tr <- Reduce(`+`, Ks_tr)
      }
      Kc_tr <- to_sq(Kc_tr, Ntr)
    } else {
      Kc_tr <- matrix(0, Ntr, Ntr)
    }
    
    if (!is.null(D_fp_full) && FP_KERNEL != "none") {
      Kf_tr <- k_fp_from(FP_KERNEL, D_fp_tr, lfp_hat)
      Kf_tr <- to_sq(Kf_tr, Ntr)
    } else {
      Kf_tr <- matrix(0, Ntr, Ntr)
    }
    
    has_cont <- (CONT_KERNEL != "none")
    has_fp   <- (FP_KERNEL   != "none")
    if (COMBINE == "additive") {
      K_tr_core <- Kc_tr + Kf_tr
    } else {
      if (has_cont && has_fp)       K_tr_core <- Kc_tr * Kf_tr
      else if (has_cont && !has_fp) K_tr_core <- Kc_tr
      else                          K_tr_core <- Kf_tr
    }
    K_tr <- (sigma_hat^2) * to_sq(K_tr_core, Ntr)
    
    Ky <- K_tr + diag(noise_hat^2, nrow = Ntr, ncol = Ntr)
    Ky <- to_sq(Ky, Ntr)
    ch <- chol_retry(Ky, base = 1e-8, tries = 8)
    L <- ch$L
    alpha <- backsolve(L, backsolve(L, y_tr_s, transpose = TRUE), transpose = FALSE)
    
    # Cross K_te
    if (!is.null(X_tr_s) && CONT_KERNEL != "none") {
      Klist_cross <- vector("list", ifelse(CONT_KERNEL=="none", 0L, P_tr))
      for (j in seq_along(Klist_cross)) {
        D2_te_j <- build_D2_cross(X_te_s[,j,drop=FALSE], X_tr_s[,j,drop=FALSE])
        if (CONT_KERNEL == "rbf")           Klist_cross[[j]] <- k_rbf_from_D2(D2_te_j, lc_hat[j])
        else if (CONT_KERNEL == "matern32") Klist_cross[[j]] <- k_m32_from_D(sqrt(pmax(0, D2_te_j)), lc_hat[j])
        else                                Klist_cross[[j]] <- k_m52_from_D(sqrt(pmax(0, D2_te_j)), lc_hat[j])
      }
      if (COMBINE == "multiplicative_all") {
        Kc_te <- Reduce(`*`, Klist_cross)
      } else if (COMBINE == "modified_multiplicative") {
        Kc_te <- Reduce(`+`, Klist_cross) / length(Klist_cross)
      } else {
        Kc_te <- Reduce(`+`, Klist_cross)
      }
      Kc_te <- to_rect(Kc_te, Nte, Ntr)
    } else {
      Kc_te <- matrix(0, nrow = Nte, ncol = Ntr)
    }
    
    if (!is.null(D_fp_full) && FP_KERNEL != "none") {
      Kf_te <- k_fp_from(FP_KERNEL, D_fp_te, lfp_hat)
      Kf_te <- to_rect(Kf_te, Nte, Ntr)
    } else {
      Kf_te <- matrix(0, nrow = Nte, ncol = Ntr)
    }
    
    if (COMBINE == "additive") {
      K_te_core <- Kc_te + Kf_te
    } else {
      if (has_cont && has_fp)       K_te_core <- Kc_te * Kf_te
      else if (has_cont && !has_fp) K_te_core <- Kc_te
      else                          K_te_core <- Kf_te
    }
    K_te <- (sigma_hat^2) * to_rect(K_te_core, Nte, Ntr)
    
    # Predictive mean (standardized)
    mu_te_s <- as.numeric(K_te %*% alpha)
    
    # Test self-covariance K_self
    if (!is.null(X_tr_s) && CONT_KERNEL != "none") {
      D2_cont_te <- build_D2_cont(X_te_s)
      Ks_self <- build_Kc_perdim_R(CONT_KERNEL, D2_cont_te, lc_hat[seq_len(ifelse(CONT_KERNEL=="none",1L,P_tr))])
      if (COMBINE == "multiplicative_all") {
        Kc_self <- Reduce(`*`, Ks_self)
      } else if (COMBINE == "modified_multiplicative") {
        Kc_self <- Reduce(`+`, Ks_self) / length(Ks_self)
      } else {
        Kc_self <- Reduce(`+`, Ks_self)
      }
      Kc_self <- to_sq(Kc_self, Nte)
    } else {
      Kc_self <- matrix(0, Nte, Nte)
    }
    if (!is.null(D_fp_full) && FP_KERNEL != "none") {
      Kf_self <- k_fp_from(FP_KERNEL, D_fp_tt, lfp_hat)
      Kf_self <- to_sq(Kf_self, Nte)
    } else {
      Kf_self <- matrix(0, Nte, Nte)
    }
    
    if (COMBINE == "additive") {
      K_self_core <- Kc_self + Kf_self
    } else {
      if (has_cont && has_fp)       K_self_core <- Kc_self * Kf_self
      else if (has_cont && !has_fp) K_self_core <- Kc_self
      else                          K_self_core <- Kf_self
    }
    K_self <- (sigma_hat^2) * to_sq(K_self_core, Nte)
    
    V <- backsolve(L, t(K_te), transpose = TRUE)
    var_te_s <- pmax(0, diag(K_self) - colSums(V^2)) + (noise_hat^2)
    
    # Back-transform
    mu_te <- mu_te_s * y_scl + y_ctr
    
    rmse <- sqrt(mean((mu_te - y_te)^2))
    mae  <- mean(abs(mu_te - y_te))
    
    if (VERBOSE) {
      cat(sprintf("[seed=%s fold=%d] RMSE=%.4f MAE=%.4f  (jitter=%.2e)\n",
                  seed, f, rmse, mae, ch$jitter))
      if (FP_KERNEL == "tanimoto") {
        cat("  NOTE: FP 'tanimoto' uses raw D_fp as kernel (non-PSD in general). HMC may be skipped.\n")
      }
    }
    
    results[[length(results)+1]] <- data.frame(
      seed = seed, fold = f, RMSE = rmse, MAE = mae,
      CONT_KERNEL = CONT_KERNEL, FP_KERNEL = FP_KERNEL, COMBINE = COMBINE,
      stringsAsFactors = FALSE
    )
    fold_info[[length(fold_info)+1]] <- list(
      seed = seed, fold = f,
      tr = tr, te = te,
      # keep training stats for this split
      ctr = if (!is.null(X_tr)) ctr else NULL,
      scl = if (!is.null(X_tr)) scl else NULL,
      y_ctr = y_ctr, y_scl = y_scl,
      # distances for FP to avoid recompute
      D_fp_tr = D_fp_tr, D_fp_te = D_fp_te, D_fp_tt = D_fp_tt,
      # store MAP hypers for warm start
      theta_map = list(sigma = sigma_hat, sigma_noise = noise_hat,
                       lc = lc_hat, l_fp = lfp_hat),
      jitter = ch$jitter,
      CONT_KERNEL = CONT_KERNEL, FP_KERNEL = FP_KERNEL, COMBINE = COMBINE
    )
    
    rm(fit, L, Ky, alpha, K_tr, Kc_tr, Kf_tr, K_te, Kc_te, Kf_te, K_self, Kc_self, Kf_self)
    gc()
  }
}

final_results <- dplyr::bind_rows(results)

if (VERBOSE) {
  cat("\nCV summary (all folds) for ", SUFFIX, ":\n", sep = "")
  print(final_results %>%
          summarise(RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
                    MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE)))
}

# Include explicit labels in the RData
labels <- list(
  CONT_KERNEL = CONT_KERNEL,
  FP_KERNEL   = FP_KERNEL,
  COMBINE     = COMBINE,
  SUFFIX      = SUFFIX
)

save(list = c("final_results","fold_info","labels",
              "CONT_KERNEL","FP_KERNEL","COMBINE",
              "DATA_CSV","TANI_PATH","SUFFIX",
              "SEEDS","K_FOLDS","MAX_SEEDS"),
     file = OUT_RDATA)


cat("\nSaved CV MAP (Stan) results to: ", OUT_RDATA, "\n")

# =======================================================================
#                     HMC REFIT PER SEED–FOLD (optional)
# =======================================================================
if (DO_HMC) {
  if (FP_CODE == 4L) {
    message("HMC disabled because FP_KERNEL == 'tanimoto' (raw distance is not PSD). Skipping HMC.")
  } else {
    # Build a work list of splits
    splits <- lapply(fold_info, function(fi) list(
      seed = fi$seed, fold = fi$fold, tr = fi$tr, te = fi$te,
      ctr = fi$ctr, scl = fi$scl, y_ctr = fi$y_ctr, y_scl = fi$y_scl,
      D_fp_tr = fi$D_fp_tr, D_fp_te = fi$D_fp_te, D_fp_tt = fi$D_fp_tt,
      theta_map = fi$theta_map
    ))
    # Cap per knob
    if (length(splits) > HMC_MAX_SPLITS) splits <- splits[seq_len(HMC_MAX_SPLITS)]
    
    # Decide which splits to plot diagnostics for
    set.seed(42)
    plot_indices <- if (PLOT_SHOW_ALL || length(splits) <= PLOT_NSPLITS) {
      seq_along(splits)
    } else {
      sort(sample.int(length(splits), PLOT_NSPLITS))
    }
    
    # Storage
    hmc_param_summaries <- list()
    hmc_pred_summaries  <- list()
    hmc_rmse_rows       <- list()
    
    # Convenience: pretty parameter names for l_cont
    pretty_lcont <- function(p_names) {
      out <- p_names
      idx <- stringr::str_match(p_names, "^l_cont\\[(\\d+)\\]$")[,2]
      for (k in which(!is.na(idx))) {
        j <- as.integer(idx[k])
        if (j >= 1 && j <= length(cont_names)) out[k] <- paste0("l_", cont_names[j])
      }
      out
    }
    
    # Open plots PDF
    pdf(PLOTS_PDF, width = 7, height = 5)
    on.exit(try(dev.off(), silent = TRUE), add = TRUE)
    
    have_bayesplot <- requireNamespace("bayesplot", quietly = TRUE)
    if (have_bayesplot) library(bayesplot)
    have_gridExtra  <- requireNamespace("gridExtra", quietly = TRUE)
    if (have_gridExtra) library(gridExtra)
    
    # HMC per split
    for (sidx in seq_along(splits)) {
      sp <- splits[[sidx]]
      tr <- sp$tr; te <- sp$te
      
      X_tr <- if (!is.null(X_cont)) X_cont[tr,,drop=FALSE] else NULL
      X_te <- if (!is.null(X_cont) && length(te)>0) X_cont[te,,drop=FALSE] else NULL
      y_tr <- y[tr]; y_te <- if (length(te)>0) y[te] else numeric(0)
      
      # standardize with stored ctr/scl (from the MAP phase), or recompute if NULL
      if (!is.null(X_tr)) {
        ctr <- if (!is.null(sp$ctr)) sp$ctr else colMeans(X_tr)
        scl <- if (!is.null(sp$scl)) sp$scl else { tmp <- apply(X_tr, 2, sd); tmp[tmp==0] <- 1; tmp }
        X_tr_s <- sweep(sweep(X_tr, 2, ctr, "-"), 2, scl, "/")
        X_te_s <- if (length(te)>0) sweep(sweep(X_te, 2, ctr, "-"), 2, scl, "/") else NULL
      } else {
        X_tr_s <- NULL; X_te_s <- NULL
      }
      y_ctr <- sp$y_ctr; y_scl <- sp$y_scl
      y_tr_s <- (y_tr - y_ctr)/y_scl
      
      # Distances
      if (!is.null(X_tr_s) && CONT_KERNEL != "none") {
        D2_cont_tr <- build_D2_cont(X_tr_s)
        P_tr <- dim(D2_cont_tr)[1]
      } else {
        P_tr <- 1L
        D2_cont_tr <- array(0, dim = c(1L, length(tr), length(tr)))
      }
      D_fp_tr <- if (!is.null(D_fp_full) && FP_KERNEL != "none") sp$D_fp_tr else matrix(0, length(tr), length(tr))
      D_fp_te <- if (!is.null(D_fp_full) && FP_KERNEL != "none" && length(te)>0) sp$D_fp_te else matrix(0, length(te), length(tr))
      D_fp_tt <- if (!is.null(D_fp_full) && FP_KERNEL != "none" && length(te)>0) sp$D_fp_tt else matrix(0, length(te), length(te))
      
      stan_data <- list(
        N = length(tr),
        P = P_tr,
        y = as.numeric(y_tr_s),
        D_fp = D_fp_tr,
        D2_cont = D2_cont_tr,
        CONT_CODE = CONT_CODE,
        FP_CODE = FP_CODE,
        COMBINE_CODE = COMBINE_CODE
      )
      
      # Warm-start at MAP for this split (if available), otherwise small positives
      mp <- sp$theta_map
      init_fun <- function() list(
        sigma       = if (!is.null(mp$sigma))       as.numeric(mp$sigma)       else 0.5,
        l_cont      = if (!is.null(mp$lc))          as.numeric(mp$lc)          else rep(1.0, P_tr),
        l_fp        = if (!is.null(mp$l_fp))        as.numeric(mp$l_fp)        else 0.5,
        sigma_noise = if (!is.null(mp$sigma_noise)) as.numeric(mp$sigma_noise) else 0.2
      )
      
      # Run HMC
      fit_hmc <- rstan::sampling(
        stan_mod,
        data    = stan_data,
        chains  = CHAINS,
        iter    = ITER,
        warmup  = WARMUP,
        control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TD),
        seed    = 1000 + sidx,
        init    = init_fun,
        init_r  = 0.02,
        pars    = c("sigma","l_cont","l_fp","sigma_noise"),
        include = TRUE
      )
      
      # -----------------------------
      # Tiny diag draws saver (for Rmd plotting later)
      # -----------------------------
      if (SAVE_DIAG_DRAWS && (SAVE_DIAG_DRAWS_ALL || (sidx %in% plot_indices))) {
        # figure out which parameter names actually exist in this fit
        aa <- rstan::extract(fit_hmc, permuted = FALSE, inc_warmup = FALSE)  # iter x chain x param
        all_pars <- dimnames(aa)[[3]]
        
        # expand DIAG_PARAMS to actual parameter names (support simple prefixes, e.g., "l_cont")
        expand_params <- function(pat, pool) {
          if (grepl("\\[", pat)) {
            intersect(pat, pool)
          } else {
            # treat as prefix (exact for scalar names like sigma; prefix for vectors like l_cont)
            grep(paste0("^", gsub("\\.", "\\\\.", pat)), pool, value = TRUE)
          }
        }
        diag_pars <- unique(unlist(lapply(DIAG_PARAMS, expand_params, pool = all_pars)))
        diag_pars <- intersect(diag_pars, all_pars)
        if (length(diag_pars) > 0) {
          a <- rstan::extract(fit_hmc, pars = diag_pars, permuted = FALSE, inc_warmup = FALSE)  # iter x chain x selected pars
          iters <- dim(a)[1]; chains <- dim(a)[2]
          
          # thin down to keep total <= DIAG_DRAWS_MAX across all chains
          total <- iters * chains
          step  <- max(1, floor(total / max(1, DIAG_DRAWS_MAX)))
          sel_iter <- seq(1, iters, by = step)
          a_small <- a[sel_iter, , , drop = FALSE]  # keep dimnames
          
          # create/append global list
          if (!exists("diag_draws", inherits = FALSE)) diag_draws <- list()
          diag_draws[[length(diag_draws) + 1]] <- list(
            meta  = list(seed = as.character(sp$seed), fold = sp$fold, split = sidx,
                         pars = diag_pars,
                         n_iter_kept = length(sel_iter), n_chains = chains,
                         step = step),
            draws = a_small  # iter x chain x param (with dimnames)
          )
        }
      }
      
      
      # Summaries of parameters
      sum_tab <- summary(
        fit_hmc,
        pars  = c("sigma","l_cont","l_fp","sigma_noise"),
        probs = c(0.025, 0.5, 0.975)
      )$summary
      sum_df <- as.data.frame(sum_tab)
      sum_df$parameter <- rownames(sum_df)
      sum_df$parameter <- pretty_lcont(sum_df$parameter)
      sum_df$seed  <- as.character(sp$seed)
      sum_df$fold  <- sp$fold
      sum_df$split <- sidx
      hmc_param_summaries[[length(hmc_param_summaries)+1]] <- sum_df
      
      # Posterior draws for predictions
      post <- rstan::extract(fit_hmc, pars = c("sigma","l_cont","l_fp","sigma_noise"), permuted = TRUE)
      S_all <- length(post$sigma)
      set.seed(1234 + sidx)
      sel <- if (S_all <= HMC_PRED_NSAMP) 1:S_all else sort(sample.int(S_all, HMC_PRED_NSAMP))
      
      # Predictive with posterior draws
      Ntr <- length(tr); Nte <- length(te)
      Pdim <- if (!is.null(X_tr_s)) ncol(X_tr_s) else 0L
      
      predict_one <- function(sid) {
        sigma       <- post$sigma[sid]
        lc          <- if (Pdim>0) post$l_cont[sid, 1:Pdim, drop = TRUE] else numeric(0)
        l_fp        <- post$l_fp[sid]
        sigma_noise <- post$sigma_noise[sid]
        
        # ----- Continuous Kc_tr (per draw) with correct combine rule -----
        if (Pdim > 0 && CONT_KERNEL != "none") {
          Klist_tr <- vector("list", Pdim)
          if (CONT_KERNEL == "rbf") {
            for (j in seq_len(Pdim)) Klist_tr[[j]] <- exp(-0.5 * D2_cont_tr[j,,] / (lc[j]^2))
          } else if (CONT_KERNEL == "matern32") {
            for (j in seq_len(Pdim)) {
              R <- sqrt(D2_cont_tr[j,,]) / lc[j]
              Klist_tr[[j]] <- (1 + sqrt(3)*R) * exp(-sqrt(3)*R)
            }
          } else { # m52
            for (j in seq_len(Pdim)) {
              R <- sqrt(D2_cont_tr[j,,]) / lc[j]
              Klist_tr[[j]] <- (1 + sqrt(5)*R + 5*R^2/3) * exp(-sqrt(5)*R)
            }
          }
          if (COMBINE == "multiplicative_all")      Kc_tr <- Reduce(`*`, Klist_tr)
          else if (COMBINE == "modified_multiplicative") Kc_tr <- Reduce(`+`, Klist_tr) / length(Klist_tr)
          else                                        Kc_tr <- Reduce(`+`, Klist_tr)
        } else {
          Kc_tr <- matrix(0, Ntr, Ntr)
        }
        
        # ----- Fingerprint Kf_tr (per draw) -----
        Kf_tr <- if (FP_KERNEL != "none") {
          if (FP_KERNEL == "rbf") {
            exp( -(D_fp_tr^2) / (2*l_fp^2) )
          } else if (FP_KERNEL == "matern32") {
            R <- D_fp_tr / l_fp; (1 + sqrt(3)*R) * exp(-sqrt(3)*R)
          } else if (FP_KERNEL == "matern52") {
            R <- D_fp_tr / l_fp; (1 + sqrt(5)*R + 5*R^2/3) * exp(-sqrt(5)*R)
          } else {
            D_fp_tr
          }
        } else matrix(0, Ntr, Ntr)
        
        # ----- Combine and factor in amplitude (train) -----
        K_core_tr <- if (COMBINE == "additive") (Kc_tr + Kf_tr) else {
          if (CONT_KERNEL != "none" && FP_KERNEL != "none") Kc_tr * Kf_tr
          else if (CONT_KERNEL != "none") Kc_tr else Kf_tr
        }
        K_tr <- (sigma^2) * K_core_tr
        Ky   <- K_tr + diag(sigma_noise^2, Ntr)
        
        L <- chol_retry(Ky, base = 1e-6, tries = 6)$L
        alpha <- backsolve(L, backsolve(L, y_tr_s, transpose = TRUE), transpose = FALSE)
        
        if (Nte == 0) return(matrix(0, 0, 2))
        
        # ----- Continuous Kc_te (per draw) with correct combine rule -----
        if (Pdim > 0 && CONT_KERNEL != "none") {
          Klist_te <- vector("list", Pdim)
          for (j in seq_len(Pdim)) {
            D2_te_j <- build_D2_cross(X_te_s[,j,drop=FALSE], X_tr_s[,j,drop=FALSE])
            if (CONT_KERNEL == "rbf") {
              Klist_te[[j]] <- exp(-0.5 * D2_te_j / (lc[j]^2))
            } else if (CONT_KERNEL == "matern32") {
              R <- sqrt(pmax(0, D2_te_j)) / lc[j]
              Klist_te[[j]] <- (1 + sqrt(3)*R) * exp(-sqrt(3)*R)
            } else {
              R <- sqrt(pmax(0, D2_te_j)) / lc[j]
              Klist_te[[j]] <- (1 + sqrt(5)*R + 5*R^2/3) * exp(-sqrt(5)*R)
            }
          }
          if (COMBINE == "multiplicative_all")      Kc_te <- Reduce(`*`, Klist_te)
          else if (COMBINE == "modified_multiplicative") Kc_te <- Reduce(`+`, Klist_te) / length(Klist_te)
          else                                        Kc_te <- Reduce(`+`, Klist_te)
        } else {
          Kc_te <- matrix(0, Nte, Ntr)
        }
        
        # ----- Fingerprint Kf_te (per draw) -----
        Kf_te <- if (FP_KERNEL != "none") {
          if (FP_KERNEL == "rbf") {
            exp( -(D_fp_te^2) / (2*l_fp^2) )
          } else if (FP_KERNEL == "matern32") {
            R <- D_fp_te / l_fp; (1 + sqrt(3)*R) * exp(-sqrt(3)*R)
          } else if (FP_KERNEL == "matern52") {
            R <- D_fp_te / l_fp; (1 + sqrt(5)*R + 5*R^2/3) * exp(-sqrt(5)*R)
          } else {
            D_fp_te
          }
        } else matrix(0, Nte, Ntr)
        
        # ----- Combine & predict mean (standardized) -----
        K_te_core <- if (COMBINE == "additive") (Kc_te + Kf_te) else {
          if (CONT_KERNEL != "none" && FP_KERNEL != "none") Kc_te * Kf_te
          else if (CONT_KERNEL != "none") Kc_te else Kf_te
        }
        K_te <- (sigma^2) * K_te_core
        mu_te_s <- as.numeric(K_te %*% alpha)
        
        # ----- Continuous Kc_self (per draw) with correct combine rule -----
        if (Pdim > 0 && CONT_KERNEL != "none") {
          D2_cont_te <- build_D2_cont(X_te_s)
          Klist_self <- vector("list", Pdim)
          if (CONT_KERNEL == "rbf") {
            for (j in seq_len(Pdim)) Klist_self[[j]] <- exp(-0.5 * D2_cont_te[j,,] / (lc[j]^2))
          } else if (CONT_KERNEL == "matern32") {
            for (j in seq_len(Pdim)) {
              R <- sqrt(D2_cont_te[j,,]) / lc[j]
              Klist_self[[j]] <- (1 + sqrt(3)*R) * exp(-sqrt(3)*R)
            }
          } else { # m52
            for (j in seq_len(Pdim)) {
              R <- sqrt(D2_cont_te[j,,]) / lc[j]
              Klist_self[[j]] <- (1 + sqrt(5)*R + 5*R^2/3) * exp(-sqrt(5)*R)
            }
          }
          if (COMBINE == "multiplicative_all")      Kc_self <- Reduce(`*`, Klist_self)
          else if (COMBINE == "modified_multiplicative") Kc_self <- Reduce(`+`, Klist_self) / length(Klist_self)
          else                                        Kc_self <- Reduce(`+`, Klist_self)
        } else {
          Kc_self <- matrix(0, Nte, Nte)
        }
        
        # ----- Fingerprint Kf_self (per draw) -----
        Kf_self <- if (FP_KERNEL != "none") {
          if (FP_KERNEL == "rbf") {
            exp( -(D_fp_tt^2) / (2*l_fp^2) )
          } else if (FP_KERNEL == "matern32") {
            R <- D_fp_tt / l_fp; (1 + sqrt(3)*R) * exp(-sqrt(3)*R)
          } else if (FP_KERNEL == "matern52") {
            R <- D_fp_tt / l_fp; (1 + sqrt(5)*R + 5*R^2/3) * exp(-sqrt(5)*R)
          } else {
            D_fp_tt
          }
        } else matrix(0, Nte, Nte)
        
        
        # ----- Combine & variance with proper self-kernel -----
        K_self_core <- if (COMBINE == "additive") (Kc_self + Kf_self) else {
          if (CONT_KERNEL != "none" && FP_KERNEL != "none") Kc_self * Kf_self
          else if (CONT_KERNEL != "none") Kc_self else Kf_self
        }
        K_self <- (sigma^2) * K_self_core
        
        v <- backsolve(L, t(K_te), transpose = TRUE)
        var_te_s <- pmax(0, diag(K_self) - colSums(v^2)) + sigma_noise^2
        
        cbind(mu_te_s, sqrt(var_te_s))
      }
      
      
      pred_draws <- lapply(sel, predict_one)
      if (length(pred_draws) > 0 && nrow(pred_draws[[1]]) > 0) {
        mu_mat <- do.call(cbind, lapply(pred_draws, function(M) M[,1]))
        sd_mat <- do.call(cbind, lapply(pred_draws, function(M) M[,2]))
        
        # Back-transform
        mu_mean <- rowMeans(mu_mat) * y_scl + y_ctr
        sd_mean <- rowMeans(sd_mat) * abs(y_scl)
        
        # Per-split posterior over RMSE/MAE via draws
        mu_draw_bt <- sweep(mu_mat, 1, y_scl, `*`)
        mu_draw_bt <- sweep(mu_draw_bt, 1, y_ctr, `+`)
        rmse_draws <- sqrt(colMeans((mu_draw_bt - matrix(y_te, nrow = length(y_te), ncol = ncol(mu_draw_bt)))^2))
        mae_draws  <- colMeans(abs(mu_draw_bt - matrix(y_te, nrow = length(y_te), ncol = ncol(mu_draw_bt))))
        
        rmse_mean <- mean(rmse_draws); rmse_sd <- sd(rmse_draws)
        rmse_lwr  <- quantile(rmse_draws, 0.025); rmse_med <- quantile(rmse_draws, 0.5)
        rmse_upr  <- quantile(rmse_draws, 0.975)
        mae_mean  <- mean(mae_draws);  mae_sd  <- sd(mae_draws)
        mae_lwr   <- quantile(mae_draws, 0.025); mae_med <- quantile(mae_draws, 0.5)
        mae_upr   <- quantile(mae_draws, 0.975)
        
        hmc_rmse_rows[[length(hmc_rmse_rows)+1]] <- data.frame(
          seed = as.character(sp$seed), fold = sp$fold, split = sidx,
          RMSE_mean = rmse_mean, RMSE_sd = rmse_sd,
          RMSE_q2.5 = as.numeric(rmse_lwr), RMSE_q50 = as.numeric(rmse_med), RMSE_q97.5 = as.numeric(rmse_upr),
          MAE_mean  = mae_mean,  MAE_sd  = mae_sd,
          MAE_q2.5  = as.numeric(mae_lwr), MAE_q50 = as.numeric(mae_med), MAE_q97.5 = as.numeric(mae_upr),
          stringsAsFactors = FALSE
        )
        
        # Store compact predictive summary (for plotting)
        hmc_pred_summaries[[length(hmc_pred_summaries)+1]] <- data.frame(
          seed = as.character(sp$seed), fold = sp$fold, split = sidx,
          obs  = te,
          true = y_te,
          mean = mu_mean,
          lower = mu_mean - 1.96*sd_mean,
          upper = mu_mean + 1.96*sd_mean,
          stringsAsFactors = FALSE
        )
        
        # Predicted mean vs true with bands for this split (diagnostic subset)
        if (sidx %in% plot_indices && length(y_te) > 0) {
          df_plot <- data.frame(true = y_te, mean = mu_mean,
                                lower = mu_mean - 1.96*sd_mean,
                                upper = mu_mean + 1.96*sd_mean)
          print(ggplot(df_plot, aes(x = true, y = mean)) +
                  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
                  geom_point(alpha = 0.5, size = 1.8, color = "firebrick") +
                  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
                  labs(title = sprintf("HMC PPI — seed=%s fold=%s", sp$seed, sp$fold),
                       x = "True y", y = "Posterior mean") +
                  theme_minimal())
        }
      }
      
      # Diagnostics plots for subset of splits
      if (have_bayesplot && (sidx %in% plot_indices)) {
        # Densities
        try({
          print(
            bayesplot::mcmc_dens(
              fit_hmc,
              pars = c("sigma","sigma_noise","l_fp"),
              regex_pars = "^l_cont\\["
            ) + ggplot2::labs(title = sprintf("Posteriors — seed=%s fold=%s", sp$seed, sp$fold))
          )
        }, silent = TRUE)
        
        # Trace
        try({
          print(
            bayesplot::mcmc_trace(
              fit_hmc,
              pars = c("sigma","sigma_noise","l_fp"),
              regex_pars = "^l_cont\\[",
              facet_args = list(ncol = 1)
            ) + ggplot2::labs(title = sprintf("Trace — seed=%s fold=%s", sp$seed, sp$fold))
          )
        }, silent = TRUE)
        
        # ACF
        try({
          print(
            bayesplot::mcmc_acf(
              fit_hmc,
              pars = c("sigma","sigma_noise","l_fp"),
              regex_pars = "^l_cont\\[",
              lags = 50
            ) + ggplot2::labs(title = sprintf("ACF — seed=%s fold=%s", sp$seed, sp$fold))
          )
        }, silent = TRUE)
      }
      
      
      hmc_param_summaries[[length(hmc_param_summaries)]]$CONT_KERNEL <- CONT_KERNEL
      hmc_param_summaries[[length(hmc_param_summaries)]]$FP_KERNEL   <- FP_KERNEL
      hmc_param_summaries[[length(hmc_param_summaries)]]$COMBINE     <- COMBINE
      
      # clean
      rm(fit_hmc, post, pred_draws)
      gc()
    } # end per-split
    
    # ---- Collate & summarize HMC results across splits ----
    if (length(hmc_param_summaries) > 0) {
      hmc_params_all <- bind_rows(hmc_param_summaries) %>%
        mutate(parameter = as.character(parameter))
      
      hmc_preds_all <- if (length(hmc_pred_summaries)>0) bind_rows(hmc_pred_summaries) else NULL
      hmc_rmse_all  <- if (length(hmc_rmse_rows)>0)      bind_rows(hmc_rmse_rows)      else NULL
      
      # Overall RMSE/MAE summary across splits
      hmc_overall <- if (!is.null(hmc_rmse_all)) {
        summarise(hmc_rmse_all,
                  RMSE_mean_of_means = mean(RMSE_mean),
                  RMSE_sd_of_means   = sd(RMSE_mean),
                  MAE_mean_of_means  = mean(MAE_mean),
                  MAE_sd_of_means    = sd(MAE_mean))
      } else NULL
      
      # -----------------------------
      # Aggregated parameter table across splits
      # -----------------------------
      agg_params <- NULL
      if (nrow(hmc_params_all) > 0) {
        agg_params <- hmc_params_all %>%
          group_by(parameter) %>%
          summarise(
            splits = dplyr::n(),
            mean_of_means    = mean(mean, na.rm = TRUE),
            sd_of_means      = sd(mean, na.rm = TRUE),
            median_of_median = median(`50%`, na.rm = TRUE),
            q2.5_of_median   = quantile(`50%`, 0.025, na.rm = TRUE),
            q97.5_of_median  = quantile(`50%`, 0.975, na.rm = TRUE),
            mean_of_q2.5     = mean(`2.5%`, na.rm = TRUE),
            mean_of_q97.5    = mean(`97.5%`, na.rm = TRUE),
            .groups = "drop"
          )
      }
      
      # -----------------------------
      # Aggregated predictive summary table across splits
      # -----------------------------
      agg_pred <- NULL
      if (!is.null(hmc_preds_all) && nrow(hmc_preds_all) > 0) {
        agg_pred <- hmc_preds_all %>%
          mutate(covered = (true >= lower & true <= upper),
                 width   = upper - lower) %>%
          summarise(
            n_points  = dplyr::n(),
            coverage  = mean(covered, na.rm = TRUE),
            avg_width = mean(width, na.rm = TRUE),
            width_sd  = sd(width, na.rm = TRUE)
          )
      }
      
      # -----------------------------
      # Global plots across all splits
      # -----------------------------
      if (!is.null(hmc_rmse_all)) {
        print(ggplot(hmc_rmse_all, aes(x = RMSE_mean)) +
                geom_histogram(bins = 30) +
                labs(title = "HMC: Distribution of per-split posterior RMSE means",
                     x = "RMSE (posterior mean per split)", y = "Count") +
                theme_minimal())
      }
      
      if (!is.null(hmc_preds_all) && nrow(hmc_preds_all) > 0) {
        print(ggplot(hmc_preds_all, aes(x = true, y = mean)) +
                geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
                geom_point(alpha = 0.35, size = 1.6) +
                geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
                labs(title = "HMC: Test Posterior Predictive (all plotted splits)",
                     x = "True y", y = "Posterior mean") +
                theme_minimal())
      }
      
      if (!is.null(agg_params) && nrow(agg_params) > 0 && have_gridExtra) {
        gridExtra::grid.table(head(agg_params, 30))
      }
      if (!is.null(agg_pred) && nrow(agg_pred) > 0 && have_gridExtra) {
        gridExtra::grid.table(agg_pred)
      }
      
      # Save HMC bundle
      objs_to_save <- c("hmc_params_all","hmc_preds_all","hmc_rmse_all","hmc_overall",
                  "agg_params","agg_pred",
                  "labels","CONT_KERNEL","FP_KERNEL","COMBINE",
                  "DATA_CSV","TANI_PATH","SUFFIX",
                  "SEEDS","K_FOLDS","MAX_SEEDS",
                  "CHAINS","ITER","WARMUP","ADAPT_DELTA","MAX_TD","HMC_MAX_SPLITS","HMC_PRED_NSAMP",
                  "PLOT_NSPLITS","PLOT_SHOW_ALL","PLOT_PARAMS")

        if (exists("diag_draws", inherits = FALSE)) objs_to_save <- c(objs_to_save, "diag_draws")

        save(list = objs_to_save, file = HMC_RDATA)

      cat("Saved HMC summaries to: ", HMC_RDATA, "\n")
      cat("Saved HMC plots to: ", PLOTS_PDF, "\n")
    } else {
      message("No HMC summaries collected (possibly HMC_MAX_SPLITS == 0).")
    }
  } # end tanimoto guard
} # end DO_HMC
