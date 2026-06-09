.libPaths("~/rlibs")

library(mvtnorm)
library(parallel)
library(optparse)

rm(list = ls())
options(scipen = 999)

source("rareF_main.R")
source("rareF_functions.R")

option_list <- list(
  make_option("--setting-id", type = "integer"),
  make_option("--output-dir", type = "character", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))
setting_id <- opt$`setting-id`
output_dir <- opt$`output-dir`

if (is.null(setting_id)) {
  stop("Please provide --setting-id.")
}

# Keep only the main-text R-Lopt(BL) setting.
settings <- list(
  list(method = "R-Lopt-BLplt", criterion = "R-Lopt")
)

method.all <- sapply(settings, function(x) x$method)
num.method <- length(method.all)

estimator.all <- c("pilot", "ssp", "cmb.union", "cmb.samples", "cmb.estimators")
num.estimator <- length(estimator.all)

safe_diag_sum <- function(A) {
  if (is.null(A) || any(is.na(A))) return(NA_real_)
  return(sum(diag(A)))
}

rareF_simu <- function(j, d_cont_setting, d_rare_setting, p_rare_setting) {
  set.seed(j + 1000)
  
  d_cont <- d_cont_setting
  d_rare <- d_rare_setting
  p_rare <- p_rare_setting
  
  beta0 <- c(0.5, rep(0.5, d_rare), rep(0.5, d_cont))
  d <- 1 + d_rare + d_cont
  
  coef.plt <- coef.ssp <- coef.cmb.union <- coef.cmb.samples <- coef.cmb.estimators <-
    array(NA_real_, dim = c(length(n.ssp), num.method, d))
  
  cov.plt <- cov.ssp <- cov.cmb.union <- cov.cmb.samples <- cov.cmb.estimators <-
    matrix(NA_real_, nrow = length(n.ssp), ncol = num.method)
  
  size.plt <- size.ssp <- size.cmb.union <- size.cmb.samples <-
    matrix(NA_real_, nrow = length(n.ssp), ncol = num.method)
  
  n.overlap <- overlap.rate <-
    matrix(NA_real_, nrow = length(n.ssp), ncol = num.method)
  
  rare.count.plt <- rare.count.ssp <- rare.count.cmb.union <-
    array(NA_real_, dim = c(length(n.ssp), num.method, d_rare))
  
  rare.counts.plt <- rare.counts.ssp <- rare.counts.cmb.union <-
    array(NA_real_, dim = c(length(n.ssp), num.method, d_rare + 1))
  
  fail.ind <- matrix(0, nrow = length(n.ssp), ncol = num.method)
  
  corr <- 0.5
  sigmax <- matrix(corr, d_cont, d_cont) + diag(1 - corr, d_cont)
  
  X.cont <- mvtnorm::rmvnorm(N, rep(0, d_cont), sigmax)
  colnames(X.cont) <- paste0("X", seq_len(d_cont))
  
  Z <- do.call(cbind, lapply(seq_along(p_rare), function(k) {
    rbinom(N, 1, p_rare[k])
  }))
  colnames(Z) <- paste0("Z", seq_len(d_rare))
  
  X <- cbind(Z, X.cont)
  
  eta <- beta0[1] + X %*% beta0[-1]
  P <- 1 / (1 + exp(-eta))
  Y <- as.integer(rbinom(N, 1, P))
  
  full.data <- as.data.frame(cbind(Y = Y, X))
  formula <- Y ~ .
  
  rareFeature.index <- 2:(1 + d_rare)
  
  for (i in seq_along(n.ssp)) {
    for (m in seq_along(settings)) {
      result <- try({
        object <- ssp.rareF.glm(
          formula = formula,
          data = full.data,
          n.plt = n.plt,
          n.ssp = n.ssp[i],
          family = "quasibinomial",
          criterion = settings[[m]]$criterion,
          sampling.method = "poisson",
          balance.X.plt = TRUE,
          balance.Y.plt = FALSE,
          balance.Y.ssp = FALSE,
          rareThreshold = 0.5,
          combine = TRUE,
          rareFeature.index = rareFeature.index
        )
        
        this.n.overlap <- length(object$index.ssp) +
          length(object$index.plt) -
          length(object$index.cmb.union)
        
        this.overlap.rate <- (this.n.overlap / length(object$index.ssp)) * 100
        
        list(
          coef.plt = object$coef.plt,
          coef.ssp = object$coef.ssp,
          coef.cmb.union = object$coef.cmb.union,
          coef.cmb.samples = object$coef.cmb.samples,
          coef.cmb.estimators = object$coef.cmb.estimators,
          cov.plt = safe_diag_sum(object$cov.plt),
          cov.ssp = safe_diag_sum(object$cov.ssp),
          cov.cmb.union = safe_diag_sum(object$cov.cmb.union),
          cov.cmb.samples = safe_diag_sum(object$cov.cmb.samples),
          cov.cmb.estimators = safe_diag_sum(object$cov.cmb.estimators),
          size.plt = length(object$index.plt),
          size.ssp = length(object$index.ssp),
          size.cmb.union = length(object$index.cmb.union),
          size.cmb.samples = length(object$index.cmb.samples),
          n.overlap = this.n.overlap,
          overlap.rate = this.overlap.rate,
          rare.count.plt = object$rare.count.plt,
          rare.count.ssp = object$rare.count.ssp,
          rare.count.cmb.union = object$rare.count.cmb.union,
          rare.counts.plt = object$rare.counts.plt,
          rare.counts.ssp = object$rare.counts.ssp,
          rare.counts.cmb.union = object$rare.counts.cmb.union,
          fail.ind = 0
        )
      }, silent = TRUE)
      
      result_is_invalid <- inherits(result, "try-error") ||
        is.null(result$coef.plt) || length(result$coef.plt) != d ||
        is.null(result$coef.ssp) || length(result$coef.ssp) != d ||
        is.null(result$coef.cmb.union) || length(result$coef.cmb.union) != d ||
        is.null(result$coef.cmb.samples) || length(result$coef.cmb.samples) != d ||
        is.null(result$coef.cmb.estimators) || length(result$coef.cmb.estimators) != d ||
        is.null(result$rare.count.plt) || length(result$rare.count.plt) != d_rare ||
        is.null(result$rare.count.ssp) || length(result$rare.count.ssp) != d_rare ||
        is.null(result$rare.count.cmb.union) || length(result$rare.count.cmb.union) != d_rare ||
        is.null(result$rare.counts.plt) || length(result$rare.counts.plt) != (d_rare + 1) ||
        is.null(result$rare.counts.ssp) || length(result$rare.counts.ssp) != (d_rare + 1) ||
        is.null(result$rare.counts.cmb.union) || length(result$rare.counts.cmb.union) != (d_rare + 1) ||
        any(is.na(result$coef.ssp)) ||
        any(is.na(result$coef.cmb.union)) ||
        any(is.na(result$coef.cmb.samples)) ||
        any(is.na(result$coef.cmb.estimators))
      
      if (result_is_invalid) {
        result <- list(
          coef.plt = rep(NA_real_, d),
          coef.ssp = rep(NA_real_, d),
          coef.cmb.union = rep(NA_real_, d),
          coef.cmb.samples = rep(NA_real_, d),
          coef.cmb.estimators = rep(NA_real_, d),
          cov.plt = NA_real_,
          cov.ssp = NA_real_,
          cov.cmb.union = NA_real_,
          cov.cmb.samples = NA_real_,
          cov.cmb.estimators = NA_real_,
          size.plt = NA_real_,
          size.ssp = NA_real_,
          size.cmb.union = NA_real_,
          size.cmb.samples = NA_real_,
          n.overlap = NA_real_,
          overlap.rate = NA_real_,
          rare.count.plt = rep(NA_real_, d_rare),
          rare.count.ssp = rep(NA_real_, d_rare),
          rare.count.cmb.union = rep(NA_real_, d_rare),
          rare.counts.plt = rep(NA_real_, d_rare + 1),
          rare.counts.ssp = rep(NA_real_, d_rare + 1),
          rare.counts.cmb.union = rep(NA_real_, d_rare + 1),
          fail.ind = 1
        )
      }
      
      coef.plt[i, m, ] <- result$coef.plt
      coef.ssp[i, m, ] <- result$coef.ssp
      coef.cmb.union[i, m, ] <- result$coef.cmb.union
      coef.cmb.samples[i, m, ] <- result$coef.cmb.samples
      coef.cmb.estimators[i, m, ] <- result$coef.cmb.estimators
      
      cov.plt[i, m] <- result$cov.plt
      cov.ssp[i, m] <- result$cov.ssp
      cov.cmb.union[i, m] <- result$cov.cmb.union
      cov.cmb.samples[i, m] <- result$cov.cmb.samples
      cov.cmb.estimators[i, m] <- result$cov.cmb.estimators
      
      size.plt[i, m] <- result$size.plt
      size.ssp[i, m] <- result$size.ssp
      size.cmb.union[i, m] <- result$size.cmb.union
      size.cmb.samples[i, m] <- result$size.cmb.samples
      
      n.overlap[i, m] <- result$n.overlap
      overlap.rate[i, m] <- result$overlap.rate
      
      rare.count.plt[i, m, ] <- result$rare.count.plt
      rare.count.ssp[i, m, ] <- result$rare.count.ssp
      rare.count.cmb.union[i, m, ] <- result$rare.count.cmb.union
      
      rare.counts.plt[i, m, ] <- result$rare.counts.plt
      rare.counts.ssp[i, m, ] <- result$rare.counts.ssp
      rare.counts.cmb.union[i, m, ] <- result$rare.counts.cmb.union
      
      fail.ind[i, m] <- result$fail.ind
    }
  }
  
  list(
    beta0 = beta0,
    coef.plt = coef.plt,
    coef.ssp = coef.ssp,
    coef.cmb.union = coef.cmb.union,
    coef.cmb.samples = coef.cmb.samples,
    coef.cmb.estimators = coef.cmb.estimators,
    cov.plt = cov.plt,
    cov.ssp = cov.ssp,
    cov.cmb.union = cov.cmb.union,
    cov.cmb.samples = cov.cmb.samples,
    cov.cmb.estimators = cov.cmb.estimators,
    size.plt = size.plt,
    size.ssp = size.ssp,
    size.cmb.union = size.cmb.union,
    size.cmb.samples = size.cmb.samples,
    n.overlap = n.overlap,
    overlap.rate = overlap.rate,
    rare.count.plt = rare.count.plt,
    rare.count.ssp = rare.count.ssp,
    rare.count.cmb.union = rare.count.cmb.union,
    rare.counts.plt = rare.counts.plt,
    rare.counts.ssp = rare.counts.ssp,
    rare.counts.cmb.union = rare.counts.cmb.union,
    fail.ind = fail.ind
  )
}

summarize_results <- function(setting.results, beta0, n.ssp, method.all,
                              d_rare,
                              remove_outliers = FALSE) {
  out <- list()
  row_id <- 1
  
  get_coef_array <- function(res, est) {
    if (est == "pilot") return(res$coef.plt)
    if (est == "ssp") return(res$coef.ssp)
    if (est == "cmb.union") return(res$coef.cmb.union)
    if (est == "cmb.samples") return(res$coef.cmb.samples)
    if (est == "cmb.estimators") return(res$coef.cmb.estimators)
    stop("Unknown estimator.")
  }
  
  get_cov_matrix <- function(res, est) {
    if (est == "pilot") return(res$cov.plt)
    if (est == "ssp") return(res$cov.ssp)
    if (est == "cmb.union") return(res$cov.cmb.union)
    if (est == "cmb.samples") return(res$cov.cmb.samples)
    if (est == "cmb.estimators") return(res$cov.cmb.estimators)
    stop("Unknown estimator.")
  }
  
  safe_mean <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
  }
  
  safe_ratio <- function(num, den) {
    if (is.na(num) || is.na(den) || den <= 0) return(NA_real_)
    num / den
  }
  
  get_bias_variance <- function(coef.keep, beta0, index = NULL) {
    if (is.null(index)) {
      index <- seq_along(beta0)
    }
    
    coef.sub <- coef.keep[, index, drop = FALSE]
    beta.sub <- beta0[index]
    
    if (nrow(coef.sub) == 0) {
      return(list(
        squared_bias = NA_real_,
        empirical_variance = NA_real_,
        empirical_MSE = NA_real_,
        bias_ratio = NA_real_,
        variance_ratio = NA_real_
      ))
    }
    
    beta.bar <- colMeans(coef.sub, na.rm = TRUE)
    
    squared.bias <- sum((beta.bar - beta.sub)^2)
    
    empirical.variance <- mean(
      rowSums(
        sweep(coef.sub, 2, beta.bar, FUN = "-")^2
      ),
      na.rm = TRUE
    )
    
    beta.mat <- matrix(
      beta.sub,
      nrow = nrow(coef.sub),
      ncol = length(beta.sub),
      byrow = TRUE
    )
    
    empirical.MSE <- mean(
      rowSums((coef.sub - beta.mat)^2),
      na.rm = TRUE
    )
    
    list(
      squared_bias = squared.bias,
      empirical_variance = empirical.variance,
      empirical_MSE = empirical.MSE,
      bias_ratio = safe_ratio(squared.bias, empirical.MSE),
      variance_ratio = safe_ratio(empirical.variance, empirical.MSE)
    )
  }
  
  d <- length(beta0)
  rare.idx <- 2:(1 + d_rare)
  nonrare.idx <- c(1, (2 + d_rare):d)
  
  for (i in seq_along(n.ssp)) {
    for (m in seq_along(method.all)) {
      for (est in estimator.all) {
        coef.mat <- do.call(rbind, lapply(setting.results, function(res) {
          get_coef_array(res, est)[i, m, ]
        }))
        
        cov.vec <- sapply(setting.results, function(res) {
          get_cov_matrix(res, est)[i, m]
        })
        
        good <- complete.cases(coef.mat)
        
        sq.err <- rep(NA_real_, length(setting.results))
        if (any(good)) {
          beta.mat <- matrix(
            beta0,
            nrow = sum(good),
            ncol = length(beta0),
            byrow = TRUE
          )
          sq.err[good] <- rowSums(
            (coef.mat[good, , drop = FALSE] - beta.mat)^2
          )
        }
        
        keep <- good
        
        if (remove_outliers && sum(good, na.rm = TRUE) >= 4) {
          outlier_iqr_mult <- 5
          q <- quantile(
            sq.err[good],
            probs = c(0.25, 0.75),
            na.rm = TRUE
          )
          iqr_val <- q[2] - q[1]
          
          lower <- q[1] - outlier_iqr_mult * iqr_val
          upper <- q[2] + outlier_iqr_mult * iqr_val
          
          keep <- good & sq.err >= lower & sq.err <= upper
        }
        
        coef.keep <- coef.mat[keep, , drop = FALSE]
        
        bv.all <- get_bias_variance(
          coef.keep = coef.keep,
          beta0 = beta0,
          index = seq_along(beta0)
        )
        
        bv.rare <- get_bias_variance(
          coef.keep = coef.keep,
          beta0 = beta0,
          index = rare.idx
        )
        
        bv.nonrare <- get_bias_variance(
          coef.keep = coef.keep,
          beta0 = beta0,
          index = nonrare.idx
        )
        
        empirical.MSE <- bv.all$empirical_MSE
        avg.formula.MSE <- safe_mean(cov.vec[keep])
        
        fail.rate <- mean(sapply(setting.results, function(res) {
          res$fail.ind[i, m]
        }), na.rm = TRUE)
        
        avg.size.plt <- mean(sapply(setting.results, function(res) {
          res$size.plt[i, m]
        }), na.rm = TRUE)
        
        avg.size.ssp <- mean(sapply(setting.results, function(res) {
          res$size.ssp[i, m]
        }), na.rm = TRUE)
        
        avg.size.union <- mean(sapply(setting.results, function(res) {
          res$size.cmb.union[i, m]
        }), na.rm = TRUE)
        
        avg.size.samples <- mean(sapply(setting.results, function(res) {
          res$size.cmb.samples[i, m]
        }), na.rm = TRUE)

        avg.n.overlap <- mean(sapply(setting.results, function(res) {
          res$n.overlap[i, m]
        }), na.rm = TRUE)
        
        avg.overlap.rate <- mean(sapply(setting.results, function(res) {
          res$overlap.rate[i, m]
        }), na.rm = TRUE)
        
        out[[row_id]] <- data.frame(
          n.ssp = n.ssp[i],
          method = method.all[m],
          estimator = est,
          empirical_MSE = empirical.MSE,
          empirical_variance = bv.all$empirical_variance,
          squared_bias = bv.all$squared_bias,
          bias_ratio = bv.all$bias_ratio,
          variance_ratio = bv.all$variance_ratio,
          empirical_MSE_rare = bv.rare$empirical_MSE,
          empirical_variance_rare = bv.rare$empirical_variance,
          squared_bias_rare = bv.rare$squared_bias,
          bias_ratio_rare = bv.rare$bias_ratio,
          variance_ratio_rare = bv.rare$variance_ratio,
          empirical_MSE_nonrare = bv.nonrare$empirical_MSE,
          empirical_variance_nonrare = bv.nonrare$empirical_variance,
          squared_bias_nonrare = bv.nonrare$squared_bias,
          bias_ratio_nonrare = bv.nonrare$bias_ratio,
          variance_ratio_nonrare = bv.nonrare$variance_ratio,
          avg_formula_MSE = avg.formula.MSE,
          formula_to_empirical_ratio = safe_ratio(avg.formula.MSE, empirical.MSE),
          fail_rate = fail.rate,
          avg_size_plt = avg.size.plt,
          avg_size_ssp = avg.size.ssp,
          avg_size_union = avg.size.union,
          avg_size_samples = avg.size.samples,
          avg_n_overlap = avg.n.overlap,
          avg_overlap_rate_percent = avg.overlap.rate
        )
        
        row_id <- row_id + 1
      }
    }
  }
  
  do.call(rbind, out)
}

rpt <- 1000
N <- 2e5

sample.rate <- c(0.005, 0.0125, 0.025, 0.05, 0.1)
n.ssp <- as.integer(sample.rate * N)

n_plt_vals  <- c(300, 1000, 10000)
d_cont_vals <- c(5)
d_rare_vals <- c(10, 5, 2)
p_rare_vals <- c(0.001, 0.005, 0.4)

combined.settings <- list()
for (d_cont in d_cont_vals) {
  for (d_rare in d_rare_vals) {
    for (p_rare in p_rare_vals) {
      for (n_plt in n_plt_vals) {
        combined.settings[[length(combined.settings) + 1]] <- list(
          d_cont = d_cont,
          d_rare = d_rare,
          p_rare = rep(p_rare, d_rare),
          n.plt  = n_plt
        )
      }
    }
  }
}

cat("setting_id:", setting_id, "\n")

set <- setting_id + 1
if (set > length(combined.settings)) {
  stop(paste("Invalid setting index:", setting_id))
}

setting <- combined.settings[[set]]

d_cont <- setting$d_cont
d_rare <- setting$d_rare
p_rare <- setting$p_rare
n.plt  <- setting$n.plt

beta0 <- c(0.5, rep(0.5, d_rare), rep(0.5, d_cont))

cat("Running setting:\n")
cat("d_cont:", d_cont, "\n")
cat("d_rare:", d_rare, "\n")
cat("p_rare:", paste(p_rare, collapse = ", "), "\n")
cat("n.plt:", n.plt, "\n")
cat("n.ssp:", paste(n.ssp, collapse = ", "), "\n")

date_tag <- format(Sys.Date(), "%m%d%Y")

if (is.null(output_dir) || output_dir == "") {
  save_dir <- file.path("results", date_tag, "combine_compare")
} else {
  save_dir <- file.path(output_dir, date_tag, "combine_compare")
}

if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
worker_cap <- as.numeric(Sys.getenv("RAREF_MAX_WORKERS", "16"))
if (is.na(worker_cap) || worker_cap < 1) {
  worker_cap <- 16
}
max_workers <- min(worker_cap, rpt)
if (is.na(n_cores) || n_cores < 1) {
  n_cores <- 1
}
if (n_cores > max_workers) {
  cat("Capping worker count from", n_cores, "to", max_workers, "to avoid worker-memory failures\n")
  n_cores <- max_workers
}
cat("Detected", n_cores, "cores for parallel execution (cap =", max_workers, ")\n")

cl <- parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

invisible(clusterEvalQ(cl, {
  .libPaths("~/rlibs")
  library(mvtnorm)
  source("rareF_main.R")
  source("rareF_functions.R")
  NULL
}))

clusterExport(
  cl = cl,
  varlist = c(
    "settings", "method.all", "num.method",
    "estimator.all", "num.estimator",
    "rareF_simu", "safe_diag_sum",
    "n.plt", "n.ssp", "N",
    "d_cont", "d_rare", "p_rare"
  ),
  envir = environment()
)

t1 <- proc.time()

setting.results <- parLapply(cl, seq_len(rpt), function(j) {
  rareF_simu(
    j = j,
    d_cont_setting = d_cont,
    d_rare_setting = d_rare,
    p_rare_setting = p_rare
  )
})

stopCluster(cl)
cl <- NULL

t2 <- proc.time()
cat("Time for setting:", t2[3] - t1[3], "seconds\n")

summary.results <- summarize_results(
  setting.results = setting.results,
  beta0 = beta0,
  n.ssp = n.ssp,
  method.all = method.all,
  d_rare = d_rare
)

print(summary.results)

save_file <- file.path(
  save_dir,
  paste0(
    "combineCompare",
    "_plt", n.plt,
    "_dcont", d_cont,
    "_drare", d_rare,
    "_prare", format(p_rare[1], scientific = FALSE),
    ".Rdata"
  )
)

save(
  setting.results,
  summary.results,
  settings,
  method.all,
  estimator.all,
  combined.settings,
  setting,
  beta0,
  n.ssp,
  rpt,
  N,
  file = save_file
)

cat("Saved results to", save_file, "\n")
