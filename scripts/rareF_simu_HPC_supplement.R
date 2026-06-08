.libPaths("~/rlibs")  # Use locally installed R packages
library(mvtnorm)
library(parallel)
library(optparse)
rm(list = ls())
options(scipen = 999)
source('rareF_main.R')
source('rareF_functions.R')

option_list <- list(
  make_option("--setting-id", type = "integer"),
  make_option("--output-dir", type = "character"),
  make_option("--rpt", type = "integer", default = NA_integer_),
  make_option("--N", type = "double", default = NA_real_),
  make_option("--sample-rate", type = "character", default = NA_character_),
  make_option("--n-plt-vals", type = "character", default = NA_character_),
  make_option("--d-cont-vals", type = "character", default = NA_character_),
  make_option("--d-rare-vals", type = "character", default = NA_character_),
  make_option("--p-rare-vals", type = "character", default = NA_character_),
  make_option("--n-cores", type = "integer", default = NA_integer_)
)
opt <- parse_args(OptionParser(option_list = option_list))
setting_id <- opt$`setting-id`
output_dir <- opt$`output-dir`

parse_numeric_csv <- function(x, as_integer = FALSE) {
  if (is.na(x) || !nzchar(x)) {
    return(NULL)
  }
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  values <- values[nzchar(values)]
  if (length(values) == 0) {
    return(NULL)
  }
  if (as_integer) {
    return(as.integer(values))
  }
  as.numeric(values)
}

make_method <- function(method, criterion, balance.X.plt) {
  list(
    method = method,
    criterion = criterion,
    balance.X.plt = balance.X.plt,
    objective.weight = "weighted",
    objective.weight.plt = "weighted"
  )
}

settings <- list(
  make_method("BL-Uni-Wobj", "BL-Uni", TRUE),
  make_method("Uni", "Uni", FALSE),
  make_method("R-Lopt-BLplt", "R-Lopt", TRUE),
  make_method("BL-Lopt-BLplt-Wpilot", "BL-Lopt", TRUE),
  make_method("Aopt-BLplt-Wpilot", "Aopt", TRUE),
  make_method("Aopt-Uplt", "Aopt", FALSE),
  make_method("Lopt-BLplt-Wpilot", "Lopt", TRUE),
  make_method("Lopt-Uplt", "Lopt", FALSE)
)

method.all <- vapply(settings, `[[`, character(1), "method")
num.method <- length(method.all)

rareF_simu <- function(j, d_cont_setting, d_rare_setting, p_rare_setting) {
  set.seed(j)
  d_cont <- d_cont_setting
  d_rare <- d_rare_setting
  p_rare <- p_rare_setting
  
  beta0 <- c(0.5, rep(0.5, d_rare), rep(0.5, d_cont))
  d <- 1 + d_rare + d_cont
  
  coef.plt <- coef.ssp <- coef.cmb <- array(NA, dim = c(length(n.ssp), num.method, d))
  cov.ssp <- cov.cmb <- array(NA, dim = c(length(n.ssp), num.method))
  SubsampleSize <- array(NA, dim = c(length(n.ssp), num.method))
  rare.count.plt <- rare.count.ssp <- rare.count.cmb <- 
    array(NA, dim = c(length(n.ssp), num.method, length(rareFeature.index)))
  rare.counts.plt <- rare.counts.ssp <- rare.counts.cmb <- 
    array(NA, dim = c(length(n.ssp), num.method, length(rareFeature.index) + 1))
  fail.ind <- matrix(0, nrow = length(n.ssp), ncol = num.method)
  
  corr <- 0.5
  sigmax  <- matrix(corr, d_cont, d_cont) + diag(1-corr, d_cont)
  X <- mvtnorm::rmvnorm(N, rep(0, d_cont), sigmax)
  Z <- do.call(cbind, lapply(seq_along(p_rare), function(i) {
    rbinom(N, 1, p_rare[i])
  }))
  
  X <- cbind(Z, X)
  
  P <- 1 - 1 / (1 + exp(beta0[1] + X %*% beta0[-1]))
  Y <- as.integer(rbinom(N, 1, P))
  full.data <- as.data.frame(cbind(Y, X))
  formula <- Y ~ .
  coef.full <- coef(glm(formula = formula,
                        data = full.data,
                        family = 'binomial'))
  for (i in seq_along(n.ssp)) {
    for (m in seq_along(settings)) {
      result <- try({
        fit <- ssp.rareF.glm(
          formula = formula, 
          data = full.data,
          n.plt = n.plt, 
          n.ssp = n.ssp[i],
          family = 'quasibinomial',
          criterion = settings[[m]]$criterion,
          sampling.method = "poisson",
          objective.weight.plt = settings[[m]]$objective.weight.plt,
          objective.weight = settings[[m]]$objective.weight,
          balance.X.plt = settings[[m]]$balance.X.plt,
          rareFeature.index = NULL
        )
        list(coef.plt = fit$coef.plt,
             coef.ssp = fit$coef.ssp,
             cov.ssp = sum(diag(fit$cov.ssp)),
             coef.cmb = fit$coef.cmb.union,
             cov.cmb = sum(diag(fit$cov.cmb.union)),
             SubsampleSize = length(fit$index.ssp),
             rare.count.plt = fit$rare.count.plt,
             rare.count.ssp = fit$rare.count.ssp,
             rare.count.cmb = fit$rare.count.cmb.union,
             rare.counts.plt = fit$rare.counts.plt,
             rare.counts.ssp = fit$rare.counts.ssp,
             rare.counts.cmb = fit$rare.counts.cmb.union,
             fail.ind = 0)
      }, silent = TRUE)
      if (inherits(result, "try-error") ||
          is.null(result$coef.ssp) || 
          any(is.na(result$coef.ssp))) {
        result <- list(coef.plt = NA,
                       coef.ssp = NA,
                       cov.ssp = NA,
                       coef.cmb = NA,
                       cov.cmb = NA,
                       SubsampleSize = NA,
                       rare.count.plt = NA,
                       rare.count.ssp = NA,
                       rare.count.cmb = NA,
                       rare.counts.plt = NA,
                       rare.counts.ssp = NA,
                       rare.counts.cmb = NA,
                       fail.ind = 1)
      }
      coef.plt[i, m, ] <- result$coef.plt
      coef.ssp[i, m, ] <- result$coef.ssp
      cov.ssp[i, m] <- result$cov.ssp
      coef.cmb[i, m, ] <- result$coef.cmb
      cov.cmb[i, m] <- result$cov.cmb
      SubsampleSize[i, m] <- result$SubsampleSize
      rare.count.plt[i, m, ] <- result$rare.count.plt
      rare.count.ssp[i, m, ] <- result$rare.count.ssp
      rare.count.cmb[i, m, ] <- result$rare.count.cmb
      rare.counts.plt[i, m, ] <- result$rare.counts.plt
      rare.counts.ssp[i, m, ] <- result$rare.counts.ssp
      rare.counts.cmb[i, m, ] <- result$rare.counts.cmb
      fail.ind[i, m] = result$fail.ind
    }
  }
  list(
    coef.plt = coef.plt,
    coef.ssp = coef.ssp, coef.cmb = coef.cmb, coef.full = coef.full,
    cov.ssp = cov.ssp, cov.cmb = cov.cmb, 
    SubsampleSize = SubsampleSize,
    rare.count.plt = rare.count.plt, rare.count.ssp = rare.count.ssp,
    rare.count.cmb = rare.count.cmb,
    rare.counts.plt = rare.counts.plt, rare.counts.ssp = rare.counts.ssp,
    rare.counts.cmb = rare.counts.cmb,
    fail.ind = fail.ind
  )
}


n_cores <- if (!is.na(opt$`n-cores`)) {
  opt$`n-cores`
} else {
  as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
}
cat("Detected", n_cores, "cores for parallel execution\n")
cl <- NULL
if (n_cores > 1) {
  cl <- parallel::makeCluster(n_cores)
}
rpt <- 1000
N <- 2e5
sample.rate <- c(0.01, 0.025, 0.05, 0.15)
n_plt_vals  <- c(300, 1000)
d_cont_vals <- c(10, 5, 2)
d_rare_vals <- c(10, 5, 2)
p_rare_vals <- c(0.001, 0.01, 0.05)

if (!is.na(opt$`rpt`)) rpt <- opt$`rpt`
if (!is.na(opt$`N`)) N <- opt$`N`
sample.rate.override <- parse_numeric_csv(opt$`sample-rate`)
n.plt.override <- parse_numeric_csv(opt$`n-plt-vals`, as_integer = TRUE)
d.cont.override <- parse_numeric_csv(opt$`d-cont-vals`, as_integer = TRUE)
d.rare.override <- parse_numeric_csv(opt$`d-rare-vals`, as_integer = TRUE)
p.rare.override <- parse_numeric_csv(opt$`p-rare-vals`)
if (!is.null(sample.rate.override)) sample.rate <- sample.rate.override
if (!is.null(n.plt.override)) n_plt_vals <- n.plt.override
if (!is.null(d.cont.override)) d_cont_vals <- d.cont.override
if (!is.null(d.rare.override)) d_rare_vals <- d.rare.override
if (!is.null(p.rare.override)) p_rare_vals <- p.rare.override

n.ssp <- sample.rate * N

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

date_tag <- format(Sys.Date(), "%m%d%Y")
save_dir <- if (!is.null(output_dir) && nzchar(output_dir)) {
  output_dir
} else {
  file.path("..", "raw_results", date_tag, "logic")
}
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

set <- setting_id + 1
if (set > length(combined.settings)) {
  stop(paste("Invalid setting index:", setting_id))
}

setting <- combined.settings[[set]]
d_cont <- setting$d_cont
d_rare <- setting$d_rare
p_rare <- setting$p_rare
n.plt  <- setting$n.plt
rareFeature.index <- 1:d_rare

cat(
  sprintf(
    "Running setting %d: n_pl=%d, d_u=%d, d_r=%d, rho=%s\n",
    setting_id,
    n.plt,
    d_cont,
    d_rare,
    format(p_rare[1], scientific = FALSE)
  )
)

if (!is.null(cl)) {
  clusterEvalQ(cl, {
    source('rareF_main.R')
    source('rareF_functions.R')
  })
  
  clusterExport(cl = cl,
                varlist = c("settings", "rareF_simu", "n.plt", "n.ssp",  "num.method",
                            "d_cont", "p_rare", "rareFeature.index", "d_rare", "N"),
                envir = environment())
}

t1 <- proc.time()
if (is.null(cl)) {
  setting.results <- lapply(1:rpt, function(j) {
    rareF_simu(j, d_cont_setting = d_cont, d_rare_setting = d_rare, p_rare_setting = p_rare)
  })
} else {
  setting.results <- parLapply(cl, 1:rpt, function(j) {
    rareF_simu(j, d_cont_setting = d_cont, d_rare_setting = d_rare, p_rare_setting = p_rare)
  })
  stopCluster(cl)
}
t2 <- proc.time()

cat("Time for setting:", t2[3] - t1[3], "seconds\n")

save_file <- file.path(save_dir,
                       paste0("plt", n.plt,
                              "_dcont", d_cont,
                              "_drare", d_rare,
                              "_prare", format(p_rare[1], scientific = FALSE),
                              ".Rdata")
)
save(setting.results, method.all, settings, file = save_file)

cat("Saved results to", save_file, "\n")
