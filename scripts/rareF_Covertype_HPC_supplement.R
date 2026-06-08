.libPaths("~/rlibs")
library(parallel)
library(optparse)

rm(list = ls())
options(scipen = 999)

source("rareF_main.R")
source("rareF_functions.R")

option_list <- list(
  make_option("--setting-id", type = "integer"),
  make_option("--output-dir", type = "character", default = NA_character_),
  make_option("--rpt-total", type = "integer", default = NA_integer_),
  make_option("--num-groups", type = "integer", default = NA_integer_),
  make_option("--n-plt-vals", type = "character", default = NA_character_),
  make_option("--sample-sizes", type = "character", default = NA_character_),
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

load("CoverType_gc_full.Rdata")
full.data <- data_model_clean
target_indices <- which(
  grepl("^gc_dummy", names(full.data)) & colMeans(full.data) < 0.05
)
rareFeature.index <- target_indices - 1
N <- nrow(full.data)
formula <- y ~ .

cat("Loaded CoverType_gc_full.Rdata, N =", N, "\n")

rpt_total <- 1000
num_groups <- 40
n_plt_vals <- c(300, 1000)
n.ssp <- c(5e3, 1e4, 5e4)

if (!is.na(opt$`rpt-total`)) rpt_total <- opt$`rpt-total`
if (!is.na(opt$`num-groups`)) num_groups <- opt$`num-groups`
n.plt.override <- parse_numeric_csv(opt$`n-plt-vals`, as_integer = TRUE)
n.ssp.override <- parse_numeric_csv(opt$`sample-sizes`, as_integer = TRUE)
if (!is.null(n.plt.override)) n_plt_vals <- n.plt.override
if (!is.null(n.ssp.override)) n.ssp <- n.ssp.override

rpt_per_group <- ceiling(rpt_total / num_groups)

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

combined.settings <- list()
for (n.plt in n_plt_vals) {
  for (g in seq_len(num_groups)) {
    combined.settings[[length(combined.settings) + 1]] <- list(
      n.plt = n.plt,
      group_id = g,
      start_idx = (g - 1) * rpt_per_group + 1,
      end_idx = min(g * rpt_per_group, rpt_total)
    )
  }
}

if (setting_id + 1 > length(combined.settings)) {
  stop(paste("Invalid setting index:", setting_id))
}
set <- combined.settings[[setting_id + 1]]
n.plt <- set$n.plt
group_id <- set$group_id
rpt_indices <- set$start_idx:set$end_idx

cat(sprintf(
  "Running setting %d (n.plt = %d, group %d, reps %d-%d)\n",
  setting_id, n.plt, group_id, min(rpt_indices), max(rpt_indices)
))

rareF_simu <- function(j) {
  set.seed(j + 1000)

  coef.plt <- coef.ssp <- coef.cmb <- array(
    NA_real_,
    dim = c(length(n.ssp), num.method, ncol(full.data))
  )
  cov.ssp <- cov.cmb <- array(NA_real_, dim = c(length(n.ssp), num.method))
  SubsampleSize <- array(NA_real_, dim = c(length(n.ssp), num.method))
  rare.count.plt <- rare.count.ssp <- rare.count.cmb <- array(
    NA_real_,
    dim = c(length(n.ssp), num.method, length(rareFeature.index))
  )
  rare.counts.plt <- rare.counts.ssp <- rare.counts.cmb <- array(
    NA_real_,
    dim = c(length(n.ssp), num.method, length(rareFeature.index) + 1)
  )
  fail.ind <- matrix(0, nrow = length(n.ssp), ncol = num.method)
  comp.time <- matrix(NA_real_, nrow = length(n.ssp), ncol = num.method)

  coef.full <- beta_full

  for (i in seq_along(n.ssp)) {
    for (m in seq_along(settings)) {
      result <- try({
        fit <- ssp.rareF.glm(
          formula = formula,
          data = full.data,
          n.plt = n.plt,
          n.ssp = n.ssp[i],
          family = "quasibinomial",
          criterion = settings[[m]]$criterion,
          sampling.method = "poisson",
          objective.weight.plt = settings[[m]]$objective.weight.plt,
          objective.weight = settings[[m]]$objective.weight,
          balance.X.plt = settings[[m]]$balance.X.plt,
          rareFeature.index = NULL,
          rareThreshold = 0.05
        )
        list(
          coef.plt = fit$coef.plt,
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
          comp.time = fit$comp.time,
          fail.ind = 0
        )
      }, silent = TRUE)

      if (inherits(result, "try-error") ||
          is.null(result$coef.ssp) ||
          any(is.na(result$coef.ssp))) {
        result <- list(
          coef.plt = NA,
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
          comp.time = NA,
          fail.ind = 1
        )
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
      comp.time[i, m] <- result$comp.time
      fail.ind[i, m] <- result$fail.ind
    }
  }

  list(
    coef.plt = coef.plt,
    coef.ssp = coef.ssp,
    coef.cmb = coef.cmb,
    coef.full = coef.full,
    cov.ssp = cov.ssp,
    cov.cmb = cov.cmb,
    SubsampleSize = SubsampleSize,
    rare.count.plt = rare.count.plt,
    rare.count.ssp = rare.count.ssp,
    rare.count.cmb = rare.count.cmb,
    rare.counts.plt = rare.counts.plt,
    rare.counts.ssp = rare.counts.ssp,
    rare.counts.cmb = rare.counts.cmb,
    fail.ind = fail.ind,
    comp.time = comp.time
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

if (!is.null(cl)) {
  clusterEvalQ(cl, {
    source("rareF_main.R")
    source("rareF_functions.R")
  })

  clusterExport(
    cl = cl,
    varlist = c(
      "settings", "rareF_simu", "n.plt", "n.ssp", "num.method",
      "rareFeature.index", "full.data", "formula", "beta_full"
    ),
    envir = environment()
  )
}

t1 <- proc.time()
if (is.null(cl)) {
  setting.results <- lapply(rpt_indices, rareF_simu)
} else {
  setting.results <- parLapply(cl, rpt_indices, rareF_simu)
  stopCluster(cl)
}
t2 <- proc.time()
cat("Total runtime:", t2[3] - t1[3], "seconds\n")

date_tag <- format(Sys.Date(), "%m%d%Y")
save_dir <- if (!is.na(output_dir) && nzchar(output_dir)) {
  output_dir
} else {
  file.path("..", "raw_results", date_tag, "CoverType")
}
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

save_file <- file.path(save_dir, paste0("plt", n.plt, "_group", group_id, ".Rdata"))
save(setting.results, method.all, settings, combined.settings, file = save_file)
cat("Saved results to", save_file, "\n")
