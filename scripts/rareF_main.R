ssp.rareF.glm <- function(formula,
                    data,
                    subset = NULL,
                    n.plt,
                    n.ssp,
                    family = 'binomial',
                    criterion = 'Lopt',
                    sampling.method = 'poisson',
                    objective.weight.plt = 'weighted',
                    objective.weight = 'weighted',
                    control = list(...),
                    contrasts = NULL,
                    balance.X.plt = FALSE,
                    balance.Y.plt = FALSE,
                    balance.Y.ssp = FALSE,
                    balance.Y.all = FALSE,
                    combine = TRUE,
                    record.stage.time = FALSE,
                    rareFeature.index = NULL,
                    rareThreshold = 0.09,
                    na.action = getOption("na.action"),
                    ...
                    ) {
  
  model.call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"),
             names(mf),
             0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  Y <- model.response(mf, "any")
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  
  X <- model.matrix(mt, mf, contrasts) 
  
  if (is.character(family)) {
    family_lookup <- c(
      binomial = "binomial",
      poisson = "poisson",
      Gamma = "Gamma",
      gaussian = "gaussian",
      Gaussian = "gaussian",
      quasibinomial = "quasibinomial",
      quasipoisson = "quasipoisson"
    )
    family <- match.arg(family, names(family_lookup))
    family <- family_lookup[[family]]
  }
  criterion <- match.arg(criterion, c('Lopt', 'Aopt', 
                                      'Uni', 'BL-Uni',
                                      'R-Lopt', 'BL-Lopt'))
  sampling.method <- match.arg(sampling.method, c('poisson'))
  objective.weight.plt <- match.arg(
    objective.weight.plt,
    c('weighted', 'unweighted')
  )
  objective.weight <- match.arg(
    objective.weight,
    c('weighted', 'unweighted')
  )
  control <- do.call("glm.control", control)
  
  ## family 
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    stop("'family' not recognized")
  }
  
  if (!(family$family %in% c("binomial", "quasibinomial"))) {
    requested_y_balancing <- balance.Y.plt || balance.Y.ssp || balance.Y.all
    balance.Y.plt <- FALSE
    balance.Y.ssp <- FALSE
    balance.Y.all <- FALSE
    if (requested_y_balancing) {
      warning("only consider balancing Y when Y is binary")
    }
  }
  
  if (criterion %in% c("Lopt", "Aopt", "R-Lopt", "BL-Lopt") && balance.Y.ssp) {
    warning(
      paste0(
        "'balance.Y.ssp' is currently only used for 'Uni' and 'BL-Uni'; ",
        "ignoring it for criterion '", criterion, "'."
      )
    )
    balance.Y.ssp <- FALSE
  }

  is.two.stage <- criterion %in% c("Lopt", "Aopt", "R-Lopt", "BL-Lopt")
  is.one.stage <- criterion %in% c("Uni", "BL-Uni")

  if (is.one.stage &&
      objective.weight == "unweighted" &&
      (balance.Y.ssp || balance.Y.all)) {
    stop(
      paste0(
        "'objective.weight = \"unweighted\"' is only valid for one-step methods ",
        "when the sampling probability does not depend on Y."
      )
    )
  }

  if (is.two.stage &&
      objective.weight.plt == "unweighted" &&
      balance.Y.plt) {
    stop(
      paste0(
        "'objective.weight.plt = \"unweighted\"' is not allowed when ",
        "'balance.Y.plt = TRUE' because the pilot sampling depends on Y."
      )
    )
  }

  if (is.two.stage && objective.weight != "weighted") {
    stop(
      paste0(
        "'objective.weight' must be 'weighted' for two-step methods ",
        "('Lopt', 'Aopt', 'R-Lopt', 'BL-Lopt')."
      )
    )
  }
  
  N <- nrow(X)
  d <- ncol(X)
  subsample.size.expect <- n.ssp
  full.Y.sum <- sum(Y, na.rm = TRUE)
  full.Y.mean <- mean(Y, na.rm = TRUE)
  stage.time <- numeric(0)
  
  ## compute balancing score
  stage_start <- proc.time()[3]
  bl.information <- balance_score(X, Y, family, 
                                  rareFeature.index = rareFeature.index,
                                  threshold = rareThreshold)
  stage.time["balance_score"] <- proc.time()[3] - stage_start
  bl <- bl.information$bl
  DN <- bl.information$DN
  full.rare.count <- bl.information$rare.count
  Y.count <- bl.information$Y.count
  binary.index <- bl.information$binary.index
  rareFeature.index <- bl.information$rareFeature.index
  Y_subgroup_prev <- bl.information$Y_subgroup_prev
  if (attr(mt, "intercept") == 1) {
    colnames(X)[1] <- "Intercept"
  }

  inputs <- list(X = X, Y = Y, N = N, d = d,
                 Y.count = Y.count,
                 n.plt = n.plt, n.ssp = n.ssp,
                 criterion = criterion, sampling.method = sampling.method,
                 objective.weight.plt = objective.weight.plt,
                 objective.weight = objective.weight,
                 family = family,
                 bl = bl, DN = DN,
                 balance.X.plt = balance.X.plt,
                 balance.Y.plt = balance.Y.plt,
                 balance.Y.ssp = balance.Y.ssp,
                 balance.Y.all = balance.Y.all,
                 combine = combine,
                 rareFeature.index = rareFeature.index, 
                 control = control
                 )
  
  if (criterion %in% c('Lopt', 'Aopt', 'R-Lopt', 'BL-Lopt')) {
    t_start <- proc.time()[3]
    ## pilot step
    stage_start <- proc.time()[3]
    plt.estimate.results <- pilot.estimate(inputs, 
                                           ...
                                           )
    stage.time["pilot.estimate"] <- proc.time()[3] - stage_start
    pi.plt <- plt.estimate.results$pi.plt
    pi.plt.N <- plt.estimate.results$pi.plt.N
    varphi.plt = plt.estimate.results$varphi.plt
    beta.plt <- plt.estimate.results$beta.plt
    ddL.plt <- plt.estimate.results$ddL.plt
    ddL.plt.design <- plt.estimate.results$ddL.plt.design
    if (is.null(ddL.plt.design)) {
      ddL.plt.design <- ddL.plt
    }
    cov.score.plt <- plt.estimate.results$cov.score.plt
    linear.predictor <- plt.estimate.results$linear.predictor # dimension: N
    index.plt <- plt.estimate.results$index.plt
    cov.plt <- plt.estimate.results$cov.plt
    n.plt <- length(index.plt)

    ## subsampling step
    stage_start <- proc.time()[3]
    ssp.results <- subsampling(inputs,
                               pi.plt = pi.plt,
                               varphi.plt = varphi.plt,
                               ddL.plt = ddL.plt.design,
                               linear.predictor = linear.predictor,
                               index.plt = index.plt
                               )
    stage.time["subsampling"] <- proc.time()[3] - stage_start
    index.ssp <- ssp.results$index.ssp
    w.ssp <- ssp.results$w.ssp
    varphi.ssp <- ssp.results$varphi.ssp
    pi.ssp <- ssp.results$pi.ssp
    pi.ssp.N <- ssp.results$pi.ssp.N
    n.ssp <- length(index.ssp)
    
    ## subsample estimating step
    stage_start <- proc.time()[3]
    ssp.estimate.results <- subsample.estimate(inputs,
                                               w.ssp = w.ssp,
                                               varphi.ssp = varphi.ssp,
                                               beta.plt = beta.plt,
                                               index.ssp = index.ssp,
                                               ...
                                               )
    stage.time["subsample.estimate"] <- proc.time()[3] - stage_start
    t_end <- proc.time()[3]
    comp.time <- t_end - t_start
    beta.ssp <- ssp.estimate.results$beta.ssp
    ddL.ssp <- ssp.estimate.results$ddL.ssp
    cov.score.ssp <- ssp.estimate.results$cov.score.ssp
    cov.ssp <- ssp.estimate.results$cov.ssp
    
    
    ## combining step: 
    stage_start <- proc.time()[3]
    combining.union.results <- combining.union(inputs,
                                   index.plt = index.plt, # dim: n.plt
                                   index.ssp = index.ssp, # dim: n.ssp
                                   pi.plt = pi.plt.N, # dim: N
                                   pi.ssp = pi.ssp.N # dim: N
                                   )
    stage.time["combining.union"] <- proc.time()[3] - stage_start
    index.cmb.union <- combining.union.results$index.cmb
    beta.cmb.union <- combining.union.results$beta.cmb
    cov.cmb.union <- combining.union.results$cov.cmb

    ## prepare for displaying results 
    names(beta.cmb.union) <- names(beta.ssp) <- names(beta.plt) <- colnames(X)
    
    
    
    stage_start <- proc.time()[3]
    rare.counts.results.plt <- rare_counts(X, 
                                           index.plt, 
                                           rareFeature.index,
                                           Y)
    rare.counts.results.ssp <- rare_counts(X, 
                                           index.ssp, 
                                           rareFeature.index,
                                           Y)
    rare.counts.results.cmb.union <- rare_counts(X, 
                                           index.cmb.union, 
                                           rareFeature.index,
                                           Y)
    stage.time["result.assembly"] <- proc.time()[3] - stage_start
    
    results <- list(model.call = model.call,
                    coef.plt = beta.plt,
                    coef.ssp = beta.ssp,
                    coef.cmb.union = beta.cmb.union,
                    cov.plt = cov.plt,
                    cov.ssp = cov.ssp,
                    cov.cmb.union = cov.cmb.union,
                    N = N,
                    family.name = family$family,
                    subsample.size.expect = subsample.size.expect,
                    subsample.size.actual = length(index.ssp),
                    full.Y.count = Y.count,
                    full.Y.sum = full.Y.sum,
                    full.Y.mean = full.Y.mean,
                    full.rare.count = full.rare.count,
                    Y_subgroup_prev = Y_subgroup_prev,
                    
                    Y.count.plt = rare.counts.results.plt$Y.count,
                    Y.proportion.plt = rare.counts.results.plt$Y.proportion,
                    rare.count.plt = rare.counts.results.plt$rare.count,
                    rare.proportion.plt = rare.counts.results.plt$rare.proportion,
                    rare.counts.plt = rare.counts.results.plt$rare.counts,
                    rows.with.rare.plt = rare.counts.results.plt$rows.with.rare,
                    
                    Y.count.ssp = rare.counts.results.ssp$Y.count,
                    Y.proportion.ssp = rare.counts.results.ssp$Y.proportion,
                    rare.count.ssp = rare.counts.results.ssp$rare.count,
                    rare.proportion.ssp = rare.counts.results.ssp$rare.proportion,
                    rare.counts.ssp = rare.counts.results.ssp$rare.counts,
                    rows.with.rare.ssp = rare.counts.results.ssp$rows.with.rare,
                    
                    Y.count.cmb.union = rare.counts.results.cmb.union$Y.count,
                    Y.proportion.cmb.union = rare.counts.results.cmb.union$Y.proportion,
                    rare.count.cmb.union = rare.counts.results.cmb.union$rare.count,
                    rare.proportion.cmb.union = rare.counts.results.cmb.union$rare.proportion,
                    rare.counts.cmb.union = rare.counts.results.cmb.union$rare.counts,
                    rows.with.rare.cmb.union = rare.counts.results.cmb.union$rows.with.rare,
                    
                    index.plt = index.plt,
                    index.ssp = index.ssp,
                    index.cmb.union = index.cmb.union,
                    rareFeature.index = rareFeature.index,
                    comp.time = comp.time,
                    stage.time = if (record.stage.time) stage.time else NULL,
                    terms = mt
                    )
    class(results) <- c("ssp.rareF.glm", "list")
    return(results)
  } else if (criterion %in% c("Uni", "BL-Uni")){
    if (criterion == "Uni") inputs$balance.X.plt = FALSE
    inputs$n.uni <- n.plt + n.ssp
    t_start <- proc.time()[3]
    
    stage_start <- proc.time()[3]
    uni.estimate.results <- uniform.estimate(inputs, ...)
    stage.time["uniform.estimate"] <- proc.time()[3] - stage_start
    
    t_end <- proc.time()[3]
    comp.time <- t_end - t_start
    pi.uni <- uni.estimate.results$pi.uni
    varphi.uni = uni.estimate.results$varphi.uni
    beta.uni <- uni.estimate.results$beta.uni
    ddL.uni <- uni.estimate.results$ddL.uni
    cov.score.uni <- uni.estimate.results$cov.score.uni
    linear.predictor <- uni.estimate.results$linear.predictor
    index.uni <- uni.estimate.results$index.uni
    cov.uni <- uni.estimate.results$cov.uni


    stage_start <- proc.time()[3]
    rare.counts.results <- rare_counts(X, index.uni, rareFeature.index, Y)
    stage.time["result.assembly"] <- proc.time()[3] - stage_start
    
    
    results <- list(model.call = model.call,
                    coef.plt = beta.uni,
                    coef.ssp = beta.uni,
                    coef.cmb.union = beta.uni,
                    cov.plt = cov.uni,
                    cov.ssp = cov.uni,
                    cov.cmb.union = cov.uni,
                    N = N,
                    family.name = family$family,
                    subsample.size.expect = inputs$n.uni,
                    subsample.size.actual = length(index.uni),
                    full.Y.count = Y.count,
                    full.Y.sum = full.Y.sum,
                    full.Y.mean = full.Y.mean,
                    full.rare.count = full.rare.count,
                    Y_subgroup_prev = Y_subgroup_prev,
                    # plt, ssp and cmb are the same.
                    Y.count.plt = rare.counts.results$Y.count,
                    Y.proportion.plt = rare.counts.results$Y.proportion,

                    rare.count.plt = rare.counts.results$rare.count,
                    rare.proportion.plt = rare.counts.results$rare.proportion,
                    rare.counts.plt = rare.counts.results$rare.counts,
                    rows.with.rare.plt = rare.counts.results$rows.with.rare,
                    
                    Y.count.ssp = rare.counts.results$Y.count,
                    Y.proportion.ssp = rare.counts.results$Y.proportion,
                    rare.count.ssp = rare.counts.results$rare.count,
                    rare.proportion.ssp = rare.counts.results$rare.proportion,
                    rare.counts.ssp = rare.counts.results$rare.counts,
                    rows.with.rare.ssp = rare.counts.results$rows.with.rare,
                    
                    Y.count.cmb.union = rare.counts.results$Y.count,
                    Y.proportion.cmb.union = rare.counts.results$Y.proportion,
                    rare.count.cmb.union = rare.counts.results$rare.count,
                    rare.proportion.cmb.union = rare.counts.results$rare.proportion,
                    rare.counts.cmb.union = rare.counts.results$rare.counts,
                    rows.with.rare.cmb.union = rare.counts.results$rows.with.rare,
                    
                    index.plt = index.uni,
                    index.ssp = index.uni,
                    index.cmb.union = index.uni,
                    rareFeature.index = rareFeature.index,
                    comp.time = comp.time,
                    stage.time = if (record.stage.time) stage.time else NULL,
                    terms = mt
    )
    class(results) <- c("ssp.rareF.glm", "list")
    return(results)
  }
}


###############################################################################
