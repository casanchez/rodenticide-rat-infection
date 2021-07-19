## utility functions

##' extract the 'conditional-model' term from a glmmTMB object;
##' otherwise, return x unchanged
collapse_cond <- function(x)
  if (is.list(x) && "cond" %in% names(x)) x[["cond"]] else x

##' Cleaned-up/adapted version of Jon Lefcheck's code from SEMfit;
##' also incorporates some stuff from MuMIn::rsquaredGLMM.
##' Computes Nakagawa/Schielzeth/Johnson analogue of R^2 for
##' GLMMs. Should work for [g]lmer(.nb), glmmTMB models ...
##'
##' @param model a fitted model
##' @return a list composed of elements "family", "link", "marginal", "conditional"
my_rsq <- function(model) {
  
  ## get basics from model (as generally as possible)
  vals <- list(
    beta=fixef(model),
    X=getME(model,"X"),
    vc=VarCorr(model),
    re=ranef(model))
  
  ## glmmTMB-safety
  if (is(model,"glmmTMB")) {
    vals <- lapply(vals,collapse_cond)
    nullEnv <- function(x) {
      environment(x) <- NULL
      return(x)
    }
    if (!identical(nullEnv(model$modelInfo$allForm$ziformula),nullEnv(~0)))
      warning("R2 ignores effects of zero-inflation")
    dform <- nullEnv(model$modelInfo$allForm$dispformula)
    if (!identical(dform,nullEnv(~1)) &&
        (!identical(dform,nullEnv(~0))))
      warning("R2 ignores effects of dispersion model")
  }
  
  ## Test for non-zero random effects
  if (any(sapply(vals$vc, function(x) any(diag(x)==0)))) {
    ## FIXME: test more generally for singularity, via theta?
    stop("Some variance components equal zero. Respecify random structure!")
  }
  
  ## set family/link info
  ret <- list()
  if (is(model,"glmmTMB") || is(model,"glmerMod")) {
    ret$family <- family(model)$family
    ret$link <- family(model)$link
  } else {
    ret$family <- "gaussian"; ret$link <- "identity"
  }
  
  ## Get variance of fixed effects: multiply coefs by design matrix
  varF <- with(vals,var(as.vector(beta %*% t(X))))
  
  ## Are random slopes present as fixed effects? Warn.
  random.slopes <- if("list" %in% class(vals$re)) {
    ## multiple RE
    unique(c(sapply(vals$re,colnames)))
  } else {
    colnames(vals$re)
  }
  if (!all(random.slopes %in% names(vals$beta))) 
    warning("Random slopes not present as fixed effects. This artificially inflates the conditional R2. Respecify fixed structure!")
  
  ## Separate observation variance from variance of random effects
  nr <- sapply(vals$re, nrow)
  not.obs.terms <- names(nr[nr != nobs(model)])
  obs.terms <- names(nr[nr==nobs(model)])
  
  ## Compute variance associated with a random-effects term
  ## (Johnson 2014)
  getVarRand <- function(terms) {
    sum(
      sapply(vals$vc[terms],
             function(Sigma) {
               Z <- vals$X[, rownames(Sigma), drop = FALSE]
               Z.m <- Z %*% Sigma
               return(sum(diag(crossprod(Z.m, Z))) / nobs(model))
             } )
    )
  }
  
  ## Variance of random effects 
  varRand <- getVarRand(not.obs.terms)
  
  if (is(model,"lmerMod") ||
      (ret$family=="gaussian" && ret$link=="identity")) {
    ## Get residual variance
    varDist <- sigma(model)^2
    varDisp <- 0
  } else {
    varDisp <- if (length(obs.terms)==0) 0 else getVarRand(obs.terms)
    
    badlink <- function(link,family) {
      warning(sprintf("Model link '%s' is not yet supported for the %s distribution",link,family))
      return(NA)
    }
    
    if(ret$family == "binomial") {
      varDist <- switch(ret$link,
                        logit=pi^2/3,
                        probit=1,
                        badlink(ret$link,ret$family))
    } else if (ret$family == "poisson" ||
               grepl("nbinom",ret$family) ||
               grepl("Negative Binomial", ret$family)) {
      ## Generate null model (intercept and random effects only, no fixed effects)
      
      ## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q4/023013.html
      ## FIXME: deparse is a *little* dangerous
      rterms <- paste0("(",sapply(findbars(formula(model)),deparse),")")
      nullform <- reformulate(rterms,response=".")
      null.model <- update(model,nullform)
      
      ## from MuMIn::rsquaredGLMM
      
      ## Get the fixed effects of the null model
      null.fixef <- unname(collapse_cond(fixef(null.model)))
      
      ## in general want log(1+var(x)/mu^2)
      logVarDist <- function(null.fixef) {
        mu <- exp(null.fixef)
        if (mu < 6)
          warning(sprintf("mu of %0.1f is too close to zero, estimate may be unreliable \n",mu))
        vv <- switch(ret$family,
                     poisson=mu,
                     nbinom1=,
                     nbinom2=family(model)$variance(mu,sigma(model)),
                     if (is(model,"merMod"))
                       mu*(1+mu/getME(model,"glmer.nb.theta"))
                     else mu*(1+mu/model$theta))
        cvsquared <- vv/mu^2
        return(log1p(cvsquared))
      }
      
      varDist <- switch(ret$link,
                        log=logVarDist(null.fixef),
                        sqrt=0.25,
                        badlink(ret$link,ret$family))
    }
  }
  ## Calculate R2 values
  ret$Marginal = varF / (varF + varRand + varDisp + varDist)
  ret$Conditional = (varF + varRand) / (varF + varRand + varDisp + varDist)
  return(ret)
}
