#dist_weight()
#' Fit a distance-weighted model of landscape effects on an environmental response
#'
#' \code{dist_weight} fits a specified model containing local and landscape variables.
#'
#' \code{dist_weight} fits the model using the function \code{optim} (Brent for 1-D
#' and L-BFGS-B for >1D problems) to find the values of the model parameters that
#' maximize the log-likelihood of the model fit to the data. It can be run with many
#' types of regression models in R: \code{lm} and \code{glm} in base R,
#' \code{lmer} and \code{glmer} in \code{lme4}, and \code{lme} and \code{gls}
#' in \code{nlme.} It adopts the model syntax of the specified regression model,
#' making it easy to use models of any type.Although \code{dist_weight} produces
#' p-values for the regression coefficients (given by the underlying \code{lm},
#' \code{glm}, \code{lmer}, \code{glmer}, \code{lme}, or \code{gls} functions,
#' these p-values are conditional on the estimate of the range parameter, and
#' consequently they will likely have inflated type I error rates. The
#' \code{dist_weight_boot}function uses a bootstrap likelihood ratio test to generate
#' a single p-value for the landscape predictor variable(s) in the model. By
#' bootstrapping, it accounts for the co-dependence of regression coefficient and
#' range parameter. Therefore, p-values reported for landscape predictor(s) should come from
#' \code{dist_weight_boot} rather than \code{dist_weight}.
#'
#' @param mod0 a "local" model object, without landscape variables.
#' @param landscape.formula formula containing the landscape predictors.
#' @param landscape.vars list of names landscape matrices (one for each landscape variable).
#' @param data a data frame with local predictors and response variables, where
#' each row as a site sorted in the same order as the original site data frame.
#' @param weight.fn the type of weighting function to use; "Gaussian" (default) or "exponential".
#' @param plot.fits specify whether to produce plots (default = TRUE).
#' @param init.range starting point distance for fitting the range parameter.
#' @param opt.range specified value range parameter value. If specified, \code{dist_weight}.
#' will fit the model using that value as the value of the range parameter rather than using.
#' \code{optim} to identify the value that maximizes the log likelihood.
#' @param optim.method specify method for optimization (default = \code{"L-BFGS-B"}).
#' See \code{optim}.
#' @param lower lower bound on the variables for the \code{"L-BFGS-B"} method. See \code{optim}.
#' @param upper upper bound on the variables for the \code{"L-BFGS-B"} method. See \code{optim}.
#' @param n.partition number of partitions to divide the log-likelihood profile, in order
#' to avoid identifying false maxima.
#'
#'
#' @return \code{dist_weight} returns an object of class \code{scalescape}. This is a list
#' containing the following:
#' \itemize{
#'    \item the estimate of the range value
#'    \item the maximum likelihood value
#'    \item  AIC value for model fit
#'    \item  BIC value for model fit
#'    \item  the number of model parameters
#'    \item  the original (local) data frame used to fit the model
#'    \item  the coefficient of the landscape variable(s)
#'    \item  the landscape model, with distance-weighted landscape effect
#'    \item  the local model, without landscape effects
#'    \item  the formula specified to fit the landscape model
#'    \item  a data frame that contains the original local variable and distance-weighted landscape variable(s)
#'    \item  the list of landscape matrices for each landscape variable
#'    \item  the specified weighting function
#'    \item  the specified maximum distance
#'    }
#'
#' @export


dist_weight <- function(mod0, landscape.formula,
                        data = NULL, landscape.vars, init.range = NULL, opt.range = NULL, weight.fn = "Gaussian",
                        plot.fits = TRUE, optim.method = "L-BFGS-B", lower = NULL, upper = NULL, n.partition = 10){

  new.vars <- names(landscape.vars)
  mod0.vars <- row.names(attr(terms(formula(mod0)), "factors"))
  mod.landscape.vars <- c(row.names(attr(terms(formula(mod0)), "factors")), new.vars)

  if(!is.null(opt.range) & (length(opt.range) != length(new.vars))) stop("If you specify opt.range, the number of values must equal the number of new spatial elements estimated.")
  if(!is.null(init.range) & (length(init.range) != length(new.vars))) stop("If you specify init.range, the number of values must equal the number of new spatial elements estimated.")

  if(is.null(data) & is.element(class(mod0)[1], c("lm","glm"))) data <- mod0$model
  if(is.null(data) & is.element(class(mod0)[1], c("lmerMod", "glmerMod"))) data <- model.frame(mod0)
  if(is.null(data) & class(mod0)[1] == "lme") data <- mod0$data
  if(is.null(data) & class(mod0)[1] == "gls"){
    if(is.null(data)) stop("For gls() you need to specify data")
    data <- data
  }
  if(class(mod0)[1] == "gls" && mod0$method == "REML"){
    warning("For gls(), fitting show be with method = ML. Therefore, we refit your model.")
    mod0 <- update(mod0, method="ML")
  }

  max.Dist <- NULL
  for(i.terms in new.vars){
    if(is.element("matrix",is(landscape.vars[[i.terms]]))) {
      Dist <- landscape.vars[[i.terms]][,1]
      max.Dist <- c(max.Dist, max(Dist))
    }
    if(is.element("list",is(landscape.vars[[i.terms]]))) {
      landscape.list <- landscape.vars[[i.terms]]
      Dist <- landscape.list[[1]][,1]
      max.Dist <- c(max.Dist, max(Dist))
    }
  }

  if(is.null(opt.range)){
    if(is.null(init.range)) {
      init.range <- rep(0.2, length(new.vars))
      opt.init.range <- init.range
      z <- weighting_funct(par=init.range, mod0=mod0, landscape.formula = landscape.formula, data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T)
      opt.logLik <- z$logLik
      for(i.range in 1:length(new.vars)){
        for(i in 1:n.partition) {
          par <- opt.init.range
          par[i.range] <- i/(n.partition+1)
          z <- weighting_funct(par=par, mod0=mod0, landscape.formula = landscape.formula,
                               data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T)

          if(z$logLik > opt.logLik){
            opt.logLik <- z$logLik
            opt.init.range <- par
          }
          # show(c(z$logLik, par))
        }
      }
      init.range <- opt.init.range
    }else{
      init.range <- init.range/max.Dist
    }
    # if(length(init.range)==1){
    # optim.method="Brent"
    # }else{
    #	optim.method="L-BFGS-B"
    # }

    if(is.element(optim.method, c("Brent", "L-BFGS-B"))){
      if(is.null(lower)) {
        lower=rep(.01, length(new.vars))
      }else{
        lower=lower/max.Dist
      }
      if(is.null(upper)) {
        upper=rep(.99, length(new.vars))
      }else{
        upper=upper/max.Dist
      }

    }else{
      lower=NULL
      upper=NULL
    }
    opt.range <- optim(par=init.range, fn=weighting_funct, method=optim.method, upper = upper, lower = lower, mod0=mod0, data = data, max.Dist = max.Dist, landscape.formula = landscape.formula, landscape.vars=landscape.vars, weight.fn = weight.fn)$par
    names(opt.range) <- new.vars
  } else {
    opt.range <- opt.range/max.Dist
  }
  sol <- weighting_funct(opt.range, mod0 = mod0, landscape.formula = landscape.formula,
                         data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T)

  npar <- length(sol$coef) + length(opt.range) + 1
  AIC <- 2*npar - 2*sol$logLik
  BIC <- 2*npar*log(nrow(data)) - 2*sol$logLik

  # if(plot.fits){
  # par(mfrow=c(length(opt.range),2), mai=c(1,1,.3,.3))
  # for(i.range in 1:length(opt.range)){
  # w <- data.frame(range=max.Dist[i.range] * .005*(1:199), logLik = NA, b = NA)
  # for(i in 1:199) {
  # opt.par <- opt.range
  # opt.par[i.range] <- .005*i
  # z <- weighting_funct(par=opt.par, mod0=mod0, landscape.formula = landscape.formula,
  # data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T)
  # if(length(z)>1){
  # w$logLik[i] <- z$logLik
  # }
  # }

  # #1
  # plot(logLik ~ range, data=w, typ="l", xlab="Distance", ylab="logLik")
  # points(max.Dist[i.range] * opt.range[i.range], w$logLik, col="red")
  # mtext(side=3, paste0("Range for ",new.vars[i.range]," = ",round(max.Dist[i.range]*opt.range[i.range], digits=3)))

  # #2
  # # Gaussian weightings
  # if(weight.fn=="Gaussian") curve(exp(-0.5*(x/(max.Dist[i.range] * opt.range[i.range]))^2),
  # from = 0, to = max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
  # abline(v=max.Dist[i.range]*opt.range[i.range], col="red", lty=2)

  # # exponential weightings
  # if(weight.fn=="exponential") curve(exp(-x/(max.Dist[i.range] * opt.range[i.range])),
  # from = 0, to = max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
  # abline(v=max.Dist[i.range]*opt.range[i.range], col="red", lty=2)
  # }
  # }
  if(plot.fits){
    par(mfrow=c(length(opt.range),2), mai=c(1,1,.3,.3))
    for(i.range in 1:length(opt.range)){
      w <- data.frame(range=max.Dist[i.range] * .005*(1:199), logLik = NA, b = NA)
      for(i in 1:199) {
        opt.par <- opt.range
        opt.par[i.range] <- .005*i
        z <- weighting_funct(par=opt.par, mod0=mod0, landscape.formula = landscape.formula,
                             data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T)
        if(length(z)>1){
          w$logLik[i] <- z$logLik
          w$b[i] <- z$coef[3]
        }
      }

      #1
      plot(logLik ~ range, data=w, typ="l", xlab="Distance", ylab="logLik")
      points(max.Dist[i.range] * opt.range[i.range], sol$logLik, col="red")
      mtext(side=3, paste0("Range for ",new.vars[i.range]," = ",round(max.Dist[i.range]*opt.range[i.range], digits=3)))

      #2
      # Gaussian weightings
      if(weight.fn=="Gaussian") curve(exp(-0.5*(x/(max.Dist[i.range] * opt.range[i.range]))^2),
                                      from = 0, to = max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
      abline(v=max.Dist[i.range]*opt.range[i.range], col="red", lty=2)

      # exponential weightings
      if(weight.fn=="exponential") curve(exp(-x/(max.Dist[i.range] * opt.range[i.range])),
                                         from = 0, to = max.Dist[i.range], ylim=c(0,1), xlab="Distance", ylab="Weighting")
      abline(v=max.Dist[i.range]*opt.range[i.range], col="red", lty=2)
    }
  }

  to.return <- list(opt.range = max.Dist * opt.range,
                    logLik = sol$logLik,
                    AIC = AIC,
                    BIC = BIC,
                    npar = npar,
                    data = data,
                    coef = sol$coef,
                    mod = sol$mod,
                    mod0 = mod0,
                    landscape.formula = landscape.formula,
                    data.with.landscape = sol$data.with.landscape,
                    landscape.vars = landscape.vars,
                    weight.fn = weight.fn,
                    max.Dist = max.Dist)
  class(to.return) <- "scalescape"
  return(to.return)
}
