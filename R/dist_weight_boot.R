#dist_weight_boot()
#' Perform a bootstrap likelihood ratio test on a scalescape object
#'
#' \code{dist_weight_boot} conducts a parametric bootstrap test of the null hypothesis that one or
#' more of the landscape predictors in a model fit using \code{dist_weight} has no
#' effect on the response variable.
#'
#' \code{dist_weight_boot} simulates datasets using the attributes and parameter values
#' from the model specified with dist_weight() and then re-fits the full and reduced models
#' to all the datasets using optim() to maximize the range parameter with each iteration.
#' It then compares them to test the null hypothesis that the reduced model fits better
#' than the full model, producing a p-value for the effect of the additional parameters in
#' the full model on the response. Although \code{dist_weight} produces p-values for the
#' regression coefficients, these p-values are conditional on the estimate of the range
#' parameter, and consequently they will likely have inflated type I error rates. By
#' bootstrapping, \code{dist_weight_boot} accounts for the co-dependence of regression
#' coefficient and range parameter. Therefore, p-values reported for landscape predictor(s)
#' should come from \code{dist_weight_boot} rather than \code{dist_weight}.
#'
#' @param mod.full a \code{scalescape} object; the landscape model fitted by \code{dist_weight}.
#' @param mod.reduced the corresponding "local" model object, without landscape variables.
#' @param nboot the number of bootstraps. To determine whether a variable is significant
#' at the alpha level 0.05, \code{nboot=2000} should be used (default).
#' @param plot.fits produce histograms of the bootstrapped range and deviance values (default = TRUE).
#' @param verbose produce output for each iteration of the bootstrap (default = FALSE)
#' @param pb.flag show a progress bar for completion of the bootstraps (default = TRUE)
#' @param n.breaks specify...
#' @param optim.method specify method for optimization (default = \code{"L-BFGS-B"}).
#' See \code{optim()}.
#' @param lower lower bound on the variables for the \code{"L-BFGS-B"} method. See \code{optim()}.
#' @param upper upper bound on the variables for the \code{"L-BFGS-B"} method. See \code{optim()}.
#' @param n.partition number of partitions to divide the log-likelihood profile, in order
#' to avoid identifying false maxima.
#' @param data a data frame with local and landscape predictors and response variables.
#' By default, this is extracted from \code{mod.full}.
#'
#' @return \code{dist_weight_boot} returns an object of class \code{scalescape.boot}. This is a list
#' containing the following:
#' \itemize{
#'   \item \code{mod.full} the full model
#'   \item \code{mod.reduced} the reduced model
#'   \item \code{dev} the deviance value of the full vs. reduced model
#'   \item \code{mean.dev} the mean deviance value of the bootstraps
#'   \item \code{sd.dev} the standard deviation of the bootstrap deviance values
#'   \item \code{P} the P value for the bootstrap likelihood ratio test
#'   \item \code{logLik.values} a data frame of log-likelihood and deviance values for each iteration of the bootstrap
#'   \item \code{coef} a data frame of coefficient values for each iteration of the bootstrap
#'   }
#'
#'   @export

dist_weight_boot <- function(mod.full, mod.reduced, nboot = 2000, plot.fits = TRUE,
                             verbose = FALSE, pb.flag = TRUE, n.breaks = NULL,
                             optim.method = "L-BFGS-B", lower = NULL, upper = NULL,
                             n.partition = 10, data = NULL){

  if(!class(mod.full)[1] == "scalescape") stop("The bootstrap must be performed on full model that is a scalescape object.")

  landscape.formula <- mod.full$landscape.formula
  landscape.vars <- mod.full$landscape.vars
  weight.fn <- mod.full$weight.fn
  max.Dist <- mod.full$max.Dist
  if(class(mod.reduced)[1] == "scalescape") {
    mod.sim <- mod.reduced$mod
    mod0 <- mod.full$mod0
  }else{
    mod.sim <- mod.reduced
    mod0 <- mod.reduced
  }
  data <- mod.full$data.with.landscape

  dev.obs <- mod.full$logLik - logLik(mod.sim)
  dat.boot <- data
  boot.df <- data.frame(i.boot=1:nboot, logLik = NA, logLik0 = NA, dev = NA)
  mod.opt.range <- mod.full$opt.range
  for(i in 1:length(mod.opt.range)) names(mod.opt.range)[i] <- paste0("range.",names(mod.opt.range)[i])
  boot.coef <- data.frame(matrix(NA, nrow=nboot, ncol=(length(mod.full$coef)+length(mod.opt.range))))
  colnames(boot.coef) <- c(names(mod.full$coef), names(mod.opt.range))

  if(pb.flag) pb <- txtProgressBar(min=0, max=nboot, initial = 0, style = 3)
  for(i.boot in 1:nboot){
    if(is.element(class(mod.sim)[1], c("lm","glm")))
      dat.boot[, which(names(data) == formula(mod.sim)[[2]])] <- as.matrix(simulate(mod.sim, nsim=1)[[1]])
    if(is.element(class(mod.sim)[1], c("lme","gls"))) {
      value <- mod.sim$modelStruct[[length(mod.sim$modelStruct)]]
      corGaus. <- Initialize(corGaus(value = value, form = ~ x + y, nugget = (length(value) > 1)), data=data)
      V <- mod0$sigma^2 * corMatrix(corGaus., corr=T)
      dat.boot[, which(names(data) == formula(mod.sim)[[2]])] <- fitted(mod.sim) + t(rmvnorm(n = 1, sigma = V))
    }

    mod0.boot <- try(update(mod0, data=dat.boot))
    if(is(mod0.boot)[1] == "try-error") {
      # show("broke here 1")
      break
    }

    init.range <- mod.full$opt.range/mod.full$max.Dist
    init.range[init.range < 0.05] <- 0.05
    opt.init.range <- init.range
    z <- weighting_funct(par=init.range, mod0=mod0, landscape.formula = landscape.formula, data = data, max.Dist = max.Dist, landscape.vars=landscape.vars, weight.fn = weight.fn, return.coef = T)
    opt.logLik <- z$logLik
    for(i.range in 1:length(init.range)){
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

    if(is.element(optim.method, c("Brent", "L-BFGS-B"))){
      if(is.null(lower)) {
        opt.lower=rep(.01, length(init.range))
      }else{
        opt.lower=lower/max.Dist
      }
      if(is.null(upper)) {
        opt.upper=rep(.99, length(init.range))
      }else{
        opt.upper=upper/max.Dist
      }

    }else{
      opt.lower=NULL
      opt.upper=NULL
    }
    opt.boot <- try(optim(par=init.range, fn=weighting_funct, method=optim.method,
                          upper = opt.upper, lower = opt.lower, mod0=mod0.boot, data = dat.boot,
                          landscape.formula = landscape.formula, landscape.vars = landscape.vars, weight.fn = weight.fn, max.Dist = mod.full$max.Dist), silent=T)

    if(is(opt.boot)[1] != "try-error" & !is.null(n.breaks)){
      for(i.break in 1:n.breaks){
        range.new <- i.break/(n.breaks+1)
        opt.boot.new <- try(optim(par=range.new, fn=weighting_funct, method=optim.method,
                                  upper = opt.upper, lower = opt.lower, mod0=mod0.boot, data = dat.boot,
                                  landscape.formula = landscape.formula, landscape.vars = landscape.vars, weight.fn = weight.fn, max.Dist = mod.full$max.Dist), silent=T)
        if((is(opt.boot.new)[1] != "try-error") & (opt.boot.new$value < opt.boot$value - 10^-5)) {
          warning("new optimal estimate found at a different range.init")
          # show(round(c(i.boot=i.boot,i.break=i.break, dif.logLik=opt.boot$value - opt.boot.new$value), digits=3))
          opt.boot <- opt.boot.new
        }
      }
    }

    if(is(opt.boot)[1] != "try-error") opt.range.boot <- opt.boot$par

    if(is(opt.boot)[1] != "try-error"){
      sol.boot <- weighting_funct(opt.range.boot, mod0=mod0.boot,
                                  landscape.formula = landscape.formula, data = dat.boot, landscape.vars=landscape.vars, weight.fn=weight.fn, max.Dist = mod.full$max.Dist, return.coef = T)
      if(length(opt.range.boot) == 1){
        boot.df$opt.range.boot[i.boot] <- mod.full$max.Dist * opt.range.boot
      }else{
        for(i.range in 1:length(opt.range.boot)) boot.df[i.boot,paste0("opt.range.boot.",i.range)] <- mod.full$max.Dist[i.range] * opt.range.boot[i.range]
      }
      boot.df$logLik[i.boot] <- sol.boot$logLik

      boot.df$logLik0[i.boot] <- logLik(mod0.boot)
      boot.df$dev[i.boot] <- boot.df$logLik[i.boot] - boot.df$logLik0[i.boot]

      boot.coef[i.boot,] <- c(sol.boot$coef, opt.range.boot)
    }else{
      boot.df$opt.range.boot[i.boot] <- NA
      boot.df$logLik[i.boot] <- NA

      boot.df$logLik0[i.boot] <- NA
      boot.df$dev[i.boot] <- NA

      boot.coef[i.boot,] <- NA
    }
    if(pb.flag) setTxtProgressBar(pb,i.boot)
    if(verbose) show(boot.df[i.boot,])
  }
  boot.df <- boot.df[!is.na(boot.df$logLik),]
  boot.mean.dev <- mean(boot.df$dev, na.rm=T)
  boot.sd.dev <- sd(boot.df$dev, na.rm=T)
  dev.TF <- (dev.obs < boot.df$dev)
  boot.P <- (sum(dev.TF, na.rm=T) + 1)/(nboot + 1)
  if(length(opt.range.boot) == 1){
    boot.opt.range <- mean(boot.df$opt.range.boot, na.rm=T)
  }else{
    boot.opt.range <- NULL
    for(i.range in 1:length(opt.range.boot)) boot.opt.range <- c(boot.opt.range, mean(boot.df[,paste0("opt.range.boot.",i.range)], na.rm=T))
  }

  # plot results
  if(plot.fits){
    par(mfcol=c(1,1 + length(mod.full$opt.range)), mai=c(1,1,.3,.3))

    #1
    if(length(mod.full$opt.range) == 1){
      hist(max.Dist * boot.df$opt.range.boot, xlab="Range", main=paste0("Bootstrap for ", names(mod.full$opt.range)))
    }else{
      for(i.range in 1:length(mod.full$opt.range)) hist(max.Dist[i.range] * boot.df[,paste0("opt.range.boot.",i.range)], xlab="Range", main=paste0("Bootstrap for ", names(mod.full$opt.range)[i.range]))
    }

    #2
    hist(boot.df$dev, xlab="Deviance", main=paste0("Bootstrap P = ", round(boot.P, digits=4)), xlim=c(-.01,max(c(boot.df$dev,dev.obs))))
    lines(c(dev.obs,dev.obs), c(0,nboot), col="red")
  }

  to.return <- list(mod.full = mod.full,
                    mod.reduced = mod.reduced,
                    dev = dev.obs,
                    mean.dev = boot.mean.dev,
                    sd.dev = boot.sd.dev,
                    P = boot.P,
                    logLik.values = boot.df,
                    coef = data.frame(boot.coef))
  class(to.return) <- c("scalescape", "scalescape.boot")
  return(to.return)
}
