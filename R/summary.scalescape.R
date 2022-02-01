#summary.scalescape
#' Summarizing scalescape model fits
#'
#' \code{summary} method for classes \code{scalescape} and \code{scalescape.boot}
#'
#' Note that although \code{summary.scalescape} displays p-values for the regression coefficients
#' (given by the underlying \code{lm()},\code{glm}, \code{lmer}, \code{glmer}, \code{lme},
#' or \code{gls} functions, p-values for landscape effets are conditional on the estimate
#' of the range parameter, and consequently they will likely have inflated type I error
#' rates. The \code{dist_weight_boot()} function uses a bootstrap likelihood ratio test to
#' generate a single p-value for the landscape predictor variable(s) in the model. By
#' bootstrapping, it accounts for the co-dependence of regression coefficient and range
#' parameter. Therefore, p-values reported for landscape predictor(s) should come from
#' \code{dist_weight_boot} rather than \code{dist_weight}.
#'
#' @param object an object of class \code{scalescape} or \code{scalescape.boot}
#' @param ... additional arguments affecting the summary produced
#'
#'
#' @return For objects of class \code{scalescape},  \code{summary.scalescape} returns a list of
#' summary statistics of the fitted model given in \code{object}, including:
#' \itemize{
#'    \item the estimated range parameter, in meters
#'    \item metrics of model fit (log-likelihood, AIC, BIC), and the number of parameters
#'    \item the call for the local (no landscape variables) model
#'    \item the call for the distance-weighted landscape model
#'    \item residuals
#'    \item coefficients, standard errors, t-statistics and corresponding p-values
#' }
#'  For objects of class \code{scalescape.boot},  \code{summary.scalescape} returns a list of
#'  summary statistics of the bootstrap likelihood ratio test comparing the full and reduced
#'  models, including:
#' \itemize{
#'   \item the number of bootstrapped datasets successfully refit
#'   \item the call for the full (e.g., landscape) model
#'   \item the call for the reduced (e.g., local) model
#'   \item the observed deviance, bootstrap deviance, standard deviation, and corresponding p-value
#'   }
#'
#' @export


summary.scalescape <- function(object, ...) {

  if(!is.element("scalescape.boot", class(object))) {

    cat("\nCall: scalescape\n")
    cat("\nEstimates of range for spatial variables\n")
    show(object$opt.range)

    fit.table <- data.frame(logLik = object$logLik, AIC = object$AIC, BIC = object$BIC, npar = object$npar)
    cat("\nOverall model fit:\n")
    print(round(fit.table, digits=3))

    cat("\n\nLocal model call: ")
    print(formula(object$mod0))

    cat("\nFitted model conditional on the estimates of the range coefficients\n")
    show(summary(object$mod))

    cat("\nWarning: the tests of the spatial coefficients are conditional on the optimal fit of the range. These P-values should not be trusted. P-values should be obtained by bootstrapping using dist_weight_boot\n", fill=TRUE)

  }else{
    cat("\nBootstrapping of coefficients involving spatial terms\n")
    num.boots <- sum(!is.na(object$bootstrap.values$opt.range.boot))
    cat("\nNumber of bootstrapped datasets successfully refit:", num.boots)
    cat("\n\nFull model:")
    print(formula(object$mod.full$mod))
    cat("\nReduced model:")
    if(class(object$mod.reduced) == "scalescape") print(object$mod.reduced$call) else print(formula(object$mod.reduced))
    bootstrap.output <- c(object$dev, object$mean.dev, object$sd.dev, object$P)
    names(bootstrap.output) <- c("Observed Dev", "Bootstrap Dev", "SD", "Pr")
    cat("\nBootstrap:\n")
    print(bootstrap.output)
  }
}
