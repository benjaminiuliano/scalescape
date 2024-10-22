# weighting_funct()

#' Internal function for distance weighting
#'
#' \code{weighting_funct} is used by \code{dist_weight} and \code{dist_weight_boot}
#' to conduct the distance weighting. It is an internal function and should not be
#' used on its own.
#'
#' @param par idk
#' @param mod0 Model object without landscape variables
#' @param landscape.vars A list of matrices generated by \code{landscape_matrix}
#'   corresponding to landcover variables of interest.
#' @param landscape.formula A character string specifying the formula containing
#'   landscape predictors
#' @param data A \code{data.frame} with local predictors and response variables,
#'   with each row as a site sorted in the same order as the site list used by
#'   \code{landscape_matrix}.
#' @param max.Dist maximum distance
#' @param weight.fn Specify the type of weighting function ("Gaussian" or "exponential").
#'   "Gausian" is the default.
#' @param return.coef Specify whether to return coeficients
#' @noRd

weighting_funct <- function(par, mod0, landscape.formula,
                            data, landscape.vars, max.Dist, weight.fn = "Gaussian", return.coef = F){

  new.vars <- names(landscape.vars)
  mod0.vars <- row.names(attr(terms(formula(mod0)), "factors"))
  mod.landscape.vars <- c(row.names(attr(terms(formula(mod0)), "factors")), new.vars)


  landscape.data <- data
  i.count <- 0
  for(i.terms in new.vars){
    if(is.element("matrix",is(landscape.vars[[i.terms]]))) {
      i.count <- i.count + 1

      Dist <- landscape.vars[[i.terms]][,1]
      if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/par[i.count]/max.Dist[i.count])^2)
      if(weight.fn=="exponential") weighting <- exp(-Dist/par[i.count]/max.Dist[i.count])

      landscape.data[as.character(i.terms)] <- t(t(weighting) %*% landscape.vars[[i.terms]][, -1]/sum(weighting))
    }
    if(is.element("list",is(landscape.vars[[i.terms]]))) {
      i.count <- i.count + 1
      landscape.list <- landscape.vars[[i.terms]]

      Dist <- landscape.list[[1]][,1]
      if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/par[i.count]/max.Dist[i.count])^2)
      if(weight.fn=="exponential") weighting <- exp(-Dist/par[i.count]/max.Dist[i.count])
      weighting <- as.matrix(weighting)/mean(weighting)

      X.landclass.sub <- NULL
      for(i.list in 1:length(landscape.list)){
        X.landclass.sub <- cbind(X.landclass.sub, as.matrix(t(t(weighting) %*% landscape.list[[i.list]][, -1]/sum(weighting)), ncol=1))
      }
      if(!is.null(names(landscape.list))) colnames(X.landclass.sub) <- names(landscape.list)
      landscape.data[as.character(i.terms)] <- X.landclass.sub

    }
  }
  formula <- formula(mod0)
  # show(par)
  mod <- try(update(mod0, landscape.formula, data=landscape.data))

  # if(is(mod)[1] == "try-error") return(-10^10)
  if(return.coef){
    if(!is.element(class(mod0)[1], c("lmerMod", "glmerMod"))) coef <- mod$coef else coef <- summary(mod)$coef[,1]
    return(list(logLik=logLik(mod), coef = coef, mod = mod, data.with.landscape = landscape.data))
  }else{
    return(-logLik(mod))
  }
}
