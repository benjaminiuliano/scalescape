#simulation_scalescape
#' Simulate a distance-weighted model of landscape effects on an environmental response
#'
#' \code{simulate} method for class \code{scalescape}.Currently only available for linear models
#'
#' @param mod.full the full (landscape) model upon which to base the simulations
#' @param mod.reduced the reduced (local) model
#' @param new.coef simulated coefficient for the landscape effect
#' @param new.range simulated range parameter (i.e. scale of effect) for the landscape effect
#' @param new.sigma simulated standard deviation. By default, this is extracted from \code{mod.reduced}
#' @param data a data frame with local and landscape predictors and response variables.
#' By default, this is extracted from \code{mod.full}
#' @param nsim number of simulations (default = 1)
#'
#'
#' @return \code{simulation_scalescape} returns a data frame of predictor and simulated response variables.
#'
#' @export

simulation_scalescape <- function(mod.full, mod.reduced, new.coef, new.range, new.sigma = NULL, data = NULL, nsim = 1){

  landscape.formula <- mod.full$landscape.formula
  landscape.vars <- mod.full$landscape.vars
  weight.fn <- mod.full$weight.fn
  max.Dist <- mod.full$max.Dist
  if(class(mod.reduced)[1] == "scalescape") {
    mod.sim <- mod.reduced$mod
    mod0 <- mod.full$mod0
    data <- mod.reduced$data.with.landscape
  }else{
    mod.sim <- mod.reduced
    mod0 <- mod.reduced
    data <- mod.full$data
  }
  if(is.null(new.sigma)) new.sigma <- sigma(mod.sim)

  mod.reduced.vars <- row.names(attr(terms(formula(mod.reduced)), "factors"))
  mod.full.vars <- row.names(attr(terms(update(formula(mod.reduced), landscape.formula)), "factors"))
  new.vars <- mod.full.vars[!is.element(mod.full.vars, mod.reduced.vars)]

  max.Dist <- mod.full$max.Dist[is.element(new.vars, mod.full.vars)]
  new.range <- new.range/max.Dist

  data.sim <- data
  sim.df <- data.frame(i.sim=1:nsim, opt.range.boot=array(NA, dim=c(nsim, length(mod.full$opt.range))))

  # construct new variables with range.new and
  i.count <- 0
  for(i.terms in new.vars){
    if(is.element("matrix",is(landscape.vars[[i.terms]]))) {
      i.count <- i.count + 1

      Dist <- landscape.vars[[i.terms]][,1]
      if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/new.range[i.count]/max.Dist[i.count])^2)
      if(weight.fn=="exponential") weighting <- exp(-Dist/new.range[i.count]/max.Dist[i.count])

      data.sim[as.character(i.terms)] <- t(t(weighting) %*% landscape.vars[[i.terms]][, -1]/sum(weighting))
    }
    if(is.element("list",is(landscape.vars[[i.terms]]))) {
      i.count <- i.count + 1
      landscape.list <- landscape.vars[[i.terms]]

      Dist <- landscape.list[[1]][,1]
      if(weight.fn=="Gaussian") weighting <- exp(-0.5*(Dist/new.range[i.count]/max.Dist[i.count])^2)
      if(weight.fn=="exponential") weighting <- exp(-Dist/new.range[i.count]/max.Dist[i.count])
      weighting <- as.matrix(weighting)/mean(weighting)

      X.landclass.sub <- NULL
      for(i.list in 1:length(landscape.list)){
        X.landclass.sub <- cbind(X.landclass.sub, as.matrix(t(t(weighting) %*% landscape.list[[i.list]][, -1]/sum(weighting)), ncol=1))
      }
      if(!is.null(names(landscape.list))) colnames(X.landclass.sub) <- names(landscape.list)
      data.sim[as.character(i.terms)] <- X.landclass.sub
    }
  }
  data.new <- data.sim[,(ncol(data.sim)-i.count+1):ncol(data.sim)]

  sim.Y <- matrix(NA, nrow(data.new), nsim)
  #simulate new data
  if(is.element(class(mod.sim)[1], c("lm"))) {
    for(i.sim in 1:nsim)
      sim.Y[,i.sim] <- fitted(mod.sim) + as.matrix(data.new) %*% new.coef + rnorm(n=nrow(data.new), sd = new.sigma)
    sim.Y <- as.data.frame(sim.Y)
    names(sim.Y) <- paste0("sim",1:nsim)
    return(sim.Y)
  }
  if(is.element(class(mod.sim)[1], c("lme","gls","glm"))) {

    stop("Sorry, this doesn't work yet for glm, lme or gls")
    value <- mod.sim$modelStruct[[length(mod.sim$modelStruct)]]
    corGaus. <- Initialize(corGaus(value = value, form = ~ x + y, nugget = (length(value) > 1)), data=data)
    V <- mod0$sigma^2 * corMatrix(corGaus.)
    data.sim[, which(names(data) == formula(mod.sim)[[2]])] <- fitted(mod.sim.update) + t(mvtnorm::rmvnorm(n = 1, sigma = V))
  }
}
