#landscape_matrix()
#' Create a matrix for the landscape surrounding each site
#'
#' \code{landscape_matrix} creates the matrix that is used by \code{dist_weight}
#'  to weight the landscape variable and fit the model.
#'
#' Each row of the \code{landscape_matrix} output corresponds to a raster pixel at a given
#' location around a site, and each column corresponds to a site. The column
#' 'dist' gives the distance of each pixel to the focal site. For binary
#' classifications, cell values in the matrix represent whether that pixel is
#' filled by a land cover type of interest (1 if yes, 0 if no), and for continuous
#' classifications they correspond to the value of the landscape variable of
#' interest. When applying \code{landscape_matrix} to a raster, you will need to specify:
#' \itemize{
#'    \item The raster associated with the landscape variable you are interested in.
#'    \item A list of sites formatted as a 2-column X-Y data frame, \code{sf} object,
#'    or \code{sp} object.
#'    \item The maximum radius you wish to evaluate.This radius should be larger than
#'    what you believe to be the maximum relevant spatial scale for your response variable
#'    (although the larger it is, the longer the computation time). All raster cells up to
#'    this distance will be included in the parameter optimization process used to define
#'    the range parameter. The value for \code{max.radius} must fall within the bounds of
#'    the raster. If \code{max.radius} around at least one site falls outside the bounds of
#'    your raster, you will receive an error message indicating this problem and should
#'    select a raster with a greater extent or choose a smaller \code{max.radius}.
#'}
#'
#' @param raster A raster object associated with the landscape variable of interest.
#' @param sites An inventory of sites formatted as a 2-column, X-Y data frame,
#'   sf object, or sp object.
#' @param max.radius The maximum radius you wish to evaluate. Must be less than the
#'   maximum extent of the raster.
#' @param is.factor Specify whether landscape variable is continuous or binary (FALSE, default)
#' or a factor with more than two levels (TRUE).
#'
#' @return A matrix used by \code{dist_weight} to do the weighting and fit the model.
#'
#' @export

landscape_matrix <- function(raster, sites, max.radius, is.factor = FALSE) {
  max.radius <- ceiling(max.radius/res(raster)[1])*res(raster)[1]
  n.rows <- 2 * max.radius/res(raster)[1]

  if (!is.factor || length(unique(raster)) == 2) {
    dat.site <- data.frame(X = 2*rep(0:(n.rows - 1), times = n.rows)/(n.rows - 1), Y = 2*rep(0:(n.rows - 1), each = n.rows)/(n.rows - 1))
    dat.site$dist <- apply(dat.site, 1, FUN = function(x) sqrt((x[1] - 1)^2 + (x[2] - 1)^2))
    dat.site$dist <- dat.site$dist/dat.site$dist[max.radius/res(raster)[1]]

    dat.site <- dat.site[dat.site$dist <= 1, ]
    landscape.matrix <- as.matrix(dat.site)
    landscape.matrix <- landscape.matrix[,-c(1,2)]
    landscape.matrix[,1] <- max.radius * landscape.matrix[,1]
    return(landscape.matrix)
  } else {
    landscape.list <- list()
    for (i.level in 2:length(unique(raster))) {
      dat.site <- data.frame(X = 2*rep(0:(n.rows - 1), times = n.rows)/(n.rows - 1), Y = 2*rep(0:(n.rows - 1), each = n.rows)/(n.rows - 1))
      dat.site$dist <- apply(dat.site, 1, FUN = function(x) sqrt((x[1] - 1)^2 + (x[2] - 1)^2))
      dat.site$dist <- dat.site$dist/dat.site$dist[max.radius/res(raster)[1]]

      if (class(sites)[1] == "data.frame") {
        for (i in 1:nrow(sites)) {
          site <- sites[i, ]
          x.min <- site[1,1] - max.radius
          x.max <- site[1,1] + max.radius
          y.min <- site[1,2] - max.radius
          y.max <- site[1,2] + max.radius
          extent.site <- extent(x.min, x.max, y.min, y.max)
          site.crop <- crop(raster, extent.site)
          site.matrix <- t(as.matrix(site.crop, ncol = n.rows))
          if (nrow(dat.site) > nrow(matrix(site.matrix, ncol=1))) stop("Maximum radius is outside the bounds of raster extent. Choose a smaller value for max.radius")
          dat.site <- cbind(dat.site, matrix(site.matrix, ncol = 1))
          names(dat.site)[i + 3] <- paste0("landclass.", i)
        }
      }
      if (class(sites)[1] == "sf") {
        for (i in 1:nrow(sites)) {
          site <- st_coordinates(sites)[i, ]
          x.min <- site[1] - max.radius
          x.max <- site[1] + max.radius
          y.min <- site[2] - max.radius
          y.max <- site[2] + max.radius
          extent.site <- extent(x.min, x.max, y.min, y.max)
          site.crop <- crop(raster, extent.site)
          site.matrix <- t(as.matrix(site.crop, ncol = n.rows))
          if (nrow(dat.site) > nrow(matrix(site.matrix, ncol=1))) stop("Maximum radius is outside the bounds of raster extent. Choose a smaller value for max.radius")
          dat.site <- cbind(dat.site, matrix(site.matrix, ncol = 1))
          names(dat.site)[i + 3] <- paste0("landclass.", i)
        }
      }
      if (class(sites)[1] %in% c("SpatialPoints", "SpatialPointsDataFrame")) {
        for (i in 1:length(sites)) {
          site <- sites[i, ]
          x.min <- site@coords[1, 1] - max.radius
          x.max <- site@coords[1, 1] + max.radius
          y.min <- site@coords[1, 2] - max.radius
          y.max <- site@coords[1, 2] + max.radius
          extent.site <- extent(x.min, x.max, y.min, y.max)
          site.crop <- crop(raster, extent.site)
          site.matrix <- t(as.matrix(site.crop, ncol = n.rows))
          if (nrow(dat.site) > nrow(matrix(site.matrix, ncol=1))) stop("Maximum radius is outside the bounds of raster extent. Choose a smaller value for max.radius")
          dat.site <- cbind(dat.site, matrix(site.matrix, ncol = 1))
          names(dat.site)[i + 3] <- paste0("landclass.", i)
        }
      }
      dat.site <- dat.site[dat.site$dist <= 1, ]
      dat.site[round(nrow(dat.site)/2) + 1, ]
      landscape.matrix <- as.matrix(dat.site)
      landscape.matrix <- landscape.matrix[,-c(1,2)]
      landscape.matrix[,1] <- max.radius * landscape.matrix[,1]
      landscape.list[[i.level-1]] <- landscape.matrix
    }
    return(landscape.list)
  }
}

