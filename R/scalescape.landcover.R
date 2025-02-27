#' Natural habitat around sampling sites
#'
#' A binary raster object showing the amount of natural habitat around  wild bee sampling locations. Example from Lowe et al. 2021.
#'
#' @docType data
#'
#' @usage data(scalescape.landcover)
#'
#' @format An object of class \code{RasterLayer} with the following characteristics:
#' \describe{
#'   \item{**dimensions**}{nrow=1539, ncol=1681, ncell=2587059}
#'   \item{**resolution**}{x=30m, y=30m}
#'   \item{**extent**}{xmin=547169, xmax=597599, ymin=389028.8, ymax=435198.8}
#'   \item{**crs**}{+proj=tmerc +lat_0=0 +lon_0=-90 +k=0.9996 +x_0=520000 +y_0=-4480000 +ellps=GRS80 +units=m +no_defs}
#'   \item{**values**}{0=non-natural habitat, 1=natural and semi-natural habitat}
#'   }
#'
#' @source https://nassgeodata.gmu.edu/CropScape/, Erin Lowe
"scalescape.landcover"
