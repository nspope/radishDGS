#' @title Map simulations
#' @author Bill Peterman
#' @description Function to map simulated demes onto conductance surface
#' @param params Object generated from `get_params` OR results from `cdpop_sim`
#' @param log (Default = TRUE) The the log transformed conductance surface will be plotted
#' @param ... Not used
#' 
#' @return Plot of conductance surface with sampled demes
#' @export
#' 
#' @usage XXX
#' 
#' @examples 
#' ## Not Run:
#' ## map_cdpop(params,
#'              log = TRUE) ##  
#' 
#' ## End (Not run)
#' 
map_cdpop <- function(params,
                      log = TRUE){

  conduct <- params$conductance_surface
  pts <- params$pts

  if(isTRUE(log)){
    raster::plot(log(conduct), main = 'Sampled demes, log conductance')
  } else {
    raster::plot(conduct, main = 'Sampled demes, log conductance')
  }
  raster::plot(pts, add = TRUE, pch = 19)
}