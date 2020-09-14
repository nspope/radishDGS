#' @author Bill Peterman
#' @title Generate correlated rater surfaces
#' @description Function to create correlated surfaces
#' 
#' @param corr Desired correlation level between fbm surfaces
#' @param FUN NLMR function to create random surface
#' @param ... Specify all necessary parameters of the function `FUN`
#' 
#' @examples 
#' ## Not Run:
#' library(raster)
#' cr <- corr_rast(corr = 0.5,
#'                 FUN = NLMR::nlm_fbm,
#'                 ncol = 50,
#'                 nrow = 50,
#'                 fract_dim = 0.5,
#'                 user_seed = 555) 
#' plot(cr) 
#' 
#' ## End (Not run)
#' @export

corr_rast <- function(corr,
                      FUN,
                      ...){
  
  sim_rast <- FUN(...)
  
  rast_corr <- 0
  
  args <- list(...)
  args[['user_seed']] <- NA
  
  while(rast_corr > corr + 0.015 | rast_corr < corr - 0.015){
    rep_sim <- do.call(FUN, args)
    
    corr_ <-  sqrt((1 / (corr^2)) - 1) 
    mat <- raster::as.matrix(sim_rast)
    
    rep_sim_ <- rep_sim * corr_
    corr_rast <- sim_rast + rep_sim_
    
    rast_corr <- raster::layerStats(raster::stack(sim_rast, corr_rast), 'pearson')$`pearson correlation coefficient`[1,2]
    cat(paste0(round(rast_corr, digits = 4), '\n \n'))
  }
  
  cat(paste0("Surface correlation is ", round(rast_corr, 4), '\n \n'))
  
  
  corr_rast[] <- (corr_rast[] - min(corr_rast[])) / (max(corr_rast[]) - min(corr_rast[]))
  
  
  return(raster::stack(list(sim_true = sim_rast,
                            corr_sim = corr_rast)))
}