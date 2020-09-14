#' @title Gaussian surface
#'
#' @description Function to create Gaussian random field raster surface for use with DGS simulation
#' 
#' @author Bill Peterman
#'
#' @param dim (Default 200). Number rows and columns in raster
#' @param autocorr_range Maximum range (raster units) of spatial autocorrelation. (Default will be 3 percent of raster 'dim')
#' @param mag_var Magnitude of variation over the entire landscape (Default 5)
#' @param nug Magnitude of variation in the scale of autocorr_range, smaller values lead to more homogeneous landscapes (Default 1)
#' @param user_seed Set random seed for the simulation
#' 
#' @export

gausDGS <- function(dim = 200,
                    autocorr_range = NULL,
                    mag_var = 5,
                    nug = 1,
                    user_seed = NULL){
  
  if(is.null(autocorr_range)){
    autocorr_range <- floor(0.03 * dim)
  }
  
  gaus <- NLMR::nlm_gaussianfield(ncol = dim,
                                  nrow = dim,
                                  autocorr_range = autocorr_range,
                                  mag_var = mag_var,
                                  nug = nug,
                                  user_seed = user_seed)
  return(gaus)
}
