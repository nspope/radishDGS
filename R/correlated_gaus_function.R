#' @author Bill Peterman
#' @title Generate correlated Gaussian surfaces
#' @description Function to create correlated Gaussian surfaces
#' 
#' @param corr Desired correlation level between fbm surfaces
#' @param dim Dimension of simulated raster 
#' @param autocorr_range Maximum range (raster units) of spatial autocorrelation. (Default = 5)
#' @param mag_var Magnitude of variation over the entire landscape. (Default = 25)
#' @param nug Magnitude of variation in the scale of autocorr_range, smaller values lead to more homogeneous landscapes.(Default = 5)
#' @param user_seed Seed can be set for reproducibility (Default = NA)
#' 
#' @export

gaus_corr <- function(corr,
                      dim,
                      autocorr_range = 5,
                      mag_var = 25,
                      nug = 5,
                      user_seed = NA){
  
  sim_rast <- NLMR::nlm_gaussianfield(ncol = dim,
                                      nrow = dim, 
                                      autocorr_range = autocorr_range,
                                      mag_var = mag_var,
                                      nug = nug,
                                      user_seed = user_seed)
  
  rast_corr <- 0
  
  while(rast_corr > corr + 0.015 | rast_corr < corr - 0.015){
    rep_sim <- NLMR::nlm_gaussianfield(ncol = dim,
                                       nrow = dim, 
                                       autocorr_range = autocorr_range,
                                       mag_var = mag_var,
                                       nug = nug, 
                                       user_seed = NA)
    
    corr_ <-  sqrt((1 / (corr^2)) - 1) 
    mat <- raster::as.matrix(sim_rast)
    
    rep_sim_ <- rep_sim * corr_
    corr_rast <- sim_rast + rep_sim_
    
    rast_corr <- raster::layerStats(stack(sim_rast, corr_rast), 'pearson')$`pearson correlation coefficient`[1,2]
    cat(paste0(round(rast_corr, digits = 4), '\n \n'))
  }
  
  cat(paste0("Surface correlation is ", round(rast_corr, 4), '\n \n'))
  
  
  corr_rast[] <- (corr_rast[] - min(corr_rast[])) / (max(corr_rast[]) - min(corr_rast[]))
  
  
  return(stack(list(sim_true = sim_rast,
                    corr_sim = corr_rast)))
}