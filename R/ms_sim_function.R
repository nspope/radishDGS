#' @author Bill Peterman
#' @title Run CDPOP simulation from R
#' @description Function to conduct landscape simulation with NLMR and genetic simulation with CDPOP
#' 
#' @param master_seed Seed value. This is important to set when wanting keep the same locations of points across simulations, but vary some other dimension(s) of the simulation (e.g., dispersal)
#' @param covariates Raster stack of input surface(s)
#' @param conductance_quantile_for_demes Threshold to set for the proportion of the landscape that is 'good' habitat where points will be distributed during simulation (Default = 0.4)
#' @param effect_size Effect size(s) of surfaces in conductance model
#' @param number_of_demes Number of demes/spatial locations to simulate (Default = 150)
#' @param sampled_proportion_of_demes Number of demes/individuals to randomly select from those simulated (Default = 0.5)
#' @param deme_size Number of individuals in each deme (Default = 10)
#' @param proportion_deme_sampled Proportion of individuals within each deme to be randomly sampled (Default = 0.5)
#' @param buffer_size Buffer to filter/exclude edge demes from the analysis (Default = 0.15). Peripheral populations will not be selected.
#' @param MeanFecundity Mean number of offspring produced by each female. Follows a Poisson process (Default = 4)
#' @param iter Can specify the iteration number of a for loop. Will be used in naming output folders. If not specified, the time stamp will be used
#' @param sim_dir Main directory where simulation results will be written
#' @param looptime Total number of steps ('generations') in CDPOP simulations (Default = 201)
#' @param output_years The year(s) that you want CDPOP to write simulation results to file (Default = 200)
#' @param loci Number of microsatellite loci to simulate (Default = 20)
#' @param alleles Number of alleles at each locus to simulate (Default = 20)
#' @param matemovethresh Dispersal threshold value (Default = 0.05). This is the quantile threshold of all resistance distances
#' @param python Optional. Only needed if Python 2.7 is not in your system path. If needed, specify the full path to the Python 2.7x program file.
#' 
#' @description  This function will compile the inputs, run CDPOP, and return the relevant data objects for future use. All generated files are also saved to the specified working directory.
#' 
#' ** NOTE: CDPOP requires Python 2.7 to run. You will need to have this version of Python installed, along with the necessary libraries (see the CDPOP manual). CDPOP itself is part if the `radishDGS` package.
#' 
#' @export
cdpop_sim <- function(master_seed,
                      covariates,
                      conductance_quantile_for_demes = 0.4,
                      effect_size,
                      number_of_demes = 150,
                      sampled_proportion_of_demes = 0.5,
                      deme_size = 10,
                      proportion_deme_sampled = 0.5,
                      buffer_size = 0.15,
                      MeanFecundity = 4,
                      iter = NULL,
                      sim_dir,
                      looptime = 201,
                      output_years = 200,
                      loci = 20,
                      alleles = 20,
                      matemovethresh = 0.05,
                      python = NULL
){
  
  dir <- try(dir.exists(sim_dir), silent = T)
  if(class(dir) == 'try-error'){
    stop('Specify a directory to save simulation results!!!')
  }
  stopifnot(class(covariates) != 'RasterStack' | class(covariates) != 'RasterLayer')
  if(class(covariates) != 'RasterStack')
    covariates <- stack(covariates)
  stopifnot(length(effect_size) == nlayers(covariates))
  stopifnot(buffer_size >= 0. & buffer_size < 0.5)
  stopifnot(number_of_demes > 0)
  stopifnot(sampled_proportion_of_demes > 0. & sampled_proportion_of_demes <= 1.)
  stopifnot(conductance_quantile_for_demes > 0.)
  
  
  set.seed(master_seed)
  
  CDPOP.py <- list.files(system.file('cdpop', package = 'radishDGS'), full.names = TRUE)[1]
  
  # Defaults ----------------------------------------------------------------
  gridformat <- 'cdpop'  
  percent_quant <- 'quant'
  n_axes <- 32
  matemoveno <- 2
  sim_name <- 'output_'
  
  # >> Create Directory -----------------------------------------------------
  if(is.null(iter)){
    iter <- as.numeric(Sys.time())
  }
  suppressWarnings(
    dir.create(paste0(sim_dir, "/Results/",
                      'iter__', iter),
               recursive = TRUE) 
  )
  
  out <- paste0(sim_dir, "/Results/",
                'iter__', iter, "/")
  
  # Create buffer -----------------------------------------------------------
  
  # create buffer surrounding study area
  buffer           <- covariates[[1]]
  buffer[]         <- 0
  buffer_size_row  <- floor(nrow(buffer) * buffer_size)
  buffer_size_col  <- floor(ncol(buffer) * buffer_size)
  if(buffer_size_row > 0.)
  {
    buffer[1:buffer_size_row,] <- 1
    buffer[nrow(buffer):(nrow(buffer)-buffer_size_row+1),] <- 1
  }
  if(buffer_size_col > 0.)
  {
    buffer[,1:buffer_size_col] <- 1
    buffer[,ncol(buffer):(ncol(buffer)-buffer_size_col+1)] <- 1
  }
  
  
  # Calculate conductance ---------------------------------------------------
  # log-linear conductance surface
  # conductance <- covariates[[1]]
  covariates_ <- covariates
  for(i in 1:nlayers(covariates)){
    covariates_[[i]][] <- scale(covariates_[[i]][])
  }
  names(covariates_) <- paste0(names(covariates_), "_s")
  
  surf_ind <- rep(1, length(effect_size))
  conductance <-  raster::stackApply(covariates_, indices = surf_ind, function(x, ...) exp(x %*% effect_size))
  
  
  # Deme locations ----------------------------------------------------------
  
  # simulate deme locations by sampling grid cells uniformly without replacement, where some condition is met
  threshold               <- quantile(conductance, 1 - conductance_quantile_for_demes, na.rm = TRUE)
  possible_cells          <- which(!is.na(conductance[]) & conductance[] >= threshold)
  all_deme_cells          <- sample(possible_cells, number_of_demes)
  all_deme_coords         <- xyFromCell(conductance, all_deme_cells, spatial=TRUE)
  
  all_ind <- rep(all_deme_cells, each = deme_size)
  all_ind_pts <- xyFromCell(conductance, all_ind, spatial = TRUE)
  # plot(conductance>threshold); plot(all_deme_coords, add = T)
  
  
  # Calculate resistance ----------------------------------------------------
  # get "true" resistance distances among demes
  # Effect sizes here are for CONDUCTANCE (the inverse of resistance)
  fmla <- as.formula(paste0("~ ", paste(names(covariates_), collapse = " + ")))
  resistance_model <- radish::conductance_surface(covariates_, 
                                                  all_deme_coords,
                                                  directions = 8)
  resistance_distance <- radish::radish_distance(matrix(effect_size, 
                                                        1, 
                                                        length(effect_size)),
                                                 fmla, 
                                                 resistance_model, 
                                                 radish::loglinear_conductance)$distance[,,1]
  geographic_distance <- radish::radish_distance(matrix(0, 1, length(effect_size)),
                                                 fmla, resistance_model, 
                                                 radish::loglinear_conductance)$distance[,,1]
  
  # standardize resistance distance to [0,1] and calculate migration matrix
  resistance_distance <- scale_to_0_1(resistance_distance)
  
  # migration_matrix <- dispersal_kernel(resistance_distance)
  # diag(migration_matrix) <- 0 #diagonal must be 0 for msprime
  
  # diploid population sizes per deme: sample single individual per deme
  pop_size <- rep(deme_size, number_of_demes)
  inds_ <- deme_size * number_of_demes
  
  ## Expand resistance matrix to individuals
  rs <- resistance_distance[rep(1:number_of_demes, each = deme_size),]
  rs <- rs[,rep(1:number_of_demes, each = deme_size)]
  
  gd_ <- geographic_distance[rep(1:number_of_demes, each = deme_size),]
  gd_ <- gd_[,rep(1:number_of_demes, each = deme_size)]
  
  # index <- seq(from = 1, 
  #              to = inds_, 
  #              by = deme_size-1)
  #  for(i in 1:(nrow(resistance_distance)-1)){
  #    deme_ind <- rep(resistance_distance[i,],  times = deme_size)
  #   res_ind[index[i]:index[i+1], ] <- deme_ind
  #   res_ind[,index[i]:index[i+1]] <- deme_ind
  # }
  
  if(percent_quant == 'quant'){
    m_thresh <- quantile(lower(resistance_distance), matemovethresh)
  } else {
    m_thresh <- max(true_res) * matemovethresh
  }
  # hist(true_res)
  # abline(v = m_thresh)
  
  
  # Run CDPOP ---------------------------------------------------------------
  cdpop_out <- cdpop(CDPOP.py = CDPOP.py,
                     sim_name = sim_name,
                     pts = all_ind_pts,
                     resist_rast = conductance, ## Conductance
                     resist_mat = rs,
                     sim_dir = out,
                     looptime = looptime,
                     output_years = output_years,
                     gridformat = gridformat,
                     loci = loci,
                     alleles = alleles,
                     matemoveno = matemoveno, ## 1 = Linear, 5 = Neg exp; 9 = custom prob matrix
                     matemovethresh = m_thresh,
                     MeanFecundity = MeanFecundity,
                     deme_size = deme_size,
                     n_demes = number_of_demes)
  
  
  cdpop_grid <- cdpop_out$grid_list[[length(cdpop_out$pop_list)]]
  pops <- cdpop_out$pop_list[[length(cdpop_out$pop_list)]]
  
  # Subsample ---------------------------------------------------------------
  # choose some number of demes inside study area (e.g. not in buffer) to sample
  number_of_sampled_demes <- floor(sampled_proportion_of_demes * number_of_demes)
  demes_not_in_buffer     <- which(!extract(buffer, all_deme_coords))
  
  if(deme_size == 1){
    n_ind_sampled <- 1
    
    
    # Sample & Calculate genetic distance(s) ---------------------------------
    
    sim_demes <- demes_not_in_buffer[demes_not_in_buffer %in% pops$ind]
    
    # >> Individual-based -----------------------------------------------------
    
    if(length(sim_demes) < number_of_sampled_demes){
      number_of_sampled_demes <- length(sim_demes)
    }
    sampled_demes <- sort(sample(sim_demes,
                                 number_of_sampled_demes,
                                 replace = F))
    
    gi_sub <- which(unique(pops$ind) %in% sampled_demes)
    gi_final <- cdpop_grid[gi_sub] 
    
    
    # ns <- ifelse(number_of_sampled_demes > length(pops$ind), 
    #              length(pops$ind), 
    #              number_of_sampled_demes)
    # ind_samp <- sort(sample(1:nInd(cdpop_grid), ns, replace = F))
    pops_ <- pops$ind[gi_sub]
    
    s_pts <- all_deme_coords[pops_]
    
    s_pops <- data.frame(pop = pops_,
                         x = all_deme_coords@coords[pops_,1],
                         y = all_deme_coords@coords[pops_,2])
    
    # PCA dist ---------------------------------------------------------------
    if(n_axes > sum(adegenet::nAll(cdpop_grid[gi_sub]))){
      n_axes <- sum(adegenet::nAll(cdpop_grid[gi_sub]))
    }
    
    Dps <- 1-adegenet::propShared(cdpop_grid[gi_sub])
    pca <- pca_dist(cdpop_grid[gi_sub], n_axes = n_axes)
    
    ## Quick plot select
    # par(mfrow = c(2,2))
    # plot(lower(Dps) ~ lower(true_res.mat[pops_,pops_]), main = "Dps ~ Resist")
    # plot(lower(pca) ~ lower(true_res.mat[pops_,pops_]), main = "PCA ~ Resist")
    # plot(lower(Dps) ~ c(dist(all_deme_coords@coords[pops_,])), main = "Dps ~ Euc dist")
    # plot(lower(pca) ~ c(dist(all_deme_coords@coords[pops_,])), main = "PCA ~ Euc dist")
    # par(mfrow = c(1,1))
    # 
    # graphics.off()
    
    
    
    # Write files -------------------------------------------------------------
    
    saveRDS(gi_final,
            paste0(out, "genind_obj.rds"))
    
    write.table(s_pops,
                paste0(out, "sampled_pops.csv"), 
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)
    
    write.csv(pca, 
              file = paste0(out, "pca_dist.csv"))
    
    write.csv(Dps, 
              file = paste0(out, "Dps_dist.csv"))
    
    writeRaster(conductance,
                paste0(out, "true_conduct.tif"),
                overwrite = T)
    
    write.csv(resistance_distance[sampled_demes,sampled_demes], 
              file = paste0(out, "true_ResistDist.csv"))
    
    write.csv(geographic_distance[sampled_demes,sampled_demes], 
              file = paste0(out, "geographic_ResistDist.csv"))
    
    writeRaster(covariates,
                filename = paste0(out,"/"),
                format = 'GTiff',
                bylayer = TRUE,
                suffix = names(covariates),
                overwrite = T)
    
    writeRaster(covariates_,
                filename = paste0(out,"/"),
                format = 'GTiff',
                bylayer = TRUE,
                suffix = names(covariates_),
                overwrite = T)
    
    results <- list(sim_genind = gi_final,
                    pts = s_pts,
                    pca = pca,
                    dps = Dps,
                    trueRes = resistance_distance[sampled_demes,sampled_demes],
                    geoD = geographic_distance[sampled_demes,sampled_demes],
                    conductance_surface = conductance,
                    covariates = covariates_,
                    covariates_Noscale = covariates,
                    out_dir = paste0(out,"/"),
                    effect_size = effect_size)
    saveRDS(results,
            paste0(out, "AllResults_list.rds")
    )
    
  } else {
    
    # >> Deme-based -----------------------------------------------------------
    
    n_ind_sampled <- floor(proportion_deme_sampled * deme_size)
    
    ## Demes not in buffer
    sim_demes <- demes_not_in_buffer[demes_not_in_buffer %in% pops$pop]
    
    if(length(sim_demes) < number_of_sampled_demes){
      number_of_sampled_demes <- length(sim_demes)
    }
    sampled_demes <- sort(sample(sim_demes, number_of_sampled_demes))
    
    # ** Sample genind -----------------------------------------------------------
    gi_sub <- which(unique(pops$pop) %in% sampled_demes)
    gi_sep <- adegenet::seppop(cdpop_grid)[gi_sub] 
    gi_sep_ <- lapply(gi_sep, function(x) x[sample(1:nrow(x$tab), n_ind_sampled)])
    
    gi_final <- adegenet::repool(gi_sep_)
    
    s_pts <- all_deme_coords[sampled_demes]
    
    s_pops <- data.frame(pop = sampled_demes,
                         x = all_deme_coords@coords[sampled_demes,1],
                         y = all_deme_coords@coords[sampled_demes,2])
    
    
    # ** Genetic distance -----------------------------------------------------
    
    Fst <- PopGenReport::pairwise.fstb(gi_final)
    
    ## Convert `genind` object to `genpop` object
    pop.gp <- adegenet::genind2genpop(gi_final)
    
    ## Calculate Dc from `genpop` object
    Dc <- as.matrix(adegenet::dist.genpop(pop.gp, 
                                          method = 2, 
                                          diag = T, 
                                          upper = T))  
    
    ##
    
    saveRDS(gi_final,
            paste0(out, "genind_obj.rds"))
    
    write.table(s_pops,
                paste0(out, "sampled_pops.csv"), 
                sep = ",",
                row.names = FALSE,
                col.names = TRUE)
    
    write.csv(Dc, 
              file = paste0(out, "Dc_dist.csv"))
    
    write.csv(Fst, 
              file = paste0(out, "Fst_dist.csv"))
    
    writeRaster(conductance,
                paste0(out, "true_conduct.tif"),
                overwrite = T)
    
    write.csv(resistance_distance[sampled_demes,sampled_demes], 
              file = paste0(out, "true_ResistDist.csv"))
    
    write.csv(geographic_distance[sampled_demes,sampled_demes], 
              file = paste0(out, "geographic_ResistDist.csv"))
    
    defaultW <- getOption("warn") 
    options(warn = -1) 
    
    suppressWarnings(
      writeRaster(covariates,
                  filename = paste0(out,"/", names(covariates), format = '.tif'),
                  # format = 'GTiff',
                  bylayer = TRUE,
                  suffix = '',
                  overwrite = T)
    )
    
    suppressWarnings(
      writeRaster(covariates_,
                  filename = paste0(out,"/", names(covariates_), format = '.tif'),
                  # format = 'GTiff',
                  bylayer = TRUE,
                  suffix = '',
                  overwrite = T)
    )
    
    options(warn = defaultW)
    
    results <- list(sim_genind = gi_final,
                    pts = s_pts,
                    dc = Dc,
                    fst = Fst,
                    trueRes = resistance_distance[sampled_demes,sampled_demes],
                    geoD = geographic_distance[sampled_demes,sampled_demes],
                    conductance_surface = conductance,
                    covariates = covariates_,
                    covariates_Noscale = covariates,
                    out_dir = paste0(out,"/"),
                    effect_size = effect_size)
    
    saveRDS(results,
            paste0(out, "AllResults_list.rds")
    )
    
    
  } ## Close Ind vs Deme
  
  ## Save plot
  png(file = paste0(out, "dist_resist_plots.png"),
      width = 10, height = 8,
      units = 'in',
      res = 300)
  plot_cdpop(results)
  dev.off()
  
  return(results)
  
} # end function
