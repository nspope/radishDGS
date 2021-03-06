#' @export
#' 
plot_one <- function(tr)
{
  library(ape)
  tr <- read.tree(text=tr)
  plot(tr, edge.col=ifelse(1:nrow(tr$edge)==which(tr$tip.label[tr$edge[,2]] == "1"),"red","black"))
}

mrca_one <- function(tr, labels)
{
  library(ape)
  library(phangorn)
  tr <- read.tree(text=tr)
  first_deme <- which(tr$tip.label %in% labels)
  phangorn::mrca.phylo(tr, node = first_deme)
}

mutate_tree <- function(tr, maf_filter = 0.)
{
  tr  <- ape::read.tree(text=tr)
  tr  <- ape::reorder.phylo(tr, "postorder")
  des <- phangorn::Descendants(tr, type="tips")
  maf <- unlist(lapply(des, length))/length(tr$tip.label)
  maf <- ifelse(maf > 0.5, 1-maf, maf)
  ok  <- maf > maf_filter
  ind <- match(1:length(des), tr$edge[,2])
  pr  <- ifelse(!ok | is.na(tr$edge.length[ind]), 0., tr$edge.length[ind])
  pr  <- pr/sum(pr)
  mut <- sample(1:length(des), 1, prob=pr)
  whe <- as.numeric(tr$tip.label[des[[mut]]])
  out <- rep(0, length(tr$tip.label))
  out[whe] <- 1
  out
}

msprime_wrapper <- 
  function(master_seed, 
           reps, 
           covariates, 
           effect_size,
           buffer_size = 0.2,
           number_of_demes = 300,
           sampled_proportion_of_demes = 0.1,
           conductance_quantile_for_demes = 0.,
           dispersal_kernel = function(dist) exp(-1./0.3 * dist^2),
           number_of_loci = 500,
           split_time = 200,
           anc_size = 500,
           deme_size = 1,
           trees_only = FALSE,
           verbose = FALSE,
           maf = 0.05)
{
  # there is something really odd going on here. The 1st population always has higher divergence, and often a negative pattern against distance.
  # this holds regardless of permutation to the migration matrix or changing the random seed.
  library(radish)

  reticulate::source_python(system.file("py/island_model.py",package="radishDGS"))

  stopifnot(reps > 0)
  stopifnot(length(effect_size) == dim(covariates)[3])
  stopifnot(buffer_size >= 0. & buffer_size < 0.5)
  stopifnot(number_of_demes > 0)
  stopifnot(sampled_proportion_of_demes > 0. & sampled_proportion_of_demes <= 1.)
  stopifnot(conductance_quantile_for_demes >= 0.)
  stopifnot(maf >= 0. & maf < 0.5)

  set.seed(master_seed)
  random_seeds <- sample.int(1e8, reps)

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

  # log-linear conductance surface
  conductance <- covariates[[1]]
  if (length(covariates) == 1)
    conductance[] <- exp(effect_size * covariates[[1]][])
  else
    conductance[] <- exp(raster::getValues(covariates) %*% effect_size)

  # simulate deme locations by sampling grid cells uniformly without replacement, where some condition is met
  threshold               <- quantile(conductance, conductance_quantile_for_demes, na.rm = TRUE)
  possible_cells          <- which(!is.na(conductance[]) & conductance[] >= threshold)
  all_deme_cells          <- sample(possible_cells, number_of_demes)
  all_deme_coords         <- xyFromCell(conductance, all_deme_cells, spatial=TRUE)

  # choose some number of demes inside study area (e.g. not in buffer) to sample
  number_of_sampled_demes <- floor(sampled_proportion_of_demes * number_of_demes)
  demes_not_in_buffer     <- which(!extract(buffer, all_deme_coords))
  sampled_demes           <- sort(sample(demes_not_in_buffer, number_of_sampled_demes)) #will fail if insufficient number

  # get "true" resistance distances among demes
  # NB: the effect sizes here are for CONDUCTANCE (the inverse of resistance)
  form <- as.formula(paste0("~", paste(names(covariates), collapse="+")))
  resistance_model        <- radish::conductance_surface(covariates, all_deme_coords, directions=8)
  resistance_distance     <- radish::radish_distance(matrix(effect_size, 1, length(effect_size)), form, resistance_model, radish::loglinear_conductance)$distance[,,1]
  geographic_distance     <- radish::radish_distance(matrix(0, 1, length(effect_size)), form, resistance_model, radish::loglinear_conductance)$distance[,,1]

  # standardize resistance distance to [0,1] and calculate migration matrix
  resistance_distance    <- scale_to_0_1(resistance_distance)
  migration_matrix       <- dispersal_kernel(resistance_distance)
  diag(migration_matrix) <- 0 #diagonal must be 0 for msprime

  output <- list()
  for (seed in random_seeds) #NB: different genetic simulations for same spatial configuration
  {
    timer <- Sys.time()
    cat("\n---\nRunning simulation with random seed", seed, "on", as.character(timer), "\n")

    try({ #this ensures that failure of a given simulation will not stop the loop

      set.seed(seed)

      # diploid population sizes per deme: sample single individual per deme
      pop_size <- rep(deme_size, number_of_demes)

      # track a diploid from each "sampled" deme
      # NB: should be less than or equal to 2*pop_size from previous step
      smp_size <- ifelse(1:number_of_demes %in% sampled_demes, 2, 0)

      # use msprime to simulate haploid genotypes
      # NB: a key point is that we only simulate haploids (e.g gametes)
      #     so diploid genotypes have to be formed from haplotypes. We also generally want to filter
      #     by minor allele frequency, as low frequency alleles are often less informative about population
      #     structure (and may violate model assumptions)
      SNPsim <- island_model(num_loci   = number_of_loci, #number of loci (e.g. SNPs)
                             split_time = split_time, #generations in past where panmixia ended
                             anc_size   = anc_size, #ancestral population size, as number of DIPLOIDS
                             migr_mat   = migration_matrix, #migration matrix. element m_ij is proportion of population i that is replaced by individuals from population j, PER GENERATION. These should be small.
                             pop_size   = pop_size, #contemporary (DIPLOID) population sizes, as number of diploids
                             smp_size   = smp_size, #number of sampled HAPLOIDS per population (can be 0, in which case populations are modelled but not sampled)
                             ran_seed   = seed,
                             maf_filter = maf,
                             keep_trees = TRUE, #if true, keep coalescent tree for each locus as newick (could be used to simulate msats)
                             keep_variants = FALSE, #!trees_only,
                             verbose = verbose,
                             use_dtwf   = TRUE) #if true, use Wright-Fisher backward time simulations (discrete generations, better when population sizes are small/migr rates high)

      # I am puzzled by weird pattern in variants, so mutating trees by hand
      #browser()
      SNPsim$genotypes <- Reduce(rbind, lapply(SNPsim$trees, mutate_tree, maf_filter = maf))

      # create diploid genotypes by pooling two haploids from same spatial location
      genotypes          <- SNPsim$genotypes[,seq(1,ncol(SNPsim$genotypes),2)] + SNPsim$genotypes[,seq(2,ncol(SNPsim$genotypes),2)]
      demes_of_genotypes <- rep(1:number_of_demes, smp_size/2) #deme id for each sampled diploid
      was_sampled        <- demes_of_genotypes %in% sampled_demes
      demes_of_samples   <- demes_of_genotypes[was_sampled]

      # filter by minor allele frequency (only applies to SNPs)
      frequency  <- rowMeans(genotypes) / 2
      maf_filter <- frequency >= maf & frequency <= 1.-maf
      if (any(!maf_filter))
        warning("MAF rejection sampling failed, check python source")
      genotypes  <- genotypes[maf_filter,]
      frequency  <- frequency[maf_filter]
      cat("Simulation produced", nrow(genotypes), "SNPs passing MAF filter\n")

      # genetic covariance (for SNPs)
      normalized_genotypes <- (genotypes - 2 * frequency)/sqrt(2 * frequency * (1-frequency)) #"normalized" genotypes 
      genetic_covariance   <- t(normalized_genotypes) %*% normalized_genotypes / (nrow(normalized_genotypes)-1)

      # store results needed to fit models/reproduce simulation (won't store anything if simulation fails)
      timestamp     <- Sys.time()
      output[[paste0("seed", seed)]] <- 
                       list(covariates    = covariates, 
                            genotypes     = genotypes[,was_sampled], 
                            coords        = all_deme_coords[demes_of_samples,], 
                            rdistance     = resistance_distance[demes_of_samples,demes_of_samples],
                            gdistance     = geographic_distance[demes_of_samples,demes_of_samples],
                            migration     = migration_matrix[demes_of_samples,demes_of_samples],
                            covariance    = genetic_covariance[was_sampled,was_sampled],
                            frequency     = frequency,
                            random_seed   = seed,
                            timestamp     = timestamp)
    })
    timer <- Sys.time() - timer
    cat("Finished in", timer, attr(timer, "units"), "\n")
  }

  # track which ones failed (will help us troubleshoot)
  cat("\nFailed for", length(random_seeds)-length(output), "of", length(random_seeds), "simulations\n")
  failures <- random_seeds[!(random_seeds %in% names(output))]
  attr(output, "failed") <- failures

  output
}

msprime_wrapper2 <- 
  function(master_seed, 
           number_of_loci, 
           reps,
           covariates, 
           effect_size,
           buffer_size = 0.2,
           number_of_demes = 300,
           sampled_proportion_of_demes = 0.1,
           conductance_quantile_for_demes = 0.,
           dispersal_kernel = function(dist) exp(-1./0.3 * dist^2),
           deme_size = 1)
{

library(radish)
library(reticulate)

# msprime via reticulate in R
msprime_src <- '

import msprime
import numpy

def island_model (reps, seed, smp_size, deme_size, migr_mat):
    migr_mat   = numpy.array(migr_mat, dtype=float)
    num_popul  = migr_mat.shape[0]
    pop_size   = numpy.array([deme_size]*num_popul, dtype=float)
    smp_size   = numpy.array(smp_size, dtype=int)
    assert(smp_size.shape[0] == num_popul)
    ran_seed   = int(seed)
    num_loci   = int(reps)
    # demes
    pop_confg  = []
    for i in range(num_popul):
      pop_confg += [msprime.PopulationConfiguration(
        sample_size=smp_size[i], initial_size=pop_size[i])]
    # simulate trees
    sim = msprime.simulate(
            population_configurations = pop_confg,
            migration_matrix = migr_mat,
            mutation_rate = 0.,
            Ne = 1.,
            num_replicates = num_loci,
            random_seed = ran_seed,
            model = "dtwf")
    trees = []
    for locus in sim:
      trees += [locus.at(0.).newick()]
    return {"trees" : trees}
'
msprime_tmp <- tempfile()
cat(msprime_src, file=msprime_tmp)
reticulate::source_python(file = msprime_tmp)

# drop neutral mutation onto a Newick-encoded tree, return genotype vector
mutate_tree <- function(newick)
{
  tree <- ape::read.tree(text=newick)
  descendants <- phangorn::Descendants(tree, type="tips")
  branch <- match(1:length(descendants), tree$edge[,2])
  branch_length <- ifelse(is.na(tree$edge.length[branch]), 0., tree$edge.length[branch])
  mutated_ancestor <- sample(1:length(descendants), 1, prob=branch_length)
  mutated_tips <- as.numeric(tree$tip.label[descendants[[mutated_ancestor]]])
  ifelse(1:length(tree$tip.label) %in% mutated_tips, 1, 0)
}

stopifnot(number_of_loci > 0)
stopifnot(length(effect_size) == dim(covariates)[3])
stopifnot(buffer_size >= 0. & buffer_size < 0.5)
stopifnot(number_of_demes > 0)
stopifnot(deme_size >= 1)
stopifnot(sampled_proportion_of_demes > 0. & sampled_proportion_of_demes <= 1.)
stopifnot(conductance_quantile_for_demes >= 0.)

#set.seed(master_seed)
random_seeds <- master_seed + 1000*(1:reps)#sample.int(1e8, reps)

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

# log-linear conductance surface
conductance <- covariates[[1]]
if (length(covariates) == 1)
  conductance[] <- exp(effect_size * covariates[[1]][])
else
  conductance[] <- exp(raster::getValues(covariates) %*% effect_size)

# simulate deme locations by sampling grid cells uniformly without replacement, where some condition is met
# NB: we need 1+number_of_demes, and drop the first, or else there are artefacts-- follow up on github/tskit-dev/msprime issues
threshold       <- quantile(conductance, conductance_quantile_for_demes, na.rm = TRUE)
possible_cells  <- which(!is.na(conductance[]) & conductance[] >= threshold)
number_of_sampled_demes <- floor(sampled_proportion_of_demes * number_of_demes)
number_of_demes <- number_of_demes + 1
all_deme_cells  <- sample(possible_cells, number_of_demes)
all_deme_coords <- xyFromCell(conductance, all_deme_cells, spatial=TRUE)

# choose some number of demes inside study area (e.g. not in buffer) to sample
demes_not_in_buffer     <- which(!extract(buffer, all_deme_coords))
if (any(demes_not_in_buffer==1)) #avoid issue described above
  demes_not_in_buffer   <- demes_not_in_buffer[-which(demes_not_in_buffer==1)]
sampled_demes           <- sort(sample(demes_not_in_buffer, number_of_sampled_demes)) #will fail if insufficient number

# get "true" resistance distances among demes
form                    <- as.formula(paste0("~", paste(names(covariates), collapse="+")))
resistance_model        <- radish::conductance_surface(covariates, all_deme_coords, directions=8)
resistance_distance     <- radish::radish_distance(matrix(effect_size, 1, length(effect_size)), form, resistance_model, radish::loglinear_conductance)$distance[,,1]
geographic_distance     <- radish::radish_distance(matrix(0, 1, length(effect_size)), form, resistance_model, radish::loglinear_conductance)$distance[,,1]

# standardize resistance distance to [0,1] and calculate migration matrix
resistance_distance    <- scale_to_0_1(resistance_distance)
migration_matrix       <- dispersal_kernel(resistance_distance)
diag(migration_matrix) <- 0 #diagonal must be 0 for msprime

output <- list()
for (seed in random_seeds) #NB: different genetic simulations for same spatial configuration
{
  timer <- Sys.time()
  cat("\n---\nRunning simulation with random seed", seed, "on", as.character(timer), "\n")

  try({ #this ensures that failure of a given simulation will not stop the loop

    # track a diploid from each "sampled" deme
    smp_size <- ifelse(1:number_of_demes %in% sampled_demes, 2, 0)

    # use msprime to simulate trees for haploids
    SNPsim <- island_model(number_of_loci, seed, smp_size, deme_size, migration_matrix)

    # simulate neutral mutations on trees
    SNPsim$genotypes <- t(sapply(SNPsim$trees, mutate_tree, USE.NAMES=FALSE))

    # create diploid genotypes by pooling two haploids from same spatial location
    genotypes <- SNPsim$genotypes[,seq(1,ncol(SNPsim$genotypes),2)] + SNPsim$genotypes[,seq(2,ncol(SNPsim$genotypes),2)]

    # genetic covariance (for SNPs)
    frequency <- rowMeans(SNPsim$genotypes)
    normalized_genotypes <- (genotypes - 2 * frequency)/sqrt(2 * frequency * (1-frequency)) #"normalized" genotypes 
    genetic_covariance   <- t(normalized_genotypes) %*% normalized_genotypes / (nrow(normalized_genotypes)-1)

    # store results needed to fit models/reproduce simulation (won't store anything if simulation fails)
    timestamp <- Sys.time()
    output[[paste0("seed", seed)]] <- 
                     list(covariates    = covariates, 
                          genotypes     = genotypes,
                          coords        = all_deme_coords[sampled_demes,], 
                          rdistance     = resistance_distance[sampled_demes,sampled_demes],
                          gdistance     = geographic_distance[sampled_demes,sampled_demes],
                          migration     = migration_matrix[sampled_demes,sampled_demes],
                          covariance    = genetic_covariance,
                          frequency     = frequency,
                          random_seed   = seed,
                          timestamp     = timestamp)
  })
  timer <- Sys.time() - timer
  cat("Finished in", timer, attr(timer, "units"), "\n")
}

output

}
