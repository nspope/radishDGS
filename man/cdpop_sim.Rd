% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms_sim_function.R
\name{cdpop_sim}
\alias{cdpop_sim}
\title{Run CDPOP simulation from R}
\usage{
cdpop_sim(
  master_seed,
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
)
}
\arguments{
\item{master_seed}{Seed value. This is important to set when wanting keep the same locations of points across simulations, but vary some other dimension(s) of the simulation (e.g., dispersal)}

\item{covariates}{Raster stack of input surface(s)}

\item{conductance_quantile_for_demes}{Threshold to set for the proportion of the landscape that is 'good' habitat where points will be distributed during simulation (Default = 0.4)}

\item{effect_size}{Effect size(s) of surfaces in conductance model}

\item{number_of_demes}{Number of demes/spatial locations to simulate (Default = 150)}

\item{sampled_proportion_of_demes}{Number of demes/individuals to randomly select from those simulated (Default = 0.5)}

\item{deme_size}{Number of individuals in each deme (Default = 10)}

\item{proportion_deme_sampled}{Proportion of individuals within each deme to be randomly sampled (Default = 0.5)}

\item{buffer_size}{Buffer to filter/exclude edge demes from the analysis (Default = 0.15). Peripheral populations will not be selected.}

\item{MeanFecundity}{Mean number of offspring produced by each female. Follows a Poisson process (Default = 4)}

\item{iter}{Can specify the iteration number of a for loop. Will be used in naming output folders. If not specified, the time stamp will be used}

\item{sim_dir}{Main directory where simulation results will be written}

\item{looptime}{Total number of steps ('generations') in CDPOP simulations (Default = 201)}

\item{output_years}{The year(s) that you want CDPOP to write simulation results to file (Default = 200)}

\item{loci}{Number of microsatellite loci to simulate (Default = 20)}

\item{alleles}{Number of alleles at each locus to simulate (Default = 20)}

\item{matemovethresh}{Dispersal threshold value (Default = 0.05). This is the quantile threshold of all resistance distances}

\item{python}{Optional. Only needed if Python 2.7 is not in your system path. If needed, specify the full path to the Python 2.7x program file.}
}
\description{
Function to conduct landscape simulation with NLMR and genetic simulation with CDPOP

This function will compile the inputs, run CDPOP, and return the relevant data objects for future use. All generated files are also saved to the specified working directory.

** NOTE: CDPOP requires Python 2.7 to run. You will need to have this version of Python installed, along with the necessary libraries (see the CDPOP manual). CDPOP itself is part if the `radishDGS` package.
}
\author{
Bill Peterman
}
