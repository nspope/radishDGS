% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_radish2sims.R
\name{fit_radish2sims}
\alias{fit_radish2sims}
\title{Fit radish models to simulated data}
\usage{
XXX
}
\arguments{
\item{Results_dir}{Full path to the top 'Results' directory where simulation results are saved. If only a single iteration of a simulation is to be analyzed, point to that directory within the 'Results' directory.}

\item{conductance_model}{Default = `loglinear_conductance`. c('loglinear_conductance', 'linear_conductance')}

\item{measurement_model}{Default = `generalized_wishart`. c('generalized_wishart', 'mlpe', 'leastsquares')}

\item{optimizer}{(Default = 'newton')}

\item{gd}{Genetic distance to use. Specify 'dps' or 'pca' if using individual-based simulations; 'dc' or 'fst' for population-based simulations.}

\item{...}{Not used}
}
\value{
List of fitted `radish` models
}
\description{
Function to fit `radish` models to CDPOP simulation data
}
\examples{
## Not Run:
## ** TO BE COMPLETED ** ##  

## End (Not run)

}
\author{
Bill Peterman
}
