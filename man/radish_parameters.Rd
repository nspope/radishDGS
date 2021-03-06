% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radish_parameters.R
\name{radish_parameters}
\alias{radish_parameters}
\title{Get parameter estimates from fitted `radish` models}
\usage{
XXX
}
\arguments{
\item{Results_dir}{Full path to the top 'Results' directory where simulation results are saved}

\item{radish_model}{(Default = 'wishart) Specify which `radish` model you want to get parameter estimates from. Should be one of c('wishart', 'mlpe', 'ls'). Assumes there is only a single fitted model object for each measurement model in the directory.}

\item{save_table}{(Default = TRUE) The parameter estimate table will be saved as a CSV file}

\item{conv}{Extra parameter to be used with other functions when checking radish convergence}

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
