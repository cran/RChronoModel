% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prog.R
\name{MultiPhaseTimeRange}
\alias{MultiPhaseTimeRange}
\title{Phase Time Range for multiple phases}
\usage{
MultiPhaseTimeRange(data, position_minimum, position_maximum = position_minimum+1,
 level = 0.95, max_decimal = 0)
}
\arguments{
  \item{data}{dataframe containing the output of the MCMC algorithm }
  \item{position_minimum}{numeric vector containing the column number corresponding to the minimum of the dates included in each phase}
  \item{position_maximum}{numeric vector containing the column number corresponding to the maximum of the dates included in each phase. By default, position_maximum = position_minimum + 1. 
}
  \item{level}{probability corresponding to the desired level of confidence}
  \item{max_decimal}{maximum number of decimal}
}
\value{
Returns a matrix of values containing the level of confidence and the endpoints of the shortest time range associated with the desired level
}
\description{
Computes the shortest interval that satisfies : P(PhaseMin < IntervalInf < IntervalSup < PhaseMax | M) = level
}
\details{
For each i, MultiPhaseTimeRange computes the time range interval for the phase defined by its minimum position_minimum[i] and its maximum position_maximum[i]. 
The default value of position_maximum corresponds to CSV files exported from ChronoModel software.
}
\author{
  Anne Philippe <Anne.Philippe@univ-nantes.fr> and

  Marie-Anne Vibet <Marie-Anne.Vibet@univ-nantes.fr>
}
\examples{
  # Data extracted from ChronoModel software
  data(Phases)
  
  # List of the name of the phases 
  names(Phases)
  
  # Stipulating position_maximum
  MultiPhaseTimeRange(Phases, c(4,2), c(5,3))
  
  # In this case, equivalent to 
  MultiPhaseTimeRange(Phases, c(4,2))
}
\keyword{ phase time range }
