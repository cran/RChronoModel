% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prog.R
\name{PhasesGap}
\alias{PhasesGap}
\title{Gap or Hiatus between two successive phases (for phases in temporal order constraint)}
\usage{
PhasesGap(Phase1Max_chain, Phase2Min_chain, level = 0.95,
  max_decimal = 0)
}
\arguments{
  \item{Phase1Max_chain}{numeric vector containing the output of the MCMC algorithm for the maximum of the dates included in the oldest phase}
  \item{Phase2Min_chain}{numeric vector containing the output of the MCMC algorithm for the minimum of the dates included in the youngest phase}
  \item{level}{probability corresponding to the level of confidence}
  \item{max_decimal}{maximum number of decimal}
}
\value{
Returns a vector of values containing the level of confidence and the endpoints of the gap between the successive phases
}
\description{
This function finds, if it exists, the gap between two successive phases. This gap or hiatus is the longest interval [IntervalInf ; IntervalSup] that satisfies : P(Phase1Max_chain < IntervalInf < IntervalSup < Phase2Min_chain | M) = level.
}
\author{
  Anne Philippe <Anne.Philippe@univ-nantes.fr> and

  Marie-Anne Vibet <Marie-Anne.Vibet@univ-nantes.fr>
}
\examples{
  data(Phases); attach(Phases)
  PhasesGap(Phase.1.beta, Phase.2.alpha, 0.95)
  PhasesGap(Phase.1.beta, Phase.2.alpha, 0.50)
}
\keyword{ gap between two phases}
\keyword{ temporal order }
\keyword{ succession of phases }