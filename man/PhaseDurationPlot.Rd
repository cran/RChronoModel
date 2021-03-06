\name{PhaseDurationPlot}
\alias{PhaseDurationPlot}
\title{Plot of the marginal posterior densities of the duration of a phase}
\usage{
PhaseDurationPlot(PhaseMin_chain, PhaseMax_chain, level=0.95, 
title = "Duration of the phase", colors = T, GridLength=1024)
}
\arguments{
  \item{PhaseMin_chain}{numeric vector containing the output of the MCMC algorithm for the minimum of the dates included in the phase}
  \item{PhaseMax_chain}{numeric vector containing the output of the MCMC algorithm for the maximum of the dates included in the phase}
  \item{level}{probability corresponding to the level of confidence used for the credible interval and the time range}
  \item{title}{title of the graph}
  \item{colors}{if TRUE  -> use of colors in the graph}
  \item{GridLength}{length of the grid used to estimate the density}
}
\value{
A plot with the marginal posterior densities of the duration of a phase and adds several summary statistics (mean, Credible interval, Time range)
}
\description{
This function draws the marginal posterior densities of the time elapsed between the minimum and the maximum of the dates included in a phase, and adds summary statistics (mean, CI)
}
\author{
  Anne Philippe <Anne.Philippe@univ-nantes.fr> and

  Marie-Anne Vibet <Marie-Anne.Vibet@univ-nantes.fr>
}
\examples{
  data(Phases); attach(Phases)

  PhaseDurationPlot(Phase.1.alpha, Phase.1.beta, 0.95, "Duration of Phase 1")
  PhaseDurationPlot(Phase.2.alpha, Phase.2.beta, 0.95, "Duration of Phase 2",colors = FALSE)
}
\keyword{ bayesian statistics }
\keyword{ highest posterior density}
\keyword{ credible interval }
\keyword{ mean }
\keyword{ marginal posterior density }
\keyword{ individual phase}
