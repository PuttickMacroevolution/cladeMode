% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyIC.R
\name{phyIC}
\alias{phyIC}
\title{independent contrasts calculation}
\usage{
phyIC(phy, yDat)
}
\arguments{
\item{phy}{dated tree in ape format}

\item{yDat}{named vector of trait values}
}
\value{
list of parameters
}
\description{
calculates independent contrasts and returns the IC of trait values, ancestral state estimates, phylogenetic mean, rate of BM evolution, root variance, the rescaled phylogeny  and the log-likelihood
}
\examples{
set.seed(30)
# simulate tree with 20 species
simTree <- rcoal(20)
simData <- cladeModeSim(n=1, phy=simTree, mode="BM")
phyIC(simTree, simData)
}
