% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggu.R
\name{work}
\alias{work}
\title{Set up a working directory (easily to convert between windows and mac)}
\usage{
work(folder = "")
}
\arguments{
\item{folder}{The path to the folder (After the path to OneDrive) that you want to set the working directory}

\item{comp}{Either "mac" or "win" for specifying in which computer are you working atm}
}
\value{
generates a working directory
}
\description{
This function creates a pie-chart on the diversity of MOTUs found grouped by the "rank"
}
\examples{
work("IMR/Data/Trawl data")
work("IMR/Data")
}
