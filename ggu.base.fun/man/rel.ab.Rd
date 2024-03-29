% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggu.R
\name{rel.ab}
\alias{rel.ab}
\title{Converts a species abundance table into relative abundance (*100\%)}
\usage{
rel.ab(input, threshold = 0)
}
\arguments{
\item{input}{Species abundances R object (species are rows x stations are columns)}

\item{threshold}{threshold=0, if you want to change the threshold where the \% will be considered 0}
}
\value{
Creates a R object with the same column & row names as the input
}
\description{
This function converts the abundances into relative abundance per station (column)
}
\examples{
output<-rel.ab(sp.abundance.data) #columns are stations and rows are species
output<-rel.ab(sp.abundance.data, threshold = 0.5)
}
