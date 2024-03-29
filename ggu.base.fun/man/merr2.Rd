% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggu.R
\name{merr2}
\alias{merr2}
\title{Joins values from another dataframe}
\usage{
merr2(input.x, input.y, input.yy)
}
\arguments{
\item{input.x}{The vector which the values will be joined to}

\item{input.y}{The bridge vector. The vector that has the same as the first vector and rows correspond to values that will be joined}

\item{input.yy}{The vector which contains the values that will be joined}
}
\value{
generates a vector with only joined values ordered for corresponding to the first vectgor
}
\description{
This function matches 2 vectors of different dataframes and joins the values attached to the second vector
}
\examples{
mm$category <- merr2(mm$samp_cat,category.transformation$samp_cat,category.transformation$category)
}
