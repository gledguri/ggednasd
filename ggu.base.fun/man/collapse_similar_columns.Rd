% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggu.R
\name{collapse_similar_columns}
\alias{collapse_similar_columns}
\title{Collapse similar columns on a dataframe}
\usage{
collapse_similar_columns(db)
}
\arguments{
\item{db}{The dataframe that columns should be collapsed}
}
\value{
A dataframe containing the collapsed dataframe. If the dataframe is not balanced then NA's will be created for those columns that have less observation
}
\description{
This function collapses similar columns or based on an identical column names. It is important to have the column names
similar for those columns that needs to be collapsed
}
\examples{
as.data.frame(extract(stanMod, par = parameter))
colnames(param) <- paste0(factor1,"_",factor2)
collapse_similar_columns(param)
}
