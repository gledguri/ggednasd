% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggu.R
\name{split_df}
\alias{split_df}
\title{Split vector to dataframe}
\usage{
split_df(input_vector, sep = "_")
}
\arguments{
\item{input_vector}{A vector with elemets to be split}
}
\value{
A dataframe with all elements splitted from the vector, each in a new column
}
\description{
This function splits a vector based on a character and the output is a dataframe
}
\examples{
split_df(blast$annotated_name,sep="\\\\ \\\\ \\\\|\\\\|\\\\ \\\\ ")
}
