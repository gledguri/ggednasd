% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggu.R
\name{get_tax_tree}
\alias{get_tax_tree}
\title{Get Taxonomic Tree}
\usage{
get_tax_tree(input)
}
\arguments{
\item{input}{Either a single species or an entire vector with species}
}
\value{
A dataframe with 7 columns (Kingdom,Phylum,Class,Order,Family,Genus,Species)
}
\description{
This function searches for the taxonomic tree of a species
}
\examples{
get_tax_tree(tax$species)
}
