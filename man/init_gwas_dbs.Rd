% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_gwas_db.R
\name{init_gwas_dbs}
\alias{init_gwas_dbs}
\title{Retrieve a current list of all queryable traits in the GWAS catalog}
\usage{
init_gwas_dbs()
}
\value{
A list of all queryable traits in the GWAS catalog
}
\description{
Wrapper function for the get_traits() function in 'gwasrapidd'. Returns a list of traits in the GWAS catalog.
}
\examples{
trait_list <- init_gwas_dbs()
}
