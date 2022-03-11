#' @title Retrieve a current list of all queryable traits in the GWAS catalog
#'
#' @description Wrapper function for the get_traits() function in 'gwasrapidd'. Returns a list of traits in the GWAS catalog.
#'
#'
#' @return A list of all queryable traits in the GWAS catalog
#' @examples
#' trait_list <- init_gwas_db()
#' @export
#' @importFrom gwasrapidd "get_traits"
#'

init_gwas_dbs <- function(){
  require(gwasrapidd)
  gwas_catalog_query_date <- Sys.Date()
  all_traits <- gwasrapidd::get_traits()
  trait_list <- list(all_traits=all_traits, query_date=gwas_catalog_query_date)
  return(trait_list)
}
