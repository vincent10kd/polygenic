#' @title Retrieve protein quantitative trait loci (pQTLs) with summary statistics
#'
#' @description Returns a formatted data.frame of variants for a given gene symbol, including summary statistics
#'
#' @param gene_symbol Expects a character vector of length 1, containing the exact gene symbol.
#' @param risk_frequency If TRUE, also includes the frequency of the risk allele in the output.
#' @param is_cis If TRUE, only includes cis pQTLs in the output.
#'
#' @return A formatted data.frame containing all pQTLs for a given gene symbol, including summary statistics.
#' @examples
#' ACE2_pqtls <- get_pQTLs('ACE2',pQTL_list, is_cis=FALSE)
#' @export

get_pQTLs <- function(gene_symbol,
                      risk_frequency=TRUE,
                      is_cis=TRUE){
  starttime <- Sys.time()
  data("pQTL_list")
  if(is_cis==TRUE){
    pQTL_list <- subset(pQTL_list, CisOrTrans=='cis')
  }
  if(!gene_symbol%in%unique(pQTL_list$Affected.Protein.Gene.name)) stop('Gene symbol not found.')
  pQTL_list <- pQTL_list[pQTL_list$Affected.Protein.Gene.name==gene_symbol,]

  SNPs_risk_alleles <- pQTL_list[c('rsName','Amin','Effect.Amin.adj',
                                   'Log10.pval.gc.cor.adj')]
  colnames(SNPs_risk_alleles) <- c('variant_id','risk_allele','beta_coefficient','pvalue')
  SNPs_risk_alleles$variant_id <- gsub(',.*','',SNPs_risk_alleles$variant_id)
  SNPs_risk_alleles$pvalue <- 10^(-1*as.numeric(SNPs_risk_alleles$pvalue))
  SNPs_risk_alleles$beta_coefficient <- as.numeric(SNPs_risk_alleles$beta_coefficient)
  SNPs_risk_alleles$standard_error <- rep(NA,nrow(SNPs_risk_alleles))
  SNPs_risk_alleles$or_per_copy_number <- rep(NA,nrow(SNPs_risk_alleles))
  SNPs_risk_alleles <- SNPs_risk_alleles[,c(1:3,5:6,4)]
  SNPs_risk_alleles$chrom <- pQTL_list$Chrom
  SNPs_risk_alleles$pos <- pQTL_list$Pos

  if(risk_frequency==TRUE){
    cat('>Obtaining risk allele frequency\n')
    SNPs_risk_alleles$risk_frequency <- as.numeric(pQTL_list[,'MAF.PC'])/100
  }

  endtime <- Sys.time()
  cat('\nRetrieving the information took',difftime(endtime,starttime,units='secs'),'seconds.\n')
  rm(starttime,endtime)
  return(invisible(SNPs_risk_alleles))
}

