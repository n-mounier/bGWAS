#' Analysis of genotypes associated with parental lifespan in the UK Biobank.
#'
#' Subset of the original dataset containing the estimated effect of SNPs on parental age of death
#'
#' @format A data frame with 400000 rows and 4 variables:
#' \describe{
#'   \item{snp}{rsid of the SNP}
#'   \item{a1}{effect allele for the SNP}
#'   \item{a0}{reference allele for the SNP}
#'   \item{beta}{estimated effect size for the SNP}
#'   \item{a1}{standard error of the estimated effect size for the SNP}
#' }
#' @source \url{https://figshare.com/articles/Plling_et_al_2017_UKB_parents_attained_age_GWAS/5439382/1}
"SmallGWAS_Pilling2017"
