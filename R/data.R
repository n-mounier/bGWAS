#' Assocation results between genotypes and parental lifespan (LifeGen Consortium).
#'
#' Subset of the original dataset containing the estimated effect of SNPs on parental survival
#'
#' @format A data frame with 100000 rows and 5 variables:
#' \describe{
#'   \item{rsid}{rsid of the SNP}
#'   \item{a1}{effect allele for the SNP}
#'   \item{a0}{reference allele for the SNP}
#'   \item{beta}{estimated effect size for the SNP}
#'   \item{se}{standard error of the estimated effect size for the SNP}
#' }
#' @source \url{https://datashare.is.ed.ac.uk/handle/10283/3209}
"SmallGWAS_Timmers2019"