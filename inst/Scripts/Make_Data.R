# to run on laptop #

set.seed(333)

library(tidyverse)

Lifespan_Data <- data.table::fread("~/Documents/SGG/Projects/LifeGen2/Data/lifegen_phase2_bothpl_alldr_2017_09_18.tsv")
head(Lifespan_Data)
Lifespan_Data$V17 <- NULL

Lifespan_Data %>%
  tidyr::drop_na() -> Lifespan_Data

# take 400,000 SNPs
SNPs <- sample(x = 1:nrow(Lifespan_Data), 400000)

Lifespan_Data %>%
  slice(SNPs) %>%
  transmute(rsid,
            a1,
            a0,
            beta=beta1,
            se) -> SmallGWAS_Timmers2019

save(SmallGWAS_Timmers2019, file="~/Documents/SGG/Projects/Packaging/bGWAS/data/SmallGWAS_Timmers2019.rda", compress='xz')


# also subset the Z-matrices files
Z_mat = "~/ZMatrices/"

# keep BMI / Years of Schooling / CAD / LDL / T2D
New_Files = c("All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz",
              "cardiogram_gwas_results.txt",
              "EDUyears_2016_sumstat.txt",
              "SBP", "DBP",
              "SSGAC_College_Rietveld2013_publicrelease.txt")

AvailableStudies <- data.table::fread(file.path(Z_mat, "AvailableStudies.tsv"))
AvailableStudies %>%
  filter(File %in% New_Files) -> New_AvailableStudies
New_AvailableStudies %>%
  pull(ID) -> IDs
New_AvailableStudies %>%
  mutate(ID=1:6) -> New_AvailableStudies
  
Full_Zmat <- data.table::fread(file.path(Z_mat, "ZMatrix_Full.csv.gz"),
                               select=c(1:5, IDs+5), showProgress = F)


Full_Zmat %>%
  filter(rs %in% SmallGWAS_Timmers2019$rsid) -> New_Full_Zmat

MR_Zmat <- data.table::fread(file.path(Z_mat, "ZMatrix_MR.csv.gz"),
                             select=c(1:5, IDs+5), showProgress = F)

MR_Zmat %>%
  filter(rs %in% SmallGWAS_Timmers2019$rsid) -> New_MR_Zmat


write.table(New_AvailableStudies, sep="\t", row.names = F, quote=F,
            file="~/Documents/SGG/Projects/Packaging/bGWAS/inst/Data/Z_Matrices/AvailableStudies.tsv")

data.table::fwrite(New_Full_Zmat,
                   file="~/Documents/SGG/Projects/Packaging/bGWAS/inst/Data/Z_Matrices/ZMatrix_Full.csv")
R.utils::gzip('~/Documents/SGG/Projects/Packaging/bGWAS/inst/Data/Z_Matrices/ZMatrix_Full.csv',
     destname='~/Documents/SGG/Projects/Packaging/bGWAS/inst/Data/Z_Matrices/ZMatrix_Full.csv.gz',
     remove=T)


data.table::fwrite(New_MR_Zmat,
                   file="~/Documents/SGG/Projects/Packaging/bGWAS/inst/Data/Z_Matrices/ZMatrix_MR.csv")
R.utils::gzip('~/Documents/SGG/Projects/Packaging/bGWAS/inst/Data/Z_Matrices/ZMatrix_MR.csv',
     destname='~/Documents/SGG/Projects/Packaging/bGWAS/inst/Data/Z_Matrices/ZMatrix_MR.csv.gz',
     remove=T)



## Rdata for tests
# last update, 2020/06/07
library(bGWAS)

data("SmallGWAS_Timmers2019")
MyStudies = select_priorGWASs(include_traits=c("Blood Pressure", "Education"),  
                              include_files=c("cardiogram_gwas_results.txt", 
                                              "All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz"))
# 6 Prior GWASs used
list_priorGWASs(MyStudies) 

A = bGWAS(name="Test_UsingSmallDataFrame",
          GWAS = SmallGWAS_Timmers2019,
          prior_studies=MyStudies,
          stepwise_threshold=0.05)



saveRDS(A, file="~/Documents/SGG/Projects/Packaging/bGWAS/inst/Data/A.RDS")
