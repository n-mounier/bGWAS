####################################################################################
#                         Lt's  make nice Z-Matrices !                             #
####################################################################################

### Make Full ZMatrix from a tidy "AvailableStudies.tsv" file
library(data.table)
library(tidyverse)

## PARAMETERS ##
# imputation quality
r2_pred = 0.8
# minor allele frequency (in Lifegen Paper : 0.5%)
maf_th = 0.005
# remove HLA region?
remove_HLA = TRUE

AvailableStudies = fread("AvailableStudies_Update.tsv")
AvailableStudies %>%
  pull(File) -> Studies 
PathToImputedStudies = "~/data/Imputed_GWASs/"


# use UK10K for alignement
UK10K <- fread("/data/sgg2/aaron/homedir/Experiments/genome-data/positions.in.TWINSUK-ALSPAC/positions.b19.TWINSUKALSPAC")
#  chr pos.b19       rs            alts   ref          aaf
#   1   52185   rs201374420           T  TTAA    0.0002644803
#   1   55249   rs200769871      CTATGG     C    0.0068764877
# ...


# start with rsid ...
# for all UK10K SNPs
UK10K %>%
  # all imputed files are swapped, so use alt=ref and ref=alts
  transmute(rs, chrm=chr, pos=pos.b19, alt=ref, ref=alts, MAF=aaf) -> Full_Z
nrow(Full_Z)
# 17,769,027 SNPs

# keep only single nucleotide polymorphism
Full_Z %>%
  filter(alt %in% c("A", "T", "C", "G") & ref %in% c("A", "T", "C", "G")) -> Full_Z
nrow(Full_Z)
# 15,609,168 SNPs

if(remove_HLA){
  Full_Z %>%
    slice(-which(chrm==6 & pos>=28.5e6 & pos<=33.5e6)) -> Full_Z
  nrow(Full_Z)
  # 15,561,505 SNPs
}

Full_Z %>%
  mutate(MAF = pmin(MAF, 1-MAF)) %>%
  filter(MAF>maf_th) %>%
  mutate(MAF = NULL) -> Full_Z



########## Create Z-Matrix ##########

# then add study at a time
# -> at the end, some SNPs might not have well-imputed data from any GWAS

for(s in Studies){
  print(paste0("Adding GWAS: ", s))

  St = fread(file.path(PathToImputedStudies, paste0(s, ".smry..fdr.uk10K")))
  
  
  ## filter on allele freq and imp quality
  print("Filtering on imputation quality...")
  nrow(St)
  # 15,471,679 SNPs
  St %>%
    filter(r2.pred>r2_pred) -> St
  nrow(St)
  # ... SNPs

  print("Adding Z-scores to the Z-matrix...")
  # St %>%
  #  slice(match(Full_Z$rs, St$SNP)) -> St # This only contains SNPs in both dataset
  # We want all SNPs from Full_Z, with NA if not well-imputed in St
  St = St[match(Full_Z$rs, St$SNP),]
  
  Full_Z %>%
    mutate({{s}} := St$z_imp) -> Full_Z
  
  
  print("Done! \n \n")
  
}

### Keep SNPs with at least 20 studies non-NA

Full_Z %>%
  mutate(count_NA = rowSums(is.na(.))) -> Full_Z

table(Full_Z$count_NA)


Full_Z %>%
  filter(count_NA<=(nrow(AvailableStudies)-20)) %>%
  mutate(count_NA=NULL) -> Full_Z

nrow(Full_Z)
# 6,811,310 SNPs

### Derive Strong Instruments ZMatrix
Full_Z %>%
  drop_na() -> MR_Z 
nrow(MR_Z)
# 5,402,848 SNPs

Zlimit = qnorm(1e-5/2, lower.tail = F)

MR_Z %>%
  filter_at(vars(-c(1:5)), 
            any_vars(abs(.) > Zlimit)) -> MR_Z
nrow(MR_Z)
# 209,840 strong instruments

AvailableStudies %>%
  mutate(N_SNPs=unlist(Full_Z %>% summarize_at(vars(-c(1:5)),
                                        function(x) sum((!is.na(x))))),
         N_Instruments=unlist(MR_Z %>% summarize_at(vars(-c(1:5)),
                                             function(x) sum((abs(x) > Zlimit))))) -> AvailableStudies
  


dir.create("ZMatrices")

write.table(AvailableStudies, sep="\t", row.names = F, quote=F,
            file="ZMatrices/AvailableStudies.tsv")

data.table::fwrite(Full_Z,
                   file="ZMatrices/ZMatrix_Full.csv")
R.utils::gzip('ZMatrices/ZMatrix_Full.csv',
              destname='ZMatrices/ZMatrix_Full.csv.gz',
              remove=T)


data.table::fwrite(MR_Z,
                   file="ZMatrices/ZMatrix_MR.csv")
R.utils::gzip('ZMatrices/ZMatrix_MR.csv',
              destname='ZMatrices/ZMatrix_MR.csv.gz',
              remove=T)


system("tar -cvzf ZMatrices.tar.gz ZMatrices")




##### Fix swapped alleles ####
# 
# MR_Z = data.table::fread("ZMatrices/ZMatrix_MR.csv.gz")
# myalt = pull(MR_Z, ref)
# myref = pull(MR_Z, alt)
# 
# MR_Z %>%
#   mutate(alt=myalt,
#          ref=myref) -> MR_Z
# 
# 
# Full_Z = data.table::fread("ZMatrices/ZMatrix_Full.csv.gz")
# myalt_full = pull(Full_Z, ref)
# myref_full = pull(Full_Z, alt)
# 
# Full_Z %>%
#   mutate(alt=myalt_full,
#          ref=myref_full) -> Full_Z
# 
# 
# 
# data.table::fwrite(Full_Z,
#                    file="ZMatrices/ZMatrix_Full.csv")
# R.utils::gzip('ZMatrices/ZMatrix_Full.csv',
#               destname='ZMatrices/ZMatrix_Full.csv.gz',
#               remove=T)
# 
# 
# data.table::fwrite(MR_Z,
#                    file="ZMatrices/ZMatrix_MR.csv")
# R.utils::gzip('ZMatrices/ZMatrix_MR.csv',
#               destname='ZMatrices/ZMatrix_MR.csv.gz',
#               remove=T)
# 
# 
# system("tar -cvzf ZMatrices.tar.gz ZMatrices")
