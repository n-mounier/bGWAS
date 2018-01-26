### Make the Z_Matrix files + SNPInfo file (really needed in the package?)

library(data.table)
library(vcfR)

# Need one file with non-imputed Z scores for identify significant studies
# + one for all SNPs to create prior

# Should I do my imputation using HRC actually ?

SNPInfo_File = FALSE

# SNPInfo_File = TRUE -> creation of the file


##### SNPInfo File #####

if(SNPInfo_File){

# Use HRC SNPs,
# that are present in our reference panel (and therefore imputed in our prior studies)
# With MAF > 0.05% ?
#MAF > 0.1%
#MAF>0.001


# Get HRC SNPs

A <- read.vcfR("/data/sgg2/ninon/data/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz",
               limit = 1e+07, cols=c(1),
               convertNA = TRUE, verbose = TRUE)
Info <- data.table(getFIX(A))
nrow(Info) # 40 405 505
# Get rid of chrm X SNPs
HRC <- Info[Info$CHROM %in% c(1:22),]
nrow(HRC) # 39 131 578
# Get rid of non-rs SNPs
HRC <- HRC[!is.na(HRC$ID),]
nrow(HRC) # 32 873 829

# Get UK10K SNPs
UK10K <- fread("/data/sgg2/aaron/homedir/Experiments/genome-data/positions.in.TWINSUK-ALSPAC/positions.b19.TWINSUKALSPAC")


table(HRC$ID %in% UK10K$rs) # 14 900 131 in common
# 17 973 698 HRC SNPs not in UK10K
table(UK10K$rs %in% HRC$ID)
# 290 0122 UK10K SNPs not in HRC

# Common SNPS
Common = HRC[HRC$ID %in% UK10K$rs, ]
Common[, c("uk10k_pos", "uk10k_alt", "uk10k_ref", "uk10k_maf", "uk10k_rs") :=
          UK10K[match(Common$ID, UK10K$rs), c("pos.b19", "alts", "ref", "aaf", "rs")]]

# Check alleles
aligned = which(Common$ALT==Common$uk10k_alt & Common$REF==Common$uk10k_ref)
length(aligned) # 14 593 655
swapped = which(Common$REF==Common$uk10k_alt & Common$ALT==Common$uk10k_ref)
length(swapped) # 118
# This means we can use 14 593 655 + 118 alleles = 14 593 773

# We would expect than the ones which are swapped are the ones with maf close to 0.5
Common$uk10k_maf[swapped]
summary(Common$uk10k_maf[swapped])
# not the case... anyway, since they are both in HRC and in UK10K, we can keep them and use UK10K for alleles

unconsistent = which(!c(1:nrow(Common)) %in% c(aligned, swapped))
length(unconsistent) # 306 358
Unconsistent = Common[unconsistent, ]
head(unconsistent)


## Check everything with an ssimp-imputed file (UK10K)
GWAS <- fread("/data/sgg2/ninon/shared/Eleonora/GWAS_Height/MetaSum.height.uk10K")
nrow(GWAS) # 15 472 006
table(HRC$ID %in% GWAS$SNP)
GWAS_NonImputed <- fread("/data/sgg2/ninon/shared/Eleonora/GWAS_Height/MetaSum.height.smry..fdr")


# SNPs in GWAS but not in UK10K
noUK10K = which(!GWAS$SNP %in% UK10K$rs)
length(noUK10K) # 44 579
# check that they all come from initial GWAS, and that they are not SNPs imputed using UK10K but
# absent from Aaron's file
# -> not this !! These are SNPs "rsXXX;rsYYY", these ones should be remove when cleaning ssimp output and this should be ok
SNPs_noUK10K = GWAS$SNP[noUK10K]
A = grep(";", SNPs_noUK10K)
table(c(1:length(SNPs_noUK10K) %in% A)) # Ok, this is the problem


# check that all the alleles from UK10K are actually imputed
notimputed = which(!UK10K$rs %in% GWAS$SNP)
length(notimputed) # 2 341 600, from UK10K Aaron's file, not imputed
notimputedC = which(!Common$uk10k_rs %in% GWAS$SNP)
length(notimputedC) # 293 866, from common HRC/UK10K, not imputed
notimputedCOK = which(!Common$uk10k_rs[c(swapped, aligned)] %in% GWAS$SNP)
length(notimputedCOK) # 58 464, from common HRC/UK10K + aligned/swapped, not imputed
head(Common[notimputedCOK,])
Common[notimputedCOK, POS] # does not seem to be a position problem

# ... for now, I will use  only the common aligned/swapped SNPs present in the imputed file
imputed = Common[which(Common$uk10k_rs[c(swapped, aligned)] %in% GWAS$SNP), c("uk10k_rs", "uk10k_maf")]
nrow(imputed[imputed$uk10k_maf>0.01,])
MySNPs <- UK10K[UK10K$rs %in% imputed$uk10k_rs,]

# Write SNPInfo
# MySNPs, should be : chrm / pos.b19 / rs / alts / ref / aaf
write.table(MySNPs, "/data/sgg2/ninon/projects/bGWAS_Packaging/SNPsInfo", quote=F, row.names=F)

} else {
  MySNPs = fread( "/data/sgg2/ninon/projects/bGWAS_Packaging/SNPsInfo")
}

# Then Create Z_Matrices
# (re-use Aaron's code)

threshold = 1e-5


## non-imputed
# Identify InstrumentsSNPs : significant SNPs in the studies
# For each study, identify the SNPs with a p-value<threshold -> InstrumentsSNPs
# The InstrumentsSNPs are identified from studies before imputation


# temporary file
#SNPs = fread("inst/Data/1e-05_Lifespan_7M.csv")
Z = fread("/data/sgg2//ninon/projects/AgingX_Tool/HugeMatrix/ZMatrix_1e-05_Lifespan_7M.csv")
nrow(Z)
Z = Z[,-c("Newlifegen2")]
Z <- Z[complete.cases(Z),]
nrow(Z)

#Z = Z[,-c("chrm", "pos", "alt", "ref")]
colnames(Z)

write.table(Z, "inst/Data/ZMatrix_NotImputed.csv", sep=",", row.names=F, quote=F)
# tar -cvzf ZMatrix_NotImputed.csv.tar.gz ZMatrix_NotImputed.csv # makes tar.gz


B = fread('inst/Data/ZMatrix_NotImputed.csv.gz')
# mac os
B = fread("zcat < inst/Data/ZMatrix_NotImputed.csv.gz")
# > B = data.table::fread(paste0("zcat < ", system.file("Data/ZMatrix_NotImputed.csv.gz", package="bGWAS")))


## imputed




## The full imputed matrix for all studies DOES NOT exist
# We have to create it
Z = fread("/data/sgg2//ninon/projects/AgingX_Tool/HugeMatrix/ZMatrix_whichSNPs_7M_Lifespan_7M.csv")
nrow(Z)
Z = Z[,-c("Newlifegen2")]

Studies=fread("SummarizeFiles_Lifespan_7M.csv", header=F)



for(i in 1:nrow(Studies)){
  print(i)
  s=Studies[i]
  if(!s %in% colnames(Z)){
    GWAS = fread(paste0("ImputedStudies/", s, ".csv"))
    # match align alleles
    GWAS = GWAS[match(Z$rs, unlist(GWAS[,"snpid"])),]
    aligned = which(GWAS[,"a1", with=F] == Z$alt &
                      GWAS[,"a2", with=F] == Z$ref)
    swapped = which(GWAS[,"a2", with=F] == Z$alt &
                      GWAS[,"a1", with=F] == Z$ref)
    weird = c(1:nrow(GWAS))[!c(1:nrow(GWAS)) %in% c(aligned, swapped)]
    GWAS[aligned, myZ:= GWAS[aligned, "z", with=F]]
    GWAS[swapped, myZ:= -GWAS[swapped, "z", with=F]]
    GWAS[weird, myZ:= NA]
    # add to Z
    Z[, as.character(s)] = GWAS$myZ
  }

}


Z2 <- Z[complete.cases(Z),]
nrow(Z2)

#Z = Z[,-c("chrm", "pos", "alt", "ref")]
colnames(Z)

write.table(Z2, "inst/Data/ZMatrix_Imputed.csv", sep=",", row.names=F, quote=F)
#system('tar -cvzf ZMatrix_NotImputed.csv.tar.gz ZMatrix_NotImputed.csv') # makes tar.gz
system("gzip ZMatrix_Imputed.csv")

B = fread('inst/Data/ZMatrix_NotImputed.csv.gz')
# mac os
B = fread("zcat < inst/Data/ZMatrix_NotImputed.csv.gz")
# > B = data.table::fread(paste0("zcat < ", system.file("Data/ZMatrix_NotImputed.csv.gz", package="bGWAS")))




