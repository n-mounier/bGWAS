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

## Create the full Z-Matrix, and then exctract the SNP with p<10^-5 in at least one study
# we don't care if imputed / non imputed

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


####################################################################################
#                         Let's make nice Z-Matrices !                             #
####################################################################################

### Make Full ZMatrix from scratch
library(data.table)

## parameters
# imputation quality
r2_pred = 0.7
# minor allele frequency (in Lifegen Paper : 0.5%)
maf_th = 0.005



# use UK10K for alignement
UK10K <- fread("/data/sgg2/aaron/homedir/Experiments/genome-data/positions.in.TWINSUK-ALSPAC/positions.b19.TWINSUKALSPAC")
#  chr pos.b19       rs            alts   ref          aaf
#   1   52185   rs201374420           T  TTAA    0.0002644803
#   1   55249   rs200769871      CTATGG     C    0.0068764877
# ...



Studies = system("ls ~/data/Imputed_GWASs/", intern=T)

# start with rsid ...
# for all UK10K SNPs
Full_Z = UK10K[,c("rs", "chr", "pos.b19", "alts", "ref")]
colnames(Full_Z) = c("rs", "chrm", "pos", "alt", "ref")

# keep only single nucleotide polymorphism
Full_Z = Full_Z[Full_Z$alt %in% c("A", "T", "C", "G"),]
Full_Z = Full_Z[Full_Z$ref %in% c("A", "T", "C", "G"),]
# 15,609,168 SNPs


########## Check imputation quality ##########
Res = data.frame(Study=character(), Total=numeric(),
                 PercR2HigherThan0.9 = numeric(),
                 PercR2HigherThan0.7 = numeric(),
                 PercR2HigherThan0.5 = numeric(),
                 PercR2HigherThan0.3 = numeric(),
                 PercMAFHigherThan0.1 = numeric(),
                 PercMAFHigherThan0.05 = numeric(),
                 PercMAFHigherThan0.01 = numeric(),
                 GoodQuality=numeric(),
                 MAF=numeric(),
                 GoodQualityAndMAF=numeric(),
                 stringsAsFactors = F)




for(s in Studies){
  print(paste0("Adding GWAS: ", s))
  St = fread(paste0("~/data/Imputed_GWASs/", s))

Res[nrow(Res)+1,] = c(as.character(s),
                      nrow(St),
                      round(nrow(St[St$r2.pred>0.9,])/nrow(St), 5),
                      round(nrow(St[St$r2.pred>0.7,])/nrow(St), 5),
                      round(nrow(St[St$r2.pred>0.5,])/nrow(St), 5),
                      round(nrow(St[St$r2.pred>0.3,])/nrow(St), 5),
                      round(nrow(St[St$maf>0.1,])/nrow(St), 5),
                      round(nrow(St[St$maf>0.05,])/nrow(St), 5),
                      round(nrow(St[St$maf>0.01,])/nrow(St), 5),
                      nrow(St[St$r2.pred>0.7,]),
                      nrow(St[St$maf>0.05,]),
                      nrow(St[St$maf>0.05 & St$r2.pred>0.7,]))
}

write.table(Res, "FileInfo.csv", sep=",", quote=F, row.names=F)

########## Check variance ##########

VarAll = c()
Var0.01WellImputed = c()
Var0.005WellImputed = c()
for(s in Studies){
  print(s)


  GWAS_Imp = fread(paste0("~/data/Imputed_GWASs/", s))

  # are SSImp file all aligned : yes.

  VarAll = c(VarAll, var(GWAS_Imp$z_imp))

  MyGWAS = GWAS_Imp[GWAS_Imp$maf > 0.01 & GWAS_Imp$r2.pred>0.7, ]
  Var0.01WellImputed = c(Var0.01WellImputed, var(MyGWAS$z_imp))

  MyGWAS2 = GWAS_Imp[GWAS_Imp$maf > 0.005 & GWAS_Imp$r2.pred>0.7, ]
  Var0.005WellImputed = c(Var0.005WellImputed, var(MyGWAS2$z_imp))

}

VAR = data.frame(GWAS=Studies, Var_All=VarAll, Var_qual0.7_1perc=Var0.01WellImputed, Var_qual0.7_0.5perc=Var0.005WellImputed)
write.csv(VAR, "VarianceNewlyImputed_GWAS.csv")



########## Create Z-Matrix ##########

# then add study at a time
# and remove SNPs without info fot this study
# -> at the end, we'll have only the SNPs with Z-score for all studies!

AllAligned=c()
# Studies_imp = Res$Study[Res$PercR2HigherThan0.7>0.5]
for(s in Studies){
  print(paste0("Adding GWAS: ", s))
  St = fread(paste0("~/data/Imputed_GWASs/", s))


  ## filter on allele freq and imp quality
  print("Filtering on maf and imputation quality...")
  St = St[St$r2.pred>r2_pred,]
  St = St[St$maf>maf_th,]


  ## reduce to UK10K SNPs and align
  print("Checking alignement with UK10 SNPs...")
  St = St[St$SNP %in% Full_Z$rs,]
  St = St[match(Full_Z$rs, St$SNP),]
  # check alignment
  aligned = which(St$a2 == Full_Z$alt &
                    St$a1 == Full_Z$ref)
  swapped = which(St$a1 == Full_Z$alt &
                    St$a2 == Full_Z$ref)
  # studies imputed using UK10K so they should not be any swapped
  St[, myZ:= numeric()]
  St[aligned, myZ:= St$z_imp[aligned]]
  St[swapped, myZ:= -St$z_imp[swapped]]
  if(length(swapped)==0) AllAligned[length(AllAligned)+1]=TRUE
  if(length(swapped)==0) print("All aligned!")
  if(length(swapped)>0) AllAligned[length(AllAligned)+1]=FALSE

  ## add z-scores to the Z-matrix
  # and remove SNPs without info
  print("Adding Z-scores to the Z-matrix...")
  Full_Z[, gsub(".smry..fdr.uk10K", "", s)] = St$myZ
  nSNPs_before = nrow(Full_Z)
  print("Removing SNPs without Z-scores...")
  Full_Z = Full_Z[complete.cases(Full_Z)]
  nSNPs_after = nrow(Full_Z)
  print(paste0(nSNPs_after, " SNPs in the full Z-Matrix..."))


  print("Done! \n \n")

}


# if some "-" in the file name, replace by "_"
colnames(Full_Z) = gsub("-", "_", colnames(Full_Z))


## removing the "low imputation quality" (less that 50 SNPs with r2.imp>0.7) GWAS -> 6,282,521 SNPs
# write.table(Full_Z, "FullZ_48GWASs.csv", sep=",", quote=F, row.names=F)
# v = apply(Full_Z[:6:ncol(Full_Z)], 2, var)

## using all GWASs -> 756,024 SNPs
write.table(Full_Z, "FullZ_58GWASs.csv", sep=",", quote=F, row.names=F)
# v = apply(Full_Z[:6:ncol(Full_Z)], 2, var)




### Derive Strong Instruments ZMatrix
SI_Z = Full_Z
Zlimit = qnorm(1e-5/2, lower.tail = F)
SNPsToKeep = apply(SI_Z[,-c(1:5)], 1, function(x) any(abs(x)>Zlimit))

# 48 studies, 750,000 SNPs -> 72,457 strong instruments
# 58 studies, 6.3M SNPs    -> 205,102 strong instruments
# tmp Non-Imputed ZMatrix  -> 103,594
SI_Z = SI_Z[SNPsToKeep,]




### Add columns to existing ZMatrix
