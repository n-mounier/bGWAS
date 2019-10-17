library(tidyverse)


#### GET ASSOCIATED TRAITS FROM GWAS CATALOG (USING SNPS IN LD) ####

if(!exists("ebicat37")) data(ebicat37, package = "gwascat")
# to use makeCurrent... need to load package before, requires data set 'si.hs.38'
library(gwascat)
if(!exists("my_ebicat37")) my_ebicat37 =  gwascat::makeCurrentGwascat(genome='GRCh37')


get_associatedTraits <- function(rsid, chr, pos, LD=0.05, distance=100000, P=5e-8, gwascatdata=NULL){
  require(gwascat)
  ### get GWAS catalog data
  if(is.null(gwascatdata)){
    #if(!exists("ebicat37")) data(ebicat37)
    #if(!exists("my_ebicat37")) my_ebicat37 =  makeCurrentGwascat(genome='GRCh37')
    gwascatdata = makeCurrentGwascat(genome='GRCh37') 
  }
   as_tibble(gwascatdata) %>%
    transmute(snp = SNPS,
              snp_with_riskAllele = STRONGEST.SNP.RISK.ALLELE,
              chrm = CHR_ID,
              posh19 = start, #CHR_POS -> GRCh38
              LD_R2 = NA_real_,
              LD_alleles = NA_character_,
              trait = DISEASE.TRAIT,
              #trait = MAPPED_TRAIT,
              p = P.VALUE,
              adjusted_p = NA_real_,
              # effect = OR.or.BETA, # useless, no direction ...
              effect = X95..CI..TEXT.,
              gene = MAPPED_GENE,
              url = LINK) -> gwascatdata_nice
  
  ### if SNP itself associated, get info
  if(rsid %in% gwascatdata_nice$snp){
    gwascatdata_nice %>%
      filter(snp==rsid,
             p<=P) %>%
      arrange(trait, p)-> my_data 
    # keep only strongest association if multiple ones reported for same trait
    my_data %>%
      filter(!duplicated(trait)) -> my_data
  }
   
  ### also test SNPs nearby
   gwascatdata_nice %>%
    filter(snp != rsid,
           chrm == chr,
           posh19 <= pos+distance,
           posh19 >= pos-distance,
           p<=P) -> snps_nearby

  # get their LD with our main SNP
  require(LDlinkR)
  get_R2 <- function (snp1, snp2) {
    d = capture.output(try(LDpair(snp1, 
                                  snp2, 
                                  pop="EUR", 
                                  token=Sys.getenv("LDLINK_TOKEN"),
                                  output="text"), silent = TRUE))
    # Tidy up data_out
    r2 = strsplit(d[22], "\\s+")[[1]]
    r2 = as.numeric(r2[3])
    return(r2)
  }
  
  get_allele <- function (snp1, snp2) {
    d = capture.output(try(LDpair(snp1, 
                                  snp2, 
                                  pop="EUR", 
                                  token=Sys.getenv("LDLINK_TOKEN"),
                                  output="text"), silent = TRUE))
    # Tidy up data_out
    alignment = strsplit(d[27], " ")[[1]]
    alignment = paste0(alignment[1], "/", alignment[6])
    return(alignment)
  }
  
  
  snps_nearby$LD_R2 = apply(snps_nearby, 1, function(x) suppressWarnings(get_R2(x[1], rsid)))
 
  # keep only the ones in LD with significant association
  snps_nearby %>%
    filter(LD_R2>LD) %>%
    mutate(adjusted_p = 2*pnorm(sqrt(LD_R2) * qnorm(0.5*p))) %>%
    filter(adjusted_p<=P) %>%
    arrange(trait, adjusted_p) -> snps_nearby
  
  # if any trait(s) in common with SNP itself -> remove
  if(exists("my_data")){
    snps_nearby %>%
      filter(!trait %in% my_data$trait) -> snps_nearby
  }
  # keep only strongest association if multiple ones reported for same trait
  snps_nearby %>%
    filter(!duplicated(trait)) -> snps_nearby
  
  # add info about alleles in LD
  snps_nearby$LD_alleles = apply(snps_nearby, 1, function(x) suppressWarnings(get_allele(x[1], rsid)))
  
  if(exists("my_data")) bind_rows(my_data, snps_nearby) -> snps_nearby

  snps_nearby %>%
    mutate(snp=snp_with_riskAllele,
           snp_with_riskAllele=NULL) -> snps_nearby
    
  return(snps_nearby)
}

# # associated with T2D, multiple studies & blood pressure
# rsid="rs34872471"
# chr=10
# pos=114754071
# get_associatedTraits(rsid, chr, pos, gwascatdata = my_ebicat37) -> A

# # not in GWAS catalog, multiple SNPs nearby associated with several traits (at 500kb)
# rsid="rs7630554"
# chr=3 
# pos=185526062
# get_associatedTraits(rsid, chr, pos,  gwascatdata = my_ebicat37) -> B




#### GET GENES NAMES USING ENSEMBL ####
# not used actually, rather use annovar

my_biomart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                           host="grch37.ensembl.org", path="/biomart/martservice", 
                           dataset="hsapiens_gene_ensembl")


get_gene <- function(chr, pos, biomart=NULL){
  require(biomaRt)
  
  ### get ensembl data
  if(is.null(biomart)){
    biomart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                      host="grch37.ensembl.org", path="/biomart/martservice", 
                      dataset="hsapiens_gene_ensembl")
  }

  # get gene(s) at this exact position
  all.genes <- biomaRt::getBM(
    attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"),
    filters=c("chromosome_name", "start", "end"),
    values=list(chromosome=chr, start=pos, end=pos),
    mart=biomart)
  
  all.genes %>% 
    filter(!hgnc_symbol=="") -> all.genes
  
  if(nrow(all.genes)==1){
    my_gene = all.genes
    
    ## what if several ones?
  } else if(nrow(all.genes)>=1){
    # get the closest?
    all.genes %>% 
      mutate(distance=pmin(abs(start_position-pos), abs(end_position-pos))) %>%
      arrange(distance)-> all.genes
    my_gene = all.genes[1,]  
    ## what if no gene exactly here?
 } else if(nrow(all.genes)==0){
   # look at the ones in a 500kb windows
   all.genes <- biomaRt::getBM(
     attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"),
     filters=c("chromosome_name", "start", "end"),
     values=list(chromosome=chr, start=pos-500000, end=pos+500000),
     mart=biomart) 
   # get the closest?
   all.genes %>% 
     filter(!hgnc_symbol=="") %>% 
     mutate(distance=pmin(abs(start_position-pos), abs(end_position-pos))) %>%
     arrange(distance)-> all.genes
   # get the most relevant one(s)
   all.genes %>% 
     filter(!hgnc_symbol=="") %>% 
     mutate(distance=pmin(abs(start_position-pos), abs(end_position-pos))) %>%
     arrange(distance)-> all.genes
   my_gene = all.genes[1,]
 }
  
  return(my_gene$hgnc_symbol)
  
}

#### GET GENES NAMES USING ANNOVAR ####


# get_geneInfo returns c(function, gene, distance)
get_geneInfo <- function(chr, pos, ref="0", alt="0"){
  
  init.dir = getwd()
  # Set needed directory (annovar needs to be installed somewhere)
  annovar.dir <- "~/bin/annovar"
  if(!dir.exists(annovar.dir)) suppressMessages(BioInstaller::install.bioinfo('annovar', annovar.dir))
  setwd(annovar.dir)
  if(!file.exists("humandb/hg19_refGene.txt")) suppressMessages(system("perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/"))
  
  write.table(data.frame(chr, pos, pos, ref, alt), "mydata", sep="\t", row.names=F, col.names=F, quote=F)
  
  
  
  #The gene-based annotation can be issued by the following command (by default, --geneanno -dbtype refGene is assumed):
  # [kaiwang@biocluster ~/]$ annotate_variation.pl -out ex1 -build hg19 example/ex1.avinput humandb/
  
  suppressMessages(system("perl annotate_variation.pl -out myres -build hg19 mydata humandb/",
                    intern = F, ignore.stdout=T, ignore.stderr=T))
  # creates myres.variant_function, myres.exonic_variant_function and myres.log / remove them afterwards
  
  data.table::fread("myres.variant_function", data.table = F) %>% 
    dplyr::select(c(1:2)) %>%
    setNames(c("Function", "Gene")) -> res
  
  # if intronic/exonic ... only gene name reported, ok
  if(res$Function %in% c("intronic", "exonic", "ncRNA_intronic")){
    res %>%
      mutate(Distance=0) -> res
  # if UTR / splicing, need to clean gene name
  } else if(res$Function %in% c("UTR3", "UTR5", "splicing")){
    my_gene = strsplit(res$Gene, "(", fixed=T)[[1]][1]
    res %>%
      mutate(Gene=my_gene,
             Distance=0) -> res
    # if downstream / upstream
  } else if(res$Function %in% c("downstream", "upstream")){
    my_gene = strsplit(res$Gene, "(", fixed=T)[[1]][1]
    my_distance = gsub('.{1}$', '', strsplit(res$Gene, "dist=", fixed=T)[[1]][2])
    res %>%
      mutate(Gene=my_gene,
        Distance=my_distance) -> res
  }else if(res$Function == "intergenic"){
    # returns both (can be more than 2?)
    my_gene1 = strsplit(res$Gene, "(", fixed=T)[[1]][1]
    my_distance1 = strsplit(gsub('.{1}$', '', strsplit(res$Gene, "dist=", fixed=T)[[1]][2]), ")", fixed=T)[[1]][1]
    part2 = strsplit(res$Gene, ",", fixed=T)[[1]][2]
    my_gene2 = strsplit(part2, "(", fixed=T)[[1]][1]
    my_distance2 = strsplit(gsub('.{1}$', '', strsplit(part2, "dist=", fixed=T)[[1]][2]), ")", fixed=T)[[1]][1]
    
    res %>%
      mutate(Gene=paste(my_gene1, my_gene2, sep="/"),
             Distance=paste(my_distance1, my_distance2, sep="/")) -> res
  } 
  
  system("rm myres.* mydata")
  setwd(init.dir)
  return(res)
}

# # exonic
# chr = 19
# pos = 45411941
# snp= "rs429358"
# ref = "T"
# alt= "C"
# get_geneInfo(chr, pos, alt, ref)
# 
# # intergenic
# chr = 6
# pos = 98322872
# snp= "rs4580876"
# ref = "G"
# alt= "A"
# get_geneInfo(chr, pos, alt, ref)
# 
# # UTR3
# chr = 16
# pos = 4013467
# snp= "rs2531995"
# ref = "C"
# alt= "T"
# get_geneInfo(chr, pos, alt, ref)
# 
# # UTR5
# chr = 1 
# pos = 948921
# snp = "rs15842"
# alt = "T"
# red = "C" 
# get_geneInfo(chr, pos, alt, ref)
# 
# # downstream
# chr = 7
# pos = 75162278
# snp = "rs62477737"
# ref = "G"
# alt = "A"
# get_geneInfo(chr, pos, alt, ref)
# 
# # intronic
# chr = 7
# pos = 75094329
# snp = "rs113160991"
# alt = "A"
# get_geneInfo(chr, pos, alt)
# 
# # splicing
# chr = 1
# pos = 5935162
# snp = "rs1287637"
# ref = "A"
# alt = "T"
# get_geneInfo(chr, pos, alt, ref)
# 
# 
# #exonic : ok
# #splicing :	ok
# #ncRNA : ?? variant overlaps a transcript without coding annotation in the gene definition (see Notes below for more explanation)	non_coding_transcript_variant (SO:0001619)
# #   ncRNA_intronic
# #UTR5	: ok
# #UTR3	: ok
# #intronic :	ok
# #upstream	: potentially two genes?
# #downstream	: potentially two genes?
# # downstream,upstream combined?
# #intergenic : more that two genes?
