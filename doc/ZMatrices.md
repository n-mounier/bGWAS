# Z-Matrix files
[//]:========================================

The Z-Matrix files are essential to create the prior. They contain summary statistics of Prior GWASs.

## Description
[//]:*******

The compressed file available contains 3 files:    
- **AvailableStudies.tsv**   
- **ZMatrix_Full.csv.gz**   
- **ZMatrix_MR.csv.gz**     
     
       
The first one, **AvailableStudies.tsv**, contains descriptive information about each prior GWAS (name of the file, trait, publication and download links). The "N_SNPs" colum contains the number of well-imputed SNPs included in the Full Z-matrix, and the "N_Instruments" contains the number of SNPs reaching p<1e-5 and included in the MR Z-Matrix for each prior GWAS.     
The second one, **ZMatrix_Full.csv.gz**, contains genome-wide summary statistics for all Prior GWASs. These summary statitics have been obtained imputing the initial summary statistics results using [SSimp](https://github.com/zkutalik/ssimp_software) version 0.1 and UK10K data as a reference panel. All SNPs with UK10K allele frequency above 0.5 \% and imputation quality higher than 0.8 for for least 20 studies have been kept. Imputed results below this threshold have been set to NA and use as 0 when calculating the prior.   
The last file, **ZMatrix_MR.csv.gz**, is actually a subset of the previous one. It contains only SNPs with complete observations (no missing value for any of the Prior GWASs) that are strongly associated (p<1e-5) with at least one of the Prior GWASs. It is used for the stepwise selection approach, to avoid loading a big file. Once the significant studies are identified, only relevant columns (Prior GWASs) of the **ZMatrix_Full.csv.gz** will be loaded to speed up prior estimation.
  
  
  
## Download
[//]:*******

You can download these files using this [link](https://drive.switch.ch/index.php/s/jvSwoIxRgCKUSI8) or following the instructions below.    

- On UNIX/MACOSX, from a terminal:    
``` bash
wget https://drive.switch.ch/index.php/s/jvSwoIxRgCKUSI8/download -O ZMatrices.tar.gz
tar xzvf ZMatrices.tar.gz
``` 
<!--- - On WINDOWS, from a terminal:   
``` bash
...
```   --->


## How to use your own data
[//]:*******

We do not provide scripts to modify these files, but the package can accept customized Z-Matrix as long as they use the same format. Please, make sure to:   
1. Modify the **AvailableStudies.tsv** file, to include your own files. The "File" column must correspond to the column names of the other two files, and ID should match the order of the columns in these files. The "Name"" and	"Trait"	 columns are needed, but you could leave the other ones empty.    
2. Create the largest Z-Matrix first, by merging all your summary statistics. The first 5 columns should be "rs", "chrm", "pos", "alt" and "ref".    
3. Subset this big matrix file to keep only strong instruments (p<1e-5 for at least one prior GWAS).    

  
## Information about the files
[//]:*******
Details about the prior GWASs used can be found in the **AvailableStudies.tsv** file or obtained directly in `R` using `**list_priorGWASs()**`.   
Please note that now, we do use 38 studies now, instead of the 58 described in McDaid et al. Removal reasons are listed below.

|            Prior GWAS            | Removal Reason(s) |    
| -------------------------------- | ----------------- |    
|   bip1.scz1.ruderfer2014..ANY    | pooled analysis of bipolar disorder/schizophrenia cases, against controls  |   
| bip1.scz1.ruderfer2014..bp_v_scz | analysis of bipolar disorder vs schizophrenia  |
|      DIAGRAMv3.2013MAY07.zip     | metabochip study, badly imputed  |
| MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt...BMI  | adjusted for BMI, could lead to bias  |
| MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt...BMI  | adjusted for BMI, could lead to bias  |
|          pgc.cross.aut           | cross-disorder analysis  |
|    pgc.cross.BIP11.2013_05.txt   | cross-disorder analysis  |
|      pgc.cross.full.2013_03      | cross-disorder analysis  |
| pgc.cross.scoring.5files..add_scz.zip | cross-disorder analysis  |
| pgc.cross.scoring.5files..aut_bip.zip | cross-disorder analysis  |
| pgc.cross.scoring.5files..bip_scz.zip | cross-disorder analysis  |
| pgc.cross.scoring.5files..mdd_scz.zip | cross-disorder analysis  |
| pgc.cross.scoring.5files..scz_bip.zip | cross-disorder analysis  |
|           pgc.cross.scz         | cross-disorder analysis  |
|         pgc.scz.2012_04         | SCZ1 file, part of the meta-analysis SCZ1 + SW that is already included  |
|   phs000124.pha002845.txt.ugz   | newest version of the study already included  |
| TransEthnic_T2D_GWAS.MegaMeta.2014OCT16.zip | non european ancestry |



## Contact
<mounier.ninon@gmail.com>

