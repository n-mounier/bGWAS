context("bGWAS parameters")


object = data.frame(a=c(1,2), b=c("aa", "bb"))
# load dataset
data("SmallGWAS_Timmers2019")
# path to existing Z-matrices
ZMats = system.file("Data/Z_Matrices", package="bGWAS")



#### name ####
test_that("name non-character", {
  expect_error( bGWAS( name=123 , verbose=F), 
                regexp = "name : non-character argument")
})



#### GWAS ####
test_that("GWAS does not exists", { # very unlikely...
  expect_error( bGWAS( name = "aaa", GWAS = non_existing_object, verbose=F),
                regexp = "* not found")
})


test_that("GWAS has uncorrect columns", {
  expect_error( bGWAS( name = "aaa", GWAS = object, verbose=F),
                regexp = "GWAS : no SNPID column")
})


test_that("GWAS file does not exists", {
  expect_error( bGWAS( name="aaa", GWAS="non_existing_file.csv", verbose=F),
                regexp = "GWAS : the file does not exist")
})



#### Z_matrices ####
test_that("ZMatrices do not exist", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = "~/Data/ZMat", verbose=F),
                regexp = "No \"ZMatrix_Full.csv.gz\" file in specified Z_matrices folder")
})


test_that("ZMatrices non-character", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = object, verbose=F),
                regexp = "Z_matrices : wrong format, should be character")
})



#### prior_studies ####
test_that("prior_studies non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       prior_studies = "bbb", verbose=F),
                regexp = "prior_studies : should be numeric if not NULL")
})


test_that("prior_studies not in the list", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       prior_studies = c(1, 2, 3, 123), verbose=F),
                regexp = "prior_studies : all the IDs provided should belong to the ones available")
})



#### MR_threshold ####
test_that("MR_threshold non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_threshold = "bbb", verbose=F),
                regexp = "MR_threshold : non-numeric argument")
})


test_that("MR_threshold out-of-range", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_threshold = 0.1, verbose=F),
                regexp = "MR_threshold : superior to the threshold limit (10^-5)")
})



#### MR_ninstruments ####
test_that("MR_ninstruments non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_ninstruments = "bbb", verbose=F),
                regexp = "MR_ninstruments : non-numeric argument")
})


test_that("MR_ninstruments out-of-range (higher)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_ninstruments = 12, verbose=F),
                regexp = "MR_ninstruments : too large, should be between 2 and 8")
})


test_that("MR_ninstruments out-of-range (smaller)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_ninstruments = 1, verbose=F),
                regexp = "MR_ninstruments : too small, should be between 2 and 8")
})



#### MR_pruning_dist ####
test_that("MR_pruning_dist non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_pruning_dist = "bbb", verbose=F),
                regexp = "MR_pruning_dist : non-numeric argument")
})


test_that("MR_pruning_dist out-of-range (higher)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_pruning_dist = 2000, verbose=F),
                regexp = "MR_pruning_dist : should be lower than 1Mb")
})

test_that("MR_pruning_dist out-of-range (smaller)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_pruning_dist = 5, verbose=F),
                regexp = "MR_pruning_dist : should be higher than 10Kb")
}) 
  


#### MR_pruning_LD ####
test_that("MR_pruning_LD non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_pruning_LD = "bbb", verbose=F),
                regexp = "MR_pruning_LD : non-numeric argument")
})


test_that("MR_pruning_LD out-of-range (higher)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_pruning_LD = 1.3, verbose=F),
                regexp = "MR_pruning_LD : should not be larger than 1")
})


test_that("MR_pruning_dist out-of-range (smaller)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_pruning_LD = -.01, verbose=F),
                regexp = "MR_pruning_LD : should be positive")
}) 
  


#### MR_shrinkage ####
test_that("MR_shrinkage non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_shrinkage = "bbb", verbose=F),
                regexp = "MR_shrinkage : non-numeric argument")
})


test_that("MR_shrinkage out-of-range (higher)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_shrinkage = 1.3, verbose=F),
                regexp = "MR_shrinkage : should not be higher than 1")
})


test_that("MR_shrinkage out-of-range (below MR threshold)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       MR_shrinkage = 1e-8, verbose=F),
                regexp = "MR_shrinkage : should be higher than the threshold used to select MR instruments")
}) 



#### stepwise_threshold ####
test_that("stepwise_threshold non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       stepwise_threshold = "bbb", verbose=F),
                regexp = "stepwise_threshold : should be numeric or NULL")
})


test_that("stepwise_threshold out-of-range (higher)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       stepwise_threshold = 0.1, verbose=F),
                regexp = "stepwise_threshold : should not be higher 0.05")
})


test_that("stepwise_threshold out-of-range (smaller)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       stepwise_threshold = 1e-8, verbose=F),
                regexp = "stepwise_threshold : should not be lower than 0.0005")
})



#### prior_shrinkage ####
test_that("prior_shrinkage non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       prior_shrinkage = "bbb", verbose=F),
                regexp = "prior_shrinkage : non-numeric argument")
})


test_that("prior_shrinkage out-of-range (higher)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       prior_shrinkage = 1.3, verbose=F),
                regexp = "prior_shrinkage : should not be higher than 1")
})


test_that("prior_shrinkage out-of-range (below MR threshold)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       prior_shrinkage = 1e-8, verbose=F),
                regexp = "prior_shrinkage : should be higher than the threshold used to select MR instruments")
})



#### sign_method ####
test_that("sign_method non accepted", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       sign_method = "bbb", verbose=F),
                regexp = "sign_method : method not accepted, should be \"p\" or \"fdr\"")
})

test_that("sign_method non-character", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       sign_method = 123, verbose=F),
                regexp = "sign_method : method not accepted, should be \"p\" or \"fdr\"")
})



#### sign_thresh ####
test_that("sign_thresh non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       sign_thresh = "bbb", verbose=F),
                regexp = "sign_thresh : non-numeric argument")
})


test_that("sign_thresh out-of-range (higher)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       sign_thresh = 1.3, verbose=F),
                regexp = "sign_thresh : a threshold higher than 1 does not make sense")
})


test_that("sign_thresh out-of-range (smaller)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       sign_thresh =-0.01, verbose=F),
                regexp = "sign_thresh : should be positive")
})



#### res_pruning_dist ####
test_that("res_pruning_dist non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       res_pruning_dist = "bbb", verbose=F),
                regexp = "res_pruning_dist : non-numeric argument")
})


test_that("res_pruning_dist out-of-range (higher)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       res_pruning_dist = 2000, verbose=F),
                regexp = "res_pruning_dist : should be lower than 1Mb")
})

test_that("res_pruning_dist out-of-range (smaller)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       res_pruning_dist = 5, verbose=F),
                regexp = "res_pruning_dist : should be higher than 10Kb")
}) 



#### res_pruning_LD ####
test_that("res_pruning_LD non-numeric", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       res_pruning_LD = "bbb", verbose=F),
                regexp = "res_pruning_LD : non-numeric argument")
})


test_that("res_pruning_LD out-of-range (higher)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       res_pruning_LD = 1.3, verbose=F),
                regexp = "res_pruning_LD : should not be larger than 1")
})


test_that("res_pruning_LD out-of-range (smaller)", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       res_pruning_LD = -.01, verbose=F),
                regexp = "res_pruning_LD : should be positive")
}) 



#### save_files ####
test_that("save_files non-logical", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       save_files = "bbb", verbose=F),
                regexp = "save_files : should be logical")
})


#### verbose ####
test_that("verbose non-logical", {
  expect_error( bGWAS( name="aaa", GWAS=SmallGWAS_Timmers2019, Z_matrices = ZMats, 
                       verbose=123),
                regexp = "verbose : should be logical")
})