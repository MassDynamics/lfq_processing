library(LFQProcessing, quietly = TRUE)
library(testthat, quietly = TRUE)
library(here)


acceptance_test<- function(tolerance = 10**-3){
  
  test_that("LFQ with fractions works", {
    
    # remove it if it exists
    if (file.exists("../data/PXD020248/output/proteinGroups_quant.txt")) {
      #Delete file if it exists
      file.remove("../data/PXD020248/output/proteinGroups_quant.txt")
    }
    
    data_folder <- file.path(here(), "tests/data/PXD020248")
    output_folder <- file.path(data_folder, "output")
    protein_quant_runner(upload_folder = "../data/PXD020248", 
                         output_folder = "../data/PXD020248/output", 
                         protein_only = TRUE)
    
    current = read.table("../data/PXD020248/output/proteinGroups_quant.txt", sep = "\t", header = TRUE)
    expected = read.table("../data/PXD020248/output/expected_proteinGroups_quant.txt", sep = "\t", header = TRUE)
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })
  
}


acceptance_test()