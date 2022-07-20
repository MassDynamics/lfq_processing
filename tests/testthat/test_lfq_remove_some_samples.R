library(LFQProcessing, quietly = TRUE)
library(testthat, quietly = TRUE)
library(here)


acceptance_test<- function(tolerance = 10**-3){
  
  test_that("LFQ Experiment where we remove some samples in the experiment design", {
    
    # remove it if it exists
    if (file.exists("../data/PXD026936/output/proteinGroups_quant_intensities.txt")) {
      #Delete file if it exists
      file.remove("../data/PXD026936/output/proteinGroups_quant_intensities.txt")
    }
    
    data_folder <- file.path(here(), "tests/data/PXD026936")
    output_folder <- file.path(data_folder, "output")
    protein_quant_runner("../data/PXD026936", output_folder, protein_only = TRUE)
    
    current = read.table("../data/PXD026936/output/proteinGroups_quant_intensities.txt", sep = "\t", header = TRUE)
    expected = read.table("../data/PXD026936/output/expected_proteinGroups_quant_intensities.txt", sep = "\t", header = TRUE)
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
    
    expect_false("OE0906_28375" %in% current$file_name) # should be out
    expect_true("OE0906_28349" %in% current$file_name) # should be in
    
  })
  
  
}

acceptance_test()
