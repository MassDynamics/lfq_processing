library(LFQProcessing, quietly = TRUE)
library(testthat, quietly = TRUE)


acceptance_test<- function(tolerance = 10**-3){
  
  test_that("LFQ Experiment where we remove some samples in the experiment design", {
    
    # remove it if it exists
    if (file.exists("../data/PXD026936/output/proteinGroups_quant.txt")) {
      #Delete file if it exists
      file.remove("../data/PXD026936/output/proteinGroups_quant.txt")
    }
    
    protein_quant_runner("../data/PXD026936", file.path("../data/PXD026936", "output"), protein_only = TRUE)
    
    current = read.table("../data/PXD026936/output/proteinGroups_quant.txt", sep = "\t", header = TRUE)
    expected = read.table("../data/PXD026936/output/expected_proteinGroups_quant.txt", sep = "\t", header = TRUE)
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })
  
}

acceptance_test()