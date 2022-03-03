library(LFQProcessing, quietly = TRUE)
library(testthat, quietly = TRUE)


acceptance_test<- function(tolerance = 10**-3){
  
  test_that("TMT experiment without mq 'experiment' label works", {
    
    # remove it if it exists
    if (file.exists("../data/PXD031774/output/proteinGroups_quant.txt")) {
      #Delete file if it exists
      file.remove("../data/PXD031774/output/proteinGroups_quant.txt")
    }
    
    protein_quant_runner("../data/PXD031774", file.path("../data/PXD031774", "output"))
    
    current = read.table("../data/PXD031774/output/proteinGroups_quant.txt", sep = "\t", header = TRUE)
    expected = read.table("../data/PXD031774/output/expected_proteinGroups_quant.txt", sep = "\t", header = TRUE)
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })
  
}

acceptance_test()