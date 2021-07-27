library(LFQProcessing)
library(testthat)

# Run Code
current_her2 =  lfq_transformer(ma_tables = example_lfq_her2_targetted_therapy_tables,
                                output_folder = "./tmp",
                                imputeStDev=0.3,
                                imputePosition=1.8)


acceptance_test<- function(current, expected,
                           tolerance = 10**-3){
  
  test_that("approximately equal", {
    
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })
  
}

load("../data/expected_her2.rda")

acceptance_test(current_her2, expected_her2, 10**-3)