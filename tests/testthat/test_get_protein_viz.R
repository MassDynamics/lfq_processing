library(LFQProcessing)
library(testthat)
library(tools)


protein_viz_is_correct <- function(){
  
  test_that("Test Files outputted as expected", {
    
    get_protein_viz(protein_viz_test_data$prot, 
                    protein_viz_test_data$prot_int, 
                    "tmp/", 
                    protein_viz_test_data$conditionComparisonMapping)
    
    expect_true(md5sum("tmp/protein_viz.json")==md5sum("../data/expected_protein_viz.json"))
    expect_true(md5sum("tmp/protein_counts_and_intensity.json")==md5sum("../data/expected_protein_counts_and_intensity.json"))
  })
  
}

load("../../data/protein_viz_test_data.rda")

protein_viz_is_correct()