library(LFQProcessing)
library(testthat)
library(tools)


protein_viz_is_correct <- function(){
  
  test_that("Test Files outputted as expected", {
    
    get_protein_viz(protein_viz_test_data$prot, 
                    protein_viz_test_data$prot_int, 
                    "tmp/", 
                    protein_viz_test_data$conditionComparisonMapping)
    
    
    current = read_json("tmp/protein_viz.json")
    expected = read_json("../data/expected_protein_viz.json")
    expect_true(all.equal(current, expected))
    
    current = read_json("tmp/protein_counts_and_intensity.json")
    expected = read_json("../data/expected_protein_counts_and_intensity.json")
    expect_true(all.equal(current, expected))

  })
  
}

load("../../data/protein_viz_test_data.rda")

protein_viz_is_correct()