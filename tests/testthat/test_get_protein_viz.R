library(LFQProcessing)
library(testthat)
library(tools)

data("protein_viz_test_data")
test_that("Test Files written as expected", {
  
  get_protein_viz(protein_viz_test_data$prot, 
                  protein_viz_test_data$prot_int, 
                  "tmp/", 
                  protein_viz_test_data$conditionComparisonMapping)
  
  
  current = read_json("tmp/protein_viz.json", simplifyVector = T)
  expected = read_json("../data/expected_protein_viz.json", simplifyVector = T)
  expect_true(all.equal(current$data[1]$FastaHeaders, expected$data[1]$FastaHeaders))
  expect_true(all.equal(current$data[1]$FoldChange, expected$data[1]$FoldChange))
  expect_true(all.equal(current$data[1]$AdjustedPValue, expected$data[1]$AdjustedPValue))
  expect_true(all.equal(current$data[1]$ProteinGroupId, expected$data[1]$ProteinGroupId))
  expect_equal(colnames(current), c("conditionComparison","up.condition","down.condition","fdrLimit","data"))
  
  current = read_json("tmp/protein_counts_and_intensity.json", simplifyVector = T)
  expected = read_json("../data/expected_protein_counts_and_intensity.json", simplifyVector = T)
  expect_true(all.equal(current, expected))
  
})
