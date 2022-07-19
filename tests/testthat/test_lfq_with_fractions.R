library(LFQProcessing, quietly = TRUE)
library(testthat, quietly = TRUE)


acceptance_test<- function(tolerance = 10**-3){
  
  test_that("LFQ with fractions works", {
    
    # remove it if it exists
    if (file.exists("../data/PXD020248/output/proteinGroups_quant.txt")) {
      #Delete file if it exists
      file.remove("../data/PXD020248/output/proteinGroups_quant.txt")
    }
    
    data_folder <- "../data/PXD020248"
    output_folder <- file.path(data_folder, "output")
    protein_quant_runner(data_folder, 
                         output_folder, 
                         protein_only = TRUE)
    
    current = read.table("../data/PXD020248/output/proteinGroups_quant.txt", sep = "\t", header = TRUE)
    expected = read.table("../data/PXD020248/output/expected_proteinGroups_quant.txt", sep = "\t", header = TRUE)
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })
  
  test_that("The only directories are the figures directories", {
    
    data_folder <- "../data/PXD020248"
    output_folder <- file.path(data_folder, "output")
    
    list_dirs_in_folder <- list.dirs(output_folder, full.names = FALSE)
    expect_true(length(list_dirs_in_folder) == 3)
    expect_true(all(c("", "figure_html", "figure_html_separate") %in% list_dirs_in_folder))
    
  })
  
}


acceptance_test()