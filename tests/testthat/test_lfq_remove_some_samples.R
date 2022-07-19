library(LFQProcessing, quietly = TRUE)
library(testthat, quietly = TRUE)


acceptance_test<- function(tolerance = 10**-3){
  
  test_that("LFQ Experiment where we remove some samples in the experiment design", {
    
    # remove it if it exists
    if (file.exists("../data/PXD026936/output/proteinGroups_quant_intensities.txt")) {
      #Delete file if it exists
      file.remove("../data/PXD026936/output/proteinGroups_quant_intensities.txt")
    }
    
    protein_quant_runner("../data/PXD026936", file.path("../data/PXD026936", "output"), protein_only = TRUE)
    
    current = read.table("../data/PXD026936/output/proteinGroups_quant_intensities.txt", sep = "\t", header = TRUE)
    expected = read.table("../data/PXD026936/output/expected_proteinGroups_quant_intensities.txt", sep = "\t", header = TRUE)
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
    
    expect_false("OE0906_28375" %in% current$file_name) # should be out
    expect_true("OE0906_28349" %in% current$file_name) # should be in
    
  })
  
  test_that("The only directories are the figures directories", {
    
    list_dirs_in_folder <- list.dirs(file.path("../data/PXD026936", "output"), full.names = FALSE)
    expect_true(length(list_dirs_in_folder) == 3)
    expect_true(all(c("", "figure_html", "figure_html_separate") %in% list_dirs_in_folder))
    
  })
  
}

acceptance_test()
