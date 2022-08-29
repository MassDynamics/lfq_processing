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

    protein_quant_runner(upload_folder = "../data/PXD020248", 
                         output_folder = "../data/PXD020248/output", 
                         protein_only = TRUE)
    
    current = read.table("../data/PXD020248/output/proteinGroups_quant.txt", sep = "\t", header = TRUE)
    expected = read.table("../data/PXD020248/output/expected_proteinGroups_quant.txt", sep = "\t", header = TRUE)
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })
  
  test_that("All QCs in LFQ Fractions are created", {
  
  expected_qcs <- file.path(output_folder, 
                            c(paste0("QC_",LFQProcessing:::get_names_qc_lfq_protein_only(), ".html"),
                              "QC_Report.html", "QC_Report.pdf"))
  qc_exists <- sapply(expected_qcs, function(x) file.exists(x))
  expect_true(all(qc_exists))
  })
}


acceptance_test()