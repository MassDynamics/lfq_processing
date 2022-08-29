library(LFQProcessing, quietly = TRUE)
library(testthat, quietly = TRUE)


acceptance_test<- function(current, expected,
                           tolerance = 10**-3){
  
  test_that("approximately equal", {
    
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })
  
}


acceptance_qcs_tmt <- function(expected_qcs){
  
  test_that("All expected QCs for TMT protein only workflow are present", {
    
    all_qcs <- sapply(expected_qcs, function(x) file.exists(x))
    expect_true(all(all_qcs))
    
  })
}

output_folder <- "./tmp"
current_tmt <- tmt_transformer(protein_groups_tmt_PXD019880,
                               des_tmt_PXD019880,
                               output_folder,
                               imputeStDev=0.3,
                               imputePosition=1.8)

expdes <- des_tmt_PXD019880
prot <- current_tmt[[1]]
prot_int <- current_tmt[[2]]
expdes <- current_tmt[[3]]
conditionComparisonMapping <- current_tmt[[4]]

# get qc report from package
file.copy(from=system.file("rmd","tmt_qc_report.Rmd", package = "LFQProcessing"),
          to=file.path(output_folder,"QC_Report.Rmd"),
          overwrite = TRUE, recursive = TRUE,
          copy.mode = TRUE)

print("\nCopy separate QCs\n")
names_tmt_protein_only <- LFQProcessing:::get_names_qc_tmt_protein_only(evidence=FALSE)

for(qc_name in names_tmt_protein_only){
  qc_report_name <- paste0("QC_", qc_name, ".Rmd")
  file.copy(from=system.file("rmd",qc_report_name, package = "LFQProcessing"),
            to=file.path(output_folder,qc_report_name),
            overwrite = TRUE, recursive = TRUE,
            copy.mode = TRUE)
  
  # I need to run it here as I need datasets loaded in the environment
  rmarkdown::render(file.path(output_folder, qc_report_name), 
                    params = list(output_figure = "figure_html_separate/"),
                    output_format = rmarkdown::html_document(
                      self_contained=FALSE, 
                      lib_dir=file.path(output_folder, "qc_report_files"))
  )
  
}

output_format = "html"
rmarkdown::render(file.path(output_folder, "QC_Report.Rmd"),
                  params = list(output_figure = "figure_html/"),
                  output_format = rmarkdown::html_document(
                    self_contained=TRUE,
                    code_folding= "hide",
                    theme="united",
                    toc = TRUE,
                    toc_float = TRUE,
                    fig_caption= TRUE,
                    df_print="paged",
                    lib_dir="qc_report_files"))

output_format = "pdf"
rmarkdown::render(file.path(output_folder, "QC_Report.Rmd"),
                  output_format=rmarkdown::pdf_document(
                    toc = TRUE,
                    fig_caption= TRUE),
                  params = list(output_figure = "figure_pdf/"))

#plot(current_tmt[[1]]$`logFC Cerebellum - Cerebrum`, -log10(current_tmt[[1]]$`P.Value Cerebellum - Cerebrum`))

load("../data/expected_tmt_PXD019880.rda")

tmt_split_qc_protein_only <- c("TMT_PCA_proteins", "TMT_PCA_screeplot_proteins", "TMT_PCA_DE_proteins", 
                               "TMT_CV_proteins",
                               "TMT_missing_by_proteins",
                               "TMT_missing_by_samples_proteins",
                               
                               "TMT_intensities_boxplots_proteins",
                               "TMT_intensities_normalised_boxplots_proteins",
                               
                               "TMT_imputed_proteins",
                               "TMT_identifications_proteins")

expected_qc_tmt <- file.path(output_folder, 
                             c(paste0("QC_",tmt_split_qc_protein_only, ".html"),
                               "QC_Report.html", "QC_Report.pdf"))
                             
acceptance_test(current_tmt[[1]], expected_tmt[[1]])
acceptance_qcs_tmt(expected_qc_tmt)

list_lib_files <- list.files(output_folder, full.names = TRUE)
unlink(list_lib_files, recursive = TRUE)
