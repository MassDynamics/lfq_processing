library(LFQProcessing)
library(testthat)

# Run Code
current_her2 =  lfq_transformer(ma_tables = example_lfq_her2_targetted_therapy_tables,
                                output_folder = "./tmp",
                                imputeStDev=0.3,
                                imputePosition=1.8,
                                protein_only=FALSE)

prot = copy(current_her2$prot)
prot_int = copy(current_her2$prot_int)
expdes <- copy(current_her2$expdes)
pept = copy(current_her2$pept)
pept_int = copy(current_her2$pept_int)
mod_pept = copy(current_her2$mod_pept)
mod_pept_int = copy(current_her2$mod_pept_int)
evidence = copy(current_her2$evidence)
msms = copy(current_her2$msms)
protein_only = FALSE

current_her2 = unname(current_her2)


# QC reports
# get qc report from package
output_folder = "./tmp"
file.copy(from=system.file("rmd","QC_Report.Rmd", package = "LFQProcessing"),
          to=output_folder,
          overwrite = TRUE, recursive = TRUE,
          copy.mode = TRUE)

print("\nCopy all separate QCs\n")
all_qcs <- LFQProcessing:::get_names_qc_lfq_all()

for(qc_name in all_qcs){
  qc_report_name <- paste0("QC_", qc_name, ".Rmd")
  cat(paste("Writing ", qc_report_name))
  file.copy(from=system.file("rmd", qc_report_name, package = "LFQProcessing"),
            to=file.path(output_folder, qc_report_name),
            overwrite = TRUE, recursive = TRUE,
            copy.mode = TRUE)
  
  # I need to run it here as I need datasets loaded in the environment
  rmarkdown::render(file.path(output_folder, qc_report_name), 
                    params = list(output_figure = "figure_html_separate/"),
                    output_format = rmarkdown::html_document(
                      self_contained=FALSE, 
                      lib_dir="qc_report_files", 
                      theme = "united",
                      fig_caption = TRUE,
                      df_print = "paged"))
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


# Tests

acceptance_test<- function(current, expected,
                           tolerance = 10**-3){
  
  test_that("approximately equal", {
    
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })
  
}

acceptance_qcs_lfq <- function(expected_qcs){
  
  test_that("All expected QCs for LFQ MaxQuant workflow are present", {
    
    all_qcs <- sapply(expected_qcs, function(x) file.exists(x))
    expect_true(all(all_qcs))
    
  })
}

load("../data/expected_her2.rda")

all_lfq_split_qc <- c("PCA_proteins", "PCA_screeplot_proteins", "PCA_DE_proteins",
                      "PCA_peptides", "PCA_screeplot_peptides", "PCA_DE_peptides",
                      "PCA_modified_peptides", "PCA_screeplot_modified_peptides", "PCA_DE_modified_peptides",
                      "CV_proteins", "CV_peptides", "CV_modified_peptides",
                      "missing_by_proteins", "missing_by_peptides", "missing_by_modified_peptides",
                      "missing_by_samples_proteins", "missing_by_samples_peptides", "missing_by_samples_modified_peptides",

                      "samples_correlations_proteins",
                      "samples_correlations_DE_proteins",
                      "samples_correlations_scatter_proteins",
                      "samples_correlations_peptides",
                      "samples_correlations_DE_peptides",
                      "samples_correlations_scatter_peptides",
                      "samples_correlations_modified_peptides",
                      "samples_correlations_DE_modified_peptides",
                      "samples_correlations_scatter_modified_peptides",

                      "identifications_proteins", "identifications_proteins_evidence",
                      "identifications_peptides", "identifications_peptides_evidence",
                      "identifications_modified_peptides", "identifications_modified_peptides_evidence",
                      "identifications_PSMs",

                      "missed_cleavages_evidence")

expected_qc_lfq <- file.path(output_folder,
                             c(paste0("QC_",all_lfq_split_qc, ".html"),
                               "QC_Report.html", "QC_Report.pdf"))

all_qcs <- sapply(expected_qc_lfq, function(x) file.exists(x))

acceptance_test(current_her2, expected_her2, 10**-3)
acceptance_qcs_lfq(expected_qc_lfq)

list_lib_files <- list.files(output_folder, full.names = TRUE)
unlink(list_lib_files, recursive = TRUE)
