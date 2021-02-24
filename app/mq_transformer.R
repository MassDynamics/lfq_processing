library(data.table)
library(foreach)
library(parallel)
library(doParallel)
library(stringr)
library(broom)
library(FNN)
library(Hmisc)
library(R.utils)
library(bit64)

mq_transfomer <- function(app_folder, mq_folder, expdes_filename, imputeStDev=0.3, imputePosition=1.8) {
  
  # DEBUG 
  # imputeStDev=0.3
  # imputePosition=1.8
  # expdes_filename = experimentalDesign_filename
  
  checkBy <- 'Experiment'
  finctions_folder <- file.path(app_folder, "app/Functions")
  sourceDirectory(finctions_folder)
  
  output_folder <- file.path(mq_folder, "results")
  dir.create(output_folder, showWarnings = FALSE)
  
  
  # read data
  ma_tables <- read_data(mq_folder, expdes_filename)
  msms <- ma_tables[[1]]
  prot <- ma_tables[[2]]
  pept <- ma_tables[[3]]
  mod_pept <- ma_tables[[4]]
  evidence <- ma_tables[[5]]
  expdes <- ma_tables[[6]]
  
  ###### ModPep ######
  mod_pep_list <- table_quant_analysis(mod_pept, expdes, id_var = "id", output_folder, "modificationSpecificPeptides_quant.txt", "modificationSpecificPeptides_quant_intensities.txt")
  mod_pept <- mod_pep_list[[1]]
  mod_pept_int <- mod_pep_list[[2]]
  
  ###### Peptides #######
  pept_list <- table_quant_analysis(pept, expdes, id_var = "id", output_folder, "peptides_quant.txt", "peptides_quant_intensities.txt")
  pept <- pept_list[[1]]
  pept_int <- pept_list[[2]]
  
  ###### Proteins ######
  prot_list <- table_quant_analysis(prot, expdes, id_var = "id", output_folder, "proteinGroups_quant.txt", "proteinGroups_quant_intensities.txt")
  prot <- prot_list[[1]]
  prot_int <- prot_list[[2]]
  
  ###### phospho ######
  # phospho_enrichment_quant(upload_folder, expdes, output_folder)
  
  # id_table <- make_id_link_table(output_folder, pept)
  
  # QC Report START
  rmarkdown::render(file.path(app_folder, "app/QC_Report.Rmd"))
  # Organise output directory
  file.copy(
    file.path(app_folder, "app/QC_Report.html"),
    output_folder
  )
  file.remove(file.path(app_folder, "app/QC_Report.html"))
  
}