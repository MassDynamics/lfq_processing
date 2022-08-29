#' Write Output of quantitative analysis
#' @param data a data.table to write
#' @param output_folder An output folder to store produced files.
#' @return Write the table to disk in the output folder
#' @example
#' write_output(prot, output_folder, "proteinGroups_quant.txt")
#' @export write_output

write_output <- function(data, output_folder, file_name){

  fwrite(data, file.path(output_folder, file_name), row.names = FALSE,
         sep = "\t", quote = FALSE,
         showProgress = T, verbose = T, nThread = detectCores() - 1)

}
