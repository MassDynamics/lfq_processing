#prot = fasta_parser(prot)

# clone prot object

# 1. determine which parser to use
# 2. apply the oarser to each column

# for row in cloned_prot:
#   t---mp_row = copy(row)
#   cloned_prot$row <- apply(row_id_parser, row)

# return cloned_prot

# disambiguate, assume no contaminants

#' @import data.table
#' @import dplyr
#' @export parse_id_columns
parse_id_columns <- function(prot){
  
  # remove extra fasta headers
  prot$`fasta headers` = sapply(strsplit(prot$`fasta headers`,";"), `[`, 1)
  
  prot$fasta_header_type <- unlist(lapply(prot$`fasta headers`, match_pattern_name))
  
  prot$Accession_id = ""
  prot$ProteinDescription = ""
  prot$Gene = ""
  
  prot <- prot %>% mutate(
    Accession_id = case_when(
      fasta_header_type=="uniprot" ~ gsub('\\|', '', str_extract(prot$`fasta headers`, '\\|(.*)\\|')),
      fasta_header_type=="accession" ~  str_extract(prot$`fasta headers`, '([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\\-?[:digit:]*'),
      fasta_header_type=="hppr" ~  prot$`fasta headers`,
      TRUE ~ sapply(strsplit(prot$`majority protein ids`, ";"), "[", 1)),
    Gene = case_when(
      fasta_header_type=="uniprot" ~ gsub("GN=", "", str_extract_all(prot$`fasta headers`, 'GN=([^\\s]+)')),
      TRUE ~ ""),
    ProteinDescription = case_when(
      fasta_header_type=="uniprot" ~ gsub("\\sOS=","",str_extract(prot$`fasta headers`, "\\s.*\\sOS=")),
      TRUE ~  prot$`fasta headers`),
  )
  prot <- prot %>% mutate(
    Accession_id = case_when(
      is.na(Accession_id) ~ prot$`fasta headers`,
      TRUE ~ prot$Accession_id
    ))
  
  prot
}

#' @export match_pattern_name
match_pattern_name <- function(fasta_header){
  
  uniprot = '(\\w{2})\\|(\\w*)\\-?[0-9]*\\|(\\w*) (.*) (OS=.*)'
  con = "(CON__.*)"
  accession = '([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\\-?[:digit:]*'
  hppr = '(HPRR\\d*)'
  missing_semi_colons = "^\\;*$"
  
  if (is.null(fasta_header)){
    return("missing")
  } else if (is.na(fasta_header)){
    return("missing")
  } else if (fasta_header == ""){
    return("missing")
  } else if (grepl(missing_semi_colons, fasta_header)){
    return("missing")
  } else if (grepl(uniprot, fasta_header)){
    return("uniprot")
  } else if (grepl(con, fasta_header)){
    return("con")
  } else if (grepl(accession, fasta_header)){
    return("accession")
  } else if (grepl(hppr, fasta_header)){
    return("hppr")
  } else {
    return ("not_found")
  }
}
