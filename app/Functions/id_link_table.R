make_id_link_table <- function(output_folder, pept) {
  
  
  ## === Create link table ================
  cl <- makeCluster(detectCores() - 1, type = "FORK")
  registerDoParallel(cl)
  
  
  id_cols <- c("id", grep("ids", colnames(pept), value = T))
  id_cols <- id_cols[!id_cols %in% c("ms/ms ids", "oxidation (m) site ids", "taxonomy ids")]
  id_cols_loop <- id_cols[id_cols != "id"]
  
  link_dt <- foreach(inx = seq(1, pept[, .N], 1), .packages = c("data.table")) %dopar% {
    # print(inx)
    DT <- pept[inx, id_cols, with=FALSE]
    one_col <- id_cols_loop[1]
    dt <-  data.table(
      id = DT[, id],
      unlist(
        str_split(DT[, one_col, with=FALSE], ";")
      )
    )
    setnames(dt, "V2", one_col)
    for (one_col in id_cols_loop[2:length(id_cols_loop)]) {
      # print(one_col)
      dt <- merge(
        dt,
        data.table(
          id = DT[, id],
          unlist(
            str_split(DT[, one_col, with=FALSE], ";")
          )
        ),
        by = "id",
        all = T,
        allow.cartesian=TRUE
      )
      setnames(dt, "V2", one_col)
      
    }
    return(dt)
  }
  
  link_dt <- rbindlist(link_dt)
  setnames(link_dt, "id", "peptides id")
  fwrite(link_dt, file.path(output_folder, 'MQ_ids_link_table.txt'), row.names = FALSE, sep = "\t", quote = FALSE, showProgress = T, verbose = T, nThread = detectCores() - 1)
  rm(link_dt)
  gc()
  stopCluster(cl)
  return(NA)
}



