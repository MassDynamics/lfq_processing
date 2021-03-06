limma_stats_fun <- function(ID_type, 
                            int_type, 
                            condition_col_name,
                            run_id_col_name,
                            rep_col_name,
                            funDT,
                            pairwise.comp = NULL,
                            all.comparisons = TRUE,
                            fix_distr = FALSE) {
  
  # # DEBUG Vars
  #
  # ID_type = "id"
  # int_type = "log2NInt"
  # condition_col_name = "experiment"
  # rep_col_name = "Replicate"
  # run_id_col_name = "run_id"
  # funDT = copy(dt_int)
  # all.comparisons = TRUE
  # fix_distr = FALSE
  # pairwise.comp = NULL # pairwisecom

  # pairwise.comp <- data.table(
  #   left = c('12AS_NMPP1', "13AS_NMPP1", "D2AS_NMPP1", "WT_NMPP1"),
  #   right = c("12AS_DMSO", "13AS_DMSO", "D2AS_DMSO", "WT_DMSO")
  # )
  # 
  
  if (all.comparisons) {
    comination_mat <- combn(x = funDT[, unique(get(condition_col_name))], 2)
    pairwise.comp <- foreach (i = 1:ncol(comination_mat)) %do% {
      return(data.table(
        left = comination_mat[1,i],
        right = comination_mat[2,i]
      ))
    }
    pairwise.comp <- rbindlist(pairwise.comp)
  }
  
  
  library(Biobase)
  library(limma)
  library(data.table)
  #library(qvalue)
  
  funDT <- funDT[str_order(get(run_id_col_name), numeric = T)]
  funDT[, ID := str_c("ID.", get(ID_type))]
  
  # Filter for IDs that are not present in at least one experiment
  filterDT <- unique(funDT[, .(get(condition_col_name), run_id)])[, .(max_count = .N), by = .(V1)]
  filterDT <- merge(funDT, filterDT, by.x = condition_col_name, by.y = "V1", all.x = T)
  filterDT <- filterDT[Imputed == 0, .(count_rep = .N, max_count = max(max_count, na.rm = T)), by = c(condition_col_name, "ID")][, repPC := count_rep/max_count]
  
  isPresent <- filterDT[repPC >= 0.5, unique(ID)]
  funDT <- funDT[ID %in% isPresent]
  
  fun.formula <- as.formula(str_c("ID ~ ", run_id_col_name))
  
  
  int_matrix <-
    dcast.data.table(funDT, fun.formula, value.var = int_type)
  int_matrix <- int_matrix[str_order(ID, numeric = T)]
  intrawnames <- int_matrix[, ID]
  intcolnames <- colnames(int_matrix[, 2:ncol(int_matrix)])
  intcolnames <- str_sort(intcolnames, numeric = TRUE)
  
  int_matrix <- as.matrix(int_matrix[, intcolnames, with = FALSE])
  rownames(int_matrix) <- intrawnames
  
  ### Imputed value mask matrix ------
  imputed_matrix <-
    dcast.data.table(funDT, fun.formula, value.var = "Imputed")
  imputed_matrix <- imputed_matrix[str_order(ID, numeric = T)]
  imp_mat_rawnames <- imputed_matrix[, ID]
  imp_mat_colnames <- colnames(imputed_matrix[, 2:ncol(imputed_matrix)])
  imp_mat_colnames <- str_sort(imp_mat_colnames, numeric = TRUE)
  
  imputed_matrix <- as.matrix(imputed_matrix[, imp_mat_colnames, with = FALSE])
  rownames(imputed_matrix) <- imp_mat_rawnames
  isZ <- imputed_matrix == 1
  
  sample_dt <- as.data.frame(unique(funDT[, .(
    run_id = get(run_id_col_name),
    condition = get(condition_col_name),
    Replicate = get(rep_col_name)
  )]))
  rownames(sample_dt) <- sample_dt$run_id
  
  # feature_dt <- as.data.frame(unique(funDT[, .(ID,
  #                                              otherinfo = NA)]))
  feature_dt <- data.frame(ID = rownames(int_matrix), otherinfo = NA)
  
  rownames(feature_dt) <- feature_dt$ID
  
  # Create ExpressionSet object
  eset <- ExpressionSet(
    assayData = int_matrix,
    phenoData = AnnotatedDataFrame(sample_dt),
    featureData = AnnotatedDataFrame(feature_dt)
  )
  
  
  design.mat <- model.matrix(~ 0 + condition,
                             data = pData(eset))
  myContrasts = NULL
  for (irow in 1:pairwise.comp[, .N]) {
    left <- pairwise.comp[irow, left]
    right <- pairwise.comp[irow, right]
    
    myContrasts = c(myContrasts,
                    str_c("condition", left, " - ", "condition", right))
  }
  
  
  contrast.matrix <- eval(as.call(c(
    as.symbol("makeContrasts"),
    as.list(myContrasts),
    levels = list(design.mat)
  )))
  
  fit <- lmFit(eset, design.mat)
  
  # exper.cond <- "UPS1"
  # fix_distr <- TRUE
  if (fix_distr) {
    for (exper.cond in pairwise.comp[, c(left, right)]) {
      Runs <- pData(eset)[pData(eset)$condition == exper.cond, "run_id"]
      i <- which(rowSums(isZ[,colnames(isZ) %in% Runs]) == length(Runs))
      eset_fix <- exprs(eset[i,])
      eset_fix[,colnames(isZ) %in% Runs] <- NA
      fit3 <- lmFit(eset_fix, design.mat)
      fit$sigma[i] <- fit3$sigma
      fit$df.residual[i] <- fit3$df.residual
    }
  }
  
  fit2 <- contrasts.fit(fit, contrasts = contrast.matrix)
  fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
  
  stats <-
    topTable(fit2,
             number = nrow(fit2),
             sort.by = "none"
    )
  stats <- as.data.table(stats)
  stats <- stats[, .(ID, AveExpr, F, adj.P.Val)]
  #ipair = 1
  #### pairwise comparisons
  for (ipair in 1:pairwise.comp[, .N]) {
    #print(ipair)
    subsecting <- funDT[get(condition_col_name) %in% pairwise.comp[ipair, c(left, right)], unique(run_id)]
    
    # Filter for IDs that are not present in at least one experiment in pairwise manner
    isPresent <- filterDT[get(condition_col_name) %in% pairwise.comp[ipair, c(left, right)] & repPC >= 0.5, unique(ID)]
    
    if (length(isPresent) > 0) {
      eset_pair <- eset[rownames(eset) %in% isPresent, colnames(eset) %in% subsecting]
      
      
      
      design.mat <- model.matrix(~ 0 + condition,
                                 data = pData(eset_pair))
      
      
      myContrasts = NULL
      left <- pairwise.comp[ipair, left]
      right <- pairwise.comp[ipair, right]
      myContrasts = c(myContrasts,
                      str_c("condition", left, " - ", "condition", right))
      contrast.matrix <- eval(as.call(c(
        as.symbol("makeContrasts"),
        as.list(myContrasts),
        levels = list(design.mat)
      )))
      
      fit <- lmFit(eset_pair, design.mat)
      fit2 <- contrasts.fit(fit, contrasts = contrast.matrix)
      fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
      
      s_dt <-
        as.data.table(topTable(
          fit2,
          number = nrow(fit2),
          sort.by = "p",
          confint = 0.95,
          coef = myContrasts
        ))
      s_dt <- s_dt[, .(ID, logFC, CI.L, CI.R, P.Value, adj.P.Val)]
      #q_vals <- qvalue(s_dt[, P.Value])
      #s_dt[, Q.Value := q_vals$qvalues]
      
      comparison <- str_replace_all(myContrasts, "condition", "")
      colnames(s_dt)[colnames(s_dt) != "ID"] <-
        str_c(colnames(s_dt)[colnames(s_dt) != "ID"], comparison, sep = " ")
      stats <- merge(stats, s_dt, by = "ID", all = T)
    }
    
   
  }
  
  
 
  
  
  stats[, ID := str_replace_all(ID, "ID.", "")]
  setnames(stats, "ID", ID_type)
  return(stats)
}