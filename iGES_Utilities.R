#Immune Gene Expression Signature Functions
# Peter Fiorica
# 2 January 2026

sig_score <- function(signature.file, GEdata){
  
  ## Ensure numeric matrix
  GEdata <- as.matrix(GEdata)
  storage.mode(GEdata) <- "numeric"
  
  sigdat <- signature.file
  
  outdat <- list()
  coverage <- numeric(nrow(sigdat))
  sig_names <- sigdat[[1]]   # first column = signature name
  
  for(i in seq_len(nrow(sigdat))){
    
    sig_genes <- as.character(
      na.omit(unlist(sigdat[i, -1]))  # drop signature name column
    )
    total_genes <- length(sig_genes)
    
    rows <- which(rownames(GEdata) %in% sig_genes)
    found_genes <- length(rows)
    
    ## Coverage
    coverage[i] <- if (total_genes > 0) {
      found_genes / total_genes * 100
    } else {
      NA_real_
    }
    
    if(found_genes > 1){
      sigscore <- apply(
        GEdata[rows, , drop = FALSE],
        2,
        median,
        na.rm = TRUE
      )
    } else if(found_genes == 1){
      sigscore <- GEdata[rows, ]
    } else {
      sigscore <- rep(NA_real_, ncol(GEdata))
    }
    
    outdat[[i]] <- sigscore
  }
  
  ## Build score data.frame
  score_df <- data.frame(
    signature = sig_names,
    do.call(rbind, outdat),
    check.names = FALSE,
    row.names = NULL
  )
  colnames(score_df)[-1] <- colnames(GEdata)
  
  ## Build coverage data.frame
  coverage_df <- data.frame(
    signature = sig_names,
    coverage = coverage,
    row.names = NULL
  )
  
  return(list(
    score = score_df,
    coverage = coverage_df
  ))
}

danaher_cell_signature <- function(signature_data, expression_data) {
  
  outdat <- NULL
  names.row <- NULL
  coverage <- numeric()
  
  unique_celltypes <- unique(signature_data$`Cell type`)
  
  for (i in seq_along(unique_celltypes)) {
    
    cell <- unique_celltypes[i]
    genes_for_cell <- signature_data %>%
      filter(`Cell type` == cell) %>%
      pull(gene_id) %>%
      as.character()
    
    total_genes <- length(genes_for_cell)
    
    # Find which genes exist in expression_data
    found_rows <- genes_for_cell[genes_for_cell %in% rownames(expression_data)]
    n_found <- length(found_rows)
    
    # Coverage calculation
    coverage[i] <- if (total_genes > 0) n_found / total_genes * 100 else NA
    
    if (n_found == 0) {
      warning(paste("Cell Type:", cell, "has no matching genes in expression data"))
      next
    }
    
    # Calculate signature score
    if (n_found > 1) {
      sigscore <- apply(expression_data[found_rows, , drop = FALSE], 2, median, na.rm = TRUE)
    } else {
      sigscore <- as.numeric(expression_data[found_rows, ])
    }
    
    outdat <- rbind(outdat, sigscore)
    names.row <- rbind(names.row, cell)
    
  }
  
  # Final output
  rownames(outdat) <- names.row
  colnames(outdat) <- colnames(expression_data)
  
  coverage_df <- data.frame(
    cell_type = unique_celltypes,
    coverage = coverage
  )
  
  return(list(
    score = outdat,
    coverage = coverage_df
  ))
}


normalize_duplicates <- function(
    df,
    gene_id_col = "Gene ID",
    drop_cols = "gene_id"
) {
  stopifnot(gene_id_col %in% colnames(df))
  ## Drop unwanted columns if present
  drop_cols <- intersect(drop_cols, colnames(df))
  df <- df %>% dplyr::filter(!is.na(.data[[gene_id_col]]) & .data[[gene_id_col]] != "")
  df <- df %>% dplyr::select(-all_of(drop_cols))
  ## Identify duplicated genes
  dup_genes <- df[[gene_id_col]][duplicated(df[[gene_id_col]])]
  ## Average duplicated genes
  averaged <- df %>%
    dplyr::filter(.data[[gene_id_col]] %in% dup_genes) %>%
    dplyr::group_by(.data[[gene_id_col]]) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), mean, na.rm = TRUE),
      .groups = "drop")
  
  ## Keep non-duplicated genes
  no_dups <- df %>%
    dplyr::filter(!.data[[gene_id_col]] %in% dup_genes)
  
  ## Combine
  out <- dplyr::bind_rows(no_dups, averaged)
  ## Set rownames and drop Gene ID column
  rownames(out) <- out[[gene_id_col]]
  out[[gene_id_col]] <- NULL
  return(out)
}

library(data.table)

clean_iges <- function(dt, start_col = 2) {
  # Ensure it's a data.table
  setDT(dt)
  
  # Columns to clean
  cols <- names(dt)[start_col:ncol(dt)]
  
  # Step 1: Replace row-wise duplicates with NA
  dt[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]  # ensure numeric
  dt[, (cols) := {
    # Convert subset of columns to matrix
    mat <- as.matrix(.SD)
    # For each row, replace duplicates with NA
    mat <- t(apply(mat, 1, function(x) {
      x[duplicated(x)] <- NA
      x
    }))
    # Convert back to data.frame list for assignment
    as.data.frame(mat)
  }, .SDcols = cols]
  
  # Step 2: Replace all 0s with NA
  dt[, (cols) := lapply(.SD, function(col) { col[col == 0] <- NA; col }), .SDcols = cols]
  
  return(dt)
}


