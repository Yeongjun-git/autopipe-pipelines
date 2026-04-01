# Create a folder at the given path if it does not exist
make_folder = function(path,name){
    folder_name = paste0(path,name)
    if (!file.exists(folder_name)) {
        dir.create(folder_name)
    }
}

# Convert an edgeR result table to a DESeq2-style column structure
TMM_to_DESeq2_format = function(TMM_DGE_df){
    colnames(TMM_DGE_df)[which(colnames(TMM_DGE_df) == "logFC")] <- "log2FoldChange"
    colnames(TMM_DGE_df)[which(colnames(TMM_DGE_df) == "PValue")] <- "pvalue"
    colnames(TMM_DGE_df)[which(colnames(TMM_DGE_df) == "FDR")] <- "padj"
    colnames(TMM_DGE_df)[which(colnames(TMM_DGE_df) == "gene_id")] <- "GENE_ID"    
    colnames(TMM_DGE_df)[which(colnames(TMM_DGE_df) == "gene_name")] <- "gene_id"
    TMM_DGE_df <- TMM_DGE_df[, c("GENE_ID", "gene_id", "logCPM", "log2FoldChange", "F", "pvalue", "padj")]
    return(TMM_DGE_df)
}

# Standardize DEG results, save to Excel, and add volcano plot columns
annotate_data = function(path, data, method, FC, pval, Comparison_name){
    result_data = data %>% filter(!is.na(padj))

    dge_description_df <- data.frame(
    Column = c("logCPM", "log2FoldChange", "F", "pvalue", "padj", "gene_id", "gene_name"),
    Meaning = c(
        "Log Normalized average expression value across all samples in two groups",
        "Log2-transformed fold change between two groups",
        "F-statistic from quasi-likelihood test from EdgeR",
        "Raw enrichment p-value with quasi-likelihood F-tests",
        "Adjusted p-value (false discovery rate (FDR))",
        "Ensembl gene id",
        "Gene Symbol"
    )
    )

    make_folder(path,"0.DGE/")

    colnames(result_data)[which(colnames(result_data) == "gene_id")] <- "gene_name"
    out_list <- list(
        "Result" = result_data,
        "Column_Description" = dge_description_df
    )
    write_xlsx(out_list, paste0(path,"0.DGE/DEG_edgeR_",Comparison_name,".xlsx"))

    colnames(result_data)[which(colnames(result_data) == "gene_name")] <- "gene_id"
    result_data = result_data[, c("gene_id", "logCPM", "log2FoldChange", "F", "pvalue", "padj")]
    result_data$diffexpressed = "NO"
    result_data$diffexpressed[(result_data$log2FoldChange >= FC) & (result_data[[method]] <= pval)] <- "UP"
    result_data$diffexpressed[(result_data$log2FoldChange <= -FC) & (result_data[[method]] <= pval)] <- "DOWN"
    result_data$delabel <- NA
    result_data$delabel[result_data$diffexpressed != "NO"] <- result_data$gene_id[result_data$diffexpressed != "NO"]
    return(result_data)
}

# Run the full downstream analysis for an individual comparison result
analysis_data = function(path, Data, method, FC, pval, Comparison_name, Sample, Keytype, norm_method="edgeR", specific_geneset = NULL, specific_genes = NULL, hallmark_geneset_gmt = NULL){
    
    if (norm_method == "edgeR"){
        dge_data = TMM_to_DESeq2_format(Data)
    }
    else if( norm_method == "DESeq2"){
        dge_data = Data
    }

    data = annotate_data(path = path,data = dge_data, method = method,FC = FC,pval = pval,Comparison_name = Comparison_name)

    draw_volcanoplot(path,data,method,FC,pval,Comparison_name)
    draw_volcanoplot_annot(path,data,method,FC,pval,Comparison_name)
    draw_volcanoplot_genes(path,data,method,FC,pval,Comparison_name, specific_genes)

    GO_DEG_analysis(path = path, data = data, Organism = Sample, Comparison_name = Comparison_name, Keytype = Keytype, Custom_geneset_gmt = specific_geneset, Hallmark_geneset_gmt = hallmark_geneset_gmt)
    GSEA_analysis(path = path, data = data, Organism = Sample, Comparison_name = Comparison_name,Keytype = Keytype, Hallmark_geneset_gmt = hallmark_geneset_gmt)

    run_specific_gsea_analysis(
        path = path,
        data = data,
        Comparison_name = Comparison_name,
        specific_geneset = specific_geneset,
        pvalue_cutoff = 1
    )
}

# Run the full pipeline from input preparation to edgeR DEG calculation
run_deg_pipeline <- function(project_settings, root_path) {
  prepare_inputs <- function(project_settings, root_path) {
    Project_name <- project_settings$Project_name
    comb_set <- project_settings$comb_set
    tx2gene_path <- project_settings$tx2gene_path
    metadata_path <- project_settings$metadata_path
    group_col <- project_settings$group_col

    if (is.null(group_col) || identical(group_col, "")) {
      stop("Please set project_settings$group_col.")
    }

    message("▶ Analyzing: ", Project_name)
    message("▶ Group column: ", group_col)
    make_folder(root_path, Project_name)
    result_path <- paste0(root_path, Project_name)

    tx2gene_data <- as.data.frame(
      read.csv(tx2gene_path, sep = "\t", col.names = c("tx_id", "gene_id", "gene_name"))
    )
    col_data <- as.data.frame(read.csv(metadata_path))

    if (!(group_col %in% colnames(col_data))) {
      stop("group_col is not present in metadata: ", group_col)
    }

    files <- setNames(col_data$file, col_data$sample)
    stopifnot(all(file.exists(files)))

    txi_salmon <- tximport(files, type = "salmon", tx2gene = tx2gene_data)
    col_data <- col_data[match(colnames(txi_salmon$counts), col_data$sample), ]

    tx2gene_unique <- tx2gene_data %>%
      dplyr::select(gene_id, gene_name) %>%
      dplyr::distinct(gene_id, .keep_all = TRUE)

    list(
      result_path = result_path,
      comb_set = comb_set,
      txi_salmon = txi_salmon,
      col_data = col_data,
      tx2gene_unique = tx2gene_unique,
      group_col = group_col
    )
  }

  save_tpm_tables <- function(txi_salmon, tx2gene_unique, result_path) {
    tpm_gene <- txi_salmon$abundance %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      left_join(tx2gene_unique, by = "gene_id") %>%
      relocate(gene_id, .before = 1) %>%
      relocate(gene_name, .after = gene_id)

    write.csv(
      tpm_gene,
      file = paste0(result_path, "TPM_gene.csv"),
      row.names = FALSE
    )

    tpm_gene_log <- tpm_gene %>%
      mutate(across(where(is.numeric), ~ log2(.x + 1)))

    write.csv(
      tpm_gene_log,
      file = paste0(result_path, "TPM_gene_log2p1.csv"),
      row.names = FALSE
    )

    invisible(NULL)
  }

  build_edger_object <- function(txi_salmon, col_data, group_col) {
    cts <- txi_salmon$counts
    normMat <- txi_salmon$length

    normMat <- normMat / exp(rowMeans(log(normMat)))
    normCts <- cts / normMat
    eff.lib <- calcNormFactors(normCts) * colSums(normCts)

    normMat <- sweep(normMat, 2, eff.lib, "*")
    normMat <- log(normMat)

    y <- DGEList(cts)
    y <- scaleOffset(y, normMat)

    group_factor <- as.factor(col_data[[group_col]])
    design <- model.matrix(~ 0 + group_factor)
    colnames(design) <- levels(group_factor)

    list(y = y, design = design)
  }

  save_tmm_tables <- function(y, tx2gene_unique, result_path) {
    tmm <- edgeR::cpm(y, offset = y$offset, log = FALSE) %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      left_join(tx2gene_unique, by = "gene_id") %>%
      relocate(gene_id, .before = 1) %>%
      relocate(gene_name, .after = gene_id)

    write.csv(
      tmm,
      file = paste0(result_path, "TMM_length_normalized_count.csv"),
      row.names = FALSE
    )

    log_tmm <- edgeR::cpm(y, offset = y$offset, log = TRUE) %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      left_join(tx2gene_unique, by = "gene_id") %>%
      relocate(gene_id, .before = 1) %>%
      relocate(gene_name, .after = gene_id)

    write.csv(
      log_tmm,
      file = paste0(result_path, "log_TMM_length_normalized_count.csv"),
      row.names = FALSE
    )

    list(tmm = tmm, log_tmm = log_tmm)
  }

  run_deg_contrasts <- function(y, design, comb_set, tx2gene_unique) {
    dge_all <- estimateDisp(y, design)
    fit_all <- glmQLFit(dge_all, design)

    deg_result <- list()
    for (pair in comb_set) {
      treat <- pair[1]
      con <- pair[2]
      key_name <- paste0(treat, "_vs_", con)

      contrast_df <- makeContrasts(
        contrasts = paste0(treat, "-", con),
        levels = design
      )

      qlf <- glmQLFTest(fit_all, contrast = contrast_df)
      result_df <- as.data.frame(topTags(qlf, n = Inf))
      result_df$gene_id <- rownames(result_df)
      result_df <- result_df[, c("gene_id", setdiff(colnames(result_df), "gene_id"))]
      result_df <- result_df %>% left_join(tx2gene_unique, by = "gene_id")

      deg_result[[key_name]] <- result_df
    }
    return(deg_result)
  }

  prep <- prepare_inputs(project_settings, root_path)
  save_tpm_tables(prep$txi_salmon, prep$tx2gene_unique, prep$result_path)

  edge <- build_edger_object(prep$txi_salmon, prep$col_data, prep$group_col)
  tmm_tables <- save_tmm_tables(edge$y, prep$tx2gene_unique, prep$result_path)

  keep <- filterByExpr(edge$y, edge$design)
  y <- edge$y[keep, ]

  tmm_filtered <- edgeR::cpm(y, offset = y$offset, log = FALSE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    left_join(prep$tx2gene_unique, by = "gene_id") %>%
    relocate(gene_id, .before = 1) %>%
    relocate(gene_name, .after = gene_id)

  write.csv(
    tmm_filtered,
    file = paste0(prep$result_path, "TMM_length_normalized_count_filtered.csv"),
    row.names = FALSE
  )

  TMM_length_normalized_count_log <- edgeR::cpm(y, offset = y$offset, log = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    left_join(prep$tx2gene_unique, by = "gene_id") %>%
    relocate(gene_id, .before = 1) %>%
    relocate(gene_name, .after = gene_id)

  plot_pca_with_condition_and_batch(
    expr_matrix = TMM_length_normalized_count_log,
    col_data = prep$col_data,
    result_path = prep$result_path,
    condition_col_2 = prep$group_col,
    gene_annotation = prep$tx2gene_unique
  )

  dge_result <- run_deg_contrasts(y, edge$design, prep$comb_set, prep$tx2gene_unique)

  list(
    Project_name = project_settings$Project_name,
    group_col = prep$group_col,
    comb_set = prep$comb_set,
    result_path = prep$result_path,
    Col_data = prep$col_data,
    tx2gene_unique = prep$tx2gene_unique,
    y = y,
    DGE_result = dge_result,
    TMM_length_normalized_count = tmm_tables$tmm,
    log_TMM_length_normalized_count = tmm_tables$log_tmm,
    TMM_length_normalized_count_log = TMM_length_normalized_count_log
  )
}