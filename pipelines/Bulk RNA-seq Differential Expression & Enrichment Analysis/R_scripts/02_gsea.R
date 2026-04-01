# Convert a SYMBOL-named ranked vector to an ENTREZID-based vector
convert_vector_to_entrez <- function(genes, orgdb, keytype) {
  # Extract names
  original_ids <- names(genes)

  # Convert IDs
  conversion <- bitr(
    original_ids,
    fromType = keytype,
    toType = "ENTREZID",
    OrgDb = orgdb
  )

  # Keep only matched ENTREZIDs
  matched <- conversion[conversion[[keytype]] %in% original_ids, ]
  
  # Map log2FC values
  matched$log2FC <- genes[matched[[keytype]]]

  # Build named vector with ENTREZID as names
  entrez_vector <- setNames(matched$log2FC, matched$ENTREZID)

  # Remove duplicates (keep first occurrence, which has highest abs value after sorting)
  entrez_vector <- entrez_vector[!duplicated(names(entrez_vector))]

  # Sort in descending order
  entrez_vector <- sort(entrez_vector, decreasing = TRUE)

  return(entrez_vector)
}

# Create an empty data frame when no GSEA result is available
empty_gsea_df <- function() {
  data.frame(
    ID = character(),
    Description = character(),
    setSize = character(),
    enrichmentScore = character(),
    NES = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    rank = character(),
    leading_edge = character(),
    core_enrichment = character(),
    stringsAsFactors = FALSE
  )
}

# Create a data frame describing the meaning of each GSEA output column
description_gsea_df <- function() {
    data.frame(
    Column = c(
        "Reference",
        "ID",
        "Description",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust",
        "qvalue",
        "rank",
        "leading_edge",
        "core_enrichment"
    ),
    Meaning = c(
        "https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html",
        "Term ID",
        "Term name",
        "Number of genes in the term",
        "Raw enrichment score (ES); reflects how concentrated gene ranks are at the top/bottom of the ranked list",
        "Normalized Enrichment Score (ES normalized by set size and permutations)",
        "The p-value of the ES is calculated using a permutation test",
        "Adjusted p-value (Benjamini-Hochberg method)",
        "Adjusted p-value (Estimated false discovery rate (FDR))",
        "Rank in the ordered gene list where the enrichment score peak occurs",
        paste0(
        "Summary of the leading edge subset contributing most to the enrichment score:\n",
        " - tags: % of gene hits before (for +ES) or after (for -ES) the peak\n",
        " - list: % of ranked genes before (for +ES) or after (for -ES) the peak\n",
        " - signal: combined enrichment strength"
        ),
        "The actual genes contributing to the enrichment signal (leading edge core genes)"
        )
    )
}

## Plot size adjustment functions

# Filter GSEA results by p-value and gene set size
filter_gsea_result <- function(df, pvalue_cutoff = 0.05, min_set_size = 10, max_set_size = 500) {
  result <- df %>%
    arrange(p.adjust) %>%
    filter(p.adjust < pvalue_cutoff, setSize >= min_set_size, setSize <= max_set_size)
  return(result)
}

# Dynamically adjust plot size based on number of items and label length
adjust_plot_size <- function(top_df) {
  # Base size
  width <- 8
  height <- 8
  
  # Adjust width based on label length
  max_label_length <- max(nchar(    as.character(top_df$Description)  ))
  counts = nrow(top_df)

  adjusted_width <- width + max_label_length * 0.1
  adjusted_height <- height - (20 - counts) * 0.3
  
  # Return adjusted width and height
  return(list(width = adjusted_width, height = adjusted_height))
}

# Select top pathways from GSEA results and reorder for visualization
select_top_GSEA_abs <- function(GSEA_result, top_n = 30) {
    Len = nrow(GSEA_result)
    GSEA_result %>%
        data.frame() %>%
        arrange(p.adjust) %>%
        head(n = min(top_n,Len)) %>%
        arrange(NES) -> top_pathway_data
    top_pathway_data$description = paste0(top_pathway_data$Description," (setsize: ",top_pathway_data$setSize,") ")
    factor(top_pathway_data$description, levels = top_pathway_data$description) -> top_pathway_data$description
    return(top_pathway_data)
}

# Save a summary GSEA bar plot as PDF
draw_gsea_summary_plot <- function(GSEA_result, prefix_suffix, path, Comparison_name, top_n = 30, pvalue_cutoff = 0.05, min_set_size = 10, max_set_size = 500) {
        set.seed(116)
        sort_df <- GSEA_result %>%
                arrange(p.adjust)
        sort_df = filter_gsea_result(sort_df, pvalue_cutoff = pvalue_cutoff, min_set_size = min_set_size, max_set_size = max_set_size)
        top_df = select_top_GSEA_abs(sort_df, top_n = top_n)

        # pdf file format
        plot_size <- adjust_plot_size(top_df)
        pdf(file= paste0(path,"3.GSEA/1.figure/",Comparison_name,"_",prefix_suffix,".pdf"), width=plot_size$width, height=plot_size$height)
        
        #draw figure
        p = ggplot(top_df ) +
        geom_bar(aes(x = description, y = NES, fill = p.adjust),
                stat='identity') +
        scale_fill_continuous(low = "red", high = "blue") +
        coord_flip(ylim = c(-3.2, 3.2)) +
        scale_y_continuous(breaks = seq(-3,3,1)) +
        theme_bw(14) +
        labs(title=paste0(Comparison_name,"_",prefix_suffix),
                x = "GO Term") +  
        theme(plot.title = element_text(hjust = 0.5),
                axis.title.y = element_blank())

        print(p)
        dev.off()
}

### GSEA functions

# Run preranked GSEA for GO BP/CC/MF
perform_go_gsea <- function(genes, keytype, orgdb, Comparison_name,path) {
  ontologies <- c("BP", "CC", "MF")
  go_results <- list()

  for (ont in ontologies) {
    ego <- gseGO(
      gene          = genes,
      OrgDb         = orgdb,
      ont           = ont,
      keyType       = keytype,
      pAdjustMethod = "BH"
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      result <- as.data.frame(ego)
      result_sorted <- result[order(result$p.adjust), ]
    } else {
      result_sorted <- empty_gsea_df()
    }
    go_results[[ont]] <- result_sorted
  }

  return(go_results)
}

# Run KEGG GSEA from a ranked gene vector
perform_kegg_gsea <- function(genes, keytype, orgdb, Org_code) {
  conv <- convert_vector_to_entrez(genes, orgdb, keytype)

  ekegg <- gseKEGG(
    gene = conv,
    organism = Org_code,
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH"
  )

  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg_readable <- setReadable(ekegg, OrgDb = orgdb, keyType = "ENTREZID")
    result <- as.data.frame(ekegg_readable)
    result <- result[order(result$p.adjust), ]
  } else {
    result <- empty_gsea_df()
  }
  return(result)
}

# Run Reactome GSEA from a ranked gene vector
perform_reactome_gsea <- function(genes, keytype, orgdb, organism) {
  conv <- convert_vector_to_entrez(genes, orgdb, keytype)

  erct <- gsePathway(
    gene          = conv,
    organism      = organism,
    pAdjustMethod = "BH"
  )
  if (!is.null(erct) && nrow(as.data.frame(erct)) > 0) {
    ekegg_readable <- setReadable(erct, OrgDb = orgdb, keyType = "ENTREZID")
    result <- as.data.frame(ekegg_readable)
    result_sorted <- result[order(result$p.adjust), ]
  } else {
    result_sorted <- empty_gsea_df()
  }

  return(result_sorted)
}

# Run preranked GSEA for a user-defined GMT gene set
perform_custom_geneset_gsea <- function(genes, Custom_geneset_gmt) {

  if (length(genes) == 0) return(empty_gsea_df())

  Custom_df <- read.gmt(Custom_geneset_gmt)

  gsea_res <- GSEA(
    geneList = genes,
    TERM2GENE = Custom_df,
    pAdjustMethod = "BH"
  )

  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    gsea_df <- as.data.frame(gsea_res)
    return(gsea_df[order(gsea_df$p.adjust), ])
  } else {
    return(empty_gsea_df())
  }
}

# Run GSEA for a specific GMT gene set and save individual plots
run_specific_gsea_analysis <- function(path, data, Comparison_name, specific_geneset, pvalue_cutoff = 1) {
    try({
            if (is.null(specific_geneset)) {
        return(invisible(NULL))
        }

        Genes <- setNames(data$log2FoldChange, data$gene_id)
        # Remove duplicate gene IDs (keep first = highest abs log2FC after sort)
        Genes <- Genes[!duplicated(names(Genes))]
        Genes <- sort(Genes, decreasing = TRUE)

        hallmark_df <- read.gmt(specific_geneset)

        gsea_res <- GSEA(
            geneList = Genes,
            TERM2GENE = hallmark_df,
            pAdjustMethod = "BH",
            pvalueCutoff = pvalue_cutoff
        )

        draw_Specific_GSEA_plot(path = path, data = gsea_res, Name = Comparison_name)
    })

}

# Run full GSEA analysis and save Excel results plus a representative figure
GSEA_analysis = function(path,data,Organism,Comparison_name,Keytype,Top_n = 20, Custom_geneset_gmt = NULL, Hallmark_geneset_gmt = NULL, prefix_suffix = "GO_BP") {
    make_folder(path,"3.GSEA/")
    make_folder(path,"3.GSEA/0.excel/")
    make_folder(path,"3.GSEA/1.figure/")
    #ready for data
    if (Organism == "human"){
        OrgDb = org.Hs.eg.db
        Org_code = "hsa"
    } else if (Organism == "mouse") {
        OrgDb = org.Mm.eg.db
        Org_code = "mmu"
    }

  run_hallmark <- !is.null(Hallmark_geneset_gmt)
  if (run_hallmark && !file.exists(Hallmark_geneset_gmt)) {
    stop("Please check Hallmark_geneset_gmt path: ", Hallmark_geneset_gmt)
  }
    
    all_results = list()

    Genes <- setNames(data$log2FoldChange, data$gene_id)
    # Remove duplicate gene IDs (keep first = highest abs log2FC after sort)
    Genes <- Genes[!duplicated(names(Genes))]
    Genes <- sort(Genes, decreasing = TRUE)
    prefix <- Comparison_name

    # Run each enrichment analysis
    go_res <- perform_go_gsea(Genes, Keytype, OrgDb)
    kegg_res <- perform_kegg_gsea(Genes, Keytype, OrgDb, Org_code)
    reactome_res <- perform_reactome_gsea(Genes, Keytype, OrgDb, Organism)
    if (run_hallmark) {
      hallmark_res <- perform_custom_geneset_gsea(Genes, Hallmark_geneset_gmt)
    }


    for (ont in names(go_res)) {
        GO_name = paste0("GO_", ont)
        all_results[[GO_name]] <- go_res[[ont]]
    }
    all_results[["KEGG"]] <- kegg_res
    all_results[["Reactome"]] <- reactome_res
    if (run_hallmark) {
      all_results[["Hallmark"]] <- hallmark_res
    }
    all_results[["Column_Description"]] = description_gsea_df()
    write_xlsx(all_results, paste0(path,"3.GSEA/0.excel/",prefix,".xlsx"))

    # Representative visualization: default - GO_BP
    GSEA_result = all_results[[prefix_suffix]]
    draw_gsea_summary_plot(
        GSEA_result = GSEA_result,
        prefix_suffix = prefix_suffix,
        path = path,
        Comparison_name = Comparison_name,
        top_n = 30
    )
}

# Save GSEA enrichment plots for each gene set one by one
draw_Specific_GSEA_plot = function(path, data, Name, P.adjust = 0.05){
    Genesets <- data$ID
    for (Geneset in Genesets){
        try({
            # setting details
            Geneset_description = data$Description[data$ID == Geneset]
            NES = format(data$NES[data$ID == Geneset], digits = 4)
            FDR = format(data$p.adjust[data$ID == Geneset], scientific = TRUE, digits = 4)
            pvalue = format(data$pvalue[data$ID == Geneset], scientific = TRUE, digits = 4)
            setSize = data$setSize[data$ID == Geneset]
            
            # P.adjust threshold
            if (data$p.adjust[data$ID == Geneset] < P.adjust){
                color = "red"
            } else {
                color = "black"
            }
            
            # make folder and save plot
            make_folder(path, paste0("3.GSEA/"))
            make_folder(path, paste0("3.GSEA/2.GeneSet/"))
            make_folder(path, paste0("3.GSEA/2.GeneSet/", Geneset_description))
            pdf(file = paste0(path, "3.GSEA/2.GeneSet/", Geneset_description, "/", Name, ".pdf"), width = 5, height = 5)
            
            p = gseaplot2(data, geneSetID = Geneset, subplots = 1:3, title = "", ES = "line")
            
            p[[1]] = p[[1]] +
                annotate("text", x = Inf, y = Inf, label = paste0("p-value: ", pvalue), hjust = 1, vjust = 2, size = 3) + 
                annotate("text", x = Inf, y = Inf, label = paste0("FDR: ", FDR), hjust = 1, vjust = 4, size = 3, color = color) +
                annotate("text", x = Inf, y = Inf, label = paste0("NES: ", NES), hjust = 1, vjust = 6, size = 3) +
                labs(title = Name, subtitle = paste0(Geneset_description, "\n (setsize: ", setSize, ") "), y = "Running Enrichment Score") + 
                theme(plot.title = element_text(hjust = 0.5, size = 14),
                      plot.subtitle = element_text(hjust = 0.5, size = 12, color = color))
            
            print(p)
            dev.off()
        })
    }
}