#!/usr/bin/env Rscript
# run_pipeline.R — entry point called by Snakemake shell rule
# All parameters are passed as environment variables.

script_dir    <- Sys.getenv("SCRIPT_DIR",    "/pipeline/R_scripts")
root_path     <- Sys.getenv("ROOT_PATH",     "/output/")
project_name  <- Sys.getenv("PROJECT_NAME",  "my_project")
group_col     <- Sys.getenv("GROUP_COL",     "Condition")
tx2gene_path  <- Sys.getenv("TX2GENE_PATH")
metadata_path <- Sys.getenv("METADATA_PATH")
comb_set_str  <- Sys.getenv("COMB_SET")       # "A,B;C,D"
method        <- Sys.getenv("METHOD",         "padj")
fc            <- as.numeric(Sys.getenv("FC",  "1"))
pval          <- as.numeric(Sys.getenv("PVAL","0.05"))
organism      <- Sys.getenv("ORGANISM",       "mouse")
keytype       <- Sys.getenv("KEYTYPE",        "SYMBOL")
norm_method   <- Sys.getenv("NORM_METHOD",    "edgeR")
specific_genes_str <- Sys.getenv("SPECIFIC_GENES", "")
hallmark_gmt  <- Sys.getenv("HALLMARK_GMT",   "")

# Source all R modules
source(file.path(script_dir, "00_libraries.R"))
source(file.path(script_dir, "01_deg_based_enrichment_analysis.R"))
source(file.path(script_dir, "02_gsea.R"))
source(file.path(script_dir, "03_visualization.R"))
source(file.path(script_dir, "04_pipeline.R"))

# Parse comb_set: "Mus_sham,Mus_UUOday2;GroupA,GroupB" -> list(c(...), c(...))
comb_set <- lapply(
    strsplit(comb_set_str, ";")[[1]],
    function(x) strsplit(trimws(x), ",")[[1]]
)

# Parse optional args
hallmark_gmt  <- if (nchar(trimws(hallmark_gmt))  == 0) NULL else hallmark_gmt
specific_genes <- if (nchar(trimws(specific_genes_str)) == 0) NULL else strsplit(trimws(specific_genes_str), ",")[[1]]

message("=== Bulk RNA-seq DE Pipeline ===")
message("Project     : ", project_name)
message("Group col   : ", group_col)
message("Comparisons : ", comb_set_str)
message("Organism    : ", organism)
message("Hallmark GMT: ", if (is.null(hallmark_gmt)) "none" else hallmark_gmt)

project_settings <- list(
    Project_name  = paste0(project_name, "/"),
    group_col     = group_col,
    comb_set      = comb_set,
    tx2gene_path  = tx2gene_path,
    metadata_path = metadata_path
)

pipeline_result <- run_deg_pipeline(
    project_settings = project_settings,
    root_path        = root_path
)

result_path <- pipeline_result$result_path
DGE_result  <- pipeline_result$DGE_result

for (comparison_name in names(DGE_result)) {
    analysis_data(
        path                 = result_path,
        Data                 = DGE_result[[comparison_name]],
        method               = method,
        FC                   = fc,
        pval                 = pval,
        Comparison_name      = comparison_name,
        Sample               = organism,
        specific_geneset     = NULL,
        specific_genes       = specific_genes,
        hallmark_geneset_gmt = hallmark_gmt,
        Keytype              = keytype,
        norm_method          = norm_method
    )
}

message("=== Pipeline completed successfully ===")