# Bulk RNA-seq Differential Expression & Enrichment Analysis Pipeline

## Overview
This pipeline performs end-to-end Bulk RNA-seq differential expression (DEG) analysis and downstream functional enrichment on **Salmon-quantified** data.

**Key steps:**
1. Import transcript-level quantifications via `tximport`
2. Build a length-aware edgeR DGE object (TMM normalisation)
3. Export TPM and TMM count matrices
4. Run quasi-likelihood F-tests for each pairwise comparison
5. For each comparison: produce DEG Excel tables, volcano plots, GO ORA, GSEA, Reactome, and optional Hallmark analyses

---

## Required Inputs

| File | Description |
|------|-------------|
| `tx2gene.tsv` | Transcript → gene mapping (3 cols: `tx_id`, `gene_id`, `gene_name`) |
| `metadata.csv` | Sample sheet with `sample`, `file` (path to `quant.sf`), and group column |
| `quant.sf` files | Salmon per-sample quantification outputs |
| GMT file *(optional)* | MSigDB Hallmark gene sets (mouse or human) |

---

## Expected Outputs

```
/output/<project_name>/
├── TPM_gene.csv
├── TPM_gene_log2p1.csv
├── TMM_length_normalized_count.csv
├── TMM_length_normalized_count_filtered.csv
├── log_TMM_length_normalized_count.csv
├── PCA_plot.png
├── 0.DGE/
│   └── DEG_edgeR_<GroupA>_vs_<GroupB>.xlsx
├── 1.GO/            ← GO ORA results & plots
├── 2.GSEA/          ← GSEA results & plots
├── 3.Reactome/      ← Reactome pathway results & plots
└── 4.Hallmark/      ← Hallmark GSEA results (if GMT provided)
```

---

## Configuration (`config.yaml`)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `project_name` | Output subdirectory name | `my_project` |
| `tx2gene_path` | Path to tx2gene TSV inside container | `/input/tx2gene.tsv` |
| `metadata_path` | Path to metadata CSV inside container | `/input/metadata.csv` |
| `group_col` | Metadata column defining groups | `Condition` |
| `comb_set` | List of `[treatment, control]` pairs | — |
| `method` | Significance column: `padj` or `pvalue` | `padj` |
| `fc_threshold` | Minimum \|log2FC\| to call DE | `1` |
| `pval_threshold` | Maximum adjusted p-value | `0.05` |
| `sample_organism` | `mouse` or `human` | `mouse` |
| `keytype` | Gene ID type: `SYMBOL`, `ENSEMBL`, `ENTREZID` | `SYMBOL` |
| `specific_genes` | Comma-separated genes to highlight on volcano | `""` |
| `hallmark_geneset_gmt` | Path to GMT file inside container (optional) | `""` |
| `threads` | Number of CPU threads | `4` |

---

## How to Run

### 1. Build the Docker image
```bash
docker build -t autopipe-bulk-rnaseq-de .
```

### 2. Run the pipeline
```bash
docker run --rm \
  -v /path/to/your/input:/input:ro \
  -v /path/to/your/output:/output \
  autopipe-bulk-rnaseq-de \
  snakemake --cores 4
```

### 3. Multi-comparison example (`config.yaml`)
```yaml
comb_set:
  - ["KO", "WT"]
  - ["KO_treated", "KO"]
  - ["KO_treated", "WT"]
```

---

## Notes
- The metadata CSV `file` column must contain **absolute paths** to `quant.sf` as they appear **inside** the container (i.e. under `/input/...`).
- Only **edgeR** normalisation is currently supported (`norm_method: "edgeR"`).
- For mouse Hallmark gene sets, download `mh.all.v*.Mm.symbols.gmt` from MSigDB and place it in your input directory.