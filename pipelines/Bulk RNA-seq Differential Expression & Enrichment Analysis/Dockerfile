FROM condaforge/miniforge3:latest

# ── Create the RNA-seq conda environment ──────────────────────────────────────
RUN conda create -n rnaseq_env r-base=4.5.1 -c conda-forge -y

# ── Install Bioconductor + CRAN packages into the environment ─────────────────
RUN conda run -n rnaseq_env conda install -y -c bioconda -c conda-forge \
    bioconductor-tximport \
    bioconductor-edger \
    bioconductor-clusterprofiler \
    bioconductor-fgsea \
    bioconductor-reactomepa \
    bioconductor-org.mm.eg.db \
    bioconductor-org.hs.eg.db \
    bioconductor-enrichplot \
    r-ggplot2 \
    r-ggrepel \
    r-dplyr \
    r-tibble \
    r-writexl \
    r-stringr \
    r-rcolorbrewer \
    && conda clean -afy

# ── AutoPipe runtime dependencies (base env) ──────────────────────────────────
RUN conda install -y -c conda-forge -c bioconda \
    snakemake-minimal \
    bash \
    && conda clean -afy

# Replace system bash with conda bash (prevents GLIBC mismatch)
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/bin/bash /bin/sh

# ── Copy pipeline files ───────────────────────────────────────────────────────
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY run_pipeline.R .
COPY R_scripts/ R_scripts/

CMD ["snakemake", "--cores", "all"]