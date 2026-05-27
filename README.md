<div align="center">

# 🧬 AMR Analyser

**A lightning-fast pseudoalignment-based pipeline for rapid detection and quantification of Antimicrobial Resistance (AMR) genes from RNA-seq data**

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?style=for-the-badge&logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-orange?style=for-the-badge)](https://github.com/yourusername/amr-analyser)


</div>

## ✨ Overview

**AMR Analyser** is a novel computational framework that bridges the gap between genomic potential and **functional gene expression** in antimicrobial resistance studies.

Traditional alignment-based RNA-seq pipelines are slow, memory-hungry, and impractical for real-time surveillance. AMR Analyser uses **pseudoalignment** powered by **k-mer indexing** and **Transcriptome de Bruijn Graph (T-DBG)**, combined with **Expectation-Maximization (EM)** probabilistic quantification — delivering **near-optimal accuracy at unprecedented speed**.

> **Process ~1.2 million reads per minute** on a standard laptop (4–8 GB RAM) — **orders of magnitude faster** than conventional methods.

---

## 🚀 Key Features

- **Ultra-Fast Pseudoalignment** using k-mer (k=8) compatibility mapping
- **Accurate Quantification** via EM algorithm (outputs TPM, RPKM, and raw counts)
- **Lightweight & Scalable** — No massive BAM files, runs efficiently on personal workstations and HPC
- **Curated Resistome Database** built from CARD + NCBI (deduplicated >95% identity)
- **Species-Agnostic** — Works on bacterial transcriptomes and host-pathogen mixed samples
- **Reproducible & User-Friendly** CLI with detailed visualizations

---

## 📊 Performance Highlights

| Metric                    | Traditional Aligners       | **AMR Analyser**              |
|--------------------------|---------------------------|-------------------------------|
| Speed                    | Hours per sample          | **1.2M reads/minute**         |
| Memory Usage             | High + large BAM files    | **Low (4-8 GB)**              |
| Index Construction       | Repeated per analysis     | **One-time reusable T-DBG**   |
| Accuracy                 | High                      | **Near-optimal**              |

**Tested successfully** on *E. coli*, *M. tuberculosis*, *A. baumannii*, *S. granuli*, *C. necator*, and human neutrophil + *S. aureus* infection models.

---

## 📁 Pipeline Architecture

```mermaid
graph LR
    A[FASTQ Reads] --> B[k-mer Decomposition]
    B --> C[T-DBG Index Query]
    C --> D[Equivalence Class Assignment]
    D --> E[EM Quantification]
    E --> F[TPM / RPKM / Counts]
    E --> G[Visualizations]

