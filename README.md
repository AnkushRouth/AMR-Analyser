<div align="center">

<!-- Hero Banner -->
<img src="https://capsule-render.vercel.app/api?type=waving&color=0:0a3d2e,50:0d6e56,100:1db88a&height=220&section=header&text=AMR%20Analyser&fontSize=64&fontColor=ffffff&fontAlignY=40&desc=A%20Pseudoalignment-Based%20Pipeline%20for%20Rapid%20Detection%20%26%20Quantification%20of%20AMR%20Genes&descSize=15&descAlignY=62&animation=fadeIn" width="100%"/>

<!-- Badges Row 1 -->
![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=flat-square&logo=python&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)
![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20HPC-FCC624?style=flat-square&logo=linux&logoColor=black)
![Status](https://img.shields.io/badge/Status-Active-brightgreen?style=flat-square)
![Thesis](https://img.shields.io/badge/MSc%20Dissertation-NIT%20Rourkela%202026-0d6e56?style=flat-square)

<!-- Badges Row 2 -->
![Database](https://img.shields.io/badge/Database-CARD%20%7C%20NCBI-E74C3C?style=flat-square)
![Method](https://img.shields.io/badge/Method-Pseudoalignment%20%7C%20k--mer-8E44AD?style=flat-square)
![Speed](https://img.shields.io/badge/Speed-1.2M%20reads%2Fmin-FF6B35?style=flat-square)
![Output](https://img.shields.io/badge/Output-TPM%20%7C%20RPKM%20%7C%20Counts-0077B5?style=flat-square)

<br/>

> **AMR Analyser** is a high-throughput, pseudoalignment-based computational framework for rapid detection and quantitative profiling of transcriptionally active antimicrobial resistance (AMR) genes directly from RNA-seq data — processing up to **1.2 million reads per minute** without full base-level alignment.

<br/>

**[📖 Read the Thesis](#-citation) · [⚙️ Installation](#-installation) · [🚀 Quick Start](#-quick-start) · [📊 Results](#-benchmarking--results) · [🧬 Methods](#-methodology)**

</div>

---

## 📌 Table of Contents

- [The Problem](#-the-problem)
- [Why AMR Analyser?](#-why-amr-analyser)
- [Methodology](#-methodology)
- [Pipeline Overview](#-pipeline-overview)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Parameters](#-parameters)
- [Benchmarking & Results](#-benchmarking--results)
- [Output Format](#-output-format)
- [Supported Databases](#-supported-databases)
- [Citation](#-citation)
- [Author](#-author)

---

## 🌍 The Problem

Antimicrobial resistance (AMR) is one of the most critical threats to global public health — responsible for an estimated **1.27 million deaths annually** (Lancet, 2022), with projections suggesting it could become the leading cause of death globally by 2050.

```
The Core Challenge:
┌──────────────────────────────────────────────────────────────────────────┐
│  Culture-based assays   → 24–72 hours · miss non-culturable organisms    │
│  Genomic WGS pipelines  → detect gene PRESENCE, not EXPRESSION           │
│  Alignment-based RNA-seq → computationally intensive · large BAM files   │
│                           · poor scalability in real-time settings        │
│                                                                          │
│  ✗ None of these distinguish between an EXPRESSED resistance gene        │
│    and a SILENT genetic element.                                          │
└──────────────────────────────────────────────────────────────────────────┘
```

> **AMR Analyser solves this** by quantifying *transcriptionally active* resistance determinants using pseudoalignment — orders of magnitude faster than traditional approaches, with no sacrifice in biological accuracy.

---

## ✅ Why AMR Analyser?

| Feature | Traditional Alignment (Bowtie2/HISAT2) | **AMR Analyser** |
|:---|:---:|:---:|
| Speed | ~100K–300K reads/min | **~1.2M reads/min** |
| Base-level alignment | ✅ Required | ❌ Not required |
| BAM intermediate files | ✅ Large (GBs) | ❌ None |
| Memory overhead | High | **Low** |
| Functional (expressed) AMR detection | ✅ Yes | **✅ Yes** |
| Multi-mapping resolution | Limited | **EM Algorithm** |
| Outputs | Counts only | **TPM, RPKM, Counts** |
| Real-time suitability | ❌ Poor | **✅ Excellent** |
| Reference | Full transcriptome | **Compact AMR index** |

---

## 🧬 Methodology

AMR Analyser is built on four core computational components:

### 1. 🗄️ Reference Database Construction

A curated AMR-specific resistome reference is built by integrating sequences from:
- **CARD** — Comprehensive Antibiotic Resistance Database
- **NCBI** — National Center for Biotechnology Information

Preprocessing pipeline:
```
Raw CARD/NCBI sequences
    │
    ▼
Redundancy filtering (≥95–100% identity threshold)
    │
    ▼
Length filtering (remove sub-threshold sequences)
    │
    ▼
FASTA standardization (consistent headers + metadata)
    │
    ▼
Curated, high-confidence resistome reference
```

### 2. 🔷 k-mer Indexing via Transcriptome de Bruijn Graph (T-DBG)

Each reference AMR gene sequence of length **L** is decomposed into overlapping k-mers of fixed length **k = 8**, yielding **(L - k + 1)** k-mers per sequence.

```
Gene sequence (length L):
ATGCTAGCTAGCTTAG...
│
▼ k-mer decomposition (k = 8)
[ATGCTAGC] [TGCTAGCT] [GCTAGCTT] [CTAGCTTA] [TAGCTTAG]...
│
▼ T-DBG Construction
Nodes  = unique (k-1)-mers
Edges  = observed k-mers connecting consecutive nodes
        + annotated with AMR gene membership sets
│
▼ Hash-based lookup table
O(1) retrieval of candidate genes during pseudoalignment
```

> The index is constructed **once** and reused across all RNA-seq datasets, making repeated analyses highly efficient.

### 3. ⚡ Pseudoalignment Algorithm

Unlike traditional aligners that compute base-by-base alignment coordinates, AMR Analyser determines **transcript compatibility** — which AMR genes are consistent with a given read — using k-mer set intersection:

```
For each sequencing read:

Read ──► Decompose into overlapping k-mers (k=8)
              │
              ▼
         Query precomputed k-mer index
              │
              ▼
         Retrieve gene sets per k-mer
              │
              ▼
         Intersect gene sets ──► Equivalence Class
         (genes consistent with all read k-mers)
              │
              ▼
         Pass to EM quantification
```

Conserved k-mers shared across many resistance gene families are treated as **non-informative** and down-weighted to prevent equivalence class inflation.

### 4. 📈 Expectation-Maximization (EM) Quantification

High sequence similarity among AMR gene families causes multi-mapping ambiguity. AMR Analyser resolves this probabilistically using an **EM algorithm**:

```
Initialize: uniform abundance estimates for all genes

┌─────────────────────────────────────────────────────────────┐
│  E-Step (Expectation):                                       │
│  Probabilistically assign each read within an equivalence   │
│  class to compatible genes, proportional to current         │
│  abundance estimates                                         │
│                                                             │
│  M-Step (Maximization):                                      │
│  Update gene abundance estimates to maximize the likelihood  │
│  of observed read assignments (fractional read counting)     │
└─────────────────────────────────────────────────────────────┘
         │
         ▼ Repeat until convergence
         (negligible change in abundance between iterations)
         │
         ▼
Final output: TPM · RPKM · Raw estimated counts
```

---

## 🔄 Pipeline Overview

```
                    ┌─────────────────────────────────┐
                    │       RNA-seq Input (FASTQ)      │
                    │   (Single-end or Paired-end)     │
                    │    75–150 bp · 10–50M reads      │
                    └────────────────┬────────────────┘
                                     │
                    ┌────────────────▼────────────────┐
                    │   Quality Control & Filtering    │
                    │  (Min quality threshold check)   │
                    └────────────────┬────────────────┘
                                     │
              ┌──────────────────────▼──────────────────────┐
              │         AMR-Specific Reference Index         │
              │     CARD + NCBI → k-mer decomposition        │
              │      → T-DBG construction → Hash lookup      │
              └──────────────────────┬──────────────────────┘
                                     │
                    ┌────────────────▼────────────────┐
                    │         Pseudoalignment          │
                    │  k-mer compatibility mapping     │
                    │  → Equivalence class assignment  │
                    └────────────────┬────────────────┘
                                     │
                    ┌────────────────▼────────────────┐
                    │     EM Algorithm (Iterative)     │
                    │  Probabilistic read assignment   │
                    │  → Fractional count estimation   │
                    └────────────────┬────────────────┘
                                     │
                    ┌────────────────▼────────────────┐
                    │       Quantified Output          │
                    │  TPM  ·  RPKM  ·  Raw Counts    │
                    │  Per AMR gene, per sample        │
                    └─────────────────────────────────┘
```

---

## 🛠️ Installation

### Prerequisites

```bash
# System requirements
OS      : Linux (Debian / Ubuntu recommended)
CPU     : Multi-core (4–8 cores minimum)
RAM     : 4–8 GB minimum
Python  : 3.8+
```

### Clone the Repository

```bash
git clone https://github.com/ankushrouth/AMR-Analyser.git
cd AMR-Analyser
```

### Install Dependencies

```bash
pip install -r requirements.txt
```

**Core dependencies:**

```
numpy>=1.21.0
pandas>=1.3.0
biopython>=1.79
```

### Verify Installation

```bash
python amr_analyser.py --version
# AMR Analyser v1.0.0 | NIT Rourkela 2026
```

---

## 🚀 Quick Start

### Step 1: Build the AMR Reference Index

```bash
python amr_analyser.py build-index \
    --db CARD \
    --fasta amr_reference.fasta \
    --kmer 8 \
    --output amr_index/
```

### Step 2: Run Pseudoalignment (Single-end)

```bash
python amr_analyser.py quantify \
    --index amr_index/ \
    --reads sample.fastq.gz \
    --output results/ \
    --threads 8
```

### Step 3: Run Pseudoalignment (Paired-end)

```bash
python amr_analyser.py quantify \
    --index amr_index/ \
    --reads-1 sample_R1.fastq.gz \
    --reads-2 sample_R2.fastq.gz \
    --output results/ \
    --threads 8
```

### Step 4: View Results

```bash
# Output files generated in results/
ls results/
# amr_abundance.tsv   ← TPM, RPKM, raw counts per gene
# summary_stats.txt   ← alignment rate, detection metrics
# equivalence_classes.txt
```

---

## ⚙️ Parameters

| Parameter | Default | Description |
|:---|:---:|:---|
| `--kmer` | `8` | k-mer length for index construction and mapping |
| `--em-iter` | `200` | Maximum EM algorithm iterations |
| `--min-quality` | `20` | Minimum Phred quality score for read filtering |
| `--threads` | `4` | Number of parallel processing threads |
| `--bootstrap` | `0` | Number of bootstrap samples for uncertainty estimation |
| `--output-format` | `tsv` | Output format: `tsv`, `csv`, or `json` |

> **Note on k-mer length:** Shorter k-mers improve tolerance to sequencing errors; longer k-mers increase specificity. k=8 balances sensitivity and specificity for the AMR gene reference space.

---

## 📊 Benchmarking & Results

AMR Analyser was benchmarked across five diverse bacterial transcriptomic datasets:

### Species-Level Performance Summary

| Parameter | *E. coli* | *S. granuli* | *M. tuberculosis* | *C. necator* | *A. baumannii* |
|:---|---:|---:|---:|---:|---:|
| **Total Reads** | 6,385,428 | 8,757,751 | 12,976,451 | 7,936,682 | 21,453,150 |
| **Aligned Reads** | 5,938 | 49,184 | 1,825 | 88,504 | 46,763 |
| **Alignment Rate** | 0.09% | 0.56% | 0.01% | 1.12% | 0.22% |
| **Detected AMR Genes** | 14 | 423 | 273 | 166 | **964** |
| **Detection Rate** | 0.23% | 6.99% | 4.51% | 2.74% | **15.93%** |
| **Mean TPM** | 71,428.57 | 2,364.07 | 3,663.00 | 6,024.10 | 1,037.34 |
| **Total Est. Counts** | 5,938.00 | 49,184.00 | 1,825.00 | 88,504.00 | 46,763.00 |

> *A. baumannii* showed the highest resistome burden (964 genes, 15.93% detection rate), consistent with its known multidrug-resistant phenotype. *E. coli* showed concentrated expression, with a single gene (hp1181, ARO:3003964) accounting for a TPM of **926,404.87**.

### Top Expressed AMR Genes in *E. coli* (Representative)

| Rank | Gene | ARO ID | TPM | Est. Counts |
|:---:|:---|:---|---:|---:|
| 1 | hp1181 | ARO:3003964 | 926,404.87 | 5,643.66 |
| 2 | CRP | ARO:3001883 | 63,085.16 | 262.75 |
| 3 | gadW | ARO:3006997 | 8,113.25 | 30.84 |
| 4 | genX | ARO: — | 1,285.32 | 6.75 |
| 5 | H-NS | ARO: — | 175.20 | 1.00 |

### Computational Efficiency

```
┌────────────────────────────────────────────────────┐
│  Throughput    : ~1.2 million reads / minute       │
│  Index size    : Compact (AMR-specific, not whole  │
│                  transcriptome)                    │
│  Intermediate  : No BAM files generated            │
│  Reusability   : Index computed once, reused       │
│                  across all datasets               │
│  Scalability   : Personal workstation → HPC cluster│
└────────────────────────────────────────────────────┘
```

---

## 📁 Output Format

### `amr_abundance.tsv`

```
gene_id         gene_name   aro_id        tpm           rpkm          est_counts
gene_001        hp1181      ARO:3003964   926404.87     481230.44     5643.66
gene_002        CRP         ARO:3001883   63085.16      32768.21      262.75
gene_003        gadW        ARO:3006997   8113.25       4214.88       30.84
...
```

### `summary_stats.txt`

```
AMR Analyser v1.0 | Run Summary
================================
Input reads          : 6,385,428
Aligned reads        : 5,938
Alignment rate       : 0.09%
AMR genes in DB      : 6,052
Detected genes       : 14
Detection rate       : 0.23%
Mean TPM             : 71,428.57
Median TPM           : 87.74
EM iterations        : 147 (converged)
Runtime              : ~5.3 min
```

---

## 🗃️ Supported Databases

| Database | Description | URL |
|:---|:---|:---|
| **CARD** | Comprehensive Antibiotic Resistance Database | [card.mcmaster.ca](https://card.mcmaster.ca) |
| **NCBI AMR** | NCBI Pathogen Detection AMR Reference Gene Catalog | [ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/) |
| **ResFinder** | Resistance gene finder database | [resfinder.cge.dtu.dk](https://resfinder.cge.dtu.dk) |

> The pipeline is database-agnostic — any FASTA-format AMR gene collection can be used to build a custom index.

---

## 🔮 Future Directions

- [ ] Integration with long-read sequencing (Oxford Nanopore / PacBio)
- [ ] Adaptive k-mer length optimization for divergent gene families
- [ ] Hybrid alignment-pseudoalignment mode for SNV resolution
- [ ] Machine learning module for resistance phenotype classification
- [ ] Multi-omics integration (proteomics + transcriptomics)
- [ ] Metatranscriptomic pipeline for mixed microbial communities
- [ ] Web interface for clinical deployment

---

## 📖 Citation

If you use **AMR Analyser** in your research, please cite:

```bibtex
@mastersthesis{rout2026amranalyser,
  author    = {Ankush Kumar Rout},
  title     = {AMR Analyser: A Pseudoalignment-Based Pipeline for Rapid
               Detection and Quantification of Antimicrobial Resistance
               Genes from RNA-seq Data},
  school    = {National Institute of Technology Rourkela},
  year      = {2026},
  month     = {May},
  department= {Department of Life Science},
  supervisor= {Dr. Akhilesh Mishra},
  note      = {M.Sc. Dissertation, Roll No. 424LS2028}
}
```

**Key references used in this work:**

> Bray NL, Pimentel H, Melsted P, Pachter L. **Near-optimal probabilistic RNA-seq quantification.** *Nat Biotechnol.* 2016;34(5):525-527. doi:10.1038/nbt.3519

> Alcock BP et al. **CARD 2023: expanded curation, support for machine learning, and resistome prediction.** *Nucleic Acids Res.* 2023;51(D1):D690-D699.

---

## 👨‍🔬 Author

<div align="center">

<img src="https://github-readme-stats.vercel.app/api?username=ankushrouth&show_icons=true&theme=merko&hide_border=true&bg_color=0d1117&title_color=0D6E56&icon_color=1a9e80&text_color=c9d1d9" height="160"/>

**Ankush Kumar Rout**
M.Sc. Life Science · NIT Rourkela (2024–2026)
*Supervisor: Dr. Akhilesh Mishra, Dept. of Life Science*

[![LinkedIn](https://img.shields.io/badge/LinkedIn-ankushrouth-0077B5?style=flat-square&logo=linkedin)](https://linkedin.com/in/ankushrouth)
[![Email](https://img.shields.io/badge/Email-ankushrouth10@gmail.com-D14836?style=flat-square&logo=gmail)](mailto:ankushrouth10@gmail.com)
[![IIT JAM](https://img.shields.io/badge/IIT%20JAM-AIR%20174-0d6e56?style=flat-square)](https://github.com/ankushrouth)

</div>

---

## 📜 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

<div align="center">

<img src="https://capsule-render.vercel.app/api?type=waving&color=0:1db88a,100:0a3d2e&height=100&section=footer&animation=fadeIn" width="100%"/>

*"Resistance may evolve. Our tools must evolve faster."*

⭐ If AMR Analyser helps your research, please consider giving it a star!

![Visitors](https://komarev.com/ghpvc/?username=ankushrouth&label=Repo+Views&color=0d6e56&style=flat-square)

</div>
