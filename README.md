# In Silico Genomic Mutation Analysis of Dengue Virus (DENV)

**Arif Rahmawan** | Independent Researcher  
| 2025

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

This repository contains all Python scripts used in the in silico analysis of dengue virus (DENV) genomic mutations, conducted as an independent research project. The analysis covers four DENV serotypes (DENV-1 through DENV-4) using 1,092 sequences from Indonesia and four regional countries (Thailand, Singapore, Philippines, Malaysia), downloaded from GISAID EpiArbo™.

The pipeline performs multi-gene analysis (envelope protein [E], non-structural protein 1 [NS1], and precursor membrane protein [prM]) to characterize mutation patterns, phylogenetic relationships, selection pressure, and genomic surveillance implications for health quarantine strengthening.

---

## Research Scope

| Parameter | Detail |
|---|---|
| **Virus** | Dengue virus (DENV) serotypes 1–4 |
| **Genes analyzed** | E (envelope), NS1 (non-structural 1), prM (precursor membrane) |
| **Total sequences** | 1,092 isolates (after QC) |
| **Countries** | Indonesia (n=376), Thailand (n=223), Singapore (n=112), Philippines (n=336), Malaysia (n=45) |
| **Collection period** | 2015–2025 |
| **Data source** | GISAID EpiArbo™ (isolates); NCBI GenBank (RefSeq) |
| **Reference sequences** | NC_001477.1 (DENV-1), NC_001474.2 (DENV-2), NC_001475.2 (DENV-3), NC_002640.1 (DENV-4) |

---

## Pipeline Overview

The analysis consists of **21 steps across 6 phases**:

```
PHASE 1 — Data Collection & Curation
  step01  Sequence curation (Indonesia + Regional)
  step02  Gene subsetting with dynamic offset detection
  step02b RefSeq nucleotide subsetting

PHASE 2 — Alignment & QC
  step03  FASTA merging & header relabeling (city metadata)
  [MAFFT via UGENE — external tool]
  validasi_alignment  Alignment QC (gap columns, consistency)

PHASE 3 — Phylogenetics & Statistics
  [IQ-TREE3 — external tool]
  [RDP5 — external tool, Windows]
  step08b Recombinant sequence identification
  step09  Fisher's Exact Test for mutation frequency

PHASE 4 — Nucleotide Mutation Analysis
  step10  GRID Profile construction
  step11  Mutation identification & hotspot detection

PHASE 5 — Protein Analysis
  step14  RefSeq protein download (NCBI)
  step15  Protein translation + BLASTp verification
  step16  Protein MSA
  step17  dN/dS ratio calculation (selection pressure)

PHASE 6 — Interpretation & Recommendations
  step18  3D structural mapping (PyMOL)
  step19  Biological interpretation
  step20  QC audit
  step21  Genomic barcode for health quarantine
```

---

## Repository Structure

```
denv-phylogenetic-analysis/
├── README.md                          # This file
├── requirements.txt                   # Python dependencies
├── LICENSE                            # MIT License
│
├── pipeline/                          # Main analysis scripts
│   ├── step01_kurasi_fasta.py         # Sequence curation & cleaning
│   ├── step02_subsetting_pipeline.py  # Gene extraction (dynamic offset)
│   ├── step02b_refseq_subset.py       # RefSeq gene subsetting
│   ├── step03_merge_relabel.py        # FASTA merge + header relabeling
│   ├── step08b_identify_recombinant.py # Recombinant sequence detection
│   ├── step17_dnds.py                 # dN/dS ratio analysis
│   └── RUNBOOK_T12_T13.txt            # Manual runbook for Steps 12–13
│
├── validation/                        # QC & validation scripts
│   ├── validasi_koordinat.py          # Gene coordinate validation
│   ├── validasi_subset_isolat.py      # Isolate subset QC
│   ├── validasi_refseq_subset.py      # RefSeq subset QC
│   ├── validasi_combined_fasta.py     # Combined FASTA QC (pre-MSA)
│   ├── validasi_alignment.py          # Alignment QC (post-MSA)
│   ├── validasi_subsetting_refseq.py  # Subsetting accuracy proof (100% identity vs NCBI)
│   ├── validasi_t12_grid_profile.py   # GRID profile validation
│   ├── validasi_t13_mutations.py      # Mutation identification validation
│   ├── validasi_t15_protein.py        # Protein translation validation
│   ├── validasi_t16_msa_protein.py    # Protein MSA validation
│   ├── validasi_t17_dnds.py           # dN/dS validation
│   ├── inventaris_revisi4.py          # Revision 4 inventory
│   ├── inventaris_revisi5.py          # Revision 5 inventory
│   └── step09_check_progress.py       # Pipeline progress checker
│
└── data/                              # NOT included (see Data Availability)
    ├── 01_raw_data/                   # GISAID sequences (requires GISAID account)
    ├── 02_subsetting/                 # Gene subsets
    ├── 03_alignment/                  # MSA outputs
    └── 07_phylogenetic/               # IQ-TREE outputs
```

---

## Installation

### Prerequisites

- Python 3.10 or higher
- macOS / Linux (tested on macOS ARM64)
- [UGENE](https://ugene.net/) with MAFFT plugin — for MSA (Phase 2)
- [IQ-TREE3](http://www.iqtree.org/) — for phylogenetics (Phase 3)
- [RDP5](http://web.cbio.uct.ac.za/~darren/rdp.html) — for recombination detection (Windows/Wine)
- [PyMOL](https://pymol.org/) — for 3D structural visualization (Phase 6)

### Python Setup

```bash
# Clone repository
git clone https://github.com/yourusername/denv-phylogenetic-analysis.git
cd denv-phylogenetic-analysis

# Create virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

---

## Usage

### Phase 1 — Data Curation & Subsetting

```bash
# Step 1: Curate raw FASTA sequences
python pipeline/step01_kurasi_fasta.py

# Step 2: Extract target genes (E, NS1, prM) using dynamic offset
python pipeline/step02_subsetting_pipeline.py

# Step 2b: Extract RefSeq nucleotide subsets
python pipeline/step02b_refseq_subset.py
```

### Phase 2 — FASTA Preparation & MSA

```bash
# Step 3: Merge sequences + relabel headers with city metadata
python pipeline/step03_merge_relabel.py

# MSA: Open UGENE → load *_combined.fasta → MAFFT → save *_alignment.fasta
# (UGENE is a GUI tool; see Methods section for parameters)
```

### Validation (run after each phase)

```bash
# Validate subsetting accuracy
python validation/validasi_subset_isolat.py
python validation/validasi_refseq_subset.py
python validation/validasi_subsetting_refseq.py

# Validate pre-MSA FASTA
python validation/validasi_combined_fasta.py

# Validate alignment quality
python validation/validasi_alignment.py
```

### Phase 3 — IQ-TREE Phylogenetics (example for DENV-1 E gene)

```bash
iqtree3 \
  -s "03_alignment/DENV1/DENV1_E_alignment.fasta" \
  -m TEST \
  -B 1000 \
  -T AUTO \
  --prefix "07_phylogenetic/DENV1/DENV1_E"
```

Repeat for all 12 gene–serotype combinations (3 genes × 4 serotypes).

---

## Methods Summary

### Gene Coordinate Extraction

Gene boundaries were obtained from NCBI RefSeq mat_peptide annotations:

| Serotype | RefSeq | prM (bp) | E (bp) | NS1 (bp) |
|---|---|---|---|---|
| DENV-1 | NC_001477.1 | 437–934 | 935–2419 | 2420–3475 |
| DENV-2 | NC_001474.2 | 439–936 | 937–2421 | 2422–3477 |
| DENV-3 | NC_001475.2 | 437–934 | 935–2413 | 2414–3469 |
| DENV-4 | NC_002640.1 | 441–938 | 939–2423 | 2424–3479 |

### Dynamic Offset Detection

Since most GISAID sequences lack UTR regions, a dynamic offset algorithm was implemented: a 30-nucleotide anchor from the RefSeq CDS start was searched within each isolate sequence, and gene coordinates were adjusted individually. Sequences without a matching anchor (n=135) used a fixed offset=0 as fallback.

### Subsetting Accuracy Validation

Subsetting accuracy was validated by translating the 12 RefSeq nucleotide subsets and comparing position-by-position against official NCBI Protein RefSeq sequences. All 12 gene–serotype combinations showed **100% amino acid identity**, confirming the validity of the extraction coordinates.

### Multiple Sequence Alignment

MSA was performed using MAFFT v7.526 (auto mode, default parameters) via UGENE v48. No trimming was required as no 100%-gap columns were detected across all 12 alignments.

### Recombination Detection

Recombination was screened using RDP5 v5.84 (Martin et al., 2021) with three methods: RDP, GENECONV, and MaxChi. A recombination event was considered significant when confirmed by ≥3 independent methods (p<0.05). One confirmed event was detected in DENV-4 E gene.

### Phylogenetic Analysis

Maximum likelihood trees were inferred using IQ-TREE3 v3.1.2 with ModelFinder for automatic model selection and 1,000 ultrafast bootstrap replicates (-B 1000).

---

## Data Availability

Raw sequence data are deposited in GISAID EpiArbo™ (accession numbers available upon request; GISAID registration required). Reference sequences were obtained from NCBI GenBank:
- DENV-1: [NC_001477.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001477.1)
- DENV-2: [NC_001474.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002640.1)
- DENV-3: [NC_001475.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_001475.2)
- DENV-4: [NC_002640.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_002640.1)

All analysis scripts are archived at Zenodo: [DOI: 10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## Software Versions

| Software | Version | Purpose |
|---|---|---|
| Python | 3.13.7 | Core scripting |
| BioPython | ≥1.81 | Sequence handling |
| NumPy | ≥1.24.0 | Numerical analysis |
| Pandas | ≥2.0.0 | Metadata processing |
| SciPy | ≥1.10.0 | Statistical tests |
| UGENE | 48.x | MSA (MAFFT interface) |
| MAFFT | 7.526 | Multiple sequence alignment |
| IQ-TREE3 | 3.1.2 | ML phylogenetics |
| RDP5 | 5.84 | Recombination detection |
| PyMOL | 3.x | 3D structural visualization |

---

## Citation

If you use these scripts in your research, please cite:

> Rahmawan, A. (2025). *In Silico Genomic Mutation Analysis of Dengue Virus as Material for Genomic Surveillance and Health Quarantine Strengthening*. Independent Research, Bandung, Indonesia. DOI: 10.5281/zenodo.XXXXXXX

Please also cite the data sources:
- **GISAID:** Shu Y & McCauley J (2017). GISAID: Global initiative on sharing all influenza data. *Euro Surveill*, 22(13):30494.
- **RDP5:** Martin DP et al. (2021). RDP5: A computer program for analyzing recombination in, and removing signals of recombination from, nucleotide sequence datasets. *Virus Evolution*, 7(1):veaa087.
- **IQ-TREE3:** Minh BQ et al. (2020). IQ-TREE 2: New models and methods for phylogenetic inference. *Mol Biol Evol*, 37(5):1530–1534.
- **MAFFT:** Katoh K & Standley DM (2013). MAFFT multiple sequence alignment software version 7. *Mol Biol Evol*, 30(4):772–780.

---

## Acknowledgements

The author acknowledges the academic environment and resources at Poltekkes Kemenkes Bandung that supported the development of this research.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

**Arif Rahmawan**  
Independent Researcher  
Bandung, West Java, Indonesia  
📧 arif.rahmawan@kemkes.go.id
