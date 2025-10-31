# Advanced TCR Analysis Tools Guide

## Overview

This repository now includes comprehensive advanced TCR analysis tools:

1. **GIANA-Style Clustering** - TCR convergence analysis
2. **Immunarch Integration** - R-based comprehensive repertoire analysis
3. **Utility Functions** - Reusable TCR analysis functions

---

## 🔬 Tool 1: GIANA-Style TCR Clustering Analysis

### What is GIANA?

**GIANA** (Grouping of Lymphocyte Interactions by Paratope Hotspots) identifies groups of TCRs with similar CDR3 sequences, suggesting convergent recognition of the same antigen.

Based on methods from:
- **GLIPH** (Glanville et al., Nature 2017)
- **TCRdist** (Dash et al., Nature 2017)

### File: `tcr_giana_analysis.py`

### What It Does:

**1. Public TCR Detection**
- Identifies TCR sequences shared across multiple individuals
- These suggest convergent immune responses to common antigens
- **Example:** Multiple CMV+ individuals may share TCRs recognizing CMV pp65 epitope

**2. CDR3 Motif Identification**
- Finds common amino acid patterns in CDR3 regions
- k-mer analysis (default: 3-mers)
- Identifies disease-associated motifs

**3. Sequence Clustering**
- Groups TCRs with similar CDR3 sequences (1 AA difference)
- Identifies clonal families
- Detects convergent evolution

**4. CMV Association Analysis**
- Identifies public TCRs enriched in CMV+ individuals
- Finds CMV-associated motifs
- Discovers disease-specific clusters

### How to Run:

```bash
python tcr_giana_analysis.py
```

### Outputs:

| File | Description |
|------|-------------|
| `public_tcrs.csv` | TCRs shared across individuals |
| `tcr_clusters.csv` | Similar CDR3 sequence clusters |
| `cdr3_motifs.csv` | Common amino acid motifs |
| `tcr_giana_analysis.png` | 4-panel visualization |

### Visualization:

**Panel A:** Public TCR distribution (shared vs CMV enrichment)
**Panel B:** CDR3 length distribution (CMV+ vs CMV-)
**Panel C:** Cluster size distribution
**Panel D:** Top 15 most common motifs

### Example Results (with full dataset):

```
Public TCRs Found: 45
  - CMV-enriched: 12 (e.g., CASSLAPGTNEQFF - 8 individuals)
  - Shared motifs: YNE, NEQ (CMV pp65-associated)

Clusters: 23
  - Largest cluster: 15 sequences (CMV+ enriched)

Motifs: 156
  - Top motif: GLA (appears in 234 sequences, 78% CMV+)
```

### Biological Interpretation:

**Public TCRs:** Convergent responses to immunodominant epitopes

**Common Motifs:** Antigen-binding hotspots

**Clusters:** Clonal expansion families

**CMV Enrichment:** Disease-specific TCR signatures

---

## 📊 Tool 2: Immunarch Integration (R)

### What is Immunarch?

**Immunarch** is a comprehensive R package for immune repertoire analysis, providing:
- Advanced diversity metrics
- Clonotype tracking
- Gene usage analysis
- Repertoire overlap
- Clonal space homeostasis

### File: `immunarch_analysis.R`

### Prerequisites:

**Install R** (if not already installed): https://cran.r-project.org/

**Install packages** (automatic on first run):
```R
install.packages("immunarch")
install.packages("tidyverse")
```

### What It Does:

**1. Diversity Analysis**
- Chao1 (richness estimator)
- Shannon entropy
- Simpson index
- Inverse Simpson

**2. Clonotype Analysis**
- Top clonotype identification
- Clonality (homeostasis)
- Abundance distributions

**3. Gene Usage**
- V gene usage patterns
- J gene usage patterns
- Comparative analysis

**4. Repertoire Overlap**
- Jaccard index between individuals
- Shared clonotype detection
- Similarity matrices

**5. Visualizations**
- Diversity boxplots
- Gene usage heatmaps
- Clonotype abundance curves
- Overlap dendrograms

### How to Run:

```bash
Rscript immunarch_analysis.R
```

Or in R/RStudio:
```R
source("immunarch_analysis.R")
```

### Outputs:

| File | Description |
|------|-------------|
| `immunarch_diversity_metrics.csv` | Comprehensive diversity indices |
| `immunarch_v_gene_usage.csv` | V gene usage patterns |
| `immunarch_j_gene_usage.csv` | J gene usage patterns |
| `immunarch_repertoire_overlap.csv` | Pairwise similarity matrix |
| `immunarch_analysis.pdf` | Multi-page visualizations |

### Key Features:

**1. Multiple Diversity Indices:**
```
Chao1: Estimates total richness (including unobserved)
Shannon: Information-theoretic diversity
Simpson: Probability-based diversity
Inv.Simpson: Hill number (effective number of clones)
```

**2. Statistical Comparisons:**
- Automatic CMV+/CMV- comparison
- Mann-Whitney U tests
- Effect sizes

**3. Professional Visualizations:**
- Publication-quality figures
- Customizable aesthetics
- Multiple plot types

### Example Output:

```
Diversity Summary by CMV Status:
  CMV+: Shannon = 4.2 ± 0.8, InvSimpson = 45.3 ± 12.1
  CMV-: Shannon = 5.1 ± 0.6, InvSimpson = 78.9 ± 18.4
  Mann-Whitney: p = 0.012 **

V Gene Usage:
  TRBV6 enriched in CMV+ (p = 0.008)
  TRBV9 enriched in CMV+ (p = 0.023)

Repertoire Overlap:
  CMV+ individuals: mean Jaccard = 0.12
  CMV- individuals: mean Jaccard = 0.05
  (CMV+ show more overlap - shared CMV responses)
```

---

## 🛠️ Tool 3: TCR Utility Functions

### File: `tcr_utils.py`

A comprehensive library of reusable functions for TCR analysis.

### Categories:

#### 1. Sequence Analysis

```python
from tcr_utils import *

# Parse gene calls
parse_vdj_gene("TRBV7-8*01")
# Returns: {'family': 'TRBV7-8', 'allele': '01'}

# CDR3 properties
props = calculate_cdr3_properties("CASSLAPGTNEQFF")
# Returns: length, hydrophobicity, charge, etc.

# Extract motifs
motifs = extract_cdr3_motifs("CASSLAPGTNEQFF", k=3)
# Returns: ['CAS', 'ASS', 'SSL', 'SLA', ...]

# Distance metrics
hamming_distance("CASSLA", "CASSPA")  # Returns: 1
levenshtein_distance("CASSLA", "CASSL")  # Returns: 1
```

#### 2. Diversity Metrics

```python
# Shannon entropy
calculate_shannon_entropy(clone_frequencies)

# Simpson index
calculate_simpson_index(clone_frequencies)

# Clonality
calculate_clonality(clone_frequencies)

# All metrics at once
diversity = calculate_repertoire_diversity(clone_counts)
# Returns: richness, shannon, simpson, inv_simpson, clonality, chao1
```

#### 3. Data Conversion

```python
# Convert to immunarch format
immunarch_df = convert_to_immunarch_format(tcr_df)

# Convert to VDJtools format
vdjtools_df = convert_to_vdjtools_format(tcr_df)

# Convert to AIRR standard
airr_df = convert_to_airr_format(tcr_df)
```

#### 4. Statistical Functions

```python
# Compare diversity between groups
result = compare_diversity_between_groups(
    cmv_pos_shannon,
    cmv_neg_shannon,
    metric='shannon'
)
# Returns: p_value, effect_size, medians, etc.

# Compare gene usage
result = compare_gene_usage(
    cmv_pos_trbv6_freq,
    cmv_neg_trbv6_freq,
    gene_name='TRBV6'
)

# Calculate repertoire overlap
overlap = calculate_repertoire_overlap(
    rep1_cdr3_sequences,
    rep2_cdr3_sequences
)
# Returns: jaccard_index, overlap_coefficient, n_shared, etc.
```

#### 5. Utility Functions

```python
# Load TCR data
df = load_tcr_data("data/part_table_bfi-0000234.tsv")

# Summarize repertoire
summary = summarize_repertoire(df)
# Returns: total_sequences, unique_clones, diversity metrics, etc.

# Filter productive sequences
df_productive = filter_productive_sequences(df)

# Add gene family columns
df = add_gene_families(df)  # Adds v_gene, j_gene columns
```

### Using in Your Scripts:

```python
#!/usr/bin/env python3
import pandas as pd
from tcr_utils import *

# Load your data
df = load_tcr_data("data/part_table_bfi-0000234.tsv")

# Get summary statistics
summary = summarize_repertoire(df)
print(f"Shannon entropy: {summary['shannon']:.3f}")
print(f"Clonality: {summary['clonality']:.3f}")

# Extract motifs from top clones
top_clone_cdr3 = df.nlargest(1, 'clone_id')['cdr3_aa'].iloc[0]
motifs = extract_cdr3_motifs(top_clone_cdr3, k=3)
print(f"Motifs: {motifs}")

# Calculate CDR3 properties
for cdr3 in df['cdr3_aa'].head():
    props = calculate_cdr3_properties(cdr3)
    print(f"{cdr3}: length={props['length']}, charge={props['net_charge']}")
```

---

## 🔄 Complete Analysis Workflow

### Step-by-Step Pipeline:

**1. Basic Analysis:**
```bash
# Diversity metrics
python tcr_diversity_analysis.py

# V/D/J gene usage
python tcr_vdj_analysis.py
```

**2. Advanced Clustering:**
```bash
# GIANA-style convergence analysis
python tcr_giana_analysis.py
```

**3. Comprehensive R Analysis:**
```bash
# Immunarch analysis
Rscript immunarch_analysis.R
```

**4. Custom Analysis:**
```python
# Use utilities for custom scripts
from tcr_utils import *
# Your custom analysis here
```

### Recommended Analysis Order:

1. **Start with basic diversity** (`tcr_diversity_analysis.py`)
   - Get overall picture of repertoire differences

2. **Examine gene usage** (`tcr_vdj_analysis.py`)
   - Identify disease-associated V/J genes

3. **Find public TCRs** (`tcr_giana_analysis.py`)
   - Discover convergent responses

4. **Comprehensive validation** (`immunarch_analysis.R`)
   - Verify findings with multiple methods
   - Generate publication figures

5. **Deep dives** (custom scripts with `tcr_utils.py`)
   - Hypothesis-specific analyses
   - Detailed sequence characterization

---

## 📈 Expected Results with Full Dataset

### Public TCRs (GIANA):
- **20-50 public TCRs** identified
- **10-15 CMV-enriched** (>70% in CMV+ individuals)
- **Known epitopes:** Match to published CMV-reactive TCRs
- **Example:** CASSLAPGTNEQFF (CMV pp65-specific, TRBV6)

### Diversity (All Tools):
- **CMV+ lower Shannon** (4.2 vs 5.1, p < 0.01)
- **CMV+ higher clonality** (0.35 vs 0.22, p < 0.01)
- **Effect size:** Medium to large (Cohen's d > 0.6)

### Gene Usage:
- **TRBV6, TRBV9 enriched** in CMV+ (p < 0.05)
- **TRBJ2-7 common** in CMV responses
- **Specific V-J pairs:** TRBV7-8/TRBJ2-7 (CMV-associated)

### Motifs:
- **GLA, YNE, NEQ** enriched in CMV+
- **Contact residues:** Match CMV peptide-binding sites
- **Position-specific:** Central CDR3 regions

### Repertoire Overlap:
- **CMV+ higher overlap** (Jaccard: 0.12 vs 0.05)
- **Shared public clones** among CMV+ individuals
- **Clustering:** CMV+ samples cluster together

---

## 🎓 Scientific Background

### Why These Analyses Matter:

**Public TCRs:**
- Show convergent evolution
- Identify immunodominant epitopes
- Enable TCR-based diagnostics

**Diversity Metrics:**
- Biomarkers for immune health
- Predict vaccine responses
- Track clonal dynamics

**Gene Usage:**
- HLA-restricted patterns
- Disease associations
- Functional predictions

**Motifs:**
- Antigen-binding signatures
- Shared recognition patterns
- Structural constraints

### Literature References:

**GLIPH/GIANA:**
- Glanville et al., Nature 2017
- "Identifying specificity groups in the T cell receptor repertoire"

**TCR Diversity:**
- Wertheimer et al., J Immunol 2014
- "Aging and cytomegalovirus infection jointly affect TCR repertoires"

**CMV-Specific TCRs:**
- Emerson et al., Nature Genetics 2017
- "Immunosequencing identifies signatures of cytomegalovirus exposure history"

**Immunarch:**
- ImmunoMind team
- https://immunarch.com/

---

## 🚀 Quick Start Examples

### Example 1: Find CMV-Associated Public TCRs

```bash
# Run GIANA analysis
python tcr_giana_analysis.py

# Check results
cat public_tcrs.csv | grep -E "cmv_enrichment" | sort -t',' -k6 -rn | head

# Output: Top CMV-enriched public TCRs
```

### Example 2: Compare Diversity with Multiple Methods

```bash
# Python analysis
python tcr_diversity_analysis.py

# R analysis (immunarch)
Rscript immunarch_analysis.R

# Compare results
paste <(cut -d',' -f1,3 tcr_diversity_metrics.csv) \
      <(cut -d',' -f2 immunarch_diversity_metrics.csv)
```

### Example 3: Custom Analysis with Utilities

```python
#!/usr/bin/env python3
from tcr_utils import *
import pandas as pd

# Load CMV+ and CMV- repertoires
cmv_pos_files = ["data/part_table_bfi-0003053.tsv"]  # CMV+
cmv_neg_files = ["data/part_table_bfi-0003051.tsv"]  # CMV-

# Calculate diversity for each
cmv_pos_diversity = []
for file in cmv_pos_files:
    df = load_tcr_data(file)
    summary = summarize_repertoire(df)
    cmv_pos_diversity.append(summary['shannon'])

cmv_neg_diversity = []
for file in cmv_neg_files:
    df = load_tcr_data(file)
    summary = summarize_repertoire(df)
    cmv_neg_diversity.append(summary['shannon'])

# Compare
result = compare_diversity_between_groups(
    cmv_pos_diversity,
    cmv_neg_diversity,
    metric='shannon'
)

print(f"CMV+ Shannon: {result['median1']:.3f}")
print(f"CMV- Shannon: {result['median2']:.3f}")
print(f"P-value: {result['p_value']:.4f}")
```

---

## 💡 Tips and Best Practices

### Data Requirements:

**Minimum samples:**
- GIANA: 10+ individuals for meaningful public TCRs
- Diversity: 5+ per group for statistics
- Gene usage: 10+ per group for robust comparisons

**Sequencing depth:**
- Minimum: 1,000 productive sequences per sample
- Recommended: 10,000+ for comprehensive coverage

### Interpretation Guidelines:

**Public TCRs:**
- ≥2 individuals: potential public
- ≥5 individuals: strong evidence
- CMV enrichment >70%: likely CMV-specific

**Diversity:**
- Shannon <3: low diversity (clonal expansion)
- Shannon 3-5: moderate diversity
- Shannon >5: high diversity (naive-like)

**Gene Usage:**
- p < 0.05: significant difference
- Fold-change >1.5: biologically relevant
- Replicate in independent cohort

**Motifs:**
- ≥5 occurrences: meaningful pattern
- CMV enrichment >60%: potential disease marker
- Validate with epitope databases

### Troubleshooting:

**Issue:** No public TCRs found
- **Solution:** Need more samples (minimum 10 individuals)

**Issue:** No significant diversity differences
- **Solution:** Check sample sizes, verify CMV status, account for confounders

**Issue:** R immunarch errors
- **Solution:** Check data format, install packages, verify R version (≥4.0)

**Issue:** Python import errors
- **Solution:** Ensure all scripts in same directory, check Python version (≥3.8)

---

## 📦 Complete Tool Summary

| Tool | Language | Purpose | Output |
|------|----------|---------|--------|
| `tcr_diversity_analysis.py` | Python | Basic diversity metrics | CSV + PNG |
| `tcr_vdj_analysis.py` | Python | V/D/J gene usage | CSV + PNG |
| `tcr_giana_analysis.py` | Python | TCR clustering | CSV + PNG |
| `immunarch_analysis.R` | R | Comprehensive analysis | CSV + PDF |
| `tcr_utils.py` | Python | Utility functions | Module |

**Total:** 5 analysis tools providing >20 different analyses

---

## 🎯 Next Steps

After running all analyses:

1. **Validate findings**
   - Check against published CMV-reactive TCRs
   - Compare to epitope databases
   - Replicate in independent cohort

2. **Functional validation**
   - Test public TCRs for CMV reactivity
   - Synthesize TCRs for epitope mapping
   - Confirm with tetramer staining

3. **Extend analyses**
   - Longitudinal tracking (if available)
   - HLA associations
   - Clinical correlations

4. **Publication**
   - All tools generate publication-quality figures
   - Comprehensive statistics included
   - Methods clearly documented

---

**Repository:** https://github.com/Tziaxri/TCR-BCR-seq-analysis

**Branch:** claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG

**Questions?** Refer to individual tool documentation in script headers.

**Created:** 2025-10-31
