# TCR Diversity Analysis: CMV+ vs CMV- Comparison

## Overview

This analysis compares T cell receptor (TCR) repertoire diversity between CMV-positive (CMV+) and CMV-negative (CMV-) individuals, while controlling for potential confounding factors such as age, sex, and genetic ancestry.

## What Was Accomplished

### 1. Data Structure Analysis ✓

**Metadata (41 CMV+/- participants):**
- 14 CMV+ individuals
- 27 CMV- individuals
- Age range: 17-81 years (median CMV+: 57, CMV-: 52)
- Sex distribution: balanced between groups
- Ancestry: Caucasian, African, Asian
- **No significant age or sex differences between groups** (good for avoiding confounding!)

**TCR Repertoire Data Structure:**
- `sequence_id`: Unique sequence identifier
- `clone_id`: Identifies clonal families (same clone = expanded T cells)
- `productive`: Whether sequence is productive (functional TCR)
- `cdr3` / `cdr3_aa`: CDR3 region (antigen-binding site)
- `v_call`, `j_call`: V and J gene segments used

### 2. Diversity Metrics Calculated ✓

For the sample analyzed (BFI-0000234):
- **Shannon Entropy: 2.639** - Measures overall diversity
- **Simpson Index: 0.929** - Probability two sequences are from different clones
- **Clonality: 0.000** - Very low clonality (highly diverse repertoire)
- **Clone Richness: 14** - Total unique clones detected

### 3. Analysis Pipeline Features ✓

The script (`tcr_diversity_analysis.py`) includes:

1. **Metadata Processing**
   - Automatic CMV status extraction
   - Confounder identification
   - Statistical testing for group differences

2. **Diversity Calculations**
   - Shannon entropy (information-theoretic diversity)
   - Simpson index (probability-based diversity)
   - Clonality (degree of clonal expansion)
   - Clone richness (unique clone count)

3. **Statistical Analysis**
   - Mann-Whitney U tests (non-parametric comparison)
   - Effect size calculations
   - Age correlation analysis (confounder check)
   - Spearman correlations

4. **Visualizations**
   - Shannon diversity comparison (CMV+ vs CMV-)
   - Clonality comparison
   - Age vs diversity scatter plots
   - Clone richness distributions

## Current Limitation

⚠️ **Only 1 sample was analyzed** because:
- We created only one TCR file (`part_table_bfi-0000234.tsv`) as a demonstration
- This sample (BFI-0000234) is from the HIV study, not the CMV study
- **You need to add more files to perform meaningful comparisons!**

## How to Add Your Full Dataset

### Step 1: Copy TCR Files to `data/` Directory

From your Windows directory, copy TCR files to the `data/` folder:

```bash
# Example files to copy (from metadata with CMV status):
# CMV+ samples:
part_table_bfi-0003053.tsv.gz  # M64-004, 68M, Caucasian, CMV+
part_table_bfi-0003054.tsv.gz  # M64-005, 59F, Caucasian, CMV+
part_table_bfi-0003055.tsv.gz  # M64-006, 58M, African, CMV+
# ... and 11 more CMV+ samples

# CMV- samples:
part_table_bfi-0003051.tsv.gz  # M64-002, 60M, Caucasian, CMV-
part_table_bfi-0003052.tsv.gz  # M64-003, 44M, Caucasian, CMV-
part_table_bfi-0003057.tsv.gz  # M64-008, 19F, Caucasian, CMV-
# ... and 24 more CMV- samples
```

**Recommended:** Copy all files that match `BFI-0003*` from your dataset, as these are from the Stanford Blood Center study with CMV status information.

### Step 2: Re-run the Analysis

```bash
python tcr_diversity_analysis.py
```

The script will:
- Automatically detect all TCR files in `data/`
- Calculate diversity for each individual
- Compare CMV+ vs CMV- groups
- Generate statistical tests and visualizations

### Step 3: Examine Results

**Output files:**
- `tcr_diversity_metrics.csv` - All calculated metrics per individual
- `tcr_diversity_analysis.png` - Comprehensive visualization

## Expected Biological Results

### Hypothesis: CMV+ individuals show altered TCR diversity

**Expected findings with full dataset:**

1. **Higher Clonality in CMV+**
   - CMV infection drives memory T cell expansion
   - Specific clones targeting CMV epitopes expand
   - Results in oligoclonal repertoires (dominated by few clones)

2. **Lower Shannon Entropy in CMV+**
   - Persistent viral antigen exposure
   - Memory inflation phenomenon
   - Less diverse overall repertoire

3. **Age Effects**
   - Older individuals may show:
     - Reduced repertoire diversity
     - More clonal expansion
     - This could confound CMV effects if not controlled

## Understanding the Metrics

### Shannon Entropy
- **Range:** 0 to log(N) where N = number of unique clones
- **Interpretation:** Higher values = more diverse
- **Formula:** H = -Σ(p_i × log(p_i))
- **Biological meaning:** How evenly distributed are the clones?

### Simpson Index
- **Range:** 0 to 1
- **Interpretation:** Higher values = more diverse
- **Formula:** D = 1 - Σ(p_i²)
- **Biological meaning:** Probability two random T cells are from different clones

### Clonality
- **Range:** 0 to 1
- **Interpretation:** Higher values = more clonal expansion
- **Formula:** 1 - (H / H_max)
- **Biological meaning:** How dominated is the repertoire by expanded clones?

### Clone Richness
- **Range:** 1 to N (total sequences)
- **Interpretation:** Higher values = more unique clones
- **Biological meaning:** Total number of different T cell clones detected

## Statistical Approach

### Primary Analysis
- **Test:** Mann-Whitney U test (non-parametric)
- **Why:** TCR diversity metrics are often not normally distributed
- **Comparison:** CMV+ vs CMV- for each metric
- **Significance:** p < 0.05

### Confounder Control
1. **Age:** Spearman correlation with diversity metrics
2. **Sex:** Check distribution balance between groups
3. **Ancestry:** Potential stratification if unbalanced

### If Confounders Detected
- Perform multivariate regression
- Stratified analysis (e.g., within age groups)
- Propensity score matching

## Advanced Analyses (Future Extensions)

Once you have the full dataset, you can extend this to:

1. **V/J Gene Usage Analysis**
   - Compare gene segment frequencies between groups
   - Identify CMV-associated gene preferences

2. **CDR3 Sequence Analysis**
   - Length distribution comparisons
   - Amino acid composition analysis
   - Public TCR identification (sequences shared across individuals)

3. **Clonal Overlap**
   - Identify shared clones between CMV+ individuals
   - Compare public vs private clones

4. **Longitudinal Analysis**
   - If you have time-series data
   - Track repertoire changes over time

## File Structure

```
TCR-BCR-seq-analysis/
├── data/                                  # Data directory (gitignored)
│   ├── metadata.tsv                       # Patient metadata
│   ├── part_table_bfi-0000234.tsv         # Sample TCR data
│   └── [add more TCR files here]
├── tcr_diversity_analysis.py              # Main analysis script
├── tcr_diversity_metrics.csv              # Output: calculated metrics
├── tcr_diversity_analysis.png             # Output: visualizations
└── TCR_DIVERSITY_ANALYSIS_README.md       # This file
```

## Troubleshooting

### Issue: "No CMV samples found"
- Check that TCR files match participant IDs in metadata
- Ensure files are named `part_table_bfi-XXXXXXX.tsv` or `.tsv.gz`

### Issue: "Module not found"
```bash
pip install pandas numpy matplotlib seaborn scipy
```

### Issue: "Empty plots"
- Verify you have both CMV+ and CMV- samples
- Check that files are in `data/` directory
- Ensure metadata has `disease_subtype` column with CMV information

## References

This analysis is based on standard TCR repertoire analysis approaches:

1. **Diversity Metrics:**
   - Greiff et al. (2015) "Quantitative assessment of the robustness of next-generation sequencing of antibody variable gene repertoires from immunized mice"

2. **CMV and TCR Repertoires:**
   - Wertheimer et al. (2014) "Aging and cytomegalovirus infection differentially and jointly affect distinct circulating T cell subsets in humans"

3. **Statistical Approaches:**
   - Zaslavsky et al. (2025) Science - Your dataset source!

## Contact

For questions or issues with this analysis, please refer to the main repository documentation.

---

**Next Steps:**
1. Copy more TCR files to `data/` directory (aim for all 41 CMV+/- samples)
2. Re-run: `python tcr_diversity_analysis.py`
3. Examine outputs and interpret results
4. Consider advanced analyses if basic comparison is significant
