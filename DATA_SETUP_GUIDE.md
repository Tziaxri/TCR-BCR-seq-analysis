# TCR Diversity Analysis - Data Setup Guide

## Overview
This guide explains where to place your TCR-seq data files for preprocessing and diversity analysis using immunarch.

## Directory Structure

```
TCR-BCR-seq-analysis/
├── data/
│   ├── raw_data/           # Place your 10 .tsv.gz files here
│   ├── processed_data/     # Preprocessed files will be saved here
│   └── metadata.tsv        # Place your metadata file here
├── scripts/
│   ├── 01_preprocess_tcr_data.R
│   └── 02_diversity_analysis.R
└── results/
    ├── figures/            # Diversity plots and visualizations
    └── tables/             # Summary statistics and tables
```

## Step 1: File Placement

### Raw TCR-seq Data (10 files for initial analysis)
**Location:** `data/raw_data/`

Place your 10 selected `.tsv.gz` files in this directory. Choose files that represent different conditions:
- Healthy controls (2-3 files)
- COVID-19 (2-3 files)
- At least 1-2 files from other conditions (HIV, lupus, T1D, or flu vaccination)

This will give you a good mix for preliminary diversity comparison.

### Metadata File
**Location:** `data/metadata.tsv`

Place your `metadata.tsv` file directly in the `data/` directory (NOT in raw_data/).

## Step 2: What the Metadata Should Contain

The metadata file should have at minimum:
- `repertoire_id` or `Sample` column (matching the .tsv.gz file identifiers)
- `condition` or `immune_state` column (COVID-19, HIV, lupus, T1D, flu, healthy)
- `age` (for controlling confounding factors)
- `geography` or `location` (if available)
- `genetic_ancestry` (if available)
- Any other relevant demographic/clinical variables

## Step 3: Current Format → Immunarch Format

### Your Current Format (Per-Read):
```
repertoire_id   Sample
productive
cdr3_aa         CDR3.aa
v_call          V.name
j_call          J.name
d_call
```

### Required Immunarch Format (Per-Clonotype):
Immunarch expects one row per unique clonotype with:
- `Clones` - Number of reads for this clonotype (abundance)
- `Proportion` - Frequency/proportion of this clonotype
- `CDR3.nt` - CDR3 nucleotide sequence
- `CDR3.aa` - CDR3 amino acid sequence
- `V.name` - V gene name
- `D.name` - D gene name (optional)
- `J.name` - J gene name

## Step 4: Preprocessing Steps

The preprocessing script will:
1. Read all .tsv.gz files from `data/raw_data/`
2. Filter for productive sequences only
3. Group by unique clonotypes (CDR3.aa + V + J combinations)
4. Count reads per clonotype
5. Calculate proportions
6. Save in immunarch-compatible format to `data/processed_data/`

## Step 5: File Selection Strategy

For your 10 test files, I recommend selecting:
1. Check the metadata.tsv for repertoire_ids
2. Select 2-3 healthy controls
3. Select 2-3 COVID-19 samples
4. Select 1-2 from each of 2 other conditions
5. Try to balance for age/geography if possible

This balanced selection will allow meaningful preliminary diversity comparisons.

## Next Steps

After placing files in the correct locations:
1. Run preprocessing script to convert format
2. Run diversity analysis comparing groups
3. Generate diversity metrics (Shannon, Simpson, Chao1, etc.)
4. Create visualizations comparing TCR diversity across conditions
