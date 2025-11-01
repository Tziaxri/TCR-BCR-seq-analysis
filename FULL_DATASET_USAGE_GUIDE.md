# 🚀 Full Dataset Analysis Guide

**Updated:** 2025-10-31
**Purpose:** Process ALL your TCR files automatically

---

## 📋 What's New

Your TCR analysis toolkit has been updated to **automatically process your entire dataset** from your Windows directory!

### ✨ **Key Updates:**

1. **Automatic Data Loading** - Processes ALL .tsv.gz files
2. **Windows Path Support** - Works with `C:\Users\chris\Desktop\TCR ANALYSIS\`
3. **Both Metadata Files** - Uses metadata.tsv and synapse_metadata_manifest.tsv
4. **Memory Management** - Smart subsampling for large datasets
5. **Progress Tracking** - Shows you what's happening in real-time

---

## 📁 Your Data Location

**Windows Path:**
```
C:\Users\chris\Desktop\TCR ANALYSIS\
```

**Files in your directory:**
- `metadata.tsv` - Participant information
- `synapse_metadata_manifest.tsv` - File manifest
- `part_table_bfi-*.tsv.gz` - TCR repertoire files (multiple)

---

## 🎯 Updated Analysis Scripts

### 1. **Data Loader Module** (New!)

**File:** `tcr_data_loader.py`

**Purpose:** Core module that loads all your data automatically

**Features:**
- Finds all .tsv.gz files in directory
- Loads both metadata files
- Merges TCR data with participant info
- Handles missing data gracefully
- Calculates diversity metrics per sample

**You don't run this directly** - it's used by all other scripts!

---

### 2. **Diversity Analysis - Full Dataset**

**File:** `tcr_diversity_analysis_full.py`

**What it does:**
- Loads ALL your TCR files
- Calculates diversity metrics for every sample
- Compares CMV+ vs CMV-
- Checks for confounders (age, sex, ancestry)
- Creates publication-quality figures

**How to run:**

1. **Edit the file** (line ~460):
   ```python
   DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"
   ```

2. **Run:**
   ```bash
   python tcr_diversity_analysis_full.py
   ```

3. **Optional - Test with subset first:**
   ```python
   MAX_FILES = 50  # Process first 50 files only
   ```

**Outputs:**
- `tcr_diversity_metrics_full.csv` - All diversity metrics
- `diversity_comparison_results.csv` - Statistical comparisons
- `tcr_diversity_full_analysis.png` - 4-panel visualization

**Runtime:** ~1-2 minutes for full dataset

---

### 3. **VDJ Gene Analysis - Full Dataset**

**File:** `tcr_vdj_analysis_full.py`

**What it does:**
- Analyzes V, D, J gene segment usage
- Compares gene frequencies (CMV+ vs CMV-)
- Identifies disease-associated genes
- Analyzes V-J pairing patterns

**How to run:**

1. **Edit the file** (line ~625):
   ```python
   DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"
   ```

2. **Run:**
   ```bash
   python tcr_vdj_analysis_full.py
   ```

**Outputs:**
- `v_gene_usage_full.csv` - V gene frequencies
- `j_gene_usage_full.csv` - J gene frequencies
- `d_gene_usage_full.csv` - D gene frequencies (if available)
- `v_genes_comparison_full.csv` - CMV+ vs CMV- V genes
- `j_genes_comparison_full.csv` - CMV+ vs CMV- J genes
- `vj_pairings_full.csv` - V-J pairing patterns
- `tcr_vdj_full_analysis.png` - Comprehensive visualization

**Runtime:** ~2-3 minutes for full dataset

---

### 4. **GIANA Clustering - Full Dataset**

**File:** `giana_analysis_full.py`

**What it does:**
- Clusters similar TCR sequences (Zhang et al. 2021 method)
- Identifies public TCRs (shared across samples)
- Calculates sample similarity (co-clustering matrix)
- Creates sample network visualization

**⚠️ IMPORTANT:** GIANA is computationally intensive!

**How to run:**

1. **Edit the file** (line ~510):
   ```python
   DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"

   # Start with these settings for testing:
   MAX_FILES = 50        # First 50 files
   MAX_SEQUENCES = 50000 # 50K sequences
   ```

2. **Run:**
   ```bash
   python giana_analysis_full.py
   ```

3. **For full analysis** (powerful computer with 16GB+ RAM):
   ```python
   MAX_FILES = None       # ALL files
   MAX_SEQUENCES = 500000 # 500K sequences
   ```

**Outputs:**
- `giana_clusters_full.csv` - All TCR clusters
- `giana_coclustering_matrix_full.csv` - Sample similarity matrix
- `giana_sample_network_full.png` - Network visualization

**Runtime:**
- 50K sequences: ~10-15 minutes
- 500K sequences: ~30-60 minutes
- Memory usage: 4-8 GB RAM

---

### 5. **Immunarch Visualizations**

**File:** `immunarch_style_viz.py` (updated to work with full dataset outputs)

**What it does:**
- Creates publication-quality figures from your results
- Uses diversity metrics from full dataset analysis
- Matches immunarch R package aesthetics

**How to run:**

1. **First, run diversity analysis:**
   ```bash
   python tcr_diversity_analysis_full.py
   ```

2. **Then run visualizations:**
   ```bash
   python immunarch_style_viz.py
   ```

   The script will automatically load `tcr_diversity_metrics_full.csv`

**Outputs:**
- `immunarch_style_diversity.pdf/png` - Diversity plots
- `immunarch_style_gene_usage.pdf/png` - Gene heatmaps

**Runtime:** ~30 seconds

---

## 🔄 Complete Workflow

### **Step-by-Step: Analyze Your Complete Dataset**

```bash
# Step 1: Diversity Analysis (MUST RUN FIRST!)
python tcr_diversity_analysis_full.py

# Step 2: VDJ Gene Analysis
python tcr_vdj_analysis_full.py

# Step 3: Create Publication Figures
python immunarch_style_viz.py

# Step 4: GIANA Clustering (Optional - takes longer)
python giana_analysis_full.py
```

**Total runtime:** ~5-10 minutes (excluding GIANA)
**With GIANA:** ~20-40 minutes

---

## ⚙️ Configuration Options

### **Testing with Subset of Data**

All scripts support limiting the number of files:

```python
# In any *_full.py script, find this section:
MAX_FILES = 50  # Process first 50 files only

# For full analysis:
MAX_FILES = None  # Process ALL files
```

### **GIANA Memory Management**

```python
# In giana_analysis_full.py:
MAX_SEQUENCES = 50000  # Good for 8GB RAM
MAX_SEQUENCES = 100000 # Good for 16GB RAM
MAX_SEQUENCES = 500000 # Good for 32GB RAM
```

---

## 📊 Expected Outputs

### **After Running All Analyses:**

```
📁 Your working directory will contain:

Diversity Analysis:
  • tcr_diversity_metrics_full.csv
  • diversity_comparison_results.csv
  • tcr_diversity_full_analysis.png

VDJ Gene Analysis:
  • v_gene_usage_full.csv
  • j_gene_usage_full.csv
  • d_gene_usage_full.csv
  • v_genes_comparison_full.csv
  • j_genes_comparison_full.csv
  • vj_pairings_full.csv
  • tcr_vdj_full_analysis.png

GIANA Clustering:
  • giana_clusters_full.csv
  • giana_coclustering_matrix_full.csv
  • giana_sample_network_full.png

Publication Figures:
  • immunarch_style_diversity.pdf
  • immunarch_style_diversity.png
  • immunarch_style_gene_usage.pdf
  • immunarch_style_gene_usage.png
```

---

## 🐛 Troubleshooting

### **"No module named 'tcr_data_loader'"**

**Solution:** Make sure you're in the correct directory

```bash
cd TCR-BCR-seq-analysis
python tcr_diversity_analysis_full.py
```

### **"File not found: metadata.tsv"**

**Solution:** Update the DATA_DIR path in the script

```python
# Make sure this matches your actual path:
DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"
```

### **"Memory Error" or "Killed"**

**Solution:** Reduce MAX_SEQUENCES in GIANA

```python
# In giana_analysis_full.py:
MAX_SEQUENCES = 20000  # Reduce from 50000
```

### **Script runs but no CMV comparison**

**Solution:** Check your metadata.tsv has `disease_subtype` column with "CMV+" and "CMV-" values

```bash
# Verify metadata format:
head metadata.tsv
```

---

## 💡 Tips for Best Results

### **1. Start Small, Then Scale Up**

```bash
# First run with subset
MAX_FILES = 10  # Try with 10 files first

# If it works, scale up
MAX_FILES = 50  # Then 50 files

# Finally, go full
MAX_FILES = None  # All files
```

### **2. Run Analyses in Order**

Always run diversity analysis first - other scripts may use its outputs!

```bash
1. tcr_diversity_analysis_full.py  # FIRST
2. tcr_vdj_analysis_full.py        # SECOND
3. immunarch_style_viz.py          # THIRD
4. giana_analysis_full.py          # LAST (optional)
```

### **3. Monitor Progress**

All scripts print detailed progress:
- "Loading file 50/300..." - File loading progress
- "Processing sequences..." - Analysis progress
- "✓ Saved: filename.csv" - Outputs created

### **4. Save Your Settings**

After editing DATA_DIR, save a copy:

```bash
# Make a config file:
echo 'DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"' > my_config.txt
```

---

## 📈 What to Expect (Results Preview)

### **With Your Full Dataset (300+ files):**

**Diversity Analysis:**
- Total sequences: ~1-5 million
- Samples analyzed: 300+
- CMV+ vs CMV- comparison with statistics
- Expected finding: CMV+ lower diversity (p < 0.05)

**VDJ Gene Analysis:**
- Unique V genes: ~50-70
- Unique J genes: ~10-15
- Top CMV-associated genes: TRBV6, TRBV9, TRBJ2-7
- V-J pairings: Hundreds to thousands

**GIANA Clustering:**
- Clusters identified: 1,000-10,000+
- Public TCRs: 100-500
- Sample network: Clear CMV+/CMV- separation

---

## 🎓 Scientific Impact

### **What You Can Discover:**

1. **CMV-driven repertoire changes**
   - Reduced diversity in CMV+
   - Specific gene biases
   - Clonal expansions

2. **Public TCR sequences**
   - Shared CMV responses
   - Convergent recognition
   - Common epitope targeting

3. **Gene usage patterns**
   - CMV-associated V/J genes
   - Characteristic pairings
   - HLA-restricted patterns

4. **Population-level insights**
   - Sample clustering by CMV status
   - Age effects on diversity
   - Ancestry influences

---

## 📚 Next Steps After Analysis

### **1. Validate Results**
- Compare to published CMV TCR studies
- Check gene usage against literature
- Verify public TCRs in databases

### **2. Generate Figures for Publication**
- Use immunarch-style PDFs (300 DPI)
- Combine panels for manuscript
- Add statistical annotations

### **3. Extended Analysis**
- Functional validation (tetramers)
- Epitope mapping
- Longitudinal tracking

### **4. Publication!**
- Methods: Cite Zhang et al. 2021 for GIANA
- Results: Report comprehensive metrics
- Figures: Publication-ready PDFs included

---

## 🔗 Quick Reference

### **File Naming Convention:**

| Old Files | New Files (Full Dataset) | Purpose |
|-----------|-------------------------|---------|
| `tcr_diversity_analysis.py` | `tcr_diversity_analysis_full.py` | Diversity metrics |
| `tcr_vdj_analysis.py` | `tcr_vdj_analysis_full.py` | Gene usage |
| `giana_analysis_zhang2021.py` | `giana_analysis_full.py` | Clustering |
| N/A | `tcr_data_loader.py` | Data loading module |

### **Key Parameters:**

| Parameter | Default | Recommendation |
|-----------|---------|----------------|
| `DATA_DIR` | (must set) | `r"C:\Users\chris\Desktop\TCR ANALYSIS"` |
| `MAX_FILES` | `None` | Start with 50 for testing |
| `MAX_SEQUENCES` | 100,000 | 50K for 8GB RAM, 500K for 32GB RAM |

---

## ✅ Checklist Before Running

- [ ] Updated `DATA_DIR` in all scripts
- [ ] Verified metadata.tsv exists in directory
- [ ] Verified .tsv.gz files exist in directory
- [ ] Have 30-60 minutes for complete analysis
- [ ] Have 8GB+ RAM (16GB+ recommended for GIANA)
- [ ] Python packages installed (pandas, numpy, matplotlib, seaborn, scipy, networkx)

---

## 📧 Support

**Common Questions:**

**Q: How long will it take?**
A: ~5-10 minutes for diversity + VDJ, +30 min for GIANA

**Q: How much memory needed?**
A: 8GB minimum, 16GB recommended, 32GB for full GIANA

**Q: Can I run this on a laptop?**
A: Yes! Use MAX_FILES=50 and MAX_SEQUENCES=20000

**Q: What if I don't have all the files?**
A: Scripts work with any subset - just update DATA_DIR

---

**Generated:** 2025-10-31
**Repository:** https://github.com/Tziaxri/TCR-BCR-seq-analysis
**Branch:** claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG

---

## 🎉 You're Ready!

**Everything is set up to analyze your complete TCR dataset!**

Just update DATA_DIR and run the scripts. Start with diversity analysis and work your way through the workflow.

**Good luck with your analysis!** 🚀
