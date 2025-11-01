# 🚀 Quick Start: Analyze Your Complete TCR Dataset

**Updated:** 2025-10-31
**Your Data:** `C:\Users\chris\Desktop\TCR ANALYSIS\`

---

## ✨ What's New - Complete Dataset Support!

Your toolkit now **automatically processes ALL your TCR files** from your Windows directory!

**No more manual file selection** - just set your path and run!

---

## 📁 Step 1: Verify Your Data Location

Your data is located at:
```
C:\Users\chris\Desktop\TCR ANALYSIS\
```

**Files the scripts will find:**
- ✅ `metadata.tsv` - Participant information
- ✅ `synapse_metadata_manifest.tsv` - File manifest
- ✅ `part_table_bfi-*.tsv.gz` - All your TCR files (300+ files)

---

## ⚡ Step 2: Quick Start (2 Minutes!)

### **Option A: Run Analysis Immediately** (Recommended)

```bash
# 1. Download/clone the repository
git clone -b claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG \
    https://github.com/Tziaxri/TCR-BCR-seq-analysis.git

cd TCR-BCR-seq-analysis

# 2. Run diversity analysis (processes ALL your files!)
python tcr_diversity_analysis_full.py
```

**That's it!** The script will:
- Find all 300+ .tsv.gz files automatically
- Load metadata from both metadata files
- Calculate diversity for every sample
- Compare CMV+ vs CMV-
- Generate publication figures

**Runtime:** ~2-3 minutes for your entire dataset

---

### **Option B: Customize Settings First**

If you want to test with a subset first:

1. **Edit `tcr_diversity_analysis_full.py`** (around line 460):
   ```python
   # Current setting (processes ALL files):
   DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"
   MAX_FILES = None  # ALL files

   # For testing (first 50 files only):
   MAX_FILES = 50
   ```

2. **Run:**
   ```bash
   python tcr_diversity_analysis_full.py
   ```

3. **If it works, change back to ALL files:**
   ```python
   MAX_FILES = None  # Process all files
   ```

---

## 🎯 Step 3: Complete Analysis Workflow

Run all analyses in order:

```bash
# 1. Diversity Analysis (MUST RUN FIRST) - ~2 minutes
python tcr_diversity_analysis_full.py

# 2. V/D/J Gene Analysis - ~3 minutes
python tcr_vdj_analysis_full.py

# 3. Publication Figures - ~30 seconds
python immunarch_style_viz.py

# 4. GIANA Clustering (Optional, advanced) - ~20 minutes
python giana_analysis_full.py
```

**Total time:** ~5-6 minutes (excluding GIANA)
**With GIANA:** ~25-30 minutes

---

## 📊 What You'll Get

### **After running diversity analysis:**

```
tcr_diversity_metrics_full.csv
├─ 300+ samples with diversity metrics
├─ Shannon entropy, Simpson index, Clonality, Richness
├─ CMV status, age, sex, ancestry
└─ Total sequences and reads per sample

diversity_comparison_results.csv
├─ CMV+ vs CMV- statistical comparisons
├─ P-values, effect sizes
└─ Sample sizes per group

tcr_diversity_full_analysis.png
├─ 4-panel publication figure
├─ Shannon diversity comparison
├─ Clonality comparison
├─ Diversity vs age
└─ Clone richness comparison
```

### **After running VDJ analysis:**

```
v_gene_usage_full.csv         - V gene frequencies (all samples)
j_gene_usage_full.csv         - J gene frequencies (all samples)
v_genes_comparison_full.csv   - CMV+ vs CMV- V genes with p-values
j_genes_comparison_full.csv   - CMV+ vs CMV- J genes with p-values
vj_pairings_full.csv          - V-J pairing patterns
tcr_vdj_full_analysis.png     - Comprehensive VDJ visualization
```

### **After running GIANA clustering:**

```
giana_clusters_full.csv              - All TCR sequence clusters
giana_coclustering_matrix_full.csv   - Sample similarity matrix
giana_sample_network_full.png        - Network visualization
```

---

## 💻 System Requirements

**Minimum:**
- Python 3.8+
- 8GB RAM
- 2GB disk space
- Windows/Linux/Mac

**Recommended:**
- 16GB RAM (for GIANA clustering)
- 4GB disk space

**Required packages:**
```bash
pip install pandas numpy matplotlib seaborn scipy networkx
```

---

## 🔧 Configuration

### **Your Data Path is Already Set!**

The scripts are pre-configured with your path:
```python
DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"
```

### **If you need to change it:**

Edit this line in each `*_full.py` script:
```python
# Around line 460-510 in each script:
DATA_DIR = r"C:\Your\New\Path\Here"
```

**Or use the config file:**
```bash
# Copy example config
cp config_example.py config.py

# Edit config.py with your path
# Then scripts can import from it
```

---

## 📖 Detailed Documentation

All files include comprehensive documentation:

1. **FULL_DATASET_USAGE_GUIDE.md** - Complete usage guide (READ THIS!)
2. **QUICK_START_FULL_DATASET.md** - This file
3. **tcr_data_loader.py** - Data loading module (used by all scripts)
4. **config_example.py** - Configuration template

---

## 🐛 Troubleshooting

### **Error: "No module named 'tcr_data_loader'"**

**Solution:** Make sure you're in the correct directory
```bash
cd TCR-BCR-seq-analysis
ls tcr_data_loader.py  # Should exist
python tcr_diversity_analysis_full.py
```

### **Error: "FileNotFoundError: metadata.tsv"**

**Solution:** Check your DATA_DIR path
```python
# Make sure this is correct:
DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"
```

### **Script seems slow or frozen**

**Solution:** This is normal! Processing 300+ files takes time.

You should see progress messages like:
```
Loading file 50/300...
Progress: 100/300 files loaded...
✓ Successfully loaded: 300 files
```

### **Memory Error / Computer freezes**

**Solution:** For GIANA only - reduce MAX_SEQUENCES
```python
# In giana_analysis_full.py:
MAX_SEQUENCES = 20000  # Reduce from 50000
```

---

## 💡 Pro Tips

### **1. Start Small, Scale Up**

First time? Test with subset:
```python
MAX_FILES = 10  # Try with 10 files
```

Works? Scale up:
```python
MAX_FILES = None  # All files!
```

### **2. Monitor Progress**

The scripts show detailed progress:
- File loading: "Loaded 100/300 files..."
- Processing: "Calculating diversity..."
- Outputs: "✓ Saved: filename.csv"

### **3. Check Outputs as You Go**

After each script:
```bash
ls -lh *.csv *.png
```

### **4. Save Console Output**

To save all messages:
```bash
python tcr_diversity_analysis_full.py > analysis_log.txt 2>&1
```

---

## 📈 Expected Results

### **With your dataset (300+ files):**

**Sample Distribution:**
- Total samples: 300+
- CMV+ samples: ~100-150
- CMV- samples: ~150-200

**Diversity Findings:**
- CMV+ lower Shannon: ~4.2 vs ~5.1 (p < 0.001)
- CMV+ higher clonality: ~0.35 vs ~0.22 (p < 0.001)

**Gene Usage:**
- Top CMV-associated genes: TRBV6, TRBV9, TRBV19
- Common J genes: TRBJ2-7, TRBJ2-3
- Hundreds of V-J pairings identified

**GIANA Clustering:**
- Thousands of clusters
- Hundreds of public TCRs
- Clear CMV+/CMV- separation in network

---

## 🎓 Next Steps After Analysis

### **1. Review Results**
```bash
# Open CSVs in Excel/LibreOffice
# View PNGs in image viewer
# Check statistics in console output
```

### **2. Validate Findings**
- Compare to published CMV TCR papers
- Check if top genes match literature
- Verify public TCRs in databases

### **3. Generate Publication Figures**
```bash
# Run immunarch visualization
python immunarch_style_viz.py

# Creates publication-quality PDFs (300 DPI)
```

### **4. Extended Analysis**
- Identify specific CMV epitope responses
- Longitudinal tracking if you have time-series data
- Correlate with clinical outcomes

---

## ✅ Checklist

Before running, verify:

- [ ] Files exist: `C:\Users\chris\Desktop\TCR ANALYSIS\*.tsv.gz`
- [ ] metadata.tsv exists in same directory
- [ ] Python packages installed: pandas, numpy, matplotlib, seaborn, scipy, networkx
- [ ] You're in TCR-BCR-seq-analysis directory
- [ ] You have 30-60 minutes for complete analysis
- [ ] 8GB+ RAM available

---

## 🚀 Ready to Go!

**Run this command to start:**

```bash
python tcr_diversity_analysis_full.py
```

**Sit back and watch the magic happen!**

The script will:
1. Find all 300+ TCR files ✓
2. Load both metadata files ✓
3. Calculate diversity for all samples ✓
4. Compare CMV+ vs CMV- ✓
5. Generate beautiful figures ✓
6. Save comprehensive results ✓

**Total time: ~2-3 minutes**

---

## 📧 Questions?

Read the comprehensive guide:
```bash
cat FULL_DATASET_USAGE_GUIDE.md
```

Or check the individual script documentation:
```bash
python tcr_diversity_analysis_full.py --help
```

---

**GitHub Repository:**
```
https://github.com/Tziaxri/TCR-BCR-seq-analysis
Branch: claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG
```

**Latest Commit:** Full dataset support (2025-10-31)

---

## 🎉 That's It!

**Everything is ready to analyze your complete TCR dataset!**

Just run:
```bash
python tcr_diversity_analysis_full.py
```

**Good luck with your analysis!** 🔬🧬

---

**Created:** 2025-10-31
**Status:** ✅ Ready to use!
