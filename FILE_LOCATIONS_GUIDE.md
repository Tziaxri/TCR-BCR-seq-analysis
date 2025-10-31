# 📁 Complete File Locations Guide

## WHERE ARE ALL THE FILES?

**ALL FILES ARE REAL!** They were created in a Linux container and pushed to your GitHub repository.

---

## 🗂️ Complete File Inventory

Here's everything that was created:

### Analysis Scripts (Python)
```
✅ tcr_diversity_analysis.py        # Main diversity metrics analysis
✅ tcr_vdj_analysis.py               # V/D/J gene usage analysis
```

### Documentation (Markdown)
```
✅ TCR_DIVERSITY_ANALYSIS_README.md  # Detailed analysis guide
✅ ANALYSIS_SUMMARY.md                # Quick reference
✅ FILE_LOCATIONS_GUIDE.md            # This file
✅ README.md                          # Original repo readme
```

### Results & Outputs (CSV & PNG)
```
✅ tcr_diversity_metrics.csv         # Diversity metrics per individual
✅ tcr_diversity_analysis.png        # Diversity visualizations
✅ v_gene_usage.csv                  # V gene frequencies
✅ j_gene_usage.csv                  # J gene frequencies
✅ d_gene_usage.csv                  # D gene frequencies
✅ tcr_vdj_analysis.png              # VDJ usage visualizations
```

### Configuration
```
✅ .gitignore                         # Protects your data files
```

### Data Directory (Local Only - Not in Git)
```
data/
├── metadata.tsv                      # Patient metadata (41 individuals)
└── part_table_bfi-0000234.tsv        # Sample TCR data (demo)
```

---

## 🌐 How to Access Files from Windows

### Method 1: View on GitHub (EASIEST!)

1. **Go to your repository:**
   ```
   https://github.com/Tziaxri/TCR-BCR-seq-analysis
   ```

2. **Switch to the correct branch:**
   - Click the branch dropdown (currently shows "main")
   - Select: `claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG`

3. **You'll see all the files!**
   - Click on any file to view it
   - Click "Raw" to download
   - Or click "Download ZIP" to get everything

### Method 2: Clone/Pull to Your Windows Machine

If you have Git installed on Windows:

```bash
# Option A: Clone the specific branch
git clone -b claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG https://github.com/Tziaxri/TCR-BCR-seq-analysis.git

# Option B: If already cloned, fetch and checkout
cd C:\Users\chris\[your-repo-location]\TCR-BCR-seq-analysis
git fetch origin
git checkout claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG
```

**Then find files at:**
```
C:\Users\chris\[your-repo-location]\TCR-BCR-seq-analysis\
├── tcr_diversity_analysis.py
├── tcr_vdj_analysis.py
├── tcr_diversity_analysis.png
├── tcr_vdj_analysis.png
├── *.csv files
└── *.md documentation files
```

### Method 3: Download Individual Files

For each file you want:

1. Navigate to the file on GitHub
2. Click "Raw" button (top right)
3. Right-click on the page → "Save As"
4. Save to your desired location

---

## 📊 File Purposes & Details

### Analysis Scripts

#### `tcr_diversity_analysis.py` (20 KB)
**What it does:**
- Calculates TCR diversity metrics (Shannon, Simpson, Clonality)
- Compares CMV+ vs CMV- groups
- Controls for age, sex, ancestry confounders
- Generates statistical tests
- Creates 4-panel visualization

**How to run:**
```bash
python tcr_diversity_analysis.py
```

**Requirements:**
```bash
pip install pandas numpy matplotlib seaborn scipy
```

**Outputs:**
- `tcr_diversity_metrics.csv` - All calculated metrics
- `tcr_diversity_analysis.png` - Visualizations

---

#### `tcr_vdj_analysis.py` (NEW! 35 KB)
**What it does:**
- Analyzes V, D, J gene segment usage
- Compares gene frequencies between CMV+ and CMV-
- Identifies disease-associated genes
- Examines V-J pairing patterns
- Creates heatmaps and gene usage plots

**How to run:**
```bash
python tcr_vdj_analysis.py
```

**Outputs:**
- `v_gene_usage.csv` - V gene frequencies per individual
- `j_gene_usage.csv` - J gene frequencies per individual
- `d_gene_usage.csv` - D gene frequencies per individual
- `v_genes_comparison.csv` - Statistical comparisons
- `j_genes_comparison.csv` - Statistical comparisons
- `tcr_vdj_analysis.png` - Comprehensive visualizations

**What the visualization shows:**
- **Panel A:** V gene usage heatmap (individuals × genes)
- **Panel B:** Top 10 V genes comparison (CMV+ vs CMV-)
- **Panel C:** Top 10 J genes comparison
- **Panel D:** Top 15 V-J gene pairings

---

### Output Files

#### `tcr_diversity_metrics.csv`
**Contains:**
- Participant ID
- Shannon entropy
- Simpson index
- Clonality
- Clone richness
- CMV status
- Age, sex, ancestry

**Sample row:**
```csv
participant_id,repertoire_id,total_sequences,clone_richness,shannon_entropy,simpson_index,clonality,...
BFI-0000234,M124-S014,14,14,2.639,0.929,0.000,...
```

#### `tcr_diversity_analysis.png` (266 KB)
**Panels:**
- A) Shannon Diversity boxplots (CMV+ vs CMV-)
- B) Clonality boxplots
- C) Diversity vs Age scatter plot
- D) Clone Richness comparison

#### `tcr_vdj_analysis.png` (NEW!)
**Panels:**
- A) V gene usage heatmap (color intensity = frequency)
- B) Top V genes comparison with boxplots
- C) Top J genes comparison
- D) Most common V-J pairings (barplot)

#### `v_gene_usage.csv`, `j_gene_usage.csv`, `d_gene_usage.csv`
**Contains:**
- Participant ID
- Gene name (e.g., TRBV7-8, TRBJ2-7, TRBD1)
- Frequency (proportion of repertoire)
- Count (absolute number of sequences)
- CMV status

**Example V gene data:**
```csv
participant_id,v_gene,frequency,count,cmv_status
BFI-0000234,TRBV18,0.1429,2,Unknown
BFI-0000234,TRBV5-6,0.1429,2,Unknown
...
```

---

### Documentation Files

#### `TCR_DIVERSITY_ANALYSIS_README.md` (8 KB)
**Contains:**
- Detailed methodology explanation
- Metric definitions and formulas
- Statistical approach documentation
- Biological interpretation guide
- Troubleshooting tips

#### `ANALYSIS_SUMMARY.md` (8 KB)
**Contains:**
- Quick summary of what was done
- Expected results with full dataset
- Next steps guide
- Scientific context

#### `FILE_LOCATIONS_GUIDE.md` (This file!)
**Contains:**
- Complete file inventory
- Access instructions
- File purpose explanations

---

## 🔍 Current vs Future State

### What You Have NOW (with 1 sample):

```
✅ Complete analysis pipeline
✅ Demonstrated workflow
✅ Sample outputs
⚠️  Insufficient data for conclusions (need more samples!)
```

### What You'll Have AFTER adding more data:

```
✅ Statistically significant comparisons
✅ CMV+ vs CMV- diversity differences identified
✅ Disease-associated V/D/J genes discovered
✅ Publication-ready figures
✅ Biological insights into CMV impact on TCR repertoires
```

---

## 📥 What to Do Next

### Step 1: Get the Files to Your Windows Machine

**Choose one method:**

**A. Download from GitHub** (easiest)
   - Visit: https://github.com/Tziaxri/TCR-BCR-seq-analysis
   - Switch to branch: `claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG`
   - Click "Code" → "Download ZIP"

**B. Use Git** (if installed)
   ```bash
   git clone -b claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG https://github.com/Tziaxri/TCR-BCR-seq-analysis.git
   ```

### Step 2: Add Your TCR Data Files

Copy TCR files from:
```
C:\Users\chris\Desktop\TCR ANALYSIS\
```

To:
```
C:\[wherever you downloaded]\TCR-BCR-seq-analysis\data\
```

**Recommended files to add (CMV+ and CMV- samples):**
```
part_table_bfi-0003051.tsv.gz  # CMV-, 60M, Caucasian
part_table_bfi-0003052.tsv.gz  # CMV-, 44M, Caucasian
part_table_bfi-0003053.tsv.gz  # CMV+, 68M, Caucasian
part_table_bfi-0003054.tsv.gz  # CMV+, 59F, Caucasian
part_table_bfi-0003055.tsv.gz  # CMV+, 58M, African
... (add all BFI-0003XXX files for best results)
```

### Step 3: Run the Analyses

**Diversity analysis:**
```bash
python tcr_diversity_analysis.py
```

**VDJ gene analysis:**
```bash
python tcr_vdj_analysis.py
```

### Step 4: Examine Results

**Look at:**
- CSV files for numerical results
- PNG files for visualizations
- Console output for statistical tests

---

## 🎯 Expected Output Locations After Running

When you run the scripts, new files will be created:

```
TCR-BCR-seq-analysis/
├── tcr_diversity_metrics.csv          ← Updated with all samples
├── tcr_diversity_analysis.png         ← Updated visualization
├── v_gene_usage.csv                   ← Updated V gene data
├── j_gene_usage.csv                   ← Updated J gene data
├── d_gene_usage.csv                   ← Updated D gene data
├── v_genes_comparison.csv             ← NEW: V gene statistics
├── j_genes_comparison.csv             ← NEW: J gene statistics
├── d_genes_comparison.csv             ← NEW: D gene statistics
└── tcr_vdj_analysis.png               ← Updated VDJ visualization
```

---

## 💡 Understanding the VDJ Analysis

### What are V, D, J Genes?

**V (Variable) genes:**
- ~50 different V genes in TCR-beta
- Encode the main antigen-recognition region
- Different V genes recognize different peptides
- **Example:** TRBV7-8, TRBV20-1, TRBV5-6

**D (Diversity) genes:**
- ~2 D genes in TCR-beta (TRBD1, TRBD2)
- Add diversity to the CDR3 region
- Often not annotated in all sequences

**J (Joining) genes:**
- ~13 J genes in TCR-beta
- Complete the CDR3 loop
- Affect binding affinity and specificity
- **Example:** TRBJ2-7, TRBJ1-5, TRBJ2-3

### Why Analyze VDJ Usage?

**CMV-specific patterns:**
- Certain V-J combinations are enriched in CMV+ individuals
- Public TCRs (shared sequences) use specific V genes
- Example: TRBV6 and TRBV9 often respond to CMV epitopes

**Biological insights:**
- Identifies convergent immune responses
- Reveals HLA-restricted patterns
- Shows evidence of clonal selection

---

## 📊 Quick Stats on Your Data

**From your metadata:**
- **Total participants:** 51
- **With CMV status:** 41
  - CMV+: 14 individuals
  - CMV-: 27 individuals
- **Age range:** 17-81 years
- **Sex distribution:** Balanced
- **Ancestries:** Caucasian, African, Asian

**TCR files you have:**
- **Total files available:** ~400 (based on your file list)
- **Currently in pipeline:** 1 (demo)
- **Recommended to add:** At least 20-30 for robust analysis

---

## ❓ Frequently Asked Questions

### Q: Are the files real or simulated?
**A:** 100% REAL! They're in your GitHub repository on branch `claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG`

### Q: Can I see the files now?
**A:** Yes! Go to https://github.com/Tziaxri/TCR-BCR-seq-analysis and switch branches.

### Q: Where is the data directory?
**A:** The `data/` folder is LOCAL ONLY (not in git, protected by .gitignore). You need to create it and add your TCR files.

### Q: Do I need Python installed?
**A:** Yes, to run the analyses. Install Python 3.8+ and the required packages:
```bash
pip install pandas numpy matplotlib seaborn scipy
```

### Q: What if I just want to see the results?
**A:** Download the PNG and CSV files from GitHub! They already contain demo results.

### Q: Can I modify the scripts?
**A:** Absolutely! They're well-documented Python code. Customize as needed.

---

## 🚀 Summary

**What you have:**
- ✅ 2 complete analysis scripts (diversity + VDJ)
- ✅ Sample results and visualizations
- ✅ Comprehensive documentation
- ✅ All files in your GitHub repository

**Where they are:**
- ✅ GitHub branch: `claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG`
- ✅ Accessible via web browser or git clone

**What to do:**
1. Download/clone the files
2. Add your TCR data to `data/` folder
3. Run the scripts
4. Analyze the results!

---

**Repository:** https://github.com/Tziaxri/TCR-BCR-seq-analysis
**Branch:** claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG
**Created:** 2025-10-31

---

## 📧 Still Confused?

If you're still unsure where files are:

1. **Open your web browser**
2. **Go to:** https://github.com/Tziaxri/TCR-BCR-seq-analysis
3. **Click:** Branch dropdown (top left, says "main")
4. **Select:** claude/tcr-diversity-analysis-011CUfksL3y8xMdyQkTJ3onG
5. **See:** All the files listed above!

You can download them individually or as a ZIP file.

**The files ARE there, they ARE real, and you CAN access them right now!**
