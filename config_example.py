#!/usr/bin/env python3
"""
TCR Analysis Configuration File
================================

Edit this file to set your data directory and analysis parameters.
Then import it in your scripts instead of editing each script individually.

Usage:
------
from config import DATA_DIR, MAX_FILES, MAX_SEQUENCES

Author: Claude Code
Date: 2025-10-31
"""

# =============================================================================
# DATA LOCATION - UPDATE THIS!
# =============================================================================

# For Windows:
DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"

# For Linux/Mac:
# DATA_DIR = "/path/to/your/data"

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

# Maximum files to process (None = all files)
# Recommended: Start with 50 for testing, then set to None for full analysis
MAX_FILES = None  # Change to 50 for testing

# Maximum sequences for GIANA clustering
# Adjust based on available RAM:
#   50,000 - Good for 8GB RAM
#  100,000 - Good for 16GB RAM
#  500,000 - Good for 32GB+ RAM
MAX_SEQUENCES = 100000

# GIANA similarity threshold (from Zhang et al. 2021)
SIMILARITY_THRESHOLD = 0.85

# Statistical significance threshold
ALPHA = 0.05

# Random seed for reproducibility
RANDOM_SEED = 42

# =============================================================================
# OUTPUT SETTINGS
# =============================================================================

# Figure DPI (resolution)
FIGURE_DPI = 300

# Figure format
FIGURE_FORMAT = 'png'  # Options: 'png', 'pdf', 'both'

# Verbose output (detailed progress messages)
VERBOSE = True

# =============================================================================
# FILE PATTERNS
# =============================================================================

# TCR file pattern (glob pattern)
TCR_FILE_PATTERN = "*.tsv.gz"

# Metadata filename
METADATA_FILE = "metadata.tsv"

# Manifest filename (optional)
MANIFEST_FILE = "synapse_metadata_manifest.tsv"

# =============================================================================
# COLUMN NAMES (adjust if your data has different column names)
# =============================================================================

# CDR3 amino acid sequence
CDR3_COLUMNS = ['amino_acid', 'cdr3_aa', 'junction_aa', 'cdr3']

# V gene
V_GENE_COLUMNS = ['v_gene', 'vMaxResolved', 'v_call', 'v_family']

# D gene
D_GENE_COLUMNS = ['d_gene', 'dMaxResolved', 'd_call']

# J gene
J_GENE_COLUMNS = ['j_gene', 'jMaxResolved', 'j_call']

# Clone count
COUNT_COLUMNS = ['duplicate_count', 'templates', 'count']

# Productive sequence indicator
PRODUCTIVE_COLUMN = 'productive'

# Sample/repertoire ID
SAMPLE_ID_COLUMNS = ['repertoire_id', 'sample_id', 'file_id', 'bfi']

# Disease status
DISEASE_COLUMN = 'disease_subtype'

# =============================================================================
# PRINT CONFIGURATION
# =============================================================================

def print_config():
    """Print current configuration settings"""
    print("="*70)
    print("CURRENT CONFIGURATION")
    print("="*70)
    print(f"\n📁 Data Settings:")
    print(f"   Data directory: {DATA_DIR}")
    print(f"   TCR file pattern: {TCR_FILE_PATTERN}")
    print(f"   Max files: {'ALL' if MAX_FILES is None else MAX_FILES}")

    print(f"\n⚙️  Analysis Settings:")
    print(f"   Max sequences (GIANA): {MAX_SEQUENCES:,}")
    print(f"   Similarity threshold: {SIMILARITY_THRESHOLD}")
    print(f"   Random seed: {RANDOM_SEED}")

    print(f"\n📊 Output Settings:")
    print(f"   Figure DPI: {FIGURE_DPI}")
    print(f"   Figure format: {FIGURE_FORMAT}")
    print(f"   Verbose: {VERBOSE}")

    print("\n" + "="*70 + "\n")


if __name__ == "__main__":
    print_config()
