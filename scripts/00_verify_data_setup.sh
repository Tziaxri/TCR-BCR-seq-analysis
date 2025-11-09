#!/bin/bash

# TCR Data Setup Verification Script
# This script checks if your data files are correctly placed before starting analysis

echo "=========================================="
echo "TCR Diversity Analysis - Data Setup Check"
echo "=========================================="
echo ""

# Check raw data directory
echo "1. Checking raw data directory..."
RAW_COUNT=$(find data/raw_data -name "*.tsv.gz" 2>/dev/null | wc -l)
echo "   Found $RAW_COUNT .tsv.gz files in data/raw_data/"

if [ $RAW_COUNT -eq 0 ]; then
    echo "   ⚠️  WARNING: No .tsv.gz files found!"
    echo "   → Please place your 10 selected .tsv.gz files in data/raw_data/"
elif [ $RAW_COUNT -lt 10 ]; then
    echo "   ⚠️  WARNING: Found only $RAW_COUNT files (expected 10)"
    echo "   → Add more .tsv.gz files to reach 10 samples"
elif [ $RAW_COUNT -eq 10 ]; then
    echo "   ✓ Perfect! Found exactly 10 files"
else
    echo "   ℹ️  Found $RAW_COUNT files (more than 10 is fine)"
fi
echo ""

# Check metadata file
echo "2. Checking metadata file..."
if [ -f "data/metadata.tsv" ]; then
    echo "   ✓ metadata.tsv found in data/"

    # Check if it's readable and show first few lines
    LINES=$(wc -l < data/metadata.tsv)
    echo "   → File has $LINES lines"

    echo ""
    echo "   First 3 lines of metadata:"
    head -n 3 data/metadata.tsv | cat -A

elif [ -f "data/metadata.csv" ]; then
    echo "   ⚠️  Found metadata.csv instead of metadata.tsv"
    echo "   → Please rename to metadata.tsv or confirm it's tab-delimited"
else
    echo "   ⚠️  WARNING: metadata.tsv not found in data/"
    echo "   → Please place metadata.tsv in the data/ directory"
fi
echo ""

# List the raw files
echo "3. List of raw data files:"
if [ $RAW_COUNT -gt 0 ]; then
    find data/raw_data -name "*.tsv.gz" -type f | sort
else
    echo "   (none found)"
fi
echo ""

# Check directory structure
echo "4. Checking directory structure..."
DIRS=("data/raw_data" "data/processed_data" "scripts" "results/figures" "results/tables")
ALL_OK=true

for dir in "${DIRS[@]}"; do
    if [ -d "$dir" ]; then
        echo "   ✓ $dir exists"
    else
        echo "   ✗ $dir missing"
        ALL_OK=false
    fi
done
echo ""

# Final summary
echo "=========================================="
echo "Setup Summary:"
echo "=========================================="
if [ $RAW_COUNT -ge 10 ] && [ -f "data/metadata.tsv" ] && [ "$ALL_OK" = true ]; then
    echo "✓ All checks passed! You're ready to run preprocessing."
    echo ""
    echo "Next step: Run the preprocessing script"
    echo "  Rscript scripts/01_preprocess_tcr_data.R"
elif [ $RAW_COUNT -ge 5 ] && [ -f "data/metadata.tsv" ]; then
    echo "⚠️  Partial setup complete. You can proceed but consider adding more samples."
else
    echo "⚠️  Setup incomplete. Please:"
    echo "   1. Place 10 .tsv.gz files in data/raw_data/"
    echo "   2. Place metadata.tsv in data/"
    echo "   3. Refer to data/DATA_SETUP_GUIDE.md for details"
fi
echo ""
