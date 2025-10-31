#!/usr/bin/env python3
"""
TCR Analysis Utility Functions
================================

Comprehensive collection of utility functions for TCR repertoire analysis.

Categories:
1. Sequence Analysis - CDR3 processing, distance calculations
2. Diversity Metrics - Shannon, Simpson, clonality
3. Data Conversion - Format conversions for different tools
4. Visualization - Plotting utilities
5. Statistical Tests - Comparison functions

Author: Claude Code
Date: 2025-10-31
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.spatial.distance import hamming
import re
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

################################################################################
# SEQUENCE ANALYSIS FUNCTIONS
################################################################################

def parse_vdj_gene(gene_call):
    """
    Parse V/D/J gene call to extract gene family and allele.

    Args:
        gene_call (str): Gene call (e.g., 'TRBV7-8*01')

    Returns:
        dict: {'family': 'TRBV7-8', 'allele': '01'}
    """
    if pd.isna(gene_call) or gene_call == '':
        return {'family': None, 'allele': None}

    parts = str(gene_call).split('*')
    return {
        'family': parts[0] if len(parts) > 0 else None,
        'allele': parts[1] if len(parts) > 1 else None
    }

def calculate_cdr3_properties(cdr3_aa):
    """
    Calculate properties of CDR3 amino acid sequence.

    Args:
        cdr3_aa (str): CDR3 amino acid sequence

    Returns:
        dict: Properties including length, hydrophobicity, charge, etc.
    """
    if pd.isna(cdr3_aa) or not isinstance(cdr3_aa, str):
        return None

    # Amino acid properties
    hydrophobic = set('AILMFVPGW')
    charged = set('DEKR')
    positive = set('KR')
    negative = set('DE')
    polar = set('STNQCY')
    aromatic = set('FWY')

    length = len(cdr3_aa)
    n_hydrophobic = sum(1 for aa in cdr3_aa if aa in hydrophobic)
    n_charged = sum(1 for aa in cdr3_aa if aa in charged)
    n_positive = sum(1 for aa in cdr3_aa if aa in positive)
    n_negative = sum(1 for aa in cdr3_aa if aa in negative)
    n_polar = sum(1 for aa in cdr3_aa if aa in polar)
    n_aromatic = sum(1 for aa in cdr3_aa if aa in aromatic)

    return {
        'length': length,
        'hydrophobic_fraction': n_hydrophobic / length if length > 0 else 0,
        'charged_fraction': n_charged / length if length > 0 else 0,
        'positive_fraction': n_positive / length if length > 0 else 0,
        'negative_fraction': n_negative / length if length > 0 else 0,
        'polar_fraction': n_polar / length if length > 0 else 0,
        'aromatic_fraction': n_aromatic / length if length > 0 else 0,
        'net_charge': n_positive - n_negative
    }

def extract_cdr3_motifs(cdr3_aa, k=3):
    """
    Extract all k-mers from CDR3 sequence.

    Args:
        cdr3_aa (str): CDR3 amino acid sequence
        k (int): Motif length (default: 3)

    Returns:
        list: All k-mers in the sequence
    """
    if pd.isna(cdr3_aa) or not isinstance(cdr3_aa, str) or len(cdr3_aa) < k:
        return []

    return [cdr3_aa[i:i+k] for i in range(len(cdr3_aa) - k + 1)]

def hamming_distance(seq1, seq2):
    """Calculate Hamming distance between two sequences of equal length."""
    if len(seq1) != len(seq2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def levenshtein_distance(seq1, seq2):
    """Calculate Levenshtein (edit) distance between two sequences."""
    if len(seq1) < len(seq2):
        return levenshtein_distance(seq2, seq1)

    if len(seq2) == 0:
        return len(seq1)

    previous_row = range(len(seq2) + 1)
    for i, c1 in enumerate(seq1):
        current_row = [i + 1]
        for j, c2 in enumerate(seq2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

def is_productive_sequence(vj_in_frame, stop_codon, productive):
    """
    Check if TCR sequence is productive.

    Args:
        vj_in_frame: V-J in-frame indicator ('t' or 'f')
        stop_codon: Stop codon indicator ('t' or 'f')
        productive: Productive indicator ('t' or 'f')

    Returns:
        bool: True if productive
    """
    return (str(productive).lower() == 't' and
            str(stop_codon).lower() == 'f' and
            str(vj_in_frame).lower() == 't')

################################################################################
# DIVERSITY METRICS
################################################################################

def calculate_shannon_entropy(frequencies):
    """
    Calculate Shannon entropy from frequency distribution.

    Args:
        frequencies (array-like): Clone frequencies (must sum to 1)

    Returns:
        float: Shannon entropy
    """
    frequencies = np.array(frequencies)
    frequencies = frequencies[frequencies > 0]  # Remove zeros

    if len(frequencies) == 0:
        return 0.0

    return -np.sum(frequencies * np.log(frequencies))

def calculate_simpson_index(frequencies):
    """
    Calculate Simpson diversity index.

    Args:
        frequencies (array-like): Clone frequencies (must sum to 1)

    Returns:
        float: Simpson index (1 - sum of squared frequencies)
    """
    frequencies = np.array(frequencies)
    return 1 - np.sum(frequencies ** 2)

def calculate_inverse_simpson(frequencies):
    """
    Calculate inverse Simpson index (Hill number).

    Args:
        frequencies (array-like): Clone frequencies

    Returns:
        float: Inverse Simpson index
    """
    frequencies = np.array(frequencies)
    sum_squared = np.sum(frequencies ** 2)
    return 1 / sum_squared if sum_squared > 0 else 0

def calculate_clonality(frequencies):
    """
    Calculate clonality (1 - normalized Shannon entropy).

    Args:
        frequencies (array-like): Clone frequencies

    Returns:
        float: Clonality (0 = maximally diverse, 1 = single clone)
    """
    shannon = calculate_shannon_entropy(frequencies)
    max_entropy = np.log(len(frequencies))

    if max_entropy == 0:
        return 0.0

    return 1 - (shannon / max_entropy)

def calculate_chao1(counts):
    """
    Calculate Chao1 richness estimator.

    Args:
        counts (array-like): Clone counts (not frequencies)

    Returns:
        float: Chao1 estimate of total richness
    """
    counts = np.array(counts)
    n = len(counts)

    f1 = np.sum(counts == 1)  # Singletons
    f2 = np.sum(counts == 2)  # Doubletons

    if f2 == 0:
        return n + (f1 * (f1 - 1)) / 2
    else:
        return n + (f1 ** 2) / (2 * f2)

def calculate_repertoire_diversity(clone_counts):
    """
    Calculate comprehensive diversity metrics from clone counts.

    Args:
        clone_counts (dict or Counter): Clone ID -> count

    Returns:
        dict: All diversity metrics
    """
    if isinstance(clone_counts, dict):
        clone_counts = Counter(clone_counts)

    counts = np.array(list(clone_counts.values()))
    total = counts.sum()
    frequencies = counts / total

    return {
        'richness': len(counts),
        'shannon': calculate_shannon_entropy(frequencies),
        'simpson': calculate_simpson_index(frequencies),
        'inv_simpson': calculate_inverse_simpson(frequencies),
        'clonality': calculate_clonality(frequencies),
        'chao1': calculate_chao1(counts),
        'top_clone_freq': frequencies.max(),
        'top10_clone_freq': frequencies[np.argsort(frequencies)[-10:]].sum() if len(frequencies) >= 10 else frequencies.sum()
    }

################################################################################
# DATA CONVERSION FUNCTIONS
################################################################################

def convert_to_immunarch_format(df):
    """
    Convert TCR dataframe to immunarch format.

    Args:
        df (DataFrame): TCR data with standard columns

    Returns:
        DataFrame: Immunarch-formatted data
    """
    # Group by unique clonotypes
    grouped = df.groupby(['cdr3_aa', 'v_call', 'j_call']).size().reset_index(name='Clones')

    # Calculate proportions
    grouped['Proportion'] = grouped['Clones'] / grouped['Clones'].sum()

    # Rename columns for immunarch
    grouped = grouped.rename(columns={
        'cdr3_aa': 'CDR3.aa',
        'v_call': 'V.name',
        'j_call': 'J.name'
    })

    # Sort by clones (most abundant first)
    grouped = grouped.sort_values('Clones', ascending=False)

    return grouped[['Clones', 'Proportion', 'CDR3.aa', 'V.name', 'J.name']]

def convert_to_vdjtools_format(df):
    """
    Convert TCR dataframe to VDJtools format.

    Args:
        df (DataFrame): TCR data

    Returns:
        DataFrame: VDJtools-formatted data
    """
    vdjtools = df.copy()

    # Required columns: count, frequency, CDR3nt, CDR3aa, V, D, J
    vdjtools = vdjtools.rename(columns={
        'sequence': 'CDR3nt',
        'cdr3_aa': 'CDR3aa',
        'v_call': 'V',
        'd_call': 'D',
        'j_call': 'J'
    })

    # Calculate counts and frequencies
    if 'clone_id' in vdjtools.columns:
        clone_counts = vdjtools.groupby('clone_id').size()
        vdjtools['count'] = vdjtools['clone_id'].map(clone_counts)
    else:
        vdjtools['count'] = 1

    total = vdjtools['count'].sum()
    vdjtools['frequency'] = vdjtools['count'] / total

    return vdjtools[['count', 'frequency', 'CDR3nt', 'CDR3aa', 'V', 'D', 'J']]

def convert_to_airr_format(df):
    """
    Convert TCR dataframe to AIRR format.

    Args:
        df (DataFrame): TCR data

    Returns:
        DataFrame: AIRR-formatted data
    """
    airr = df.copy()

    # AIRR standard column names
    column_mapping = {
        'sequence': 'sequence',
        'sequence_id': 'sequence_id',
        'productive': 'productive',
        'v_call': 'v_call',
        'd_call': 'd_call',
        'j_call': 'j_call',
        'junction': 'junction',
        'junction_aa': 'junction_aa'
    }

    airr = airr.rename(columns={k: v for k, v in column_mapping.items() if k in airr.columns})

    return airr

################################################################################
# STATISTICAL FUNCTIONS
################################################################################

def compare_diversity_between_groups(group1_metrics, group2_metrics, metric='shannon'):
    """
    Compare diversity metric between two groups.

    Args:
        group1_metrics (list): Diversity values for group 1
        group2_metrics (list): Diversity values for group 2
        metric (str): Metric name (for reporting)

    Returns:
        dict: Statistical test results
    """
    group1 = np.array([x for x in group1_metrics if not np.isnan(x)])
    group2 = np.array([x for x in group2_metrics if not np.isnan(x)])

    if len(group1) == 0 or len(group2) == 0:
        return {
            'metric': metric,
            'n1': len(group1),
            'n2': len(group2),
            'mean1': np.nan,
            'mean2': np.nan,
            'p_value': np.nan,
            'significant': False
        }

    # Mann-Whitney U test (non-parametric)
    statistic, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')

    # Effect size (rank-biserial correlation)
    n1, n2 = len(group1), len(group2)
    effect_size = 1 - (2 * statistic) / (n1 * n2)

    return {
        'metric': metric,
        'n1': n1,
        'n2': n2,
        'mean1': np.mean(group1),
        'mean2': np.mean(group2),
        'median1': np.median(group1),
        'median2': np.median(group2),
        'statistic': statistic,
        'p_value': p_value,
        'effect_size': effect_size,
        'significant': p_value < 0.05
    }

def compare_gene_usage(group1_freqs, group2_freqs, gene_name):
    """
    Compare gene usage frequency between two groups.

    Args:
        group1_freqs (list): Gene frequencies in group 1
        group2_freqs (list): Gene frequencies in group 2
        gene_name (str): Gene name

    Returns:
        dict: Comparison results
    """
    group1 = np.array([x for x in group1_freqs if not np.isnan(x)])
    group2 = np.array([x for x in group2_freqs if not np.isnan(x)])

    if len(group1) == 0 or len(group2) == 0:
        return None

    statistic, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')

    mean1 = np.mean(group1)
    mean2 = np.mean(group2)
    fold_change = mean1 / mean2 if mean2 > 0 else np.inf

    return {
        'gene': gene_name,
        'mean_group1': mean1,
        'mean_group2': mean2,
        'fold_change': fold_change,
        'p_value': p_value,
        'significant': p_value < 0.05,
        'direction': 'enriched_g1' if fold_change > 1 else 'enriched_g2'
    }

def calculate_repertoire_overlap(rep1_cdr3s, rep2_cdr3s):
    """
    Calculate overlap between two repertoires.

    Args:
        rep1_cdr3s (set or list): CDR3 sequences from repertoire 1
        rep2_cdr3s (set or list): CDR3 sequences from repertoire 2

    Returns:
        dict: Overlap metrics
    """
    set1 = set(rep1_cdr3s)
    set2 = set(rep2_cdr3s)

    intersection = set1 & set2
    union = set1 | set2

    jaccard = len(intersection) / len(union) if len(union) > 0 else 0
    overlap_coef = len(intersection) / min(len(set1), len(set2)) if min(len(set1), len(set2)) > 0 else 0

    return {
        'n_unique_rep1': len(set1),
        'n_unique_rep2': len(set2),
        'n_shared': len(intersection),
        'jaccard_index': jaccard,
        'overlap_coefficient': overlap_coef,
        'percent_shared_rep1': len(intersection) / len(set1) * 100 if len(set1) > 0 else 0,
        'percent_shared_rep2': len(intersection) / len(set2) * 100 if len(set2) > 0 else 0
    }

################################################################################
# UTILITY FUNCTIONS
################################################################################

def load_tcr_data(file_path, productive_only=True):
    """
    Load TCR data from file with automatic format detection.

    Args:
        file_path (str): Path to TCR data file
        productive_only (bool): Filter to productive sequences only

    Returns:
        DataFrame: Loaded TCR data
    """
    df = pd.read_csv(file_path, sep='\t')

    if productive_only and 'productive' in df.columns:
        df = df[df['productive'] == 't']

    return df

def summarize_repertoire(df):
    """
    Generate summary statistics for a TCR repertoire.

    Args:
        df (DataFrame): TCR data

    Returns:
        dict: Summary statistics
    """
    summary = {
        'total_sequences': len(df),
        'unique_cdr3': df['cdr3_aa'].nunique() if 'cdr3_aa' in df.columns else 0,
        'unique_clones': df['clone_id'].nunique() if 'clone_id' in df.columns else 0
    }

    if 'cdr3_aa' in df.columns:
        cdr3_lengths = df['cdr3_aa'].dropna().apply(len)
        summary['mean_cdr3_length'] = cdr3_lengths.mean()
        summary['median_cdr3_length'] = cdr3_lengths.median()

    if 'v_call' in df.columns:
        summary['unique_v_genes'] = df['v_call'].str.split('*').str[0].nunique()

    if 'j_call' in df.columns:
        summary['unique_j_genes'] = df['j_call'].str.split('*').str[0].nunique()

    if 'clone_id' in df.columns:
        clone_counts = df['clone_id'].value_counts()
        summary.update(calculate_repertoire_diversity(clone_counts))

    return summary

def filter_productive_sequences(df):
    """Filter dataframe to productive sequences only."""
    if 'productive' in df.columns:
        return df[df['productive'] == 't'].copy()
    return df.copy()

def add_gene_families(df):
    """Add gene family columns (without alleles)."""
    if 'v_call' in df.columns:
        df['v_gene'] = df['v_call'].str.split('*').str[0]
    if 'd_call' in df.columns:
        df['d_gene'] = df['d_call'].str.split('*').str[0]
    if 'j_call' in df.columns:
        df['j_gene'] = df['j_call'].str.split('*').str[0]
    return df

################################################################################
# EXAMPLE USAGE
################################################################################

if __name__ == '__main__':
    print("TCR Analysis Utilities")
    print("=" * 70)
    print("\nAvailable functions:")
    print("\nSequence Analysis:")
    print("  - parse_vdj_gene()")
    print("  - calculate_cdr3_properties()")
    print("  - extract_cdr3_motifs()")
    print("  - hamming_distance(), levenshtein_distance()")

    print("\nDiversity Metrics:")
    print("  - calculate_shannon_entropy()")
    print("  - calculate_simpson_index()")
    print("  - calculate_clonality()")
    print("  - calculate_repertoire_diversity()")

    print("\nData Conversion:")
    print("  - convert_to_immunarch_format()")
    print("  - convert_to_vdjtools_format()")
    print("  - convert_to_airr_format()")

    print("\nStatistical Functions:")
    print("  - compare_diversity_between_groups()")
    print("  - compare_gene_usage()")
    print("  - calculate_repertoire_overlap()")

    print("\nUtility Functions:")
    print("  - load_tcr_data()")
    print("  - summarize_repertoire()")
    print("  - filter_productive_sequences()")
    print("=" * 70)

    # Example usage
    print("\nExample: Calculate diversity from clone counts")
    example_counts = Counter({1: 100, 2: 50, 3: 25, 4: 10, 5: 5})
    diversity = calculate_repertoire_diversity(example_counts)
    print(f"\nRichness: {diversity['richness']}")
    print(f"Shannon: {diversity['shannon']:.3f}")
    print(f"Simpson: {diversity['simpson']:.3f}")
    print(f"Clonality: {diversity['clonality']:.3f}")
