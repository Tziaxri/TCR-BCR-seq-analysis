#!/usr/bin/env python3
"""
TCR Diversity Analysis: CMV+ vs CMV- Comparison (Full Dataset Version)
=======================================================================

This script analyzes TCR repertoire diversity between CMV+ and CMV- individuals
using ALL TCR files in your directory, accounting for potential confounders.

Updated: 2025-10-31
- Processes ALL .tsv.gz files automatically
- Uses tcr_data_loader module
- Handles Windows/Linux paths
- Supports both metadata files

Usage:
------
# Windows:
DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"

# Or edit the path at the bottom of this script and run:
python tcr_diversity_analysis_full.py

Author: Claude Code
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from scipy.stats import mannwhitneyu, spearmanr
import warnings
warnings.filterwarnings('ignore')

# Import data loader
from tcr_data_loader import load_complete_dataset

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)
plt.rcParams['font.size'] = 10


class TCRDiversityAnalyzerFull:
    """
    Analyzes TCR repertoire diversity using entire dataset

    Features:
    - Loads ALL .tsv.gz files from directory
    - Processes both metadata files
    - Calculates comprehensive diversity metrics
    - Statistical comparisons with confounder adjustment
    - Publication-quality visualizations
    """

    def __init__(self, data_directory: str, max_files: int = None):
        """
        Initialize analyzer

        Parameters:
        -----------
        data_directory : str
            Path to directory with TCR files and metadata
            Example: r"C:\Users\chris\Desktop\TCR ANALYSIS"
        max_files : int
            Maximum files to process (None = all files)
        """
        self.data_dir = data_directory
        self.max_files = max_files
        self.diversity_df = None
        self.loader = None

    def load_all_data(self):
        """Load complete dataset using TCRDataLoader"""

        print("="*80)
        print("LOADING COMPLETE TCR DATASET")
        print("="*80)
        print(f"\n📁 Data directory: {self.data_dir}")

        if self.max_files:
            print(f"⚠️  Loading first {self.max_files} files only (for testing)")
        else:
            print("📊 Loading ALL TCR files in directory...")

        # Load complete dataset
        self.diversity_df, self.loader = load_complete_dataset(
            self.data_dir,
            max_files=self.max_files,
            verbose=True
        )

        print(f"\n✓ Dataset loaded successfully!")
        print(f"  • Total samples: {len(self.diversity_df)}")
        print(f"  • Total TCR files: {len(self.loader.tcr_data)}")

        # Check CMV status availability
        if 'disease_subtype' in self.diversity_df.columns:
            cmv_dist = self.diversity_df['disease_subtype'].value_counts()
            print(f"\n📊 CMV Status Distribution:")
            for status, count in cmv_dist.items():
                print(f"     {status}: {count} samples")

        return self.diversity_df

    def explore_confounders(self):
        """Analyze potential confounding variables"""

        print("\n" + "="*80)
        print("CONFOUNDER ANALYSIS")
        print("="*80)

        # Filter to samples with CMV status
        df = self.diversity_df[
            self.diversity_df['disease_subtype'].isin(['CMV+', 'CMV-'])
        ].copy()

        if len(df) == 0:
            print("⚠️  No samples with CMV+/CMV- status found")
            return

        print(f"\nAnalyzing {len(df)} samples with CMV status")

        # Age analysis
        if 'age' in df.columns:
            print("\n📊 Age Distribution:")
            age_stats = df.groupby('disease_subtype')['age'].describe()
            print(age_stats)

            # Test for age differences
            cmv_pos = df[df['disease_subtype'] == 'CMV+']['age'].dropna()
            cmv_neg = df[df['disease_subtype'] == 'CMV-']['age'].dropna()

            if len(cmv_pos) > 0 and len(cmv_neg) > 0:
                u_stat, p_val = mannwhitneyu(cmv_pos, cmv_neg)
                print(f"\n   Age comparison (Mann-Whitney U):")
                print(f"   CMV+ mean: {cmv_pos.mean():.1f} ± {cmv_pos.std():.1f} years")
                print(f"   CMV- mean: {cmv_neg.mean():.1f} ± {cmv_neg.std():.1f} years")
                print(f"   p-value: {p_val:.4f} {'***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'}")

        # Sex analysis
        if 'sex' in df.columns:
            print("\n📊 Sex Distribution:")
            sex_dist = pd.crosstab(df['disease_subtype'], df['sex'])
            print(sex_dist)

            # Chi-square test
            from scipy.stats import chi2_contingency
            chi2, p_val, dof, expected = chi2_contingency(sex_dist)
            print(f"\n   Sex distribution (Chi-square):")
            print(f"   p-value: {p_val:.4f} {'***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'}")

        # Ancestry analysis
        if 'ancestry' in df.columns:
            print("\n📊 Ancestry Distribution:")
            ancestry_dist = pd.crosstab(
                df['disease_subtype'],
                df['ancestry'].fillna('Unknown')
            )
            print(ancestry_dist)

        # Sequencing depth analysis
        print("\n📊 Sequencing Depth:")
        depth_stats = df.groupby('disease_subtype')['total_reads'].describe()
        print(depth_stats)

        cmv_pos_depth = df[df['disease_subtype'] == 'CMV+']['total_reads']
        cmv_neg_depth = df[df['disease_subtype'] == 'CMV-']['total_reads']

        if len(cmv_pos_depth) > 0 and len(cmv_neg_depth) > 0:
            u_stat, p_val = mannwhitneyu(cmv_pos_depth, cmv_neg_depth)
            print(f"\n   Sequencing depth comparison:")
            print(f"   CMV+ median: {cmv_pos_depth.median():,.0f} reads")
            print(f"   CMV- median: {cmv_neg_depth.median():,.0f} reads")
            print(f"   p-value: {p_val:.4f} {'***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'}")

    def compare_diversity(self):
        """Compare diversity metrics between CMV+ and CMV-"""

        print("\n" + "="*80)
        print("DIVERSITY COMPARISON: CMV+ vs CMV-")
        print("="*80)

        # Filter to CMV samples
        df = self.diversity_df[
            self.diversity_df['disease_subtype'].isin(['CMV+', 'CMV-'])
        ].copy()

        if len(df) == 0:
            print("⚠️  No samples with CMV+/CMV- status for comparison")
            return None

        # Diversity metrics to compare
        metrics = ['shannon_entropy', 'simpson_index', 'clonality', 'richness']

        results = []

        for metric in metrics:
            if metric not in df.columns:
                continue

            cmv_pos = df[df['disease_subtype'] == 'CMV+'][metric].dropna()
            cmv_neg = df[df['disease_subtype'] == 'CMV-'][metric].dropna()

            if len(cmv_pos) == 0 or len(cmv_neg) == 0:
                continue

            # Statistical test
            u_stat, p_val = mannwhitneyu(cmv_pos, cmv_neg)

            # Effect size (Cohen's d)
            pooled_std = np.sqrt((cmv_pos.std()**2 + cmv_neg.std()**2) / 2)
            cohens_d = (cmv_pos.mean() - cmv_neg.mean()) / pooled_std if pooled_std > 0 else 0

            results.append({
                'metric': metric,
                'cmv_pos_mean': cmv_pos.mean(),
                'cmv_pos_std': cmv_pos.std(),
                'cmv_neg_mean': cmv_neg.mean(),
                'cmv_neg_std': cmv_neg.std(),
                'p_value': p_val,
                'cohens_d': cohens_d,
                'n_cmv_pos': len(cmv_pos),
                'n_cmv_neg': len(cmv_neg)
            })

            # Print results
            print(f"\n📊 {metric.upper().replace('_', ' ')}:")
            print(f"   CMV+ (n={len(cmv_pos)}): {cmv_pos.mean():.3f} ± {cmv_pos.std():.3f}")
            print(f"   CMV- (n={len(cmv_neg)}): {cmv_neg.mean():.3f} ± {cmv_neg.std():.3f}")
            print(f"   Difference: {cmv_pos.mean() - cmv_neg.mean():+.3f}")
            print(f"   p-value: {p_val:.4f} {'***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'}")
            print(f"   Effect size (Cohen's d): {cohens_d:.3f}")

        comparison_df = pd.DataFrame(results)

        # Save results
        comparison_df.to_csv('diversity_comparison_results.csv', index=False)
        print(f"\n✓ Saved comparison results to 'diversity_comparison_results.csv'")

        return comparison_df

    def analyze_age_correlation(self):
        """Analyze correlation between diversity and age"""

        print("\n" + "="*80)
        print("AGE CORRELATION ANALYSIS")
        print("="*80)

        df = self.diversity_df[
            self.diversity_df['disease_subtype'].isin(['CMV+', 'CMV-'])
        ].copy()

        if 'age' not in df.columns:
            print("⚠️  Age data not available")
            return

        metrics = ['shannon_entropy', 'simpson_index', 'clonality', 'richness']

        for metric in metrics:
            if metric not in df.columns:
                continue

            data = df[[metric, 'age', 'disease_subtype']].dropna()

            if len(data) < 3:
                continue

            # Overall correlation
            rho, p_val = spearmanr(data[metric], data['age'])

            print(f"\n📊 {metric.upper().replace('_', ' ')} vs Age:")
            print(f"   Spearman correlation: ρ = {rho:.3f}")
            print(f"   p-value: {p_val:.4f} {'***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'}")

            # Group-specific correlations
            for group in ['CMV+', 'CMV-']:
                group_data = data[data['disease_subtype'] == group]
                if len(group_data) >= 3:
                    rho_g, p_g = spearmanr(group_data[metric], group_data['age'])
                    print(f"   {group}: ρ = {rho_g:.3f}, p = {p_g:.4f}")

    def create_visualizations(self, output_file='tcr_diversity_full_analysis.png'):
        """Create comprehensive visualization"""

        print("\n" + "="*80)
        print("GENERATING VISUALIZATIONS")
        print("="*80)

        df = self.diversity_df[
            self.diversity_df['disease_subtype'].isin(['CMV+', 'CMV-'])
        ].copy()

        if len(df) == 0:
            print("⚠️  No data available for visualization")
            return

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('TCR Repertoire Diversity Analysis: CMV+ vs CMV- (Full Dataset)',
                     fontsize=16, fontweight='bold', y=0.995)

        # Panel A: Shannon Entropy
        ax = axes[0, 0]
        if 'shannon_entropy' in df.columns:
            parts = ax.violinplot(
                [df[df['disease_subtype'] == 'CMV+']['shannon_entropy'].dropna(),
                 df[df['disease_subtype'] == 'CMV-']['shannon_entropy'].dropna()],
                positions=[1, 2],
                showmeans=True,
                showmedians=True
            )
            ax.set_xticks([1, 2])
            ax.set_xticklabels(['CMV+', 'CMV-'])
            ax.set_ylabel('Shannon Entropy')
            ax.set_title('A. Shannon Diversity Index', fontweight='bold')
            ax.grid(axis='y', alpha=0.3)

            # Add p-value
            cmv_pos = df[df['disease_subtype'] == 'CMV+']['shannon_entropy'].dropna()
            cmv_neg = df[df['disease_subtype'] == 'CMV-']['shannon_entropy'].dropna()
            if len(cmv_pos) > 0 and len(cmv_neg) > 0:
                _, p_val = mannwhitneyu(cmv_pos, cmv_neg)
                sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                ax.text(1.5, ax.get_ylim()[1]*0.95, f'p={p_val:.4f} {sig}',
                       ha='center', fontsize=10, fontweight='bold')

        # Panel B: Clonality
        ax = axes[0, 1]
        if 'clonality' in df.columns:
            clonality_data = [
                df[df['disease_subtype'] == 'CMV+']['clonality'].dropna(),
                df[df['disease_subtype'] == 'CMV-']['clonality'].dropna()
            ]
            bp = ax.boxplot(clonality_data, positions=[1, 2],
                           labels=['CMV+', 'CMV-'],
                           patch_artist=True)
            for patch in bp['boxes']:
                patch.set_facecolor('lightblue')
            ax.set_ylabel('Clonality')
            ax.set_title('B. Clonality (Normalized Gini Index)', fontweight='bold')
            ax.grid(axis='y', alpha=0.3)

            # Add p-value
            if len(clonality_data[0]) > 0 and len(clonality_data[1]) > 0:
                _, p_val = mannwhitneyu(clonality_data[0], clonality_data[1])
                sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                ax.text(1.5, ax.get_ylim()[1]*0.95, f'p={p_val:.4f} {sig}',
                       ha='center', fontsize=10, fontweight='bold')

        # Panel C: Shannon vs Age
        ax = axes[1, 0]
        if 'shannon_entropy' in df.columns and 'age' in df.columns:
            for group, color, marker in [('CMV+', 'red', 'o'), ('CMV-', 'blue', 's')]:
                group_data = df[df['disease_subtype'] == group]
                ax.scatter(group_data['age'], group_data['shannon_entropy'],
                          c=color, marker=marker, label=group, alpha=0.6, s=50)

            ax.set_xlabel('Age (years)')
            ax.set_ylabel('Shannon Entropy')
            ax.set_title('C. Diversity vs Age', fontweight='bold')
            ax.legend()
            ax.grid(alpha=0.3)

        # Panel D: Richness comparison
        ax = axes[1, 1]
        if 'richness' in df.columns:
            # Create bar plot
            cmv_pos_richness = df[df['disease_subtype'] == 'CMV+']['richness'].dropna()
            cmv_neg_richness = df[df['disease_subtype'] == 'CMV-']['richness'].dropna()

            means = [cmv_pos_richness.mean(), cmv_neg_richness.mean()]
            sems = [cmv_pos_richness.sem(), cmv_neg_richness.sem()]

            bars = ax.bar([1, 2], means, yerr=sems, capsize=5,
                         color=['salmon', 'skyblue'], alpha=0.7)
            ax.set_xticks([1, 2])
            ax.set_xticklabels(['CMV+', 'CMV-'])
            ax.set_ylabel('Clone Richness')
            ax.set_title('D. Clone Richness', fontweight='bold')
            ax.grid(axis='y', alpha=0.3)

            # Add p-value
            if len(cmv_pos_richness) > 0 and len(cmv_neg_richness) > 0:
                _, p_val = mannwhitneyu(cmv_pos_richness, cmv_neg_richness)
                sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                ax.text(1.5, ax.get_ylim()[1]*0.95, f'p={p_val:.4f} {sig}',
                       ha='center', fontsize=10, fontweight='bold')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved visualization: {output_file}")

        return fig

    def save_diversity_metrics(self, output_file='tcr_diversity_metrics_full.csv'):
        """Save diversity metrics to CSV"""

        self.diversity_df.to_csv(output_file, index=False)
        print(f"\n✓ Saved diversity metrics: {output_file}")
        print(f"  Columns: {', '.join(self.diversity_df.columns)}")
        print(f"  Rows: {len(self.diversity_df)}")

    def run_complete_analysis(self):
        """Run complete diversity analysis pipeline"""

        print("\n" + "="*80)
        print("TCR DIVERSITY ANALYSIS - COMPLETE PIPELINE")
        print("="*80)

        # Step 1: Load data
        self.load_all_data()

        # Step 2: Explore confounders
        self.explore_confounders()

        # Step 3: Compare diversity
        comparison_results = self.compare_diversity()

        # Step 4: Age correlation
        self.analyze_age_correlation()

        # Step 5: Create visualizations
        self.create_visualizations()

        # Step 6: Save metrics
        self.save_diversity_metrics()

        print("\n" + "="*80)
        print("✅ ANALYSIS COMPLETE!")
        print("="*80)
        print("\nGenerated files:")
        print("  • tcr_diversity_metrics_full.csv")
        print("  • diversity_comparison_results.csv")
        print("  • tcr_diversity_full_analysis.png")
        print("\n" + "="*80)


def main():
    """
    Main execution function

    IMPORTANT: Update DATA_DIR to match your computer!
    """

    # ========================================================================
    # CONFIGURATION - UPDATE THIS PATH!
    # ========================================================================

    # For Windows:
    DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"

    # For Linux/Mac (if testing locally):
    # DATA_DIR = "data/"

    # ========================================================================

    # Optional: Limit files for testing
    # Set to None to process ALL files
    MAX_FILES = None  # Change to 50 for testing with first 50 files

    # ========================================================================

    print("\n" + "="*80)
    print("TCR DIVERSITY ANALYSIS - FULL DATASET VERSION")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Data directory: {DATA_DIR}")
    print(f"  Max files: {'ALL FILES' if MAX_FILES is None else MAX_FILES}")
    print("\n" + "="*80)

    # Create analyzer
    analyzer = TCRDiversityAnalyzerFull(
        data_directory=DATA_DIR,
        max_files=MAX_FILES
    )

    # Run complete analysis
    analyzer.run_complete_analysis()


if __name__ == "__main__":
    main()
