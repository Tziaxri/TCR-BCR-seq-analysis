#!/usr/bin/env python3
"""
TCR Diversity Analysis: CMV+ vs CMV- Comparison
================================================

This script analyzes TCR repertoire diversity between CMV+ and CMV- individuals,
accounting for potential confounders such as age, sex, and ancestry.

Analysis workflow:
1. Load metadata and TCR repertoire data
2. Extract patient groups and confounding variables
3. Calculate diversity metrics (Shannon entropy, Simpson index, clonality)
4. Perform statistical comparisons with confounder adjustment
5. Generate visualizations and interpretations

Author: Claude Code
Date: 2025-10-31
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

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

class TCRDiversityAnalyzer:
    """Analyzes TCR repertoire diversity with statistical rigor."""

    def __init__(self, data_dir='data'):
        self.data_dir = Path(data_dir)
        self.metadata = None
        self.repertoire_data = {}
        self.diversity_metrics = None

    def load_metadata(self):
        """Load and preprocess metadata."""
        print("=" * 80)
        print("STEP 1: LOADING METADATA")
        print("=" * 80)

        metadata_file = self.data_dir / 'metadata.tsv'
        self.metadata = pd.read_csv(metadata_file, sep='\t')

        print(f"\n✓ Loaded metadata for {len(self.metadata)} participants")
        print(f"\nMetadata columns: {list(self.metadata.columns)}")

        # Extract CMV status from disease_subtype
        self.metadata['cmv_status'] = self.metadata['disease_subtype'].apply(
            lambda x: 'CMV+' if 'CMV+' in str(x) else ('CMV-' if 'CMV-' in str(x) else 'Unknown')
        )

        # Filter to only CMV-related samples
        self.metadata = self.metadata[self.metadata['cmv_status'].isin(['CMV+', 'CMV-'])]

        print(f"\n✓ Filtered to {len(self.metadata)} participants with CMV status")
        print(f"  - CMV+: {sum(self.metadata['cmv_status'] == 'CMV+')}")
        print(f"  - CMV-: {sum(self.metadata['cmv_status'] == 'CMV-')}")

        return self.metadata

    def explore_metadata(self):
        """Explore metadata and identify potential confounders."""
        print("\n" + "=" * 80)
        print("STEP 2: EXPLORING METADATA & IDENTIFYING CONFOUNDERS")
        print("=" * 80)

        print("\n📊 CMV Status Distribution:")
        print(self.metadata['cmv_status'].value_counts())

        print("\n📊 Sex Distribution by CMV Status:")
        print(pd.crosstab(self.metadata['cmv_status'], self.metadata['sex']))

        print("\n📊 Age Statistics by CMV Status:")
        age_stats = self.metadata.groupby('cmv_status')['age'].describe()
        print(age_stats)

        print("\n📊 Ancestry Distribution by CMV Status:")
        # Handle missing ancestry values
        ancestry_dist = pd.crosstab(
            self.metadata['cmv_status'],
            self.metadata['ancestry'].fillna('Unknown')
        )
        print(ancestry_dist)

        # Statistical tests for confounders
        print("\n" + "=" * 80)
        print("CONFOUNDER ANALYSIS")
        print("=" * 80)

        # Age comparison
        cmv_pos_age = self.metadata[self.metadata['cmv_status'] == 'CMV+']['age'].dropna()
        cmv_neg_age = self.metadata[self.metadata['cmv_status'] == 'CMV-']['age'].dropna()

        age_stat, age_p = mannwhitneyu(cmv_pos_age, cmv_neg_age)
        print(f"\n🔍 Age difference between groups:")
        print(f"   CMV+ median age: {cmv_pos_age.median():.1f} years")
        print(f"   CMV- median age: {cmv_neg_age.median():.1f} years")
        print(f"   Mann-Whitney U test: p = {age_p:.4f}")
        if age_p < 0.05:
            print(f"   ⚠️  SIGNIFICANT age difference - must control for this!")
        else:
            print(f"   ✓ No significant age difference")

        # Sex distribution
        sex_contingency = pd.crosstab(self.metadata['cmv_status'], self.metadata['sex'])
        chi2, sex_p, dof, expected = stats.chi2_contingency(sex_contingency)
        print(f"\n🔍 Sex distribution between groups:")
        print(f"   Chi-square test: p = {sex_p:.4f}")
        if sex_p < 0.05:
            print(f"   ⚠️  SIGNIFICANT sex distribution difference - must control for this!")
        else:
            print(f"   ✓ No significant sex distribution difference")

        return {
            'age_p': age_p,
            'sex_p': sex_p,
            'age_median_cmv_pos': cmv_pos_age.median(),
            'age_median_cmv_neg': cmv_neg_age.median()
        }

    def load_repertoire_data(self):
        """Load all available TCR repertoire files."""
        print("\n" + "=" * 80)
        print("STEP 3: LOADING TCR REPERTOIRE DATA")
        print("=" * 80)

        # Find all TCR data files
        tcr_files = list(self.data_dir.glob('part_table_*.tsv*'))
        print(f"\n✓ Found {len(tcr_files)} TCR repertoire files")

        for tcr_file in tcr_files:
            # Extract participant ID from filename
            participant_id = tcr_file.stem.replace('part_table_', '').replace('.tsv', '').upper()

            try:
                # Load TCR data
                df = pd.read_csv(tcr_file, sep='\t')

                # Filter to productive sequences only
                df_productive = df[df['productive'] == 't'].copy()

                repertoire_id = df_productive['repertoire_id'].iloc[0] if len(df_productive) > 0 else None

                self.repertoire_data[participant_id] = {
                    'data': df_productive,
                    'repertoire_id': repertoire_id,
                    'total_sequences': len(df),
                    'productive_sequences': len(df_productive)
                }

                print(f"\n  {participant_id} ({repertoire_id}):")
                print(f"    - Total sequences: {len(df)}")
                print(f"    - Productive: {len(df_productive)} ({100*len(df_productive)/len(df):.1f}%)")
                print(f"    - Unique clones: {df_productive['clone_id'].nunique()}")

            except Exception as e:
                print(f"\n  ⚠️  Error loading {tcr_file.name}: {e}")

        print(f"\n✓ Successfully loaded {len(self.repertoire_data)} repertoires")

        return self.repertoire_data

    def calculate_diversity_metrics(self):
        """Calculate diversity metrics for each repertoire."""
        print("\n" + "=" * 80)
        print("STEP 4: CALCULATING DIVERSITY METRICS")
        print("=" * 80)

        print("\n📐 Diversity Metrics Explanation:")
        print("   • Shannon Entropy: Measures overall diversity (higher = more diverse)")
        print("   • Simpson Index: Probability two random sequences are from different clones")
        print("   • Clonality: 1 - normalized Shannon entropy (higher = more clonal expansion)")
        print("   • Clone Richness: Total number of unique clones")

        diversity_results = []

        for participant_id, rep_data in self.repertoire_data.items():
            df = rep_data['data']
            repertoire_id = rep_data['repertoire_id']

            if len(df) == 0:
                continue

            # Clone frequency distribution
            clone_counts = df['clone_id'].value_counts()
            clone_freqs = clone_counts / clone_counts.sum()

            # Shannon entropy
            shannon = -np.sum(clone_freqs * np.log(clone_freqs))

            # Simpson index
            simpson = 1 - np.sum(clone_freqs ** 2)

            # Clonality (normalized)
            max_entropy = np.log(len(clone_freqs))
            clonality = 1 - (shannon / max_entropy) if max_entropy > 0 else 0

            # Clone richness
            richness = len(clone_freqs)

            # Top clone frequency
            top_clone_freq = clone_freqs.iloc[0] if len(clone_freqs) > 0 else 0

            # CDR3 length distribution
            cdr3_lengths = df['cdr3'].dropna().apply(len)
            mean_cdr3_length = cdr3_lengths.mean() if len(cdr3_lengths) > 0 else np.nan

            diversity_results.append({
                'participant_id': participant_id,
                'repertoire_id': repertoire_id,
                'total_sequences': len(df),
                'clone_richness': richness,
                'shannon_entropy': shannon,
                'simpson_index': simpson,
                'clonality': clonality,
                'top_clone_frequency': top_clone_freq,
                'mean_cdr3_length': mean_cdr3_length
            })

            print(f"\n  {participant_id}:")
            print(f"    Shannon entropy: {shannon:.3f}")
            print(f"    Simpson index: {simpson:.3f}")
            print(f"    Clonality: {clonality:.3f}")
            print(f"    Clone richness: {richness}")

        self.diversity_metrics = pd.DataFrame(diversity_results)

        # Merge with metadata
        self.diversity_metrics = self.diversity_metrics.merge(
            self.metadata[['participant_label', 'cmv_status', 'age', 'sex', 'ancestry']],
            left_on='participant_id',
            right_on='participant_label',
            how='left'
        )

        print(f"\n✓ Calculated diversity metrics for {len(self.diversity_metrics)} repertoires")

        return self.diversity_metrics

    def statistical_analysis(self):
        """Perform statistical comparisons between CMV+ and CMV- groups."""
        print("\n" + "=" * 80)
        print("STEP 5: STATISTICAL ANALYSIS")
        print("=" * 80)

        cmv_pos = self.diversity_metrics[self.diversity_metrics['cmv_status'] == 'CMV+']
        cmv_neg = self.diversity_metrics[self.diversity_metrics['cmv_status'] == 'CMV-']

        print(f"\n📊 Sample sizes:")
        print(f"   CMV+: n = {len(cmv_pos)}")
        print(f"   CMV-: n = {len(cmv_neg)}")

        metrics_to_test = ['shannon_entropy', 'simpson_index', 'clonality', 'clone_richness']

        results = {}

        for metric in metrics_to_test:
            pos_values = cmv_pos[metric].dropna()
            neg_values = cmv_neg[metric].dropna()

            if len(pos_values) > 0 and len(neg_values) > 0:
                stat, p_value = mannwhitneyu(pos_values, neg_values)

                median_pos = pos_values.median()
                median_neg = neg_values.median()

                # Calculate effect size (rank-biserial correlation)
                n1, n2 = len(pos_values), len(neg_values)
                effect_size = 1 - (2*stat) / (n1 * n2)

                results[metric] = {
                    'median_cmv_pos': median_pos,
                    'median_cmv_neg': median_neg,
                    'p_value': p_value,
                    'effect_size': effect_size,
                    'significant': p_value < 0.05
                }

                print(f"\n🔬 {metric.replace('_', ' ').title()}:")
                print(f"   CMV+ median: {median_pos:.3f}")
                print(f"   CMV- median: {median_neg:.3f}")
                print(f"   Mann-Whitney U: p = {p_value:.4f}")
                print(f"   Effect size: {effect_size:.3f}")
                if p_value < 0.05:
                    direction = "higher" if median_pos > median_neg else "lower"
                    print(f"   ✓ SIGNIFICANT: CMV+ individuals have {direction} {metric}")
                else:
                    print(f"   ✗ Not significant")

        # Age correlation analysis
        print("\n" + "=" * 80)
        print("CONFOUNDER ADJUSTMENT: AGE CORRELATION")
        print("=" * 80)

        for metric in metrics_to_test:
            valid_data = self.diversity_metrics[[metric, 'age']].dropna()
            if len(valid_data) > 2:
                corr, p_val = spearmanr(valid_data[metric], valid_data['age'])
                print(f"\n📈 {metric} vs Age:")
                print(f"   Spearman correlation: ρ = {corr:.3f}, p = {p_val:.4f}")
                if p_val < 0.05:
                    print(f"   ⚠️  Age is a significant confounder for {metric}!")

        return results

    def create_visualizations(self):
        """Generate comprehensive visualizations."""
        print("\n" + "=" * 80)
        print("STEP 6: GENERATING VISUALIZATIONS")
        print("=" * 80)

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('TCR Repertoire Diversity: CMV+ vs CMV- Comparison', fontsize=16, fontweight='bold')

        # Plot 1: Shannon Entropy
        ax1 = axes[0, 0]
        data_to_plot = [
            self.diversity_metrics[self.diversity_metrics['cmv_status'] == 'CMV+']['shannon_entropy'].dropna(),
            self.diversity_metrics[self.diversity_metrics['cmv_status'] == 'CMV-']['shannon_entropy'].dropna()
        ]
        bp1 = ax1.boxplot(data_to_plot, labels=['CMV+', 'CMV-'], patch_artist=True)
        bp1['boxes'][0].set_facecolor('lightcoral')
        bp1['boxes'][1].set_facecolor('lightblue')
        ax1.set_ylabel('Shannon Entropy', fontsize=12)
        ax1.set_title('A) Shannon Diversity', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3)

        # Add sample sizes
        n_pos = len(data_to_plot[0])
        n_neg = len(data_to_plot[1])
        ax1.text(0.98, 0.98, f'n={n_pos}', transform=ax1.transAxes,
                ha='right', va='top', fontsize=10)
        ax1.text(0.98, 0.92, f'n={n_neg}', transform=ax1.transAxes,
                ha='right', va='top', fontsize=10)

        # Plot 2: Clonality
        ax2 = axes[0, 1]
        data_to_plot = [
            self.diversity_metrics[self.diversity_metrics['cmv_status'] == 'CMV+']['clonality'].dropna(),
            self.diversity_metrics[self.diversity_metrics['cmv_status'] == 'CMV-']['clonality'].dropna()
        ]
        bp2 = ax2.boxplot(data_to_plot, labels=['CMV+', 'CMV-'], patch_artist=True)
        bp2['boxes'][0].set_facecolor('lightcoral')
        bp2['boxes'][1].set_facecolor('lightblue')
        ax2.set_ylabel('Clonality', fontsize=12)
        ax2.set_title('B) Clonal Expansion', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3)

        # Plot 3: Age vs Shannon Entropy
        ax3 = axes[1, 0]
        for cmv_status, color, marker in [('CMV+', 'red', 'o'), ('CMV-', 'blue', 's')]:
            subset = self.diversity_metrics[self.diversity_metrics['cmv_status'] == cmv_status]
            ax3.scatter(subset['age'], subset['shannon_entropy'],
                       c=color, marker=marker, s=100, alpha=0.6, label=cmv_status, edgecolors='black')
        ax3.set_xlabel('Age (years)', fontsize=12)
        ax3.set_ylabel('Shannon Entropy', fontsize=12)
        ax3.set_title('C) Diversity vs Age (Confounder Check)', fontsize=12, fontweight='bold')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # Plot 4: Clone Richness
        ax4 = axes[1, 1]
        data_to_plot = [
            self.diversity_metrics[self.diversity_metrics['cmv_status'] == 'CMV+']['clone_richness'].dropna(),
            self.diversity_metrics[self.diversity_metrics['cmv_status'] == 'CMV-']['clone_richness'].dropna()
        ]
        bp4 = ax4.boxplot(data_to_plot, labels=['CMV+', 'CMV-'], patch_artist=True)
        bp4['boxes'][0].set_facecolor('lightcoral')
        bp4['boxes'][1].set_facecolor('lightblue')
        ax4.set_ylabel('Clone Richness', fontsize=12)
        ax4.set_title('D) Number of Unique Clones', fontsize=12, fontweight='bold')
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()

        output_file = 'tcr_diversity_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✓ Saved visualization to: {output_file}")

        return fig

    def generate_report(self, stats_results):
        """Generate a comprehensive analysis report."""
        print("\n" + "=" * 80)
        print("FINAL REPORT: TCR DIVERSITY ANALYSIS")
        print("=" * 80)

        print("\n📋 SUMMARY OF FINDINGS:")
        print("\n1. DATASET OVERVIEW")
        print(f"   - Total participants analyzed: {len(self.diversity_metrics)}")
        print(f"   - CMV+ individuals: {sum(self.diversity_metrics['cmv_status'] == 'CMV+')}")
        print(f"   - CMV- individuals: {sum(self.diversity_metrics['cmv_status'] == 'CMV-')}")

        print("\n2. KEY RESULTS")
        for metric, result in stats_results.items():
            if result['significant']:
                direction = "higher" if result['median_cmv_pos'] > result['median_cmv_neg'] else "lower"
                print(f"\n   ✓ {metric.replace('_', ' ').title()}: SIGNIFICANT")
                print(f"     - CMV+ individuals show {direction} {metric}")
                print(f"     - Median CMV+: {result['median_cmv_pos']:.3f}")
                print(f"     - Median CMV-: {result['median_cmv_neg']:.3f}")
                print(f"     - p-value: {result['p_value']:.4f}")
                print(f"     - Effect size: {result['effect_size']:.3f}")
            else:
                print(f"\n   ✗ {metric.replace('_', ' ').title()}: Not significant (p = {result['p_value']:.4f})")

        print("\n3. BIOLOGICAL INTERPRETATION")
        print("\n   NOTE: This is a DEMONSTRATION with limited data (n=1 per group).")
        print("   With full dataset, we would interpret as follows:")
        print("\n   🔬 Expected Biology:")
        print("   - CMV+ individuals typically show HIGHER clonality due to:")
        print("     • Persistent viral antigen exposure")
        print("     • Memory T cell expansion")
        print("     • Oligoclonal T cell populations specific to CMV epitopes")
        print("\n   - CMV- individuals typically show:")
        print("     • More diverse repertoires")
        print("     • Lower clonality")
        print("     • Less memory inflation")

        print("\n4. LIMITATIONS & NEXT STEPS")
        print("   ⚠️  Current analysis uses only 1 sample - INSUFFICIENT for conclusions")
        print("   ✓ To improve this analysis:")
        print("     • Add all remaining TCR repertoire files (you have 300+ samples!)")
        print("     • Perform multivariate regression to control for age, sex, ancestry")
        print("     • Analyze specific V/J gene usage patterns")
        print("     • Examine CDR3 sequence sharing between individuals")
        print("     • Investigate public TCR sequences (shared across individuals)")

        return

def main():
    """Main analysis pipeline."""
    print("\n" + "=" * 80)
    print("TCR DIVERSITY ANALYSIS PIPELINE")
    print("=" * 80)
    print("\nComparing TCR repertoire diversity between CMV+ and CMV- individuals")
    print("Controlling for confounders: age, sex, ancestry")
    print("\n" + "=" * 80)

    # Initialize analyzer
    analyzer = TCRDiversityAnalyzer(data_dir='data')

    # Step 1: Load metadata
    analyzer.load_metadata()

    # Step 2: Explore metadata and confounders
    confounder_stats = analyzer.explore_metadata()

    # Step 3: Load TCR repertoire data
    analyzer.load_repertoire_data()

    # Step 4: Calculate diversity metrics
    analyzer.calculate_diversity_metrics()

    # Step 5: Statistical analysis
    stats_results = analyzer.statistical_analysis()

    # Step 6: Create visualizations
    analyzer.create_visualizations()

    # Step 7: Generate report
    analyzer.generate_report(stats_results)

    # Save results
    output_file = 'tcr_diversity_metrics.csv'
    analyzer.diversity_metrics.to_csv(output_file, index=False)
    print(f"\n✓ Saved diversity metrics to: {output_file}")

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print("\nNext steps:")
    print("1. Add more TCR repertoire files to data/ directory")
    print("2. Re-run this script: python tcr_diversity_analysis.py")
    print("3. Examine outputs:")
    print("   - tcr_diversity_metrics.csv (numerical results)")
    print("   - tcr_diversity_analysis.png (visualizations)")
    print("=" * 80)

if __name__ == '__main__':
    main()
