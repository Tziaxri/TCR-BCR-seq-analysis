#!/usr/bin/env python3
"""
TCR V/D/J Gene Usage Analysis (Full Dataset Version)
====================================================

Analyzes V, D, and J gene segment usage patterns in CMV+ vs CMV- individuals
using ALL TCR files in your directory.

Features:
- Processes all .tsv.gz files automatically
- V/D/J gene frequency analysis
- V-J pairing patterns
- Statistical comparisons (CMV+ vs CMV-)
- Gene enrichment analysis
- Publication-quality visualizations

Updated: 2025-10-31

Usage:
------
# Edit DATA_DIR at the bottom of this script, then run:
python tcr_vdj_analysis_full.py

Author: Claude Code
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu, fisher_exact
import warnings
warnings.filterwarnings('ignore')

# Import data loader
from tcr_data_loader import TCRDataLoader

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (16, 12)
plt.rcParams['font.size'] = 10


class VDJAnalyzerFull:
    """
    Comprehensive V/D/J gene usage analyzer for full dataset

    Analyzes:
    - V gene segment frequencies
    - D gene segment frequencies
    - J gene segment frequencies
    - V-J pairing patterns
    - CMV-associated gene enrichment
    """

    def __init__(self, data_directory: str, max_files: int = None):
        """
        Initialize VDJ analyzer

        Parameters:
        -----------
        data_directory : str
            Path to directory with TCR files
        max_files : int
            Maximum files to process (None = all)
        """
        self.data_dir = data_directory
        self.max_files = max_files
        self.loader = None
        self.combined_data = None
        self.v_gene_usage = None
        self.j_gene_usage = None
        self.d_gene_usage = None

    def load_data(self):
        """Load all TCR files"""

        print("="*80)
        print("LOADING TCR DATA FOR VDJ ANALYSIS")
        print("="*80)

        # Create loader
        self.loader = TCRDataLoader(self.data_dir)

        # Load metadata
        try:
            self.loader.load_metadata()
        except FileNotFoundError:
            print("⚠️  metadata.tsv not found")

        # Load TCR files
        self.loader.load_all_tcr_files(
            pattern="*.tsv.gz",
            max_files=self.max_files,
            verbose=True
        )

        # Merge all data
        print("\n🔗 Combining all TCR data...")
        self.combined_data = self.loader.merge_with_metadata()

        # Filter to productive sequences only
        if 'productive' in self.combined_data.columns:
            before_count = len(self.combined_data)
            self.combined_data = self.combined_data[
                self.combined_data['productive'] == True
            ].copy()
            print(f"   ✓ Filtered to productive sequences:")
            print(f"     Before: {before_count:,} sequences")
            print(f"     After: {len(self.combined_data):,} sequences")

        # Print summary
        self.loader.print_summary()

        return self.combined_data

    def analyze_v_gene_usage(self):
        """Analyze V gene usage patterns"""

        print("\n" + "="*80)
        print("V GENE USAGE ANALYSIS")
        print("="*80)

        # Find V gene column
        v_col = None
        for possible_col in ['v_gene', 'vMaxResolved', 'v_call', 'v_family']:
            if possible_col in self.combined_data.columns:
                v_col = possible_col
                break

        if v_col is None:
            print("⚠️  No V gene column found")
            return None

        print(f"\n📊 Using column: '{v_col}'")

        # Calculate V gene frequencies per sample
        sample_v_usage = []

        for sample_id in self.combined_data['repertoire_id'].unique():
            sample_data = self.combined_data[
                self.combined_data['repertoire_id'] == sample_id
            ].copy()

            if len(sample_data) == 0:
                continue

            # Get V gene counts
            v_counts = sample_data[v_col].value_counts()
            total = v_counts.sum()

            # Get metadata for this sample
            metadata_cols = ['disease_subtype', 'age', 'sex', 'ancestry']
            sample_metadata = {}
            for col in metadata_cols:
                if col in sample_data.columns:
                    sample_metadata[col] = sample_data[col].iloc[0]

            # Create row for each V gene
            for v_gene, count in v_counts.items():
                if pd.isna(v_gene) or v_gene == '':
                    continue

                row = {
                    'repertoire_id': sample_id,
                    'v_gene': str(v_gene),
                    'count': count,
                    'frequency': count / total,
                    'total_sequences': total
                }
                row.update(sample_metadata)
                sample_v_usage.append(row)

        self.v_gene_usage = pd.DataFrame(sample_v_usage)

        print(f"\n✓ Calculated V gene usage:")
        print(f"  • Samples: {self.v_gene_usage['repertoire_id'].nunique()}")
        print(f"  • Unique V genes: {self.v_gene_usage['v_gene'].nunique()}")
        print(f"  • Total observations: {len(self.v_gene_usage)}")

        # Top V genes overall
        top_v = self.v_gene_usage.groupby('v_gene')['count'].sum().sort_values(ascending=False).head(10)
        print(f"\n📊 Top 10 most abundant V genes:")
        for i, (gene, count) in enumerate(top_v.items(), 1):
            print(f"  {i}. {gene}: {count:,} occurrences")

        # Save to CSV
        self.v_gene_usage.to_csv('v_gene_usage_full.csv', index=False)
        print(f"\n✓ Saved: v_gene_usage_full.csv")

        return self.v_gene_usage

    def analyze_j_gene_usage(self):
        """Analyze J gene usage patterns"""

        print("\n" + "="*80)
        print("J GENE USAGE ANALYSIS")
        print("="*80)

        # Find J gene column
        j_col = None
        for possible_col in ['j_gene', 'jMaxResolved', 'j_call']:
            if possible_col in self.combined_data.columns:
                j_col = possible_col
                break

        if j_col is None:
            print("⚠️  No J gene column found")
            return None

        print(f"\n📊 Using column: '{j_col}'")

        # Calculate J gene frequencies per sample
        sample_j_usage = []

        for sample_id in self.combined_data['repertoire_id'].unique():
            sample_data = self.combined_data[
                self.combined_data['repertoire_id'] == sample_id
            ].copy()

            if len(sample_data) == 0:
                continue

            # Get J gene counts
            j_counts = sample_data[j_col].value_counts()
            total = j_counts.sum()

            # Get metadata
            metadata_cols = ['disease_subtype', 'age', 'sex', 'ancestry']
            sample_metadata = {}
            for col in metadata_cols:
                if col in sample_data.columns:
                    sample_metadata[col] = sample_data[col].iloc[0]

            # Create row for each J gene
            for j_gene, count in j_counts.items():
                if pd.isna(j_gene) or j_gene == '':
                    continue

                row = {
                    'repertoire_id': sample_id,
                    'j_gene': str(j_gene),
                    'count': count,
                    'frequency': count / total,
                    'total_sequences': total
                }
                row.update(sample_metadata)
                sample_j_usage.append(row)

        self.j_gene_usage = pd.DataFrame(sample_j_usage)

        print(f"\n✓ Calculated J gene usage:")
        print(f"  • Samples: {self.j_gene_usage['repertoire_id'].nunique()}")
        print(f"  • Unique J genes: {self.j_gene_usage['j_gene'].nunique()}")

        # Top J genes
        top_j = self.j_gene_usage.groupby('j_gene')['count'].sum().sort_values(ascending=False).head(10)
        print(f"\n📊 Top 10 most abundant J genes:")
        for i, (gene, count) in enumerate(top_j.items(), 1):
            print(f"  {i}. {gene}: {count:,} occurrences")

        # Save to CSV
        self.j_gene_usage.to_csv('j_gene_usage_full.csv', index=False)
        print(f"\n✓ Saved: j_gene_usage_full.csv")

        return self.j_gene_usage

    def analyze_d_gene_usage(self):
        """Analyze D gene usage patterns (if available)"""

        print("\n" + "="*80)
        print("D GENE USAGE ANALYSIS")
        print("="*80)

        # Find D gene column
        d_col = None
        for possible_col in ['d_gene', 'dMaxResolved', 'd_call']:
            if possible_col in self.combined_data.columns:
                d_col = possible_col
                break

        if d_col is None:
            print("⚠️  No D gene column found (this is normal for many datasets)")
            return None

        print(f"\n📊 Using column: '{d_col}'")

        # Calculate D gene frequencies per sample
        sample_d_usage = []

        for sample_id in self.combined_data['repertoire_id'].unique():
            sample_data = self.combined_data[
                self.combined_data['repertoire_id'] == sample_id
            ].copy()

            if len(sample_data) == 0:
                continue

            # Get D gene counts
            d_counts = sample_data[d_col].value_counts()
            total = d_counts.sum()

            # Get metadata
            metadata_cols = ['disease_subtype', 'age', 'sex', 'ancestry']
            sample_metadata = {}
            for col in metadata_cols:
                if col in sample_data.columns:
                    sample_metadata[col] = sample_data[col].iloc[0]

            # Create row for each D gene
            for d_gene, count in d_counts.items():
                if pd.isna(d_gene) or d_gene == '' or d_gene == 'unresolved':
                    continue

                row = {
                    'repertoire_id': sample_id,
                    'd_gene': str(d_gene),
                    'count': count,
                    'frequency': count / total,
                    'total_sequences': total
                }
                row.update(sample_metadata)
                sample_d_usage.append(row)

        if len(sample_d_usage) > 0:
            self.d_gene_usage = pd.DataFrame(sample_d_usage)

            print(f"\n✓ Calculated D gene usage:")
            print(f"  • Samples: {self.d_gene_usage['repertoire_id'].nunique()}")
            print(f"  • Unique D genes: {self.d_gene_usage['d_gene'].nunique()}")

            # Save to CSV
            self.d_gene_usage.to_csv('d_gene_usage_full.csv', index=False)
            print(f"\n✓ Saved: d_gene_usage_full.csv")

        return self.d_gene_usage

    def compare_gene_usage_cmv(self, gene_type='v'):
        """
        Compare gene usage between CMV+ and CMV-

        Parameters:
        -----------
        gene_type : str
            'v', 'd', or 'j'
        """

        print(f"\n" + "="*80)
        print(f"{gene_type.upper()} GENE COMPARISON: CMV+ vs CMV-")
        print("="*80)

        # Select appropriate dataset
        if gene_type == 'v':
            usage_df = self.v_gene_usage
            gene_col = 'v_gene'
        elif gene_type == 'j':
            usage_df = self.j_gene_usage
            gene_col = 'j_gene'
        elif gene_type == 'd':
            usage_df = self.d_gene_usage
            gene_col = 'd_gene'
        else:
            raise ValueError("gene_type must be 'v', 'd', or 'j'")

        if usage_df is None or len(usage_df) == 0:
            print(f"⚠️  No {gene_type.upper()} gene data available")
            return None

        # Filter to CMV samples
        cmv_data = usage_df[usage_df['disease_subtype'].isin(['CMV+', 'CMV-'])].copy()

        if len(cmv_data) == 0:
            print("⚠️  No CMV+/CMV- samples found")
            return None

        # Get top genes to compare
        top_genes = usage_df.groupby(gene_col)['count'].sum().sort_values(ascending=False).head(20).index

        comparison_results = []

        for gene in top_genes:
            gene_data = cmv_data[cmv_data[gene_col] == gene]

            if len(gene_data) == 0:
                continue

            # Calculate mean frequency per group
            cmv_pos = gene_data[gene_data['disease_subtype'] == 'CMV+']['frequency']
            cmv_neg = gene_data[gene_data['disease_subtype'] == 'CMV-']['frequency']

            if len(cmv_pos) == 0 or len(cmv_neg) == 0:
                continue

            # Statistical test
            u_stat, p_val = mannwhitneyu(cmv_pos, cmv_neg, alternative='two-sided')

            # Fold change
            fold_change = (cmv_pos.mean() + 1e-10) / (cmv_neg.mean() + 1e-10)

            comparison_results.append({
                'gene': gene,
                'cmv_pos_mean_freq': cmv_pos.mean(),
                'cmv_pos_std': cmv_pos.std(),
                'cmv_pos_n': len(cmv_pos),
                'cmv_neg_mean_freq': cmv_neg.mean(),
                'cmv_neg_std': cmv_neg.std(),
                'cmv_neg_n': len(cmv_neg),
                'fold_change': fold_change,
                'p_value': p_val
            })

        comparison_df = pd.DataFrame(comparison_results)

        if len(comparison_df) > 0:
            # Sort by p-value
            comparison_df = comparison_df.sort_values('p_value')

            # Save results
            output_file = f'{gene_type}_genes_comparison_full.csv'
            comparison_df.to_csv(output_file, index=False)
            print(f"\n✓ Saved comparison results: {output_file}")

            # Print top results
            print(f"\n📊 Top {gene_type.upper()} gene differences (by p-value):")
            for i, row in comparison_df.head(10).iterrows():
                sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else 'ns'
                print(f"  {row['gene']}:")
                print(f"    CMV+: {row['cmv_pos_mean_freq']:.4f} ± {row['cmv_pos_std']:.4f}")
                print(f"    CMV-: {row['cmv_neg_mean_freq']:.4f} ± {row['cmv_neg_std']:.4f}")
                print(f"    Fold change: {row['fold_change']:.2f}x")
                print(f"    p={row['p_value']:.4f} {sig}")

        return comparison_df

    def analyze_vj_pairing(self):
        """Analyze V-J pairing patterns"""

        print("\n" + "="*80)
        print("V-J PAIRING ANALYSIS")
        print("="*80)

        # Find V and J columns
        v_col = None
        j_col = None

        for col in ['v_gene', 'vMaxResolved', 'v_call']:
            if col in self.combined_data.columns:
                v_col = col
                break

        for col in ['j_gene', 'jMaxResolved', 'j_call']:
            if col in self.combined_data.columns:
                j_col = col
                break

        if v_col is None or j_col is None:
            print("⚠️  V or J gene columns not found")
            return None

        # Count V-J pairs
        vj_pairs = self.combined_data.groupby([v_col, j_col, 'disease_subtype']).size().reset_index(name='count')
        vj_pairs.columns = ['v_gene', 'j_gene', 'disease_subtype', 'count']

        # Filter to top pairs
        top_pairs = vj_pairs.groupby(['v_gene', 'j_gene'])['count'].sum().sort_values(ascending=False).head(20)

        print(f"\n📊 Top 20 V-J pairings:")
        for i, ((v, j), count) in enumerate(top_pairs.items(), 1):
            print(f"  {i}. {v} + {j}: {count:,} occurrences")

        # Save V-J pairing data
        vj_pairs.to_csv('vj_pairings_full.csv', index=False)
        print(f"\n✓ Saved: vj_pairings_full.csv")

        return vj_pairs

    def create_visualizations(self, output_file='tcr_vdj_full_analysis.png'):
        """Create comprehensive VDJ visualization"""

        print("\n" + "="*80)
        print("GENERATING VISUALIZATIONS")
        print("="*80)

        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

        fig.suptitle('TCR V/D/J Gene Usage Analysis (Full Dataset)',
                     fontsize=18, fontweight='bold')

        # Panel A: V gene heatmap (top genes)
        if self.v_gene_usage is not None:
            ax = fig.add_subplot(gs[0, :])

            # Get top 15 V genes
            top_v_genes = self.v_gene_usage.groupby('v_gene')['count'].sum().sort_values(ascending=False).head(15).index

            # Create matrix for heatmap
            v_matrix = self.v_gene_usage[
                (self.v_gene_usage['v_gene'].isin(top_v_genes)) &
                (self.v_gene_usage['disease_subtype'].isin(['CMV+', 'CMV-']))
            ].pivot_table(
                index='repertoire_id',
                columns='v_gene',
                values='frequency',
                fill_value=0
            )

            # Sort by CMV status if available
            if hasattr(self.loader, 'metadata') and self.loader.metadata is not None:
                # Try to add CMV status
                pass

            sns.heatmap(v_matrix[top_v_genes], cmap='YlOrRd', ax=ax,
                       cbar_kws={'label': 'Frequency'}, linewidths=0.5)
            ax.set_title('A. V Gene Usage Heatmap (Top 15 Genes)', fontweight='bold', fontsize=12)
            ax.set_xlabel('V Gene', fontweight='bold')
            ax.set_ylabel('Sample ID', fontweight='bold')

        # Panel B: V gene comparison (CMV+ vs CMV-)
        if self.v_gene_usage is not None:
            ax = fig.add_subplot(gs[1, 0])

            cmv_data = self.v_gene_usage[self.v_gene_usage['disease_subtype'].isin(['CMV+', 'CMV-'])]
            top_10_v = self.v_gene_usage.groupby('v_gene')['count'].sum().sort_values(ascending=False).head(10).index

            plot_data = cmv_data[cmv_data['v_gene'].isin(top_10_v)]

            if len(plot_data) > 0:
                sns.boxplot(data=plot_data, x='v_gene', y='frequency',
                           hue='disease_subtype', ax=ax, palette=['salmon', 'skyblue'])
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
                ax.set_title('B. Top 10 V Genes (CMV+ vs CMV-)', fontweight='bold')
                ax.set_ylabel('Frequency', fontweight='bold')
                ax.set_xlabel('')
                ax.legend(title='CMV Status')

        # Panel C: J gene comparison
        if self.j_gene_usage is not None:
            ax = fig.add_subplot(gs[1, 1])

            cmv_data = self.j_gene_usage[self.j_gene_usage['disease_subtype'].isin(['CMV+', 'CMV-'])]
            top_10_j = self.j_gene_usage.groupby('j_gene')['count'].sum().sort_values(ascending=False).head(10).index

            plot_data = cmv_data[cmv_data['j_gene'].isin(top_10_j)]

            if len(plot_data) > 0:
                sns.boxplot(data=plot_data, x='j_gene', y='frequency',
                           hue='disease_subtype', ax=ax, palette=['salmon', 'skyblue'])
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
                ax.set_title('C. Top 10 J Genes (CMV+ vs CMV-)', fontweight='bold')
                ax.set_ylabel('Frequency', fontweight='bold')
                ax.set_xlabel('')
                ax.legend(title='CMV Status')

        # Panel D: V-J pairing (top pairs)
        ax = fig.add_subplot(gs[2, :])

        # Load VJ pairing data
        try:
            vj_data = pd.read_csv('vj_pairings_full.csv')
            top_pairs = vj_data.groupby(['v_gene', 'j_gene'])['count'].sum().sort_values(ascending=False).head(15)

            pairs_list = [f"{v}\n{j}" for v, j in top_pairs.index]
            counts_list = top_pairs.values

            bars = ax.barh(range(len(pairs_list)), counts_list, color='steelblue', alpha=0.7)
            ax.set_yticks(range(len(pairs_list)))
            ax.set_yticklabels(pairs_list, fontsize=8)
            ax.set_xlabel('Count', fontweight='bold')
            ax.set_title('D. Top 15 V-J Pairings', fontweight='bold')
            ax.grid(axis='x', alpha=0.3)

        except:
            ax.text(0.5, 0.5, 'V-J pairing data not available',
                   ha='center', va='center', fontsize=12)

        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✓ Saved visualization: {output_file}")

        return fig

    def run_complete_analysis(self):
        """Run complete VDJ analysis pipeline"""

        print("\n" + "="*80)
        print("TCR V/D/J GENE USAGE ANALYSIS - COMPLETE PIPELINE")
        print("="*80)

        # Step 1: Load data
        self.load_data()

        # Step 2: Analyze V genes
        self.analyze_v_gene_usage()

        # Step 3: Analyze J genes
        self.analyze_j_gene_usage()

        # Step 4: Analyze D genes
        self.analyze_d_gene_usage()

        # Step 5: Compare CMV+ vs CMV-
        self.compare_gene_usage_cmv('v')
        self.compare_gene_usage_cmv('j')
        if self.d_gene_usage is not None:
            self.compare_gene_usage_cmv('d')

        # Step 6: V-J pairing
        self.analyze_vj_pairing()

        # Step 7: Visualize
        self.create_visualizations()

        print("\n" + "="*80)
        print("✅ VDJ ANALYSIS COMPLETE!")
        print("="*80)
        print("\nGenerated files:")
        print("  • v_gene_usage_full.csv")
        print("  • j_gene_usage_full.csv")
        print("  • d_gene_usage_full.csv (if D genes available)")
        print("  • v_genes_comparison_full.csv")
        print("  • j_genes_comparison_full.csv")
        print("  • vj_pairings_full.csv")
        print("  • tcr_vdj_full_analysis.png")
        print("\n" + "="*80)


def main():
    """Main execution function"""

    # ========================================================================
    # CONFIGURATION - UPDATE THIS PATH!
    # ========================================================================

    # For Windows:
    DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"

    # For Linux/Mac:
    # DATA_DIR = "data/"

    # ========================================================================

    # Optional: Limit files for testing
    MAX_FILES = None  # Set to 50 for testing with first 50 files

    # ========================================================================

    print("\n" + "="*80)
    print("TCR V/D/J GENE USAGE ANALYSIS - FULL DATASET VERSION")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Data directory: {DATA_DIR}")
    print(f"  Max files: {'ALL FILES' if MAX_FILES is None else MAX_FILES}")
    print("\n" + "="*80)

    # Create analyzer
    analyzer = VDJAnalyzerFull(
        data_directory=DATA_DIR,
        max_files=MAX_FILES
    )

    # Run analysis
    analyzer.run_complete_analysis()


if __name__ == "__main__":
    main()
