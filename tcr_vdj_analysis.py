#!/usr/bin/env python3
"""
TCR V/D/J Gene Usage Analysis: CMV+ vs CMV- Comparison
=======================================================

This script analyzes V, D, and J gene segment usage patterns in TCR repertoires,
comparing CMV+ and CMV- individuals to identify disease-associated gene preferences.

V/D/J genes encode the variable, diversity, and joining regions of TCRs.
Different pathogens can drive preferential usage of specific gene segments.

Analysis workflow:
1. Load TCR repertoire data with V/D/J annotations
2. Calculate gene usage frequencies per individual
3. Compare gene usage between CMV+ and CMV- groups
4. Identify significantly enriched/depleted genes
5. Visualize gene usage patterns with heatmaps

Author: Claude Code
Date: 2025-10-31
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from scipy.stats import mannwhitneyu, fisher_exact
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (16, 10)
plt.rcParams['font.size'] = 10

class VDJAnalyzer:
    """Analyzes V/D/J gene usage patterns in TCR repertoires."""

    def __init__(self, data_dir='data'):
        self.data_dir = Path(data_dir)
        self.metadata = None
        self.repertoire_data = {}
        self.v_usage = None
        self.j_usage = None
        self.d_usage = None

    def load_data(self):
        """Load metadata and TCR repertoire data."""
        print("=" * 80)
        print("VDJ GENE USAGE ANALYSIS")
        print("=" * 80)

        # Load metadata
        metadata_file = self.data_dir / 'metadata.tsv'
        self.metadata = pd.read_csv(metadata_file, sep='\t')

        # Extract CMV status
        self.metadata['cmv_status'] = self.metadata['disease_subtype'].apply(
            lambda x: 'CMV+' if 'CMV+' in str(x) else ('CMV-' if 'CMV-' in str(x) else 'Unknown')
        )
        self.metadata = self.metadata[self.metadata['cmv_status'].isin(['CMV+', 'CMV-'])]

        print(f"\n✓ Loaded metadata for {len(self.metadata)} participants")
        print(f"  - CMV+: {sum(self.metadata['cmv_status'] == 'CMV+')}")
        print(f"  - CMV-: {sum(self.metadata['cmv_status'] == 'CMV-')}")

        # Load TCR data
        tcr_files = list(self.data_dir.glob('part_table_*.tsv*'))
        print(f"\n✓ Found {len(tcr_files)} TCR repertoire files\n")

        for tcr_file in tcr_files:
            participant_id = tcr_file.stem.replace('part_table_', '').replace('.tsv', '').upper()

            try:
                df = pd.read_csv(tcr_file, sep='\t')
                df_productive = df[df['productive'] == 't'].copy()

                if len(df_productive) > 0:
                    repertoire_id = df_productive['repertoire_id'].iloc[0]
                    self.repertoire_data[participant_id] = {
                        'data': df_productive,
                        'repertoire_id': repertoire_id
                    }
                    print(f"  Loaded {participant_id}: {len(df_productive)} productive sequences")
            except Exception as e:
                print(f"  ⚠️  Error loading {tcr_file.name}: {e}")

        print(f"\n✓ Successfully loaded {len(self.repertoire_data)} repertoires")

    def calculate_v_gene_usage(self):
        """Calculate V gene usage frequencies for each individual."""
        print("\n" + "=" * 80)
        print("STEP 1: V GENE USAGE ANALYSIS")
        print("=" * 80)

        print("\n📊 V genes encode the variable region of TCRs")
        print("   Different V genes recognize different antigen peptides")
        print("   CMV infection may drive expansion of specific V genes\n")

        v_usage_data = []

        for participant_id, rep_data in self.repertoire_data.items():
            df = rep_data['data']

            # Extract V gene calls (remove allele info, keep gene family)
            df['v_gene'] = df['v_call'].str.split('*').str[0]

            # Calculate frequency of each V gene
            v_counts = df['v_gene'].value_counts()
            v_freqs = v_counts / v_counts.sum()

            # Store for this individual
            for v_gene, freq in v_freqs.items():
                v_usage_data.append({
                    'participant_id': participant_id,
                    'v_gene': v_gene,
                    'frequency': freq,
                    'count': v_counts[v_gene]
                })

        self.v_usage = pd.DataFrame(v_usage_data)

        # Merge with metadata
        self.v_usage = self.v_usage.merge(
            self.metadata[['participant_label', 'cmv_status']],
            left_on='participant_id',
            right_on='participant_label',
            how='left'
        )

        print(f"✓ Calculated V gene usage for {len(self.repertoire_data)} individuals")
        print(f"✓ Detected {self.v_usage['v_gene'].nunique()} unique V genes")

        # Show top V genes
        top_v = self.v_usage.groupby('v_gene')['frequency'].mean().sort_values(ascending=False).head(10)
        print("\n📊 Top 10 most common V genes (average frequency):")
        for gene, freq in top_v.items():
            print(f"   {gene:20s} {freq*100:6.2f}%")

        return self.v_usage

    def calculate_j_gene_usage(self):
        """Calculate J gene usage frequencies for each individual."""
        print("\n" + "=" * 80)
        print("STEP 2: J GENE USAGE ANALYSIS")
        print("=" * 80)

        print("\n📊 J genes encode the joining region of TCRs")
        print("   J gene selection affects CDR3 loop structure")
        print("   Certain J genes may pair preferentially with V genes\n")

        j_usage_data = []

        for participant_id, rep_data in self.repertoire_data.items():
            df = rep_data['data']

            # Extract J gene calls
            df['j_gene'] = df['j_call'].str.split('*').str[0]

            # Calculate frequency
            j_counts = df['j_gene'].value_counts()
            j_freqs = j_counts / j_counts.sum()

            for j_gene, freq in j_freqs.items():
                j_usage_data.append({
                    'participant_id': participant_id,
                    'j_gene': j_gene,
                    'frequency': freq,
                    'count': j_counts[j_gene]
                })

        self.j_usage = pd.DataFrame(j_usage_data)

        # Merge with metadata
        self.j_usage = self.j_usage.merge(
            self.metadata[['participant_label', 'cmv_status']],
            left_on='participant_id',
            right_on='participant_label',
            how='left'
        )

        print(f"✓ Calculated J gene usage for {len(self.repertoire_data)} individuals")
        print(f"✓ Detected {self.j_usage['j_gene'].nunique()} unique J genes")

        # Show top J genes
        top_j = self.j_usage.groupby('j_gene')['frequency'].mean().sort_values(ascending=False).head(10)
        print("\n📊 Top 10 most common J genes (average frequency):")
        for gene, freq in top_j.items():
            print(f"   {gene:20s} {freq*100:6.2f}%")

        return self.j_usage

    def calculate_d_gene_usage(self):
        """Calculate D gene usage frequencies for each individual."""
        print("\n" + "=" * 80)
        print("STEP 3: D GENE USAGE ANALYSIS")
        print("=" * 80)

        print("\n📊 D genes encode the diversity region of TCRs (TCR-beta only)")
        print("   D genes contribute to CDR3 sequence diversity")
        print("   Note: Not all sequences have D gene annotations\n")

        d_usage_data = []

        for participant_id, rep_data in self.repertoire_data.items():
            df = rep_data['data']

            # Filter to sequences with D gene calls
            df_with_d = df[df['d_call'].notna() & (df['d_call'] != '')].copy()

            if len(df_with_d) == 0:
                continue

            # Extract D gene calls
            df_with_d['d_gene'] = df_with_d['d_call'].str.split('*').str[0]

            # Calculate frequency
            d_counts = df_with_d['d_gene'].value_counts()
            d_freqs = d_counts / d_counts.sum()

            for d_gene, freq in d_freqs.items():
                d_usage_data.append({
                    'participant_id': participant_id,
                    'd_gene': d_gene,
                    'frequency': freq,
                    'count': d_counts[d_gene]
                })

        if len(d_usage_data) > 0:
            self.d_usage = pd.DataFrame(d_usage_data)

            # Merge with metadata
            self.d_usage = self.d_usage.merge(
                self.metadata[['participant_label', 'cmv_status']],
                left_on='participant_id',
                right_on='participant_label',
                how='left'
            )

            print(f"✓ Calculated D gene usage for {len(d_usage_data)} individuals")
            print(f"✓ Detected {self.d_usage['d_gene'].nunique()} unique D genes")

            # Show all D genes (usually small number)
            top_d = self.d_usage.groupby('d_gene')['frequency'].mean().sort_values(ascending=False)
            print("\n📊 D gene frequencies (average):")
            for gene, freq in top_d.items():
                print(f"   {gene:20s} {freq*100:6.2f}%")
        else:
            print("⚠️  No D gene annotations found in data")
            self.d_usage = pd.DataFrame()

        return self.d_usage

    def compare_gene_usage(self):
        """Compare gene usage between CMV+ and CMV- groups."""
        print("\n" + "=" * 80)
        print("STEP 4: STATISTICAL COMPARISON (CMV+ vs CMV-)")
        print("=" * 80)

        results = {'v_genes': [], 'j_genes': [], 'd_genes': []}

        # V gene comparison
        print("\n🔬 V GENE COMPARISON:\n")
        if self.v_usage is not None and len(self.v_usage) > 0:
            # Get genes present in at least 50% of samples
            gene_presence = self.v_usage.groupby('v_gene')['participant_id'].nunique()
            common_genes = gene_presence[gene_presence >= max(2, len(self.repertoire_data) * 0.3)].index

            for v_gene in common_genes:
                gene_data = self.v_usage[self.v_usage['v_gene'] == v_gene]

                cmv_pos = gene_data[gene_data['cmv_status'] == 'CMV+']['frequency']
                cmv_neg = gene_data[gene_data['cmv_status'] == 'CMV-']['frequency']

                if len(cmv_pos) > 0 and len(cmv_neg) > 0:
                    # Mann-Whitney U test
                    stat, p_value = mannwhitneyu(cmv_pos, cmv_neg, alternative='two-sided')

                    mean_pos = cmv_pos.mean()
                    mean_neg = cmv_neg.mean()
                    fold_change = mean_pos / mean_neg if mean_neg > 0 else np.inf

                    results['v_genes'].append({
                        'gene': v_gene,
                        'mean_cmv_pos': mean_pos,
                        'mean_cmv_neg': mean_neg,
                        'fold_change': fold_change,
                        'p_value': p_value,
                        'n_cmv_pos': len(cmv_pos),
                        'n_cmv_neg': len(cmv_neg)
                    })

                    if p_value < 0.05:
                        direction = "ENRICHED" if fold_change > 1 else "DEPLETED"
                        print(f"   ✓ {v_gene:20s} {direction} in CMV+ (FC={fold_change:.2f}, p={p_value:.4f})")

            if not any(r['p_value'] < 0.05 for r in results['v_genes']):
                print("   No significant V gene differences detected")

        # J gene comparison
        print("\n🔬 J GENE COMPARISON:\n")
        if self.j_usage is not None and len(self.j_usage) > 0:
            common_genes = self.j_usage.groupby('j_gene')['participant_id'].nunique()
            common_genes = common_genes[common_genes >= max(2, len(self.repertoire_data) * 0.3)].index

            for j_gene in common_genes:
                gene_data = self.j_usage[self.j_usage['j_gene'] == j_gene]

                cmv_pos = gene_data[gene_data['cmv_status'] == 'CMV+']['frequency']
                cmv_neg = gene_data[gene_data['cmv_status'] == 'CMV-']['frequency']

                if len(cmv_pos) > 0 and len(cmv_neg) > 0:
                    stat, p_value = mannwhitneyu(cmv_pos, cmv_neg, alternative='two-sided')

                    mean_pos = cmv_pos.mean()
                    mean_neg = cmv_neg.mean()
                    fold_change = mean_pos / mean_neg if mean_neg > 0 else np.inf

                    results['j_genes'].append({
                        'gene': j_gene,
                        'mean_cmv_pos': mean_pos,
                        'mean_cmv_neg': mean_neg,
                        'fold_change': fold_change,
                        'p_value': p_value,
                        'n_cmv_pos': len(cmv_pos),
                        'n_cmv_neg': len(cmv_neg)
                    })

                    if p_value < 0.05:
                        direction = "ENRICHED" if fold_change > 1 else "DEPLETED"
                        print(f"   ✓ {j_gene:20s} {direction} in CMV+ (FC={fold_change:.2f}, p={p_value:.4f})")

            if not any(r['p_value'] < 0.05 for r in results['j_genes']):
                print("   No significant J gene differences detected")

        # D gene comparison
        print("\n🔬 D GENE COMPARISON:\n")
        if self.d_usage is not None and len(self.d_usage) > 0:
            for d_gene in self.d_usage['d_gene'].unique():
                gene_data = self.d_usage[self.d_usage['d_gene'] == d_gene]

                cmv_pos = gene_data[gene_data['cmv_status'] == 'CMV+']['frequency']
                cmv_neg = gene_data[gene_data['cmv_status'] == 'CMV-']['frequency']

                if len(cmv_pos) > 0 and len(cmv_neg) > 0:
                    stat, p_value = mannwhitneyu(cmv_pos, cmv_neg, alternative='two-sided')

                    mean_pos = cmv_pos.mean()
                    mean_neg = cmv_neg.mean()
                    fold_change = mean_pos / mean_neg if mean_neg > 0 else np.inf

                    results['d_genes'].append({
                        'gene': d_gene,
                        'mean_cmv_pos': mean_pos,
                        'mean_cmv_neg': mean_neg,
                        'fold_change': fold_change,
                        'p_value': p_value,
                        'n_cmv_pos': len(cmv_pos),
                        'n_cmv_neg': len(cmv_neg)
                    })

                    if p_value < 0.05:
                        direction = "ENRICHED" if fold_change > 1 else "DEPLETED"
                        print(f"   ✓ {d_gene:20s} {direction} in CMV+ (FC={fold_change:.2f}, p={p_value:.4f})")

            if len(results['d_genes']) == 0 or not any(r['p_value'] < 0.05 for r in results['d_genes']):
                print("   No significant D gene differences detected")
        else:
            print("   No D gene data available for comparison")

        return results

    def create_visualizations(self):
        """Generate comprehensive VDJ usage visualizations."""
        print("\n" + "=" * 80)
        print("STEP 5: GENERATING VISUALIZATIONS")
        print("=" * 80)

        fig = plt.figure(figsize=(18, 12))
        gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

        # Plot 1: V gene usage heatmap
        ax1 = fig.add_subplot(gs[0, :])
        if self.v_usage is not None and len(self.v_usage) > 0:
            # Create matrix: individuals x V genes
            v_pivot = self.v_usage.pivot_table(
                index='participant_id',
                columns='v_gene',
                values='frequency',
                fill_value=0
            )

            # Select top 20 most variable genes
            v_var = v_pivot.var().sort_values(ascending=False).head(20)
            v_pivot_top = v_pivot[v_var.index]

            # Add CMV status as row colors
            cmv_colors = []
            for pid in v_pivot_top.index:
                if pid in self.metadata['participant_label'].values:
                    status = self.metadata[self.metadata['participant_label'] == pid]['cmv_status'].iloc[0]
                    cmv_colors.append('red' if status == 'CMV+' else 'blue')
                else:
                    cmv_colors.append('gray')

            sns.heatmap(v_pivot_top, cmap='YlOrRd', ax=ax1, cbar_kws={'label': 'Frequency'})
            ax1.set_title('A) V Gene Usage Heatmap (Top 20 Most Variable)', fontsize=14, fontweight='bold')
            ax1.set_xlabel('V Gene', fontsize=12)
            ax1.set_ylabel('Individual', fontsize=12)

            # Add color bar for CMV status
            for i, color in enumerate(cmv_colors):
                ax1.add_patch(plt.Rectangle((-0.5, i), 0.3, 1, color=color, transform=ax1.get_yaxis_transform(), clip_on=False))

        # Plot 2: V gene usage comparison (top genes)
        ax2 = fig.add_subplot(gs[1, 0])
        if self.v_usage is not None and len(self.v_usage) > 0:
            # Get top 10 V genes
            top_v = self.v_usage.groupby('v_gene')['frequency'].mean().sort_values(ascending=False).head(10)

            plot_data = []
            for gene in top_v.index:
                gene_data = self.v_usage[self.v_usage['v_gene'] == gene]
                for _, row in gene_data.iterrows():
                    plot_data.append({
                        'gene': gene,
                        'frequency': row['frequency'],
                        'cmv_status': row['cmv_status']
                    })

            plot_df = pd.DataFrame(plot_data)

            # Box plot
            for i, gene in enumerate(top_v.index):
                gene_subset = plot_df[plot_df['gene'] == gene]
                cmv_pos = gene_subset[gene_subset['cmv_status'] == 'CMV+']['frequency']
                cmv_neg = gene_subset[gene_subset['cmv_status'] == 'CMV-']['frequency']

                positions = [i*2, i*2+0.8]
                bp = ax2.boxplot([cmv_pos, cmv_neg] if len(cmv_pos) > 0 and len(cmv_neg) > 0 else [[0], [0]],
                                 positions=positions, widths=0.6, patch_artist=True,
                                 showfliers=False)
                bp['boxes'][0].set_facecolor('lightcoral')
                bp['boxes'][1].set_facecolor('lightblue')

            ax2.set_xticks([i*2 + 0.4 for i in range(len(top_v))])
            ax2.set_xticklabels(top_v.index, rotation=45, ha='right')
            ax2.set_ylabel('Frequency', fontsize=12)
            ax2.set_title('B) Top 10 V Genes: CMV+ vs CMV-', fontsize=12, fontweight='bold')
            ax2.legend([plt.Rectangle((0,0),1,1,fc='lightcoral'), plt.Rectangle((0,0),1,1,fc='lightblue')],
                      ['CMV+', 'CMV-'], loc='upper right')
            ax2.grid(True, alpha=0.3)

        # Plot 3: J gene usage comparison
        ax3 = fig.add_subplot(gs[1, 1])
        if self.j_usage is not None and len(self.j_usage) > 0:
            # Get top 10 J genes
            top_j = self.j_usage.groupby('j_gene')['frequency'].mean().sort_values(ascending=False).head(10)

            plot_data = []
            for gene in top_j.index:
                gene_data = self.j_usage[self.j_usage['j_gene'] == gene]
                for _, row in gene_data.iterrows():
                    plot_data.append({
                        'gene': gene,
                        'frequency': row['frequency'],
                        'cmv_status': row['cmv_status']
                    })

            plot_df = pd.DataFrame(plot_data)

            for i, gene in enumerate(top_j.index):
                gene_subset = plot_df[plot_df['gene'] == gene]
                cmv_pos = gene_subset[gene_subset['cmv_status'] == 'CMV+']['frequency']
                cmv_neg = gene_subset[gene_subset['cmv_status'] == 'CMV-']['frequency']

                positions = [i*2, i*2+0.8]
                bp = ax3.boxplot([cmv_pos, cmv_neg] if len(cmv_pos) > 0 and len(cmv_neg) > 0 else [[0], [0]],
                                 positions=positions, widths=0.6, patch_artist=True,
                                 showfliers=False)
                bp['boxes'][0].set_facecolor('lightcoral')
                bp['boxes'][1].set_facecolor('lightblue')

            ax3.set_xticks([i*2 + 0.4 for i in range(len(top_j))])
            ax3.set_xticklabels(top_j.index, rotation=45, ha='right')
            ax3.set_ylabel('Frequency', fontsize=12)
            ax3.set_title('C) Top 10 J Genes: CMV+ vs CMV-', fontsize=12, fontweight='bold')
            ax3.legend([plt.Rectangle((0,0),1,1,fc='lightcoral'), plt.Rectangle((0,0),1,1,fc='lightblue')],
                      ['CMV+', 'CMV-'], loc='upper right')
            ax3.grid(True, alpha=0.3)

        # Plot 4: VJ pairing analysis
        ax4 = fig.add_subplot(gs[2, :])
        # Combine V and J data for pairing analysis
        vj_pairs = []
        for participant_id, rep_data in self.repertoire_data.items():
            df = rep_data['data']
            df['v_gene'] = df['v_call'].str.split('*').str[0]
            df['j_gene'] = df['j_call'].str.split('*').str[0]

            # Count VJ pairs
            pair_counts = df.groupby(['v_gene', 'j_gene']).size()
            total = len(df)

            for (v, j), count in pair_counts.items():
                vj_pairs.append({
                    'participant_id': participant_id,
                    'v_gene': v,
                    'j_gene': j,
                    'count': count,
                    'frequency': count / total
                })

        if len(vj_pairs) > 0:
            vj_df = pd.DataFrame(vj_pairs)

            # Get top VJ pairs
            top_pairs = vj_df.groupby(['v_gene', 'j_gene'])['frequency'].mean().sort_values(ascending=False).head(15)

            pair_labels = [f"{v}\n{j}" for v, j in top_pairs.index]

            ax4.barh(range(len(top_pairs)), top_pairs.values, color='steelblue')
            ax4.set_yticks(range(len(top_pairs)))
            ax4.set_yticklabels(pair_labels, fontsize=9)
            ax4.set_xlabel('Average Frequency', fontsize=12)
            ax4.set_title('D) Top 15 V-J Gene Pairings', fontsize=12, fontweight='bold')
            ax4.grid(True, alpha=0.3, axis='x')

        plt.suptitle('TCR V/D/J Gene Usage Analysis: CMV+ vs CMV-', fontsize=16, fontweight='bold', y=0.995)

        output_file = 'tcr_vdj_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✓ Saved visualization to: {output_file}")

        return fig

    def save_results(self, comparison_results):
        """Save detailed results to CSV files."""
        print("\n" + "=" * 80)
        print("STEP 6: SAVING RESULTS")
        print("=" * 80)

        # Save V gene usage
        if self.v_usage is not None and len(self.v_usage) > 0:
            self.v_usage.to_csv('v_gene_usage.csv', index=False)
            print("✓ Saved V gene usage to: v_gene_usage.csv")

        # Save J gene usage
        if self.j_usage is not None and len(self.j_usage) > 0:
            self.j_usage.to_csv('j_gene_usage.csv', index=False)
            print("✓ Saved J gene usage to: j_gene_usage.csv")

        # Save D gene usage
        if self.d_usage is not None and len(self.d_usage) > 0:
            self.d_usage.to_csv('d_gene_usage.csv', index=False)
            print("✓ Saved D gene usage to: d_gene_usage.csv")

        # Save comparison results
        for gene_type, results in comparison_results.items():
            if len(results) > 0:
                results_df = pd.DataFrame(results)
                results_df = results_df.sort_values('p_value')
                filename = f'{gene_type}_comparison.csv'
                results_df.to_csv(filename, index=False)
                print(f"✓ Saved {gene_type} comparison to: {filename}")

        print("\n" + "=" * 80)
        print("VDJ ANALYSIS COMPLETE!")
        print("=" * 80)

def main():
    """Main analysis pipeline."""
    print("\n" + "=" * 80)
    print("TCR V/D/J GENE USAGE ANALYSIS PIPELINE")
    print("=" * 80)
    print("\nAnalyzing V, D, and J gene segment usage patterns")
    print("Comparing CMV+ vs CMV- individuals")
    print("=" * 80 + "\n")

    # Initialize analyzer
    analyzer = VDJAnalyzer(data_dir='data')

    # Load data
    analyzer.load_data()

    # Calculate gene usage
    analyzer.calculate_v_gene_usage()
    analyzer.calculate_j_gene_usage()
    analyzer.calculate_d_gene_usage()

    # Compare between groups
    comparison_results = analyzer.compare_gene_usage()

    # Visualize
    analyzer.create_visualizations()

    # Save results
    analyzer.save_results(comparison_results)

    print("\n📊 BIOLOGICAL INTERPRETATION:")
    print("\nWith full dataset, expect to see:")
    print("  • CMV+ individuals may show enrichment of specific V genes")
    print("    (e.g., TRBV6, TRBV9 often associated with CMV epitopes)")
    print("  • Preferential VJ pairings for CMV-reactive TCRs")
    print("  • Public TCR sequences (shared V-J combinations)")
    print("\nThese patterns reflect:")
    print("  • Convergent immune responses to common CMV epitopes")
    print("  • HLA-restricted antigen presentation")
    print("  • Memory T cell expansion\n")

    print("=" * 80)
    print("Add more samples to data/ directory for robust comparisons!")
    print("=" * 80)

if __name__ == '__main__':
    main()
