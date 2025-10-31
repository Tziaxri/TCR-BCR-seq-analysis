#!/usr/bin/env python3
"""
GIANA-Style TCR Clustering Analysis
====================================

GIANA (Grouping of Lymphocyte Interactions by Paratope Hotspots) identifies
TCR clusters with similar CDR3 sequences, suggesting convergent antigen recognition.

This implementation provides:
1. CDR3 sequence similarity clustering
2. Motif identification in CDR3 regions
3. Public TCR detection (shared sequences across individuals)
4. CMV-associated cluster identification

Based on methods from:
- GLIPH (Glanville et al., Nature 2017)
- TCRdist (Dash et al., Nature 2017)

Author: Claude Code
Date: 2025-10-31
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (16, 10)

class GIANAAnalyzer:
    """GIANA-style TCR clustering and convergence analysis."""

    def __init__(self, data_dir='data'):
        self.data_dir = Path(data_dir)
        self.metadata = None
        self.repertoire_data = {}
        self.all_cdr3s = []
        self.clusters = []
        self.motifs = {}
        self.public_tcrs = {}

    def load_data(self):
        """Load metadata and TCR repertoire data."""
        print("=" * 80)
        print("GIANA-STYLE TCR CLUSTERING ANALYSIS")
        print("=" * 80)
        print("\n🔬 GIANA identifies TCR clusters with similar CDR3 sequences")
        print("   suggesting convergent immune responses to common antigens.\n")

        # Load metadata
        metadata_file = self.data_dir / 'metadata.tsv'
        self.metadata = pd.read_csv(metadata_file, sep='\t')

        # Extract CMV status
        self.metadata['cmv_status'] = self.metadata['disease_subtype'].apply(
            lambda x: 'CMV+' if 'CMV+' in str(x) else ('CMV-' if 'CMV-' in str(x) else 'Unknown')
        )
        self.metadata = self.metadata[self.metadata['cmv_status'].isin(['CMV+', 'CMV-'])]

        print(f"✓ Loaded metadata for {len(self.metadata)} participants")

        # Load TCR data
        tcr_files = list(self.data_dir.glob('part_table_*.tsv*'))
        print(f"✓ Found {len(tcr_files)} TCR repertoire files\n")

        for tcr_file in tcr_files:
            participant_id = tcr_file.stem.replace('part_table_', '').replace('.tsv', '').upper()

            try:
                df = pd.read_csv(tcr_file, sep='\t')
                df_productive = df[df['productive'] == 't'].copy()

                if len(df_productive) > 0:
                    # Extract relevant fields
                    df_productive['participant_id'] = participant_id
                    df_productive['v_gene'] = df_productive['v_call'].str.split('*').str[0]
                    df_productive['j_gene'] = df_productive['j_call'].str.split('*').str[0]

                    # Get CMV status
                    cmv_status = self.metadata[
                        self.metadata['participant_label'] == participant_id
                    ]['cmv_status'].iloc[0] if participant_id in self.metadata['participant_label'].values else 'Unknown'

                    df_productive['cmv_status'] = cmv_status

                    self.repertoire_data[participant_id] = df_productive
                    self.all_cdr3s.extend(df_productive[['cdr3_aa', 'v_gene', 'j_gene',
                                                          'participant_id', 'cmv_status']].to_dict('records'))

                    print(f"  {participant_id}: {len(df_productive)} sequences, CMV status: {cmv_status}")
            except Exception as e:
                print(f"  ⚠️  Error loading {tcr_file.name}: {e}")

        print(f"\n✓ Loaded {len(self.repertoire_data)} repertoires")
        print(f"✓ Total CDR3 sequences: {len(self.all_cdr3s)}")

    def calculate_hamming_distance(self, seq1, seq2):
        """Calculate Hamming distance between two sequences."""
        if len(seq1) != len(seq2):
            return float('inf')
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

    def calculate_levenshtein_distance(self, seq1, seq2):
        """Calculate Levenshtein (edit) distance between two sequences."""
        if len(seq1) < len(seq2):
            return self.calculate_levenshtein_distance(seq2, seq1)

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

    def find_public_tcrs(self, min_individuals=2):
        """Identify public TCRs (sequences shared across multiple individuals)."""
        print("\n" + "=" * 80)
        print("STEP 1: IDENTIFYING PUBLIC TCRs")
        print("=" * 80)
        print(f"\n📊 Public TCRs are sequences shared by ≥{min_individuals} individuals")
        print("   These suggest convergent responses to common antigens.\n")

        # Group by CDR3 amino acid sequence
        cdr3_dict = defaultdict(list)
        for cdr3_data in self.all_cdr3s:
            if cdr3_data['cdr3_aa'] and isinstance(cdr3_data['cdr3_aa'], str):
                cdr3_dict[cdr3_data['cdr3_aa']].append(cdr3_data)

        # Find sequences present in multiple individuals
        public_tcrs = {}
        for cdr3_aa, records in cdr3_dict.items():
            unique_individuals = set(r['participant_id'] for r in records)
            if len(unique_individuals) >= min_individuals:
                cmv_pos = sum(1 for r in records if r['cmv_status'] == 'CMV+')
                cmv_neg = sum(1 for r in records if r['cmv_status'] == 'CMV-')

                public_tcrs[cdr3_aa] = {
                    'cdr3_aa': cdr3_aa,
                    'n_individuals': len(unique_individuals),
                    'n_sequences': len(records),
                    'individuals': list(unique_individuals),
                    'v_genes': [r['v_gene'] for r in records],
                    'j_genes': [r['j_gene'] for r in records],
                    'cmv_pos_count': cmv_pos,
                    'cmv_neg_count': cmv_neg,
                    'cmv_enrichment': cmv_pos / (cmv_pos + cmv_neg) if (cmv_pos + cmv_neg) > 0 else 0
                }

        self.public_tcrs = public_tcrs

        print(f"✓ Found {len(public_tcrs)} public TCR sequences")

        # Show top public TCRs
        if len(public_tcrs) > 0:
            sorted_public = sorted(public_tcrs.items(), key=lambda x: x[1]['n_individuals'], reverse=True)
            print(f"\n📊 Top 10 Most Shared TCRs:")
            print(f"{'CDR3 Sequence':<20} {'# Indiv':<10} {'CMV+':<8} {'CMV-':<8} {'V/J Genes':<30}")
            print("-" * 80)

            for cdr3, data in sorted_public[:10]:
                v_gene = Counter(data['v_genes']).most_common(1)[0][0] if data['v_genes'] else 'N/A'
                j_gene = Counter(data['j_genes']).most_common(1)[0][0] if data['j_genes'] else 'N/A'
                print(f"{cdr3:<20} {data['n_individuals']:<10} {data['cmv_pos_count']:<8} "
                      f"{data['cmv_neg_count']:<8} {v_gene}/{j_gene}")

        return public_tcrs

    def find_cdr3_motifs(self, min_length=3, min_occurrences=3):
        """Identify common motifs in CDR3 sequences."""
        print("\n" + "=" * 80)
        print("STEP 2: IDENTIFYING CDR3 MOTIFS")
        print("=" * 80)
        print(f"\n📊 Searching for {min_length}-mer motifs appearing ≥{min_occurrences} times\n")

        motif_counts = defaultdict(lambda: {'count': 0, 'cmv_pos': 0, 'cmv_neg': 0, 'sequences': []})

        for cdr3_data in self.all_cdr3s:
            cdr3_aa = cdr3_data.get('cdr3_aa')
            if not cdr3_aa or not isinstance(cdr3_aa, str):
                continue

            # Extract all k-mers
            for i in range(len(cdr3_aa) - min_length + 1):
                motif = cdr3_aa[i:i+min_length]
                motif_counts[motif]['count'] += 1
                motif_counts[motif]['sequences'].append(cdr3_aa)

                if cdr3_data['cmv_status'] == 'CMV+':
                    motif_counts[motif]['cmv_pos'] += 1
                elif cdr3_data['cmv_status'] == 'CMV-':
                    motif_counts[motif]['cmv_neg'] += 1

        # Filter to common motifs
        common_motifs = {
            motif: data for motif, data in motif_counts.items()
            if data['count'] >= min_occurrences
        }

        self.motifs = common_motifs

        print(f"✓ Found {len(common_motifs)} common motifs")

        # Show top motifs
        if len(common_motifs) > 0:
            sorted_motifs = sorted(common_motifs.items(), key=lambda x: x[1]['count'], reverse=True)
            print(f"\n📊 Top 15 Most Common Motifs:")
            print(f"{'Motif':<10} {'Count':<10} {'CMV+':<8} {'CMV-':<8} {'CMV+ Ratio':<12}")
            print("-" * 60)

            for motif, data in sorted_motifs[:15]:
                total = data['cmv_pos'] + data['cmv_neg']
                ratio = data['cmv_pos'] / total if total > 0 else 0
                print(f"{motif:<10} {data['count']:<10} {data['cmv_pos']:<8} {data['cmv_neg']:<8} {ratio:.3f}")

        return common_motifs

    def cluster_similar_cdr3s(self, max_distance=1):
        """Cluster CDR3 sequences by similarity (same length, 1 AA difference)."""
        print("\n" + "=" * 80)
        print("STEP 3: CLUSTERING SIMILAR CDR3 SEQUENCES")
        print("=" * 80)
        print(f"\n📊 Clustering CDR3s with ≤{max_distance} amino acid difference (same length)\n")

        # Group by length first (Hamming distance requires same length)
        length_groups = defaultdict(list)
        for cdr3_data in self.all_cdr3s:
            cdr3_aa = cdr3_data.get('cdr3_aa')
            if cdr3_aa and isinstance(cdr3_aa, str):
                length_groups[len(cdr3_aa)].append(cdr3_data)

        all_clusters = []
        cluster_id = 0

        for length, sequences in length_groups.items():
            if len(sequences) < 2:
                continue

            # Create similarity matrix
            n = len(sequences)
            visited = [False] * n

            for i in range(n):
                if visited[i]:
                    continue

                cluster = [sequences[i]]
                visited[i] = True

                for j in range(i + 1, n):
                    if visited[j]:
                        continue

                    # Check if similar to any sequence in cluster
                    for seq_in_cluster in cluster:
                        dist = self.calculate_hamming_distance(
                            sequences[j]['cdr3_aa'],
                            seq_in_cluster['cdr3_aa']
                        )
                        if dist <= max_distance:
                            cluster.append(sequences[j])
                            visited[j] = True
                            break

                # Only keep clusters with multiple sequences
                if len(cluster) > 1:
                    cmv_pos = sum(1 for s in cluster if s['cmv_status'] == 'CMV+')
                    cmv_neg = sum(1 for s in cluster if s['cmv_status'] == 'CMV-')

                    all_clusters.append({
                        'cluster_id': cluster_id,
                        'size': len(cluster),
                        'cdr3_length': length,
                        'sequences': [s['cdr3_aa'] for s in cluster],
                        'v_genes': [s['v_gene'] for s in cluster],
                        'j_genes': [s['j_gene'] for s in cluster],
                        'cmv_pos_count': cmv_pos,
                        'cmv_neg_count': cmv_neg
                    })
                    cluster_id += 1

        self.clusters = all_clusters

        print(f"✓ Found {len(all_clusters)} CDR3 clusters")

        # Show top clusters
        if len(all_clusters) > 0:
            sorted_clusters = sorted(all_clusters, key=lambda x: x['size'], reverse=True)
            print(f"\n📊 Top 10 Largest Clusters:")
            print(f"{'Cluster ID':<12} {'Size':<8} {'Length':<10} {'CMV+':<8} {'CMV-':<8} {'Example CDR3':<25}")
            print("-" * 80)

            for cluster in sorted_clusters[:10]:
                example = cluster['sequences'][0]
                print(f"{cluster['cluster_id']:<12} {cluster['size']:<8} {cluster['cdr3_length']:<10} "
                      f"{cluster['cmv_pos_count']:<8} {cluster['cmv_neg_count']:<8} {example:<25}")

        return all_clusters

    def identify_cmv_associated_features(self):
        """Identify TCR features associated with CMV status."""
        print("\n" + "=" * 80)
        print("STEP 4: CMV-ASSOCIATED TCR FEATURES")
        print("=" * 80)

        # Public TCRs enriched in CMV+
        print("\n📊 Public TCRs Enriched in CMV+ Individuals:")
        cmv_public = [
            (cdr3, data) for cdr3, data in self.public_tcrs.items()
            if data['cmv_enrichment'] > 0.7 and data['n_individuals'] >= 2
        ]

        if len(cmv_public) > 0:
            print(f"\n{'CDR3 Sequence':<20} {'# Indiv':<10} {'CMV+ Ratio':<12} {'V/J':<25}")
            print("-" * 70)
            for cdr3, data in sorted(cmv_public, key=lambda x: x[1]['cmv_enrichment'], reverse=True)[:10]:
                v_gene = Counter(data['v_genes']).most_common(1)[0][0] if data['v_genes'] else 'N/A'
                j_gene = Counter(data['j_genes']).most_common(1)[0][0] if data['j_genes'] else 'N/A'
                print(f"{cdr3:<20} {data['n_individuals']:<10} {data['cmv_enrichment']:.3f}        "
                      f"{v_gene}/{j_gene}")
        else:
            print("  No CMV-enriched public TCRs found (need more samples)")

        # Motifs enriched in CMV+
        print("\n📊 Motifs Enriched in CMV+ Individuals:")
        cmv_motifs = [
            (motif, data) for motif, data in self.motifs.items()
            if data['cmv_pos'] + data['cmv_neg'] >= 3
        ]

        if len(cmv_motifs) > 0:
            cmv_motifs_sorted = sorted(
                cmv_motifs,
                key=lambda x: x[1]['cmv_pos'] / (x[1]['cmv_pos'] + x[1]['cmv_neg']),
                reverse=True
            )
            print(f"\n{'Motif':<10} {'CMV+':<8} {'CMV-':<8} {'CMV+ Ratio':<12}")
            print("-" * 40)
            for motif, data in cmv_motifs_sorted[:10]:
                total = data['cmv_pos'] + data['cmv_neg']
                ratio = data['cmv_pos'] / total if total > 0 else 0
                print(f"{motif:<10} {data['cmv_pos']:<8} {data['cmv_neg']:<8} {ratio:.3f}")
        else:
            print("  Insufficient data for motif enrichment analysis")

        # Clusters enriched in CMV+
        print("\n📊 Clusters Enriched in CMV+ Individuals:")
        cmv_clusters = [
            cluster for cluster in self.clusters
            if cluster['cmv_pos_count'] + cluster['cmv_neg_count'] >= 2
        ]

        if len(cmv_clusters) > 0:
            cmv_clusters_sorted = sorted(
                cmv_clusters,
                key=lambda x: x['cmv_pos_count'] / (x['cmv_pos_count'] + x['cmv_neg_count']),
                reverse=True
            )
            print(f"\n{'Cluster ID':<12} {'Size':<8} {'CMV+':<8} {'CMV-':<8} {'CMV+ Ratio':<12}")
            print("-" * 50)
            for cluster in cmv_clusters_sorted[:10]:
                total = cluster['cmv_pos_count'] + cluster['cmv_neg_count']
                ratio = cluster['cmv_pos_count'] / total if total > 0 else 0
                print(f"{cluster['cluster_id']:<12} {cluster['size']:<8} "
                      f"{cluster['cmv_pos_count']:<8} {cluster['cmv_neg_count']:<8} {ratio:.3f}")
        else:
            print("  Insufficient data for cluster enrichment analysis")

    def create_visualizations(self):
        """Generate comprehensive visualizations."""
        print("\n" + "=" * 80)
        print("STEP 5: GENERATING VISUALIZATIONS")
        print("=" * 80)

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('GIANA-Style TCR Clustering Analysis', fontsize=16, fontweight='bold')

        # Plot 1: Public TCR frequency
        ax1 = axes[0, 0]
        if len(self.public_tcrs) > 0:
            public_df = pd.DataFrame([
                {
                    'n_individuals': data['n_individuals'],
                    'cmv_enrichment': data['cmv_enrichment']
                }
                for data in self.public_tcrs.values()
            ])

            ax1.scatter(public_df['n_individuals'], public_df['cmv_enrichment'],
                       alpha=0.6, s=100, c=public_df['cmv_enrichment'],
                       cmap='RdBu_r', edgecolors='black')
            ax1.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
            ax1.set_xlabel('Number of Individuals Sharing TCR', fontsize=12)
            ax1.set_ylabel('CMV+ Fraction', fontsize=12)
            ax1.set_title('A) Public TCR Distribution', fontsize=12, fontweight='bold')
            ax1.grid(True, alpha=0.3)

        # Plot 2: CDR3 length distribution
        ax2 = axes[0, 1]
        cdr3_lengths = defaultdict(lambda: {'CMV+': 0, 'CMV-': 0})
        for cdr3_data in self.all_cdr3s:
            if cdr3_data.get('cdr3_aa') and isinstance(cdr3_data['cdr3_aa'], str):
                length = len(cdr3_data['cdr3_aa'])
                status = cdr3_data['cmv_status']
                if status in ['CMV+', 'CMV-']:
                    cdr3_lengths[length][status] += 1

        lengths = sorted(cdr3_lengths.keys())
        cmv_pos_counts = [cdr3_lengths[l]['CMV+'] for l in lengths]
        cmv_neg_counts = [cdr3_lengths[l]['CMV-'] for l in lengths]

        x = np.arange(len(lengths))
        width = 0.35
        ax2.bar(x - width/2, cmv_pos_counts, width, label='CMV+', color='lightcoral', alpha=0.8)
        ax2.bar(x + width/2, cmv_neg_counts, width, label='CMV-', color='lightblue', alpha=0.8)
        ax2.set_xlabel('CDR3 Length (amino acids)', fontsize=12)
        ax2.set_ylabel('Count', fontsize=12)
        ax2.set_title('B) CDR3 Length Distribution', fontsize=12, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(lengths)
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # Plot 3: Cluster size distribution
        ax3 = axes[1, 0]
        if len(self.clusters) > 0:
            cluster_sizes = [c['size'] for c in self.clusters]
            ax3.hist(cluster_sizes, bins=20, color='steelblue', alpha=0.7, edgecolor='black')
            ax3.set_xlabel('Cluster Size', fontsize=12)
            ax3.set_ylabel('Number of Clusters', fontsize=12)
            ax3.set_title('C) TCR Cluster Size Distribution', fontsize=12, fontweight='bold')
            ax3.grid(True, alpha=0.3)

        # Plot 4: Top motifs
        ax4 = axes[1, 1]
        if len(self.motifs) > 0:
            sorted_motifs = sorted(self.motifs.items(), key=lambda x: x[1]['count'], reverse=True)[:15]
            motif_names = [m[0] for m in sorted_motifs]
            motif_counts = [m[1]['count'] for m in sorted_motifs]

            ax4.barh(range(len(motif_names)), motif_counts, color='darkorange', alpha=0.7)
            ax4.set_yticks(range(len(motif_names)))
            ax4.set_yticklabels(motif_names, fontsize=10)
            ax4.set_xlabel('Frequency', fontsize=12)
            ax4.set_title('D) Top 15 CDR3 Motifs', fontsize=12, fontweight='bold')
            ax4.grid(True, alpha=0.3, axis='x')

        plt.tight_layout()
        output_file = 'tcr_giana_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✓ Saved visualization to: {output_file}")

    def save_results(self):
        """Save results to CSV files."""
        print("\n" + "=" * 80)
        print("STEP 6: SAVING RESULTS")
        print("=" * 80)

        # Save public TCRs
        if len(self.public_tcrs) > 0:
            public_df = pd.DataFrame([
                {
                    'cdr3_aa': cdr3,
                    'n_individuals': data['n_individuals'],
                    'n_sequences': data['n_sequences'],
                    'cmv_pos_count': data['cmv_pos_count'],
                    'cmv_neg_count': data['cmv_neg_count'],
                    'cmv_enrichment': data['cmv_enrichment'],
                    'most_common_v': Counter(data['v_genes']).most_common(1)[0][0] if data['v_genes'] else '',
                    'most_common_j': Counter(data['j_genes']).most_common(1)[0][0] if data['j_genes'] else ''
                }
                for cdr3, data in self.public_tcrs.items()
            ])
            public_df = public_df.sort_values('n_individuals', ascending=False)
            public_df.to_csv('public_tcrs.csv', index=False)
            print("✓ Saved public TCRs to: public_tcrs.csv")

        # Save clusters
        if len(self.clusters) > 0:
            cluster_df = pd.DataFrame([
                {
                    'cluster_id': c['cluster_id'],
                    'size': c['size'],
                    'cdr3_length': c['cdr3_length'],
                    'cmv_pos_count': c['cmv_pos_count'],
                    'cmv_neg_count': c['cmv_neg_count'],
                    'example_cdr3': c['sequences'][0] if c['sequences'] else ''
                }
                for c in self.clusters
            ])
            cluster_df = cluster_df.sort_values('size', ascending=False)
            cluster_df.to_csv('tcr_clusters.csv', index=False)
            print("✓ Saved TCR clusters to: tcr_clusters.csv")

        # Save motifs
        if len(self.motifs) > 0:
            motif_df = pd.DataFrame([
                {
                    'motif': motif,
                    'count': data['count'],
                    'cmv_pos': data['cmv_pos'],
                    'cmv_neg': data['cmv_neg'],
                    'cmv_enrichment': data['cmv_pos'] / (data['cmv_pos'] + data['cmv_neg'])
                                      if (data['cmv_pos'] + data['cmv_neg']) > 0 else 0
                }
                for motif, data in self.motifs.items()
            ])
            motif_df = motif_df.sort_values('count', ascending=False)
            motif_df.to_csv('cdr3_motifs.csv', index=False)
            print("✓ Saved CDR3 motifs to: cdr3_motifs.csv")

        print("\n" + "=" * 80)
        print("GIANA ANALYSIS COMPLETE!")
        print("=" * 80)

def main():
    """Main analysis pipeline."""
    print("\n" + "=" * 80)
    print("GIANA-STYLE TCR CLUSTERING PIPELINE")
    print("=" * 80)
    print("\nIdentifying convergent TCR responses through sequence clustering")
    print("=" * 80 + "\n")

    analyzer = GIANAAnalyzer(data_dir='data')

    # Load data
    analyzer.load_data()

    # Find public TCRs
    analyzer.find_public_tcrs(min_individuals=2)

    # Find motifs
    analyzer.find_cdr3_motifs(min_length=3, min_occurrences=3)

    # Cluster sequences
    analyzer.cluster_similar_cdr3s(max_distance=1)

    # Identify CMV-associated features
    analyzer.identify_cmv_associated_features()

    # Visualize
    analyzer.create_visualizations()

    # Save results
    analyzer.save_results()

    print("\n📊 BIOLOGICAL INTERPRETATION:")
    print("\nPublic TCRs and clustered sequences suggest:")
    print("  • Convergent immune responses to common epitopes")
    print("  • CMV-specific TCR signatures")
    print("  • HLA-restricted antigen recognition patterns")
    print("\nWith full dataset, expect to see:")
    print("  • CMV-associated public TCRs (shared across CMV+ individuals)")
    print("  • Specific motifs in CMV-reactive CDR3 regions")
    print("  • Larger clusters in CMV+ individuals (clonal expansion)")
    print("\n" + "=" * 80)

if __name__ == '__main__':
    main()
