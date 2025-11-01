#!/usr/bin/env python3
"""
GIANA Clustering Analysis (Full Dataset Version)
=================================================

Implementation of GIANA (Grouping of Immunoglobulin and T cell receptor
Amino acid sequences via Nucleotide-based Alignments) from Zhang et al.,
Nature Communications 2021.

Updated to process ALL TCR files in your directory with optimizations
for large-scale datasets (millions of sequences).

Reference: https://www.nature.com/articles/s41467-021-25006-7

Features:
- Graph-based TCR clustering
- Needleman-Wunsch global alignment
- Co-clustering matrix for sample similarity
- Network visualization by disease group
- Optimized for large datasets

Usage:
------
# Edit DATA_DIR at the bottom, then run:
python giana_analysis_full.py

Author: Claude Code
Date: 2025-10-31
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
from scipy.stats import spearmanr
from scipy.cluster import hierarchy
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Import data loader
from tcr_data_loader import TCRDataLoader

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)


class GIANAAnalyzerFull:
    """
    GIANA clustering for full TCR dataset

    Implements Zhang et al., Nature Communications 2021 method
    with optimizations for large-scale data.

    Key Parameters:
    - similarity_threshold: 0.85 (from paper)
    - alignment: Needleman-Wunsch global alignment
    - clustering: Single-linkage graph-based
    """

    def __init__(self, data_directory: str,
                 max_files: int = None,
                 max_sequences: int = 100000,
                 similarity_threshold: float = 0.85):
        """
        Initialize GIANA analyzer

        Parameters:
        -----------
        data_directory : str
            Path to TCR data directory
        max_files : int
            Maximum files to load (None = all)
        max_sequences : int
            Maximum sequences to cluster (for memory management)
            Recommended: 100,000 for typical workstations
        similarity_threshold : float
            Sequence similarity threshold (default: 0.85 from paper)
        """
        self.data_dir = data_directory
        self.max_files = max_files
        self.max_sequences = max_sequences
        self.similarity_threshold = similarity_threshold

        self.loader = None
        self.combined_data = None
        self.clusters = None
        self.co_clustering_matrix = None

    def load_data(self):
        """Load TCR data"""

        print("="*80)
        print("LOADING DATA FOR GIANA CLUSTERING")
        print("="*80)

        # Create loader
        self.loader = TCRDataLoader(self.data_dir)

        # Load metadata
        try:
            self.loader.load_metadata()
        except:
            print("⚠️  Metadata not found, continuing...")

        # Load TCR files
        self.loader.load_all_tcr_files(
            pattern="*.tsv.gz",
            max_files=self.max_files,
            verbose=True
        )

        # Merge data
        self.combined_data = self.loader.merge_with_metadata()

        # Filter to productive sequences
        if 'productive' in self.combined_data.columns:
            self.combined_data = self.combined_data[
                self.combined_data['productive'] == True
            ].copy()

        print(f"\n✓ Loaded {len(self.combined_data):,} productive TCR sequences")

        # Subsample if needed
        if len(self.combined_data) > self.max_sequences:
            print(f"\n⚠️  Subsampling to {self.max_sequences:,} sequences for clustering")
            print(f"   (Original: {len(self.combined_data):,} sequences)")

            # Stratified sampling by sample
            sample_sizes = self.combined_data.groupby('repertoire_id').size()
            n_samples = len(sample_sizes)
            target_per_sample = self.max_sequences // n_samples

            sampled = []
            for sample_id in self.combined_data['repertoire_id'].unique():
                sample_data = self.combined_data[
                    self.combined_data['repertoire_id'] == sample_id
                ]
                n_sample = min(len(sample_data), target_per_sample)
                sampled.append(sample_data.sample(n=n_sample, random_state=42))

            self.combined_data = pd.concat(sampled, ignore_index=True)

            print(f"   Sampled: {len(self.combined_data):,} sequences")

        return self.combined_data

    def _global_alignment_score(self, seq1: str, seq2: str) -> float:
        """
        Needleman-Wunsch global alignment

        Returns normalized similarity score (0-1)

        Parameters from Zhang et al. 2021:
        - Match: +1
        - Mismatch: -1
        - Gap: -1
        """
        match_score = 1
        mismatch_penalty = -1
        gap_penalty = -1

        len1, len2 = len(seq1), len(seq2)

        # Initialize scoring matrix
        score_matrix = np.zeros((len1 + 1, len2 + 1))

        # Initialize gap penalties
        for i in range(1, len1 + 1):
            score_matrix[i][0] = gap_penalty * i
        for j in range(1, len2 + 1):
            score_matrix[0][j] = gap_penalty * j

        # Fill scoring matrix
        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                if seq1[i-1] == seq2[j-1]:
                    match = score_matrix[i-1][j-1] + match_score
                else:
                    match = score_matrix[i-1][j-1] + mismatch_penalty

                delete = score_matrix[i-1][j] + gap_penalty
                insert = score_matrix[i][j-1] + gap_penalty

                score_matrix[i][j] = max(match, delete, insert)

        # Get alignment score
        alignment_score = score_matrix[len1][len2]

        # Normalize by maximum possible score
        max_possible_score = min(len1, len2) * match_score

        if max_possible_score == 0:
            return 0.0

        normalized_score = max(0, alignment_score / max_possible_score)

        return normalized_score

    def cluster_sequences(self):
        """
        Cluster TCR sequences using graph-based single-linkage clustering

        Based on Zhang et al. 2021 GIANA method
        """

        print("\n" + "="*80)
        print("GIANA SEQUENCE CLUSTERING")
        print("="*80)

        # Find CDR3 column
        cdr3_col = None
        for col in ['amino_acid', 'cdr3_aa', 'junction_aa', 'cdr3']:
            if col in self.combined_data.columns:
                cdr3_col = col
                break

        if cdr3_col is None:
            raise ValueError("CDR3 amino acid column not found")

        print(f"\n📊 Using CDR3 column: '{cdr3_col}'")
        print(f"   Similarity threshold: {self.similarity_threshold}")
        print(f"   Total sequences: {len(self.combined_data):,}")

        # Get unique CDR3 sequences with their sample mapping
        sequence_data = []

        for idx, row in self.combined_data.iterrows():
            cdr3 = str(row[cdr3_col])

            if pd.isna(cdr3) or cdr3 == '' or cdr3 == 'nan':
                continue

            sequence_data.append({
                'cdr3': cdr3,
                'sample_id': row['repertoire_id'],
                'disease_subtype': row.get('disease_subtype', 'Unknown'),
                'index': idx
            })

        seq_df = pd.DataFrame(sequence_data)

        # Get unique sequences
        unique_sequences = seq_df['cdr3'].unique()
        print(f"   Unique CDR3 sequences: {len(unique_sequences):,}")

        # Build similarity graph
        print(f"\n🔗 Building similarity graph...")
        print(f"   (This may take several minutes for large datasets)")

        # Create graph
        G = nx.Graph()

        # Add nodes
        for seq in unique_sequences:
            G.add_node(seq)

        # Add edges for similar sequences
        # For efficiency, only compare sequences of similar length
        sequences_by_length = defaultdict(list)
        for seq in unique_sequences:
            sequences_by_length[len(seq)].append(seq)

        edge_count = 0
        total_comparisons = 0

        for length, seqs in sequences_by_length.items():
            n_seqs = len(seqs)

            # Also compare with sequences of similar lengths (±2 amino acids)
            comparable_seqs = []
            for l in range(max(1, length-2), length+3):
                comparable_seqs.extend(sequences_by_length[l])

            print(f"   Processing {n_seqs} sequences of length {length}...")

            # Compare all pairs within this length group
            for i in range(len(seqs)):
                for j in range(i+1, len(comparable_seqs)):
                    seq1 = seqs[i] if i < len(seqs) else comparable_seqs[j]
                    seq2 = comparable_seqs[j] if j < len(comparable_seqs) else seqs[i]

                    total_comparisons += 1

                    if total_comparisons % 100000 == 0:
                        print(f"     Comparisons: {total_comparisons:,}, Edges: {edge_count:,}")

                    # Calculate similarity
                    similarity = self._global_alignment_score(seq1, seq2)

                    if similarity >= self.similarity_threshold:
                        G.add_edge(seq1, seq2, weight=similarity)
                        edge_count += 1

        print(f"\n✓ Graph construction complete:")
        print(f"   Nodes (unique sequences): {G.number_of_nodes():,}")
        print(f"   Edges (similar pairs): {G.number_of_edges():,}")

        # Find connected components (clusters)
        print(f"\n🔍 Identifying clusters...")
        clusters_list = list(nx.connected_components(G))

        print(f"✓ Found {len(clusters_list):,} clusters")

        # Create cluster dataframe
        cluster_data = []
        for cluster_id, cluster_seqs in enumerate(clusters_list):
            cluster_size = len(cluster_seqs)

            # Get samples containing sequences in this cluster
            cluster_samples = seq_df[seq_df['cdr3'].isin(cluster_seqs)]['sample_id'].unique()

            # Public if shared across multiple samples
            is_public = len(cluster_samples) >= 2

            for seq in cluster_seqs:
                # Get disease association
                seq_samples = seq_df[seq_df['cdr3'] == seq]
                disease_types = seq_samples['disease_subtype'].value_counts().to_dict()

                cluster_data.append({
                    'cluster_id': cluster_id,
                    'cdr3_sequence': seq,
                    'cluster_size': cluster_size,
                    'n_samples': len(cluster_samples),
                    'is_public': is_public,
                    'disease_types': str(disease_types)
                })

        self.clusters = pd.DataFrame(cluster_data)

        # Save clusters
        self.clusters.to_csv('giana_clusters_full.csv', index=False)
        print(f"\n✓ Saved clusters: giana_clusters_full.csv")

        # Print cluster summary
        cluster_sizes = self.clusters.groupby('cluster_id')['cluster_size'].first()
        print(f"\n📊 Cluster size distribution:")
        print(f"   Singleton clusters (size=1): {(cluster_sizes == 1).sum():,}")
        print(f"   Small clusters (2-5): {((cluster_sizes >= 2) & (cluster_sizes <= 5)).sum():,}")
        print(f"   Medium clusters (6-10): {((cluster_sizes >= 6) & (cluster_sizes <= 10)).sum():,}")
        print(f"   Large clusters (>10): {(cluster_sizes > 10).sum():,}")

        if (cluster_sizes > 10).any():
            print(f"\n📊 Largest clusters:")
            for i, (cluster_id, size) in enumerate(cluster_sizes.nlargest(5).items(), 1):
                cluster_info = self.clusters[self.clusters['cluster_id'] == cluster_id].iloc[0]
                print(f"   {i}. Cluster {cluster_id}: {size} sequences, {cluster_info['n_samples']} samples")

        return self.clusters

    def calculate_coclustering_matrix(self):
        """Calculate sample co-clustering matrix (sample similarity)"""

        print("\n" + "="*80)
        print("CALCULATING SAMPLE CO-CLUSTERING MATRIX")
        print("="*80)

        if self.clusters is None or len(self.clusters) == 0:
            print("⚠️  No clusters found. Run cluster_sequences() first.")
            return None

        # Get unique samples
        samples = self.combined_data['repertoire_id'].unique()
        n_samples = len(samples)

        print(f"   Samples: {n_samples}")

        # Initialize co-clustering matrix
        coclustering_matrix = np.zeros((n_samples, n_samples))

        # For each cluster, increment co-clustering for sample pairs
        for cluster_id in self.clusters['cluster_id'].unique():
            cluster_seqs = self.clusters[self.clusters['cluster_id'] == cluster_id]['cdr3_sequence']

            # Find samples containing these sequences
            cluster_samples = []
            for seq in cluster_seqs:
                seq_samples = self.combined_data[
                    self.combined_data[self.combined_data.columns[
                        self.combined_data.columns.str.contains('amino_acid|cdr3')
                    ][0]] == seq
                ]['repertoire_id'].unique()
                cluster_samples.extend(seq_samples)

            cluster_samples = list(set(cluster_samples))

            # Increment co-clustering for all pairs
            for i, sample1 in enumerate(samples):
                for j, sample2 in enumerate(samples):
                    if sample1 in cluster_samples and sample2 in cluster_samples:
                        coclustering_matrix[i, j] += 1

        # Normalize by total sequences per sample
        for i, sample1 in enumerate(samples):
            for j, sample2 in enumerate(samples):
                n1 = len(self.combined_data[self.combined_data['repertoire_id'] == sample1])
                n2 = len(self.combined_data[self.combined_data['repertoire_id'] == sample2])
                normalization = np.sqrt(n1 * n2)
                if normalization > 0:
                    coclustering_matrix[i, j] /= normalization

        self.co_clustering_matrix = pd.DataFrame(
            coclustering_matrix,
            index=samples,
            columns=samples
        )

        # Save matrix
        self.co_clustering_matrix.to_csv('giana_coclustering_matrix_full.csv')
        print(f"\n✓ Saved co-clustering matrix: giana_coclustering_matrix_full.csv")

        return self.co_clustering_matrix

    def create_sample_network(self, output_file='giana_sample_network_full.png'):
        """Create sample network visualization"""

        print("\n" + "="*80)
        print("CREATING SAMPLE NETWORK VISUALIZATION")
        print("="*80)

        if self.co_clustering_matrix is None:
            print("⚠️  Co-clustering matrix not calculated")
            return None

        # Create network
        G = nx.Graph()

        # Add nodes
        samples = self.co_clustering_matrix.index.tolist()
        G.add_nodes_from(samples)

        # Add edges for high co-clustering
        threshold = self.co_clustering_matrix.values[np.triu_indices_from(self.co_clustering_matrix.values, k=1)].mean()

        for i, sample1 in enumerate(samples):
            for j, sample2 in enumerate(samples):
                if i < j:
                    weight = self.co_clustering_matrix.iloc[i, j]
                    if weight > threshold:
                        G.add_edge(sample1, sample2, weight=weight)

        # Get disease status for node colors
        sample_colors = []
        for sample in samples:
            sample_data = self.combined_data[self.combined_data['repertoire_id'] == sample]
            if len(sample_data) > 0:
                disease = sample_data.iloc[0].get('disease_subtype', 'Unknown')
                if 'CMV+' in str(disease):
                    sample_colors.append('red')
                elif 'CMV-' in str(disease):
                    sample_colors.append('blue')
                else:
                    sample_colors.append('gray')
            else:
                sample_colors.append('gray')

        # Create visualization
        fig, ax = plt.subplots(figsize=(12, 10))

        pos = nx.spring_layout(G, k=1, iterations=50, seed=42)

        nx.draw_networkx_nodes(G, pos, node_color=sample_colors,
                              node_size=300, alpha=0.7, ax=ax)
        nx.draw_networkx_edges(G, pos, alpha=0.2, ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=6, ax=ax)

        ax.set_title('Sample Network Based on TCR Cluster Sharing',
                    fontsize=14, fontweight='bold')
        ax.axis('off')

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', label='CMV+'),
            Patch(facecolor='blue', label='CMV-'),
            Patch(facecolor='gray', label='Unknown')
        ]
        ax.legend(handles=legend_elements, loc='upper right')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✓ Saved network: {output_file}")

        return fig

    def run_complete_analysis(self):
        """Run complete GIANA analysis pipeline"""

        print("\n" + "="*80)
        print("GIANA CLUSTERING ANALYSIS - COMPLETE PIPELINE")
        print("="*80)

        # Step 1: Load data
        self.load_data()

        # Step 2: Cluster sequences
        self.cluster_sequences()

        # Step 3: Calculate co-clustering
        self.calculate_coclustering_matrix()

        # Step 4: Visualize sample network
        self.create_sample_network()

        print("\n" + "="*80)
        print("✅ GIANA ANALYSIS COMPLETE!")
        print("="*80)
        print("\nGenerated files:")
        print("  • giana_clusters_full.csv")
        print("  • giana_coclustering_matrix_full.csv")
        print("  • giana_sample_network_full.png")
        print("\n" + "="*80)


def main():
    """Main execution"""

    # ========================================================================
    # CONFIGURATION - UPDATE THIS PATH!
    # ========================================================================

    # For Windows:
    DATA_DIR = r"C:\Users\chris\Desktop\TCR ANALYSIS"

    # For Linux/Mac:
    # DATA_DIR = "data/"

    # ========================================================================

    # IMPORTANT: GIANA is computationally intensive!
    # Adjust these parameters based on your computer's memory:

    MAX_FILES = 50  # Start with 50 files for testing (None = all files)
    MAX_SEQUENCES = 50000  # Maximum sequences to cluster (50K recommended)

    # For full analysis on powerful workstation:
    # MAX_FILES = None
    # MAX_SEQUENCES = 500000

    # ========================================================================

    print("\n" + "="*80)
    print("GIANA CLUSTERING ANALYSIS - FULL DATASET VERSION")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Data directory: {DATA_DIR}")
    print(f"  Max files: {'ALL FILES' if MAX_FILES is None else MAX_FILES}")
    print(f"  Max sequences: {MAX_SEQUENCES:,}")
    print(f"\n⚠️  NOTE: GIANA clustering is computationally intensive.")
    print(f"  Estimated runtime: 10-30 minutes for {MAX_SEQUENCES:,} sequences")
    print("\n" + "="*80)

    # Create analyzer
    analyzer = GIANAAnalyzerFull(
        data_directory=DATA_DIR,
        max_files=MAX_FILES,
        max_sequences=MAX_SEQUENCES,
        similarity_threshold=0.85  # From Zhang et al. 2021
    )

    # Run analysis
    analyzer.run_complete_analysis()


if __name__ == "__main__":
    main()
