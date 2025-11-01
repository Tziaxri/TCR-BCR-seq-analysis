#!/usr/bin/env python3
"""
TCR Data Loader Module
Loads all TCR repertoire files and metadata from user's directory

Author: Claude Code
Date: 2025-10-31
"""

import pandas as pd
import numpy as np
import os
import glob
from pathlib import Path
import gzip
from typing import Dict, List, Tuple

class TCRDataLoader:
    """
    Comprehensive data loader for TCR repertoire analysis

    Handles:
    - Multiple .tsv.gz TCR repertoire files
    - metadata.tsv (participant information)
    - synapse_metadata_manifest.tsv (file manifest)
    - Cross-platform path handling (Windows/Linux)
    """

    def __init__(self, data_directory: str):
        """
        Initialize data loader

        Parameters:
        -----------
        data_directory : str
            Path to directory containing TCR files and metadata
            Example: "C:\\Users\\chris\\Desktop\\TCR ANALYSIS\\"
                 or: "data/" (for Linux)
        """
        self.data_dir = Path(data_directory)
        self.metadata = None
        self.manifest = None
        self.tcr_data = {}  # Dictionary: {sample_id: DataFrame}

        print(f"📁 TCR Data Loader initialized")
        print(f"   Directory: {self.data_dir}")

    def load_metadata(self, metadata_file: str = "metadata.tsv") -> pd.DataFrame:
        """
        Load participant metadata

        Expected columns:
        - participant_id
        - repertoire_id
        - disease_subtype (CMV+, CMV-, Unknown)
        - age, sex, ancestry, etc.
        """
        metadata_path = self.data_dir / metadata_file

        print(f"\n📋 Loading metadata...")
        print(f"   File: {metadata_path}")

        if not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

        self.metadata = pd.read_csv(metadata_path, sep='\t')

        print(f"   ✓ Loaded {len(self.metadata)} samples")
        print(f"   Columns: {', '.join(self.metadata.columns[:5])}...")

        # Display CMV status distribution
        if 'disease_subtype' in self.metadata.columns:
            cmv_counts = self.metadata['disease_subtype'].value_counts()
            print(f"\n   CMV Status Distribution:")
            for status, count in cmv_counts.items():
                print(f"     • {status}: {count} samples")

        return self.metadata

    def load_manifest(self, manifest_file: str = "synapse_metadata_manifest.tsv") -> pd.DataFrame:
        """
        Load Synapse metadata manifest (file-level metadata)

        This file typically contains:
        - File IDs
        - File paths/names
        - Sample mappings
        - Additional file-level annotations
        """
        manifest_path = self.data_dir / manifest_file

        print(f"\n📄 Loading manifest...")
        print(f"   File: {manifest_path}")

        if not manifest_path.exists():
            print(f"   ⚠️  Manifest file not found: {manifest_path}")
            print(f"   Continuing without manifest...")
            return None

        self.manifest = pd.read_csv(manifest_path, sep='\t')

        print(f"   ✓ Loaded {len(self.manifest)} file records")
        print(f"   Columns: {', '.join(self.manifest.columns[:5])}...")

        return self.manifest

    def find_tcr_files(self, pattern: str = "*.tsv.gz") -> List[Path]:
        """
        Find all TCR repertoire files matching pattern

        Parameters:
        -----------
        pattern : str
            Glob pattern for TCR files (default: "*.tsv.gz")
            Examples: "part_table_*.tsv.gz", "*.tsv", "*.tsv.gz"

        Returns:
        --------
        List[Path] : List of file paths
        """
        print(f"\n🔍 Searching for TCR files...")
        print(f"   Directory: {self.data_dir}")
        print(f"   Pattern: {pattern}")

        # Find all matching files
        files = list(self.data_dir.glob(pattern))

        # Exclude metadata files
        excluded_patterns = ['metadata', 'manifest', 'synapse']
        files = [f for f in files if not any(pattern in f.name.lower() for pattern in excluded_patterns)]

        # Sort by filename for consistent processing
        files = sorted(files)

        print(f"   ✓ Found {len(files)} TCR files")

        if len(files) > 0:
            print(f"\n   Sample files:")
            for f in files[:5]:
                print(f"     • {f.name}")
            if len(files) > 5:
                print(f"     • ... and {len(files) - 5} more files")

        return files

    def load_single_tcr_file(self, filepath: Path) -> pd.DataFrame:
        """
        Load a single TCR repertoire file

        Handles both .tsv and .tsv.gz formats
        """
        try:
            # Determine if file is gzipped
            if filepath.suffix == '.gz':
                with gzip.open(filepath, 'rt') as f:
                    df = pd.read_csv(f, sep='\t')
            else:
                df = pd.read_csv(filepath, sep='\t')

            # Add source filename
            df['source_file'] = filepath.name

            return df

        except Exception as e:
            print(f"   ⚠️  Error loading {filepath.name}: {e}")
            return None

    def load_all_tcr_files(self, pattern: str = "*.tsv.gz",
                          max_files: int = None,
                          verbose: bool = True) -> Dict[str, pd.DataFrame]:
        """
        Load all TCR repertoire files

        Parameters:
        -----------
        pattern : str
            Glob pattern for files
        max_files : int
            Maximum number of files to load (None = all)
        verbose : bool
            Print progress messages

        Returns:
        --------
        Dict[str, pd.DataFrame] : {filename: DataFrame}
        """
        files = self.find_tcr_files(pattern)

        if max_files:
            files = files[:max_files]
            print(f"\n   ⚠️  Limiting to first {max_files} files")

        print(f"\n📊 Loading {len(files)} TCR repertoire files...")

        self.tcr_data = {}
        loaded_count = 0
        error_count = 0

        for i, filepath in enumerate(files, 1):
            if verbose and i % 10 == 0:
                print(f"   Progress: {i}/{len(files)} files loaded...")

            df = self.load_single_tcr_file(filepath)

            if df is not None and len(df) > 0:
                self.tcr_data[filepath.name] = df
                loaded_count += 1
            else:
                error_count += 1

        print(f"\n   ✓ Successfully loaded: {loaded_count} files")
        if error_count > 0:
            print(f"   ⚠️  Errors: {error_count} files")

        return self.tcr_data

    def extract_sample_id_from_filename(self, filename: str) -> str:
        """
        Extract sample/repertoire ID from filename

        Examples:
        - "part_table_bfi-0000234.tsv.gz" -> "0000234"
        - "part_table_bfi-0003053.tsv.gz" -> "0003053"
        """
        # Try to extract numeric ID
        import re
        match = re.search(r'(\d{7})', filename)
        if match:
            return match.group(1)

        # Otherwise use filename without extension
        return filename.replace('.tsv.gz', '').replace('.tsv', '')

    def merge_with_metadata(self) -> pd.DataFrame:
        """
        Merge TCR data with metadata

        Returns:
        --------
        pd.DataFrame : Combined dataframe with TCR data + metadata
        """
        if self.metadata is None:
            raise ValueError("Metadata not loaded. Call load_metadata() first.")

        if not self.tcr_data:
            raise ValueError("TCR data not loaded. Call load_all_tcr_files() first.")

        print(f"\n🔗 Merging TCR data with metadata...")

        combined_data = []

        for filename, tcr_df in self.tcr_data.items():
            # Extract sample ID from filename
            sample_id = self.extract_sample_id_from_filename(filename)

            # Find matching metadata
            # Try multiple possible column names for repertoire ID
            possible_id_columns = ['repertoire_id', 'sample_id', 'file_id', 'bfi']

            metadata_row = None
            for col in possible_id_columns:
                if col in self.metadata.columns:
                    # Try exact match
                    matches = self.metadata[self.metadata[col].astype(str).str.contains(sample_id, na=False)]
                    if len(matches) > 0:
                        metadata_row = matches.iloc[0]
                        break

            if metadata_row is None:
                # Try matching by index
                print(f"   ⚠️  No metadata match for {filename}, using defaults")
                # Create default metadata
                metadata_row = pd.Series({
                    'repertoire_id': sample_id,
                    'disease_subtype': 'Unknown',
                    'source_file': filename
                })

            # Add metadata columns to TCR data
            for col in metadata_row.index:
                if col not in tcr_df.columns:
                    tcr_df[col] = metadata_row[col]

            combined_data.append(tcr_df)

        # Combine all dataframes
        combined_df = pd.concat(combined_data, ignore_index=True)

        print(f"   ✓ Combined {len(self.tcr_data)} files")
        print(f"   ✓ Total TCR sequences: {len(combined_df):,}")

        return combined_df

    def calculate_diversity_per_sample(self) -> pd.DataFrame:
        """
        Calculate diversity metrics for each sample

        Returns:
        --------
        pd.DataFrame : Diversity metrics per sample
        """
        from tcr_utils import calculate_diversity_metrics

        if not self.tcr_data:
            raise ValueError("TCR data not loaded. Call load_all_tcr_files() first.")

        print(f"\n📊 Calculating diversity metrics per sample...")

        diversity_results = []

        for filename, tcr_df in self.tcr_data.items():
            sample_id = self.extract_sample_id_from_filename(filename)

            # Filter to productive sequences
            if 'productive' in tcr_df.columns:
                tcr_df = tcr_df[tcr_df['productive'] == True].copy()

            # Get clone counts
            if 'duplicate_count' in tcr_df.columns:
                clone_counts = tcr_df['duplicate_count'].values
            elif 'templates' in tcr_df.columns:
                clone_counts = tcr_df['templates'].values
            else:
                # Assume each row is a unique clone with count=1
                clone_counts = np.ones(len(tcr_df))

            # Calculate diversity
            metrics = calculate_diversity_metrics(clone_counts)
            metrics['sample_id'] = sample_id
            metrics['repertoire_id'] = sample_id
            metrics['source_file'] = filename
            metrics['total_sequences'] = len(tcr_df)
            metrics['total_reads'] = clone_counts.sum()

            diversity_results.append(metrics)

        diversity_df = pd.DataFrame(diversity_results)

        # Merge with metadata
        if self.metadata is not None:
            # Try to merge on multiple possible columns
            for id_col in ['repertoire_id', 'sample_id', 'bfi']:
                if id_col in self.metadata.columns:
                    diversity_df = diversity_df.merge(
                        self.metadata,
                        left_on='repertoire_id',
                        right_on=id_col,
                        how='left'
                    )
                    break

        print(f"   ✓ Calculated metrics for {len(diversity_df)} samples")

        return diversity_df

    def get_summary_statistics(self) -> Dict:
        """
        Get summary statistics for the entire dataset
        """
        if not self.tcr_data:
            return {}

        total_files = len(self.tcr_data)
        total_sequences = sum(len(df) for df in self.tcr_data.values())

        # Get productive sequences
        total_productive = 0
        for df in self.tcr_data.values():
            if 'productive' in df.columns:
                total_productive += df['productive'].sum()
            else:
                total_productive += len(df)

        # Get unique CDR3 sequences
        all_cdr3 = []
        for df in self.tcr_data.values():
            if 'amino_acid' in df.columns:
                all_cdr3.extend(df['amino_acid'].dropna().tolist())
            elif 'cdr3_aa' in df.columns:
                all_cdr3.extend(df['cdr3_aa'].dropna().tolist())

        unique_cdr3 = len(set(all_cdr3))

        stats = {
            'total_files': total_files,
            'total_sequences': total_sequences,
            'total_productive': total_productive,
            'unique_cdr3': unique_cdr3,
            'avg_sequences_per_file': total_sequences / total_files if total_files > 0 else 0
        }

        # Add metadata stats
        if self.metadata is not None:
            if 'disease_subtype' in self.metadata.columns:
                stats['cmv_distribution'] = self.metadata['disease_subtype'].value_counts().to_dict()

        return stats

    def print_summary(self):
        """
        Print comprehensive dataset summary
        """
        stats = self.get_summary_statistics()

        print("\n" + "="*70)
        print("📊 DATASET SUMMARY")
        print("="*70)

        print(f"\n📁 Data Files:")
        print(f"   • Total TCR files: {stats.get('total_files', 0):,}")
        print(f"   • Total sequences: {stats.get('total_sequences', 0):,}")
        print(f"   • Productive sequences: {stats.get('total_productive', 0):,}")
        print(f"   • Unique CDR3 sequences: {stats.get('unique_cdr3', 0):,}")
        print(f"   • Avg sequences/file: {stats.get('avg_sequences_per_file', 0):,.1f}")

        if 'cmv_distribution' in stats:
            print(f"\n🔬 Sample Distribution:")
            for status, count in stats['cmv_distribution'].items():
                print(f"   • {status}: {count} samples")

        if self.metadata is not None:
            print(f"\n👥 Participant Metadata:")
            print(f"   • Total participants: {len(self.metadata)}")
            if 'age' in self.metadata.columns:
                print(f"   • Age range: {self.metadata['age'].min():.0f}-{self.metadata['age'].max():.0f} years")
            if 'sex' in self.metadata.columns:
                sex_dist = self.metadata['sex'].value_counts()
                print(f"   • Sex distribution: {dict(sex_dist)}")

        print("\n" + "="*70 + "\n")


def load_complete_dataset(data_directory: str,
                          max_files: int = None,
                          verbose: bool = True) -> Tuple[pd.DataFrame, TCRDataLoader]:
    """
    Convenience function to load complete dataset

    Parameters:
    -----------
    data_directory : str
        Path to data directory
        Example: "C:\\Users\\chris\\Desktop\\TCR ANALYSIS\\"
    max_files : int
        Maximum files to load (None = all)
    verbose : bool
        Print progress

    Returns:
    --------
    Tuple[pd.DataFrame, TCRDataLoader]
        - diversity_df: Diversity metrics per sample
        - loader: TCRDataLoader object for additional analysis

    Example:
    --------
    >>> diversity_df, loader = load_complete_dataset(
    ...     "C:\\Users\\chris\\Desktop\\TCR ANALYSIS\\",
    ...     max_files=50  # Load first 50 files for testing
    ... )
    >>> loader.print_summary()
    """
    # Create loader
    loader = TCRDataLoader(data_directory)

    # Load metadata
    try:
        loader.load_metadata("metadata.tsv")
    except FileNotFoundError:
        print("⚠️  metadata.tsv not found, continuing without it...")

    # Load manifest (optional)
    try:
        loader.load_manifest("synapse_metadata_manifest.tsv")
    except:
        pass

    # Load all TCR files
    loader.load_all_tcr_files(pattern="*.tsv.gz", max_files=max_files, verbose=verbose)

    # Calculate diversity metrics
    diversity_df = loader.calculate_diversity_per_sample()

    # Print summary
    if verbose:
        loader.print_summary()

    return diversity_df, loader


if __name__ == "__main__":
    """
    Example usage and testing
    """
    print("="*70)
    print("TCR DATA LOADER - Example Usage")
    print("="*70)

    # Example 1: Load from Windows directory
    print("\n📌 Example: Loading from Windows directory")
    print("   data_dir = r'C:\\Users\\chris\\Desktop\\TCR ANALYSIS\\'")

    # Example 2: Load from Linux directory (for testing)
    print("\n📌 Example: Loading from Linux directory")
    data_dir = "data/"

    if os.path.exists(data_dir):
        diversity_df, loader = load_complete_dataset(
            data_dir,
            max_files=5,  # Load only 5 files for testing
            verbose=True
        )

        print("\n✓ Diversity metrics calculated:")
        print(diversity_df.head())
    else:
        print(f"⚠️  Directory '{data_dir}' not found")
        print("\n💡 To use this module:")
        print("   1. Update data_directory path to your location")
        print("   2. Run: python tcr_data_loader.py")
