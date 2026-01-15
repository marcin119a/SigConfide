#!/usr/bin/env python3
"""
Main script for running SigConfide on test data.

This script demonstrates how to use SigConfide to analyze mutational signatures
on sample data from the tests/data directory.

Usage:
    # Run with default examples
    python main.py

    # Run with custom parameters
    python main.py --samples tests/data/reduced_data.dat --output output/custom --signatures 3.4
"""

import os
import sys
import argparse
from pathlib import Path

# Add the project root to the path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from sigconfide.modelselection.analyzer import fit


def run_single_analysis(samples_file, output_dir, signatures, threshold=0.01, 
                       mutation_count=None, R=100, significance_level=0.01, 
                       drop_zeros_columns=False):
    """Run a single analysis with specified parameters."""
    print(f"Input file: {samples_file}")
    print(f"Output directory: {output_dir}")
    
    if isinstance(signatures, (int, float)):
        print(f"Signatures: COSMIC v{signatures}")
    else:
        print(f"Signatures file: {signatures}")
    
    print(f"Parameters: threshold={threshold}, R={R}, significance_level={significance_level}")
    if mutation_count is not None:
        print(f"  mutation_count={mutation_count}")
    print(f"  drop_zeros_columns={drop_zeros_columns}")
    print()
    
    try:
        fit(
            str(samples_file),
            str(output_dir),
            signatures=signatures,
            threshold=threshold,
            mutation_count=mutation_count,
            R=R,
            significance_level=significance_level,
            drop_zeros_columns=drop_zeros_columns
        )
        print(f"\n✓ Analysis completed successfully!")
        print(f"  Results saved to: {Path(output_dir) / 'Assignment_Solution_Activities.csv'}")
        return True
    except Exception as e:
        print(f"\n✗ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return False


def run_examples():
    """Run example analyses with test data."""
    # Get the directory of this script
    script_dir = Path(__file__).parent
    test_data_dir = script_dir / "tests" / "data"
    output_base_dir = script_dir / "output"
    
    print("=" * 70)
    print("SigConfide - Analysis of Mutational Signatures")
    print("Running example analyses...")
    print("=" * 70)
    print()
    
    # Example 1: Using COSMIC signatures version 2.0
    print("Example 1: Analyzing with COSMIC signatures v2.0")
    print("-" * 70)
    samples_file = test_data_dir / "reduced_data.dat"
    output_dir = output_base_dir / "cosmic_v2"
    
    if samples_file.exists():
        run_single_analysis(
            samples_file, output_dir, signatures=2.0,
            threshold=0.01, R=100, significance_level=0.01,
            drop_zeros_columns=False
        )
    else:
        print(f"✗ Test data file not found: {samples_file}")
    
    print()
    print()
    
    # Example 2: Using custom Breast Signatures
    print("Example 2: Analyzing with custom Breast Signatures")
    print("-" * 70)
    samples_file = test_data_dir / "reduced_data.dat"
    signatures_file = test_data_dir / "Breast_Signatures.csv"
    output_dir = output_base_dir / "breast_signatures"
    
    if samples_file.exists() and signatures_file.exists():
        run_single_analysis(
            samples_file, output_dir, signatures=str(signatures_file),
            threshold=0.01, R=100, significance_level=0.01,
            drop_zeros_columns=True
        )
    else:
        missing = []
        if not samples_file.exists():
            missing.append(f"  - {samples_file}")
        if not signatures_file.exists():
            missing.append(f"  - {signatures_file}")
        print(f"✗ Test data files not found:")
        for f in missing:
            print(f)
    
    print()
    print()
    
    # Example 3: Using COSMIC signatures version 3.4 (default)
    print("Example 3: Analyzing with COSMIC signatures v3.4 (default)")
    print("-" * 70)
    samples_file = test_data_dir / "reduced_data.dat"
    output_dir = output_base_dir / "cosmic_v3.4"
    
    if samples_file.exists():
        run_single_analysis(
            samples_file, output_dir, signatures=3.4,
            threshold=0.01, R=100, significance_level=0.01,
            mutation_count=None, drop_zeros_columns=True
        )
    else:
        print(f"✗ Test data file not found: {samples_file}")
    
    print()
    print("=" * 70)
    print("All analyses completed!")
    print("=" * 70)


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Run SigConfide analysis on mutational signature data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default examples
  python main.py

  # Run with custom parameters
  python main.py --samples tests/data/reduced_data.dat --output output/custom --signatures 3.4

  # Use custom signatures file
  python main.py --samples tests/data/reduced_data.dat --output output/custom \\
                 --signatures tests/data/Breast_Signatures.csv

  # Specify all parameters
  python main.py --samples tests/data/reduced_data.dat --output output/custom \\
                 --signatures 2.0 --threshold 0.02 --R 50 --significance 0.05
        """
    )
    
    parser.add_argument(
        '--samples',
        type=str,
        help='Path to samples file (mutational matrix)'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        help='Path to output directory'
    )
    
    parser.add_argument(
        '--signatures',
        type=str,
        help='COSMIC version (1.0, 2.0, 3.0, 3.1, 3.4) or path to custom signatures file'
    )
    
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.01,
        help='Threshold for significant exposures (default: 0.01)'
    )
    
    parser.add_argument(
        '--mutation-count',
        type=int,
        default=None,
        help='Total mutation count (default: None, computed from data)'
    )
    
    parser.add_argument(
        '--R',
        type=int,
        default=100,
        help='Number of bootstrap replicates (default: 100)'
    )
    
    parser.add_argument(
        '--significance',
        type=float,
        default=0.01,
        dest='significance_level',
        help='Statistical significance level (default: 0.01)'
    )
    
    parser.add_argument(
        '--drop-zeros',
        action='store_true',
        help='Drop columns with all zero values from output'
    )
    
    return parser.parse_args()


def main():
    """Main function to run SigConfide analysis."""
    args = parse_arguments()
    
    # If no arguments provided, run examples
    if args.samples is None:
        run_examples()
        return
    
    # Otherwise, run with provided arguments
    if args.output is None:
        print("Error: --output is required when --samples is provided")
        sys.exit(1)
    
    if args.signatures is None:
        print("Error: --signatures is required when --samples is provided")
        sys.exit(1)
    
    # Convert signatures to float if it's a numeric version
    signatures = args.signatures
    try:
        version = float(signatures)
        if version in [1.0, 2.0, 3.0, 3.1, 3.4]:
            signatures = version
    except ValueError:
        # Not a version number, treat as file path
        pass
    
    samples_path = Path(args.samples)
    if not samples_path.exists():
        print(f"Error: Samples file not found: {samples_path}")
        sys.exit(1)
    
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("SigConfide - Analysis of Mutational Signatures")
    print("=" * 70)
    print()
    
    success = run_single_analysis(
        samples_path,
        output_path,
        signatures=signatures,
        threshold=args.threshold,
        mutation_count=args.mutation_count,
        R=args.R,
        significance_level=args.significance_level,
        drop_zeros_columns=args.drop_zeros
    )
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
