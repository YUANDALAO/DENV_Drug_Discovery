"""
QSAR Distribution Analysis Tool for REINVENT4 Results

Usage:
    python scripts/utils/QSAR_distribution.py <results_dir>
    
Example:
    python scripts/utils/QSAR_distribution.py experiments/runs/run13b
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import re


def clean_smiles_for_rdkit(smiles):
    """
    Clean SMILES string for RDKit parsing
    Remove LibInvent/LinkInvent attachment points and decorations
    """
    if pd.isna(smiles):
        return None
    
    # Remove everything after | (attachment point decorations)
    smiles = smiles.split('|')[0]
    
    # Replace [*] with [H] for valid parsing
    smiles = smiles.replace('[*]', '[H]')
    
    return smiles


def find_csv_file(results_dir):
    """Find the main results CSV file in the directory"""
    csv_files = list(results_dir.glob("*.csv"))
    
    # Priority order: look for common patterns
    patterns = ["scaffold", "results", "output", "summary", "memory"]
    
    for pattern in patterns:
        matching = [f for f in csv_files if pattern in f.name.lower()]
        if matching:
            return matching[0]
    
    # If no pattern match, return the largest CSV file
    if csv_files:
        return max(csv_files, key=lambda f: f.stat().st_size)
    
    return None


def analyze_qsar_distribution(results_dir):
    """
    Analyze QSAR score distribution and training progress
    
    Args:
        results_dir: Path to the results directory containing CSV files
    """
    results_dir = Path(results_dir)
    
    if not results_dir.exists():
        print(f"❌ Error: Directory not found: {results_dir}")
        sys.exit(1)
    
    # Find CSV file
    csv_file = find_csv_file(results_dir)
    
    if csv_file is None:
        print(f"❌ Error: No CSV files found in {results_dir}")
        sys.exit(1)
    
    print(f"📊 Analyzing: {csv_file.name}")
    
    # Read data
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"❌ Error reading CSV: {e}")
        sys.exit(1)
    
    print(f"✓ Loaded {len(df):,} rows")
    print(f"✓ Columns: {df.columns.tolist()}")
    
    # Identify QSAR column (flexible naming)
    qsar_col = None
    possible_names = [
        'DENV_Activity (raw)', 'QSAR Score', 
        'QSAR_Score', 'pIC50', 'Activity'
    ]
    for name in possible_names:
        if name in df.columns:
            qsar_col = name
            break
    
    if qsar_col is None:
        print(f"❌ Error: No QSAR score column found. Available columns: {df.columns.tolist()}")
        sys.exit(1)
    
    print(f"✓ Found QSAR column: '{qsar_col}'")
    
    # Identify other key columns
    step_col = 'step' if 'step' in df.columns else None
    smiles_col = 'SMILES' if 'SMILES' in df.columns else None
    score_col = 'Score' if 'Score' in df.columns else ('total_score' if 'total_score' in df.columns else None)
    
    # Calculate statistics
    print("\n" + "="*70)
    print("QSAR SCORE STATISTICS")
    print("="*70)
    qsar_values = df[qsar_col].dropna()
    print(f"Total rows:     {len(df):,}")
    print(f"Valid values:   {len(qsar_values):,}")
    print(f"Missing:        {df[qsar_col].isna().sum():,}")
    print(f"\nMin:            {qsar_values.min():.4f}")
    print(f"Q1 (25%):       {qsar_values.quantile(0.25):.4f}")
    print(f"Median (50%):   {qsar_values.median():.4f}")
    print(f"Q3 (75%):       {qsar_values.quantile(0.75):.4f}")
    print(f"Max:            {qsar_values.max():.4f}")
    print(f"Mean:           {qsar_values.mean():.4f}")
    print(f"Std:            {qsar_values.std():.4f}")
    
    # Count molecules in different pIC50 ranges
    print(f"\n--- Activity Distribution ---")
    ranges = [(0, 5), (5, 6), (6, 7), (7, 8), (8, 100)]
    for low, high in ranges:
        count = ((qsar_values >= low) & (qsar_values < high)).sum()
        pct = count / len(qsar_values) * 100
        print(f"pIC50 [{low}-{high}):  {count:6,}  ({pct:5.1f}%)")
    
    # Total Score statistics if available
    if score_col:
        print("\n" + "="*70)
        print("TOTAL SCORE STATISTICS")
        print("="*70)
        score_values = df[score_col].dropna()
        print(f"Min:            {score_values.min():.4f}")
        print(f"Q1 (25%):       {score_values.quantile(0.25):.4f}")
        print(f"Median (50%):   {score_values.median():.4f}")
        print(f"Q3 (75%):       {score_values.quantile(0.75):.4f}")
        print(f"Max:            {score_values.max():.4f}")
        print(f"Mean:           {score_values.mean():.4f}")
        print(f"Std:            {score_values.std():.4f}")
    
    # Scaffold diversity analysis if SMILES available
    if smiles_col:
        print("\n" + "="*70)
        print("SCAFFOLD DIVERSITY ANALYSIS")
        print("="*70)
        print("(Cleaning SMILES for RDKit parsing...)")
        
        scaffolds = []
        valid_smiles = 0
        parse_errors = 0
        
        for smi in df[smiles_col]:
            try:
                cleaned_smi = clean_smiles_for_rdkit(smi)
                if cleaned_smi:
                    mol = Chem.MolFromSmiles(cleaned_smi)
                    if mol:
                        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                        scaffolds.append(Chem.MolToSmiles(scaffold))
                        valid_smiles += 1
                    else:
                        parse_errors += 1
            except:
                parse_errors += 1
                continue
        
        unique_scaffolds = len(set(scaffolds))
        diversity_ratio = unique_scaffolds / valid_smiles if valid_smiles > 0 else 0
        
        print(f"Total SMILES:       {len(df):,}")
        print(f"Valid parsed:       {valid_smiles:,}")
        print(f"Parse errors:       {parse_errors:,}")
        print(f"Unique scaffolds:   {unique_scaffolds:,}")
        print(f"Diversity ratio:    {diversity_ratio:.2%}")
        
        # Mode collapse warning
        if diversity_ratio < 0.1:
            print("\n⚠️  WARNING: Low scaffold diversity (<10%)!")
            print("   This indicates potential mode collapse.")
        elif diversity_ratio < 0.3:
            print("\n⚠️  CAUTION: Moderate scaffold diversity (10-30%)")
            print("   Consider monitoring for mode collapse.")
        else:
            print("\n✓ Good scaffold diversity (>30%)")
    
    # Create visualization
    fig = plt.figure(figsize=(16, 5))
    
    # Plot 1: QSAR Score Distribution
    ax1 = plt.subplot(131)
    n, bins, patches = ax1.hist(qsar_values, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.axvline(qsar_values.mean(), color='red', linestyle='--', 
                linewidth=2, label=f'Mean: {qsar_values.mean():.2f}')
    ax1.axvline(qsar_values.median(), color='green', linestyle='--', 
                linewidth=2, label=f'Median: {qsar_values.median():.2f}')
    ax1.set_xlabel(f'{qsar_col}', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax1.set_title(f'{qsar_col} Distribution', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(alpha=0.3, linestyle='--')
    
    # Plot 2: QSAR Score over Training
    if step_col:
        ax2 = plt.subplot(132)
        ax2.scatter(df[step_col], df[qsar_col], alpha=0.3, s=1, color='steelblue')
        ax2.set_xlabel('Step', fontsize=12, fontweight='bold')
        ax2.set_ylabel(f'{qsar_col}', fontsize=12, fontweight='bold')
        ax2.set_title(f'{qsar_col} over Training', fontsize=14, fontweight='bold')
        ax2.grid(alpha=0.3, linestyle='--')
        
        # Add moving average
        window = min(100, len(df) // 20)
        if window > 1:
            moving_avg = df[qsar_col].rolling(window=window).mean()
            ax2.plot(df[step_col], moving_avg, color='red', linewidth=2, 
                    label=f'Moving Avg ({window})', alpha=0.8)
            ax2.legend(fontsize=10)
    
    # Plot 3: Total Score over Training
    if step_col and score_col:
        ax3 = plt.subplot(133)
        ax3.scatter(df[step_col], df[score_col], alpha=0.3, s=1, color='steelblue')
        ax3.set_xlabel('Step', fontsize=12, fontweight='bold')
        ax3.set_ylabel('Total Score', fontsize=12, fontweight='bold')
        ax3.set_title('Total Score over Training', fontsize=14, fontweight='bold')
        ax3.grid(alpha=0.3, linestyle='--')
        
        # Add moving average
        if window > 1:
            moving_avg = df[score_col].rolling(window=window).mean()
            ax3.plot(df[step_col], moving_avg, color='red', linewidth=2, 
                    label=f'Moving Avg ({window})', alpha=0.8)
            ax3.legend(fontsize=10)
    
    plt.tight_layout()
    
    # Save figure in the results directory
    output_file = results_dir / f"{results_dir.name}_qsar_analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved: {output_file}")
    
    # Save statistics to text file
    stats_file = results_dir / f"{results_dir.name}_qsar_statistics.txt"
    with open(stats_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write(f"QSAR Analysis for: {csv_file.name}\n")
        f.write(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*70 + "\n\n")
        
        f.write("QSAR SCORE STATISTICS\n")
        f.write("-"*70 + "\n")
        f.write(f"Total rows:     {len(df):,}\n")
        f.write(f"Valid values:   {len(qsar_values):,}\n")
        f.write(f"Missing:        {df[qsar_col].isna().sum():,}\n\n")
        f.write(f"Min:            {qsar_values.min():.4f}\n")
        f.write(f"Q1 (25%):       {qsar_values.quantile(0.25):.4f}\n")
        f.write(f"Median (50%):   {qsar_values.median():.4f}\n")
        f.write(f"Q3 (75%):       {qsar_values.quantile(0.75):.4f}\n")
        f.write(f"Max:            {qsar_values.max():.4f}\n")
        f.write(f"Mean:           {qsar_values.mean():.4f}\n")
        f.write(f"Std:            {qsar_values.std():.4f}\n\n")
        
        f.write("Activity Distribution:\n")
        for low, high in ranges:
            count = ((qsar_values >= low) & (qsar_values < high)).sum()
            pct = count / len(qsar_values) * 100
            f.write(f"  pIC50 [{low}-{high}): {count:6,}  ({pct:5.1f}%)\n")
        
        if score_col:
            f.write("\n" + "="*70 + "\n")
            f.write("TOTAL SCORE STATISTICS\n")
            f.write("-"*70 + "\n")
            score_values = df[score_col].dropna()
            f.write(f"Min:            {score_values.min():.4f}\n")
            f.write(f"Q1 (25%):       {score_values.quantile(0.25):.4f}\n")
            f.write(f"Median (50%):   {score_values.median():.4f}\n")
            f.write(f"Q3 (75%):       {score_values.quantile(0.75):.4f}\n")
            f.write(f"Max:            {score_values.max():.4f}\n")
            f.write(f"Mean:           {score_values.mean():.4f}\n")
            f.write(f"Std:            {score_values.std():.4f}\n")
        
        if smiles_col:
            f.write("\n" + "="*70 + "\n")
            f.write("SCAFFOLD DIVERSITY\n")
            f.write("-"*70 + "\n")
            f.write(f"Total SMILES:       {len(df):,}\n")
            f.write(f"Valid parsed:       {valid_smiles:,}\n")
            f.write(f"Parse errors:       {parse_errors:,}\n")
            f.write(f"Unique scaffolds:   {unique_scaffolds:,}\n")
            f.write(f"Diversity ratio:    {diversity_ratio:.2%}\n")
            
            if diversity_ratio < 0.1:
                f.write("\n⚠️  WARNING: Low scaffold diversity! Potential mode collapse.\n")
            elif diversity_ratio < 0.3:
                f.write("\n⚠️  CAUTION: Moderate scaffold diversity.\n")
            else:
                f.write("\n✓ Good scaffold diversity.\n")
    
    print(f"✓ Statistics saved: {stats_file}")
    
    # Analysis summary
    print("\n" + "="*70)
    print("📊 ANALYSIS SUMMARY")
    print("="*70)
    
    # Check for early convergence
    if step_col and len(df) > 100:
        first_100 = df.head(100)[qsar_col].mean()
        last_100 = df.tail(100)[qsar_col].mean()
        improvement = ((last_100 - first_100) / first_100) * 100
        
        print(f"First 100 steps mean: {first_100:.4f}")
        print(f"Last 100 steps mean:  {last_100:.4f}")
        print(f"Improvement:          {improvement:+.2f}%")
        
        if abs(improvement) < 5:
            print("\n⚠️  WARNING: Limited improvement detected!")
            print("   Model may have converged early or hit a plateau.")
            print("   Consider:")
            print("   - Reducing learning rate")
            print("   - Adjusting sigma parameter")
            print("   - Relaxing diversity filter")
    
    plt.show()


def main():
    if len(sys.argv) != 2:
        print("Usage: python scripts/utils/QSAR_distribution.py <results_dir>")
        print("\nExample:")
        print("  python scripts/utils/QSAR_distribution.py experiments/runs/run13b")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    analyze_qsar_distribution(results_dir)


if __name__ == "__main__":
    main()