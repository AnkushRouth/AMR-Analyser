import sys
import argparse
from Bio import SeqIO
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def build_index(fasta_file, k=8):
    """Build k-mer index from FASTA file and keep full headers."""
    print(f"Building k-mer index (k={k})...")
    index = defaultdict(set)
    gene_lengths = {}
    header_map = {}
    gene_count = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id                       # short identifier (first token)
        full_header = record.description.strip()  # entire FASTA header
        header_map[gene_id] = full_header
        
        seq = str(record.seq).upper()
        gene_lengths[gene_id] = len(seq)
        gene_count += 1
        
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if 'N' not in kmer:  # skip ambiguous bases
                index[kmer].add(gene_id)
    
    print(f"  Indexed {gene_count} genes with {len(index)} unique k-mers")
    return index, gene_lengths, header_map


def pseudoalign_reads(fastq_file, index, k=8):
    """Pseudoalign reads to genes using the k-mer index."""
    print("Pseudoaligning reads...")
    eq_classes = Counter()
    total_reads = 0
    aligned_reads = 0
    
    for record in SeqIO.parse(fastq_file, "fastq"):
        total_reads += 1
        seq = str(record.seq).upper()
        
        if len(seq) < k:
            continue
        
        gene_sets = [
            index[seq[i:i + k]]
            for i in range(len(seq) - k + 1)
            if seq[i:i + k] in index and 'N' not in seq[i:i + k]
        ]
        
        if not gene_sets:
            continue
        
        compatible_genes = set.intersection(*gene_sets)
        if compatible_genes:
            eq_classes[frozenset(compatible_genes)] += 1
            aligned_reads += 1
        
        if total_reads % 100000 == 0:
            print(f"  Processed {total_reads} reads, aligned {aligned_reads}")
    
    alignment_pct = 100 * aligned_reads / total_reads if total_reads else 0.0
    print(f"  Total reads: {total_reads}, Aligned: {aligned_reads} ({alignment_pct:.2f}%)")
    return eq_classes, total_reads, aligned_reads


def quantify(eq_classes, gene_lengths, gene_list, max_iter=1000, tol=1e-6):
    """Quantify gene abundances using the EM algorithm."""
    print("Quantifying gene abundances (EM algorithm)...")
    N = len(gene_list)
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}
    lengths = np.array([gene_lengths[g] for g in gene_list])
    
    alphas = np.ones(N) / N
    prev_alphas = np.zeros(N)
    
    for iteration in range(max_iter):
        y = np.zeros(N)
        
        for S, c in eq_classes.items():
            if not S:
                continue
            S_idx = np.array([gene_to_idx[g] for g in S])
            rho_S = alphas[S_idx] * lengths[S_idx]
            sum_rho = rho_S.sum()
            
            if sum_rho == 0:
                continue
            
            y[S_idx] += c * (rho_S / sum_rho)
        
        if y.sum() == 0:
            break
        
        alphas = y / y.sum()
        
        if np.max(np.abs(alphas - prev_alphas)) < tol:
            print(f"  Converged after {iteration + 1} iterations")
            break
        
        prev_alphas = alphas.copy()
    else:
        print(f"  Reached maximum iterations ({max_iter})")
    
    est_counts = y
    tpm = (est_counts / lengths) / ((est_counts / lengths).sum() + 1e-10) * 1e6
    total_counts = est_counts.sum()
    rpkm = (est_counts * 1e9) / (lengths * total_counts + 1e-10)
    
    return est_counts, tpm, rpkm


def create_visualizations(df, output_prefix):
    """Create multiple visualizations of the results."""
    print("Creating visualizations...")
    
    df_filtered = df[df['est_count'] > 0].copy()
    if len(df_filtered) == 0:
        print("  Warning: No genes detected, skipping visualizations")
        return
    
    df_filtered = df_filtered.sort_values('tpm', ascending=False)
    
    fig = plt.figure(figsize=(20, 12))
    top_n = min(20, len(df_filtered))
    top_genes = df_filtered.head(top_n)
    colors = plt.cm.viridis(np.linspace(0, 1, top_n))
    
    ax1 = plt.subplot(2, 3, 1)
    ax1.barh(range(top_n), top_genes['tpm'], color=colors)
    ax1.set_yticks(range(top_n))
    ax1.set_yticklabels(top_genes['gene_id'], fontsize=8)
    ax1.set_xlabel('TPM', fontsize=10)
    ax1.set_title(f'Top {top_n} AMR Genes by TPM', fontsize=12, fontweight='bold')
    ax1.invert_yaxis()
    
    ax2 = plt.subplot(2, 3, 2)
    ax2.barh(range(top_n), top_genes['est_count'], color=colors)
    ax2.set_yticks(range(top_n))
    ax2.set_yticklabels(top_genes['gene_id'], fontsize=8)
    ax2.set_xlabel('Estimated Count', fontsize=10)
    ax2.set_title(f'Top {top_n} AMR Genes by Count', fontsize=12, fontweight='bold')
    ax2.invert_yaxis()
    
    ax3 = plt.subplot(2, 3, 3)
    ax3.barh(range(top_n), top_genes['rpkm'], color=colors)
    ax3.set_yticks(range(top_n))
    ax3.set_yticklabels(top_genes['gene_id'], fontsize=8)
    ax3.set_xlabel('RPKM', fontsize=10)
    ax3.set_title(f'Top {top_n} AMR Genes by RPKM', fontsize=12, fontweight='bold')
    ax3.invert_yaxis()
    
    ax4 = plt.subplot(2, 3, 4)
    heatmap_data = top_genes[['est_count', 'tpm', 'rpkm']].T
    heatmap_data_norm = (heatmap_data - heatmap_data.min(axis=1).values.reshape(-1, 1)) / \
                        (heatmap_data.max(axis=1).values.reshape(-1, 1) -
                         heatmap_data.min(axis=1).values.reshape(-1, 1) + 1e-10)
    sns.heatmap(heatmap_data_norm,
                xticklabels=top_genes['gene_id'],
                yticklabels=['Est. Count', 'TPM', 'RPKM'],
                cmap='YlOrRd',
                cbar_kws={'label': 'Normalized Value'},
                ax=ax4,
                linewidths=0.5)
    ax4.set_title('Heatmap: Top Genes (Normalized)', fontsize=12, fontweight='bold')
    plt.setp(ax4.get_xticklabels(), rotation=45, ha='right', fontsize=7)
    
    ax5 = plt.subplot(2, 3, 5)
    ax5.hist(np.log10(df_filtered['tpm'] + 1), bins=30, color='skyblue', edgecolor='black')
    ax5.set_xlabel('log10(TPM + 1)', fontsize=10)
    ax5.set_ylabel('Frequency', fontsize=10)
    ax5.set_title('Distribution of TPM Values', fontsize=12, fontweight='bold')
    ax5.grid(alpha=0.3)
    
    ax6 = plt.subplot(2, 3, 6)
    cumsum = df_filtered['tpm'].cumsum() / df_filtered['tpm'].sum() * 100
    ax6.plot(range(len(cumsum)), cumsum, linewidth=2, color='darkblue')
    ax6.set_xlabel('Number of Genes', fontsize=10)
    ax6.set_ylabel('Cumulative % of Total TPM', fontsize=10)
    ax6.set_title('Cumulative Gene Expression', fontsize=12, fontweight='bold')
    ax6.grid(alpha=0.3)
    ax6.axhline(y=90, color='r', linestyle='--', alpha=0.5, label='90%')
    ax6.legend()
    
    plt.tight_layout()
    viz_file = f"{output_prefix}_visualization.png"
    plt.savefig(viz_file, dpi=300, bbox_inches='tight')
    print(f"  Saved visualization to {viz_file}")
    plt.close()
    
    if len(df_filtered) > 10:
        create_detailed_heatmap(df_filtered, output_prefix)


def create_detailed_heatmap(df, output_prefix):
    """Create a detailed heatmap for all detected genes."""
    fig, ax = plt.subplots(figsize=(12, max(8, len(df) * 0.3)))
    
    heatmap_data = df[['est_count', 'tpm', 'rpkm']].T
    heatmap_data_norm = (heatmap_data - heatmap_data.min(axis=1).values.reshape(-1, 1)) / \
                        (heatmap_data.max(axis=1).values.reshape(-1, 1) -
                         heatmap_data.min(axis=1).values.reshape(-1, 1) + 1e-10)
    
    sns.heatmap(heatmap_data_norm,
                xticklabels=df['gene_id'],
                yticklabels=['Est. Count', 'TPM', 'RPKM'],
                cmap='RdYlGn',
                cbar_kws={'label': 'Normalized Value'},
                ax=ax,
                linewidths=0.5,
                annot=False)
    
    ax.set_title('Detailed Heatmap: All Detected AMR Genes (Normalized)',
                 fontsize=14, fontweight='bold', pad=20)
    plt.setp(ax.get_xticklabels(), rotation=90, ha='right', fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=10)
    
    plt.tight_layout()
    heatmap_file = f"{output_prefix}_detailed_heatmap.png"
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    print(f"  Saved detailed heatmap to {heatmap_file}")
    plt.close()


def write_summary_stats(df, total_reads, aligned_reads, output_prefix):
    """Write summary statistics to a file (full FASTA headers included)."""
    summary_file = f"{output_prefix}_summary.txt"
    
    detected_genes = len(df[df['est_count'] > 0])
    total_genes = len(df)
    
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("AMR ANALYSIS SUMMARY\n")
        f.write("=" * 60 + "\n")
        
        f.write("READ STATISTICS:\n")
        f.write(f"  Total reads: {total_reads:,}\n")
        alignment_pct = 100 * aligned_reads / total_reads if total_reads else 0.0
        f.write(f"  Aligned reads: {aligned_reads:,} ({alignment_pct:.2f}%)\n")
        
        f.write("GENE DETECTION:\n")
        f.write(f"  Total AMR genes in database: {total_genes}\n")
        f.write(f"  Detected genes (count > 0): {detected_genes}\n")
        detection_pct = 100 * detected_genes / total_genes if total_genes else 0.0
        f.write(f"  Detection rate: {detection_pct:.2f}%\n")
        
        if detected_genes > 0:
            df_detected = df[df['est_count'] > 0].sort_values('tpm', ascending=False)
            
            f.write("TOP 10 DETECTED GENES:\n")
            for _, row in df_detected.head(10).iterrows():
                header = row['full_header'][:11]
                tpm_value = row['tpm']
                count_value = row['est_count']
                f.write(f"  {header}: TPM={tpm_value:.2f}, Count={count_value:.2f}\n")
            
            f.write("\nSTATISTICS FOR DETECTED GENES:\n")
            f.write(f"  Mean TPM: {df_detected['tpm'].mean():.2f}\n")
            f.write(f"  Median TPM: {df_detected['tpm'].median():.2f}\n")
            f.write(f"  Total estimated counts: {df_detected['est_count'].sum():.2f}\n")
        else:
            f.write("No genes detected in this sample.\n")
    
    print(f"  Saved summary statistics to {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description='AMR Gene Quantification and Visualization Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python amr_analyzer.py amr.fasta reads.fastq output
  python amr_analyzer.py amr.fasta reads.fastq output -k 25
  python amr_analyzer.py amr.fasta reads.fastq output --no-viz
        """
    )
    
    parser.add_argument('amr_fasta', help='AMR gene database (FASTA format)')
    parser.add_argument('reads_fastq', help='Sequencing reads (FASTQ format)')
    parser.add_argument('output_prefix', help='Output file prefix')
    parser.add_argument('-k', '--kmer-size', type=int, default=8,
                        help='K-mer size for pseudoalignment (default: 8)')
    parser.add_argument('--max-iter', type=int, default=200,
                        help='Maximum EM iterations (default: 200)')
    parser.add_argument('--no-viz', action='store_true',
                        help='Skip visualization generation')
    
    args = parser.parse_args()
    
    print("\n" + "=" * 60)
    print("AMR GENE QUANTIFICATION PIPELINE")
    print("=" * 60 + "\n")
    
    index, gene_lengths, header_map = build_index(args.amr_fasta, args.kmer_size)
    gene_list = list(gene_lengths.keys())
    
    eq_classes, total_reads, aligned_reads = pseudoalign_reads(
        args.reads_fastq, index, args.kmer_size
    )
    
    est_counts, tpm, rpkm = quantify(
        eq_classes, gene_lengths, gene_list, max_iter=args.max_iter
    )
    
    df = pd.DataFrame({
        'gene_id'    : gene_list,
        'full_header': [header_map[g] for g in gene_list],
        'gene_length': [gene_lengths[g] for g in gene_list],
        'est_count'  : est_counts,
        'tpm'        : tpm,
        'rpkm'       : rpkm
    }).sort_values('tpm', ascending=False)
    
    csv_file = f"{args.output_prefix}.csv"
    df.to_csv(csv_file, index=False, float_format='%.4f')
    print(f"\nWriting output to {csv_file}...")
    print(f"  Wrote {len(df)} genes")
    
    if not args.no_viz:
        try:
            create_visualizations(df, args.output_prefix)
        except Exception as e:
            print(f"  Warning: Could not create visualizations: {e}")
    
    write_summary_stats(df, total_reads, aligned_reads, args.output_prefix)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
