#!/usr/bin/env python3

import sys
import os

def filter_tissues(matrix_file, tissues_file, ignore_tissues, output_dir):
    """
    Filter matrix and tissues files by removing specified tissues and any rows
    that have no remaining tissues after filtering.
    
    Args:
        matrix_file (str): Path to the input matrix TSV file
        tissues_file (str): Path to the input tissues TSV file
        ignore_tissues (list): List of tissue names to remove
        output_dir (str): Directory to write output files
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read tissues file and filter out ignored tissues
    tissues_data = {}
    with open(tissues_file) as f:
        header = next(f)  # Skip header
        for line in f:
            asv, tissues = line.strip().split('\t')
            tissue_list = tissues.split(',')
            # Remove ignored tissues
            remaining_tissues = [t for t in tissue_list if t not in ignore_tissues]
            if remaining_tissues:  # Only keep ASVs with remaining tissues
                tissues_data[asv] = remaining_tissues
    
    # Write filtered tissues file
    with open(os.path.join(output_dir, 'tissues.tsv'), 'w') as f:
        f.write('group_name\ttissues\n')
        for asv in tissues_data:
            line = asv + '\t' + ','.join(tissues_data[asv]) + '\n'
            f.write(line)
    
    # Read and filter matrix file
    with open(matrix_file) as f:
        matrix_header = next(f)
        matrix_data = {}
        for line in f:
            parts = line.strip().split('\t')
            asv = parts[0]
            if asv in tissues_data:  # Only keep ASVs that have remaining tissues
                matrix_data[asv] = parts[1:]
    
    # Write filtered matrix file
    with open(os.path.join(output_dir, 'matrix.tsv'), 'w') as f:
        f.write(matrix_header)
        for asv in matrix_data:
            line = asv + '\t' + '\t'.join(matrix_data[asv]) + '\n'
            f.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python filter_for_tissue_scenarios.py <matrix_file> <tissues_file> <ignore_tissues> <output_dir>")
        sys.exit(1)
        
    matrix_file = sys.argv[1]
    tissues_file = sys.argv[2]
    ignore_tissues = sys.argv[3].split(',')
    output_dir = sys.argv[4]

    filter_tissues(matrix_file, tissues_file, ignore_tissues, output_dir) 