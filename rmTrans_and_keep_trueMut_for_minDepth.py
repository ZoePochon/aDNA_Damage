#!/usr/bin/env python3

import argparse
import pysam
from collections import defaultdict
from tqdm import tqdm

"""
Authors:
    Zoé Pochon: zoe.pochon@gmail.com

Description:
    This script processes a BAM file to remove C->T and G->A transitions, which are common damage patterns in ancient DNA,
    while retaining transitions at positions likely representing true mutations based on specific thresholds.

    The script uses two adjustable thresholds:
    - Transition frequency threshold (default 0.9): Positions with at least 90% transition frequency are considered true mutations.
    - Minimum depth threshold (default 3): Only positions with a minimum read depth of 3 are considered for mutation analysis, the other ones have transition removed.

    Usage:
        python rmTrans_and_keep_trueMut_for_minDepth.py -b input.bam -o output.bam -r reference.fasta -t 0.9 -d 2
"""

# Parse arguments
parser = argparse.ArgumentParser(description="Process BAM file to remove C->T and G->A transitions unless they appear as true mutations.")
parser.add_argument("-b", "--bamfile", required=True, help="Path to the input BAM file.")
parser.add_argument("-o", "--output", required=True, help="Path to the output BAM file.")
parser.add_argument("-r", "--reference", required=True, help="Path to the reference FASTA file.")
parser.add_argument("-t", "--threshold", type=float, default=0.9, help="Threshold for transition percentage to consider as a true mutation.")
parser.add_argument("-d", "--min_depth", type=int, default=3, help="Minimum read depth required to consider a position a true mutation.")

args = parser.parse_args()

# Load the reference FASTA file
reference = pysam.FastaFile(args.reference)

# Initialize a dictionary to store counts for each position
transition_counts = defaultdict(lambda: {"total": 0, "transitions": 0})

# Function to process each read and collect transition statistics
def count_transitions(bamfile):
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in tqdm(bam.fetch(), desc="Counting transitions"):
            if read.mapping_quality == 0 or read.is_unmapped:
                continue

            read_seq = read.query_sequence
            ref_pos = read.reference_start
            ref_name = read.reference_name

            for base_index in range(len(read_seq)):
                ref_base = reference.fetch(ref_name, ref_pos + base_index, ref_pos + base_index + 1).upper()
                read_base = read_seq[base_index]
                current_position = f"{ref_name}:{ref_pos + base_index}"

                transition_counts[current_position]["total"] += 1

                if (ref_base == "C" and read_base == "T") or (ref_base == "G" and read_base == "A"):
                    transition_counts[current_position]["transitions"] += 1

# Identify positions where transitions meet the threshold and minimum depth
def identify_true_mutations(threshold, min_depth):
    true_mutations = set()
    for position, counts in transition_counts.items():
        if counts["total"] >= min_depth and counts["transitions"] / counts["total"] >= threshold:
            true_mutations.add(position)
    return true_mutations

# Open input BAM and output BAM with header preserved
with pysam.AlignmentFile(args.bamfile, "rb") as bam, pysam.AlignmentFile(args.output, "wb", header=bam.header) as out_bam:
    count_transitions(args.bamfile)
    true_mutations = identify_true_mutations(args.threshold, args.min_depth)
    print(f"Identified {len(true_mutations)} true mutation positions with depth >= {args.min_depth}.")

    for read in tqdm(bam.fetch(), desc="Processing reads"):
        if read.mapping_quality == 0 or read.is_unmapped:
            continue

        read_seq = read.query_sequence
        read_qualities = read.query_qualities
        ref_pos = read.reference_start
        ref_name = read.reference_name
        updated_seq = list(read_seq)

        for base_index in range(len(read_seq)):
            current_position = f"{ref_name}:{ref_pos + base_index}"

            if current_position in true_mutations:
                continue

            ref_base = reference.fetch(ref_name, ref_pos + base_index, ref_pos + base_index + 1).upper()
            read_base = read_seq[base_index]

            if ref_base == "C" and read_base == "T":
                updated_seq[base_index] = "N"
            elif ref_base == "G" and read_base == "A":
                updated_seq[base_index] = "N"

        read.query_sequence = "".join(updated_seq)
        read.query_qualities = read_qualities
        out_bam.write(read)
