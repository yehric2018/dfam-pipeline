#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
dfam_pipeline.py - Given a set of consensus sequences, compute 10
score thresholds for each sequence for the given genome, one for each
GC tuned matrix. Also generates alignments for the consensus sequence
against the given genome.

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import argparse
import os

from sequence_util import consensusSize, genomeSize
from generate_alignments import generateAlignments, splitConsensus
from score_thresholds import scoreThresholds

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("consensus",
            help="either a fa file containing one or more consensus" +
                " sequences or a directory containing one or more " +
                "fa files, each with exactly one consensus sequence.")
    """
    parser.add_argument("genome",
            help="a fa file containing a single genome, which the" +
                " given consensus sequences are associated with.")
    parser.add_argument("benchmark",
            help="a fa file containing a benchmark for the given " +
                "genome.")
    """
    args = parser.parse_args()

    consensus_dir = args.consensus
    # Check if the consensus sequence is a file or directory
    if os.path.isfile(args.consensus):
        # if a file, must split into new directory with several consensus sequences
        print("consensus sequence file")
        consensus_dir = splitConsensus(args.consensus)
    elif not os.path.isdir(args.consensus):
        print("no consensus sequence found")
        # if a directory, can use it as it is

    n = 3209286105

    filelist = os.listdir(consensus_dir)
    for f in filelist:
        name = f[:-3]
        m = consensusSize(os.path.join(consensus_dir, f))

        print("Generating genomic alignments for " + name)
        generateAlignments(os.path.join(consensus_dir, f),
            "../data/hg38bins/dfamseq_bins", "../results/genomic_hits/")
        print("Generating benchmark alignments for " + name)
        generateAlignments(os.path.join(consensus_dir, f),
            "../data/hg38bins/benchmark_bins", "../results/benchmark_hits/")

        print("Calculating score thresholds for " + name)
        scoreThresholds(os.path.join("../results/genomic_hits/", name),
                        os.path.join("../results/benchmark_hits/", name),
                        query_size=m, subject_size=n)
    # For now, we will assume all consensus sequences passed in are from the human genomes.
    # If we are able to do several genomes as well, we must add the following steps:
        # Check if the genome already has binned batches produced for it, same for benchmark
        # if not yet produced, run bin_genome on both the genome and the benchmark
    # Generate alignments for each consensus sequence against the given genome and benchmark
    # Finally, for each set of benchmark/genomic alignments, calculate the score threshold and output

if __name__ == '__main__':
    main()
