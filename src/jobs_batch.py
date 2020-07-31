#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
jobs_batch.py - For each of the given fa files in args, generate
alignments against hg38 and compute a score threshold per GC tuned
matrix. For one instance of this program, all the consensus sequences
will be run sequentially.

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
    parser.add_argument("dirname",
            help="dir with consensus sequence fa files whose score " +
                "thresholds will be computed.")
    args = parser.parse_args()

    n = 3209286105 # size of hg38

    for fpath in [os.path.join(args.dirname, f)
                for f in os.listdir(args.dirname)]:
        name = fpath.split("/")[-1][:-3]
        m = consensusSize(fpath)

        print("Generating genomic alignments for " + name)
        generateAlignments(fpath,
            "../data/hg38bins/dfamseq_bins", "../results/genomic_hits/")
        print("Generating benchmark alignments for " + name)
        generateAlignments(fpath,
            "../data/hg38bins/benchmark_bins", "../results/benchmark_hits/")

        print("Calculating score thresholds for " + name)
        scoreThresholds(os.path.join("../results/genomic_hits/", name),
                        os.path.join("../results/benchmark_hits/", name),
                        query_size=m, subject_size=n)

if __name__ == '__main__':
    main()
