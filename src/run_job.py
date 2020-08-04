#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
run_job.py - For the given consensus sequence fa file, generate
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
    parser.add_argument("consensus",
            help="consensus sequence fa file whose score " +
                "thresholds will be computed.")
    args = parser.parse_args()

    fpath = args.consensus

    n = 3209286105 # size of hg38

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
