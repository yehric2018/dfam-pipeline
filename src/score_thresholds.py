#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
score_thresholds.py: Calculate score thresholds for a consensus
sequence such that the false discovery rate is below 0.2%. Uses
empirical and theoretical FDR calculations and chooses the more
conservative score threshold value for each GC background.

Usage:
Args:

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import re

def generateScoreThreshold(genome_file, benchmark_file):
    """
    generateScoreThreshold(genome_file, benchmark_file) -
    Take in the file names of two alignment files produced from
    RMBlast, one produced from genome bins and one produced from
    benchmark bins. Returns a conservative score threshold for
    this consensus sequence.

    Args:
        genome_file - alignment file against genome bins.
        benchmark_file - alignment file against benchmark bins.

    Returns: score threshold for this consensus sequence computed
        from the given alignments.
    """
    genomic_hits = []
    benchmark_hits = []
    gf = open(genome_file, "r")
    bf = open(benchmark_file, "r")
    regex = re.compile(r"^\s*(\d+)\s+\d+\.\d+\s+\d+\.\d+")
    count = 0
    line = gf.readline()
    while line != "":
        mo = regex.search(line)
        if mo:
            count += 1
            genomic_hits.append(int(mo.group(1)))
        line = gf.readline()
    print(genomic_hits)

if __name__ == '__main__':
    generateScoreThreshold("../results/test_alignments/dfamseq/DF0000001_25p35g.sc",
                            "../results/test_alignments/benchmark/DF0000001_25p35g.sc")

