#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
score_thresholds.py: Calculate score thresholds for a consensus
sequence such that the false discovery rate is below 0.2%. Uses
empirical and theoretical FDR calculations and chooses the more
conservative score threshold value for each GC background.

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import re

def readScoresFromFile(sc_file):
    """
    readScoresFromFile(sc_file) - Reads the given alignment file and
    produces a list containing for each alignment the score and the
    lengths of the overlapping sequences, sorted in decreasing order
    by score.

    Args:
        sc_file - path to alignment file produced from RMBlast.

    Returns: list of Hit objects produced from sc_file sorted
        in decreasing order by score.
    """
    f = open(sc_file, "r")
    scoreRegex = re.compile(r"^\s*(\d+)\s+\d+\.\d+\s+\d+\.\d+")
    rangeRegex1 = re.compile(r"\s*(\d+)\s+(\d+)\s+\(\d+\)\s+\S+\s+(\d+)\s+(\d+)\s+\(\d+\)$")
    rangeRegex2 = re.compile(r"\s*(\d+)\s+(\d+)\s+\(\d+\)\s+C\s+\S+\s+\(\d+\)\s+(\d+)\s+(\d+)$")
    hits = []
    line = f.readline()
    while line != "":
        mo = scoreRegex.search(line)
        if mo:
            hits.append(int(mo.group(1)))
        mo = rangeRegex1.search(line)
        if mo:
            print(mo)
            print(int(mo.group(2)) - int(mo.group(1)))
        mo = rangeRegex2.search(line)
        if mo:
            print(mo)
            print(int(mo.group(3)) - int(mo.group(4)))
        line = f.readline()
    hits.sort(reverse=True)
    print(hits)
    return hits

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
    genomic_hits = readScoresFromFile(genome_file)
    benchmark_hits = readScoresFromFile(benchmark_file)

if __name__ == '__main__':
    generateScoreThreshold("../results/test_alignments/dfamseq/DF0000001_25p35g.sc",
                            "../results/test_alignments/benchmark/DF0000001_25p35g.sc")

