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
import sqlite3
import re
import math

GUMBEL = {}

class Hit:
    """
    Stores information on hits needed to compute E-values.

    Fields:
        score - score of alignment
        m - length of query sequence
        n - length of subject sequence
    """
    def __init__(self, score, m, n):
        """
        Initialize this hit object with the given params.
        """
        self.score = score
        self.m = m
        self.n = n

    def __repr__(self):
        return str(self.score)

def computeEValue(hit, matrix):
    """
    computeEValues(hit) - Returns the E-value of the given hit. The
    E-value is the number of expected hits of similar score that
    could be found just by chance.

    The E-value is computed via the following formula:
        E = K * m * n * e^(-lambda * S)

    lambda and K are parameters of a fitted Gumbel distribution, S is
    the raw score of the alignment and "m" and "n" are the lengths of
    the aligned query and subject sequences.
    """
    return (GUMBEL[matrix]["k"] * hit.m * hit.n *
                math.exp(-GUMBEL[matrix]["lambda"] * hit.score))

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
    matrix = sc_file[-9:-3]
    print(matrix)
    f = open(sc_file, "r")
    scoreRegex = re.compile(r"^\s*(\d+)\s+\d+\.\d+\s+\d+\.\d+")
    rangeRegex1 = re.compile(r"\s*(\d+)\s+(\d+)\s+\(\d+\)\s+\S+\s+(\d+)\s+(\d+)\s+\(\d+\)$")
    rangeRegex2 = re.compile(r"\s*(\d+)\s+(\d+)\s+\(\d+\)\s+C\s+\S+\s+\(\d+\)\s+(\d+)\s+(\d+)$")
    hits = []
    line = f.readline()
    while line != "":
        mo = scoreRegex.search(line)
        if mo:
            score = int(mo.group(1))
            mo1 = rangeRegex1.search(line)
            mo2 = rangeRegex2.search(line)
            if mo1:
                n = int(mo1.group(2)) - int(mo1.group(1))
                m = int(mo1.group(4)) - int(mo1.group(3))
                hits.append(Hit(score, m, n))
            elif mo2:
                n = int(mo2.group(2)) - int(mo2.group(1))
                m = int(mo2.group(3)) - int(mo2.group(4))
                hits.append(Hit(score, m, n))
        line = f.readline()
    hits.sort(key=lambda hit: hit.score, reverse=True)
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
    if not GUMBEL:
        conn = sqlite3.connect("../data/gumbel.db")
        c = conn.cursor()
        for row in c.execute("SELECT * FROM params"):
            matrix = row[0]
            GUMBEL[matrix] = {}
            GUMBEL[matrix]["lambda"] = row[3]
            GUMBEL[matrix]["k"] = row[4]
        conn.close()
    genomic_hits = readScoresFromFile(genome_file)
    genomic_e = [computeEValue(h, matrix) for h in genomic_hits]
    print(genomic_e)
    #benchmark_hits = readScoresFromFile(benchmark_file)

if __name__ == '__main__':
    generateScoreThreshold("../results/test_alignments/dfamseq/DF0000001_25p35g.sc",
                            "../results/test_alignments/benchmark/DF0000001_25p35g.sc")

