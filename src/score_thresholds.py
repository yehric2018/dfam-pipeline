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
m, n = 0
FDR_THRESHOLD = 0.002
FDR_THEORY_TARGET = 0.01
TEMP_GENOME_SIZE = 3099000000
TEMP_CONSENSUS_SIZE = 262

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
    return (GUMBEL[matrix]["k"] * TEMP_CONSENSUS_SIZE * TEMP_GENOME_SIZE *
                math.exp(-GUMBEL[matrix]["lambda"] * hit))

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
    hits = []
    line = f.readline()
    while line != "":
        mo = scoreRegex.search(line)
        if mo:
            score = int(mo.group(1))
            hits.append(score)
        line = f.readline()
    hits.sort(reverse=True)
    return hits

def empiricalFDRCalculation(genomic_hits, benchmark_hits):
    """
    empiricalFDRCalculation(genomic_hits, benchmark_hits) - Uses
    empirical FDR calculation to compute a score threshold for this
    consensus sequence for a certain GC background. There must be
    sufficient alignments for a family to a benchmark (FP) sequence.

    Searches genomic/benchmark hits sorted (high to low) by raw score,
    stops when FP_hits / cumulative_genomic > FDR_target. The
    threshold is chosen to be 0.05 + score of hit that exceeded the
    threshold.

    Args:
        genomic_hits: List of Hit objects obtained from alignment
            of consensus against genomic sequence.
        benchmark_hits: List of Hit objects obtained from alignment
            of consensus against bnechmark sequence.

    Returns: score threshold that keeps the false discovery rate
        below 0.2%.
    """
    fp_hits = 0
    hits = 0
    i = 0
    j = 0
    if genomic_hits[i]< benchmark_hits[j]:
        return genomic_hits[i]+ 0.05
    i += 1
    hits += 1
    while fp_hits / hits < FDR_THRESHOLD:
        if genomic_hits[i]< benchmark_hits[j]:
            fp_hits += 1
            j += 1
        else:
            hits += 1
            i += 1
    print(genomic_hits[i]+ 0.05)
    print(i)
    return genomic_hits[i]+ 0.05

def theoreticalFDRCalculation(genomic_hits, benchmark_hits, matrix):
    """
    theoreticalFDRCalculation(genomic_hits, benchmark_hits) - Uses
    theoretical FDR calculation to compute a score threshold for this
    consensus sequence for a certain GC background.

    The score threshold is required to be at least as high as the
    score corresponding to an E-value of min(1000, (genomic_hits -
    benchmark_hits) * FDR_threshold). This allows us to avoid
    unreasonably liberal scores for families with a large number of
    genomic hits and few, if any, benchmark hits.

    Args:
        genomic_hits: List of Hit objects obtained from alignment
            of consensus against genomic sequence.
        benchmark_hits: List of Hit objects obtained from alignment
            of consensus against bnechmark sequence.

    Returns: theoretical score threshold that should keep the false
        discovery rate below 0.2%.
    """
    print("Theoretical calculation:")
    print(len(genomic_hits))
    print(len(benchmark_hits))
    tp_estimate = len(genomic_hits) - len(benchmark_hits)
    if tp_estimate <= 0 or tp_estimate > 1000:
        tp_estimate = 1000
    target = tp_estimate * FDR_THEORY_TARGET

    # Convert e-value target into score threshold
    in_log = GUMBEL[matrix]['k'] * TEMP_CONSENSUS_SIZE * TEMP_GENOME_SIZE
    raw = (math.log(in_log) - math.log(target)) / GUMBEL[matrix]['lambda']
    print(raw)
    return raw

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
    conn = sqlite3.connect("../data/consensus.db")
    c = conn.cursor()
    #m = c.execute("SELECT size FROM consensus WHERE sequence = ")[0][0]
    if not GUMBEL:
        for row in c.execute("SELECT * FROM gumbel"):
            matrix = row[0]
            GUMBEL[matrix] = {}
            GUMBEL[matrix]["lambda"] = row[3]
            GUMBEL[matrix]["k"] = row[4]
    conn.close()
    genomic_hits = readScoresFromFile(genome_file)
    genomic_e = [computeEValue(h, matrix) for h in genomic_hits]
    benchmark_hits = readScoresFromFile(benchmark_file)
    print("first benchmark: " + str(benchmark_hits[0]))
    benchmark_e = [computeEValue(h, matrix) for h in benchmark_hits]

    print(max(empiricalFDRCalculation(genomic_hits, benchmark_hits), theoreticalFDRCalculation(genomic_hits, benchmark_hits, matrix)))

if __name__ == '__main__':
    generateScoreThreshold("../results/test_alignments/dfamseq/DF0000001_25p37g.sc",
                            "../results/test_alignments/benchmark/DF0000001_25p37g.sc")

