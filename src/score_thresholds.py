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
import argparse
import re
import math
import os

from sequence_util import consensusSize

GUMBEL = {'25p53g': {'lambda': 0.109152, 'k': 0.111427},
        '14p51g': {'lambda': 0.126273, 'k': 0.271705},
        '25p41g': {'lambda': 0.113153, 'k': 0.142018},
        '14p39g': {'lambda': 0.120541, 'k': 0.267775},
        '14p35g': {'lambda': 0.126797, 'k': 0.284773},
        '18p43g': {'lambda': 0.121582, 'k': 0.234784},
        '25p51g': {'lambda': 0.114421, 'k': 0.13735},
        '18p35g': {'lambda': 0.124283, 'k': 0.246276},
        '25p35g': {'lambda': 0.108087, 'k': 0.126783},
        '25p45g': {'lambda': 0.108282, 'k': 0.126305},
        '14p37g': {'lambda': 0.126554, 'k': 0.284202},
        '14p43g': {'lambda': 0.123083, 'k': 0.279322},
        '18p41g': {'lambda': 0.114068, 'k': 0.215272},
        '25p39g': {'lambda': 0.110804, 'k': 0.133886},
        '14p41g': {'lambda': 0.123274, 'k': 0.278371},
        '20p43g': {'lambda': 0.116161, 'k': 0.207456},
        '18p45g': {'lambda': 0.121927, 'k': 0.241789},
        '20p51g': {'lambda': 0.11795, 'k': 0.183554},
        '14p47g': {'lambda': 0.119249, 'k': 0.258223},
        '14p53g': {'lambda': 0.123762, 'k': 0.258946},
        '20p49g': {'lambda': 0.121086, 'k': 0.190207},
        '18p49g': {'lambda': 0.127971, 'k': 0.257942},
        '18p51g': {'lambda': 0.115054, 'k': 0.191605},
        '18p39g': {'lambda': 0.126286, 'k': 0.261353},
        '25p43g': {'lambda': 0.123854, 'k': 0.176985},
        '18p53g': {'lambda': 0.113257, 'k': 0.190965},
        '14p49g': {'lambda': 0.119393, 'k': 0.255202},
        '25p37g': {'lambda': 0.109651, 'k': 0.129716},
        '20p35g': {'lambda': 0.117394, 'k': 0.195042},
        '20p37g': {'lambda': 0.118879, 'k': 0.203742},
        '25p47g': {'lambda': 0.108085, 'k': 0.127639},
        '25p49g': {'lambda': 0.117029, 'k': 0.151864},
        '20p41g': {'lambda': 0.114931, 'k': 0.197257},
        '20p47g': {'lambda': 0.124492, 'k': 0.226433},
        '20p53g': {'lambda': 0.12917, 'k': 0.213967},
        '20p39g': {'lambda': 0.127701, 'k': 0.234063},
        '18p37g': {'lambda': 0.124933, 'k': 0.246589},
        '18p47g': {'lambda': 0.11989, 'k': 0.233897},
        '14p45g': {'lambda': 0.130133, 'k': 0.302126},
        '20p45g': {'lambda': 0.116422, 'k': 0.189823}}
m = n = 0
FDR_THRESHOLD = 0.002
FDR_THEORY_TARGET = 0.01
MAX_E_TARGET = 1000
THRESHOLDS_TABLE = open("../results/thresholds.txt", "a")
TEMP_GENOME_SIZE = 3209286105
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
    tp_estimate = len(genomic_hits) - len(benchmark_hits)
    if tp_estimate < 0:
        tp_estimate = 0
    target = tp_estimate * FDR_THEORY_TARGET
    if target == 0 or target > MAX_E_TARGET:
        target = MAX_E_TARGET

    # Convert e-value target into score threshold
    in_log = GUMBEL[matrix]['k'] * TEMP_CONSENSUS_SIZE * TEMP_GENOME_SIZE
    raw = (math.log(in_log) - math.log(target)) / GUMBEL[matrix]['lambda']
    return raw

def generateScoreThreshold(genome_file, benchmark_file, thresholds_table):
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
    consensus = genome_file.split("/")[-1].split("_")[0]
    matrix = genome_file[-9:-3]
    #print("computing score threshold for " + consensus + " with matrix " + matrix)
    genomic_hits = readScoresFromFile(genome_file)
    genomic_e = [computeEValue(h, matrix) for h in genomic_hits]
    benchmark_hits = readScoresFromFile(benchmark_file)
    benchmark_e = [computeEValue(h, matrix) for h in benchmark_hits]

    empirical = empiricalFDRCalculation(genomic_hits, benchmark_hits)
    theoretical = theoreticalFDRCalculation(genomic_hits, benchmark_hits, matrix)
    #print("empirical score: " + str(empirical))
    #print("theoretical score: " + str(theoretical))
    #print("final score threshold: " + str(max(empirical, theoretical)))
    line = (consensus + "\t" + matrix + "\t" + str(empirical) + "\t" +
            str(theoretical) + "\t" + str((max(empirical, theoretical))))
    print(line)
    thresholds_table.write(line + "\n")

def scoreThresholds(genome_dir, benchmark_dir,
        query_size=TEMP_CONSENSUS_SIZE, subject_size=TEMP_GENOME_SIZE):
    """
    scoreThresholds(genome_dir, benchmark_dir) -
    Creates a table for all the alignment files found in the given
    genomic and benchmark directories for a single consensus
    sequence, entries tab-separated, in the following format:

    sequence_name\tmatrix\tempirical\ttheoretical\tfinal

    where sequence_name is the name of the consensus sequence, matrix
    is the matrix used for this score calculation, empirical is the
    threshold generated using the empirical method, theoretical is
    the threshold generated using the theoretical method, and final
    is the more conservative of the empirical and theoretical score
    thresholds found.

    Args:
        genome_dir: Path to the directory containing genomic
            alignments for a consensus sequence.
        benchmark_dir: Path to the directory containing benchmark
            alignments for the same consensus sequence used in
            genome_dir, and with the same file names.
        query_size: number of bps in consensus sequence.
        subject_size: number of bps in subject sequence/genome.
    """
    #thresholds_table = open("../results/thresholds.txt", "a")
    consensus = genome_dir[:-1].split("/")[-1].split("_")[0]
    thresholds_table = open("../results/thresholds/" + consensus +
                    ".thresh", "w")
    print("../results/thresholds/" + consensus + ".thresh")
    genome_list = os.listdir(genome_dir)
    benchmark_list = os.listdir(benchmark_dir)
    m = query_size
    n = subject_size
    for f in genome_list:
        generateScoreThreshold(os.path.join(genome_dir, f),
            os.path.join(benchmark_dir, f), thresholds_table)
    thresholds_table.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("genomic_hits",
            help="path to the genomic alignments for consensus " +
                "sequence of interest.")
    parser.add_argument("benchmark_hits",
            help="path to the benchmark alignments for consensus " +
                "sequence of interest.")
    parser.add_argument("--m",type=int, default=TEMP_CONSENSUS_SIZE,
            help="size of given consensus sequence (query)")
    parser.add_argument("--n",type=int, default=TEMP_GENOME_SIZE,
            help="size of genome searched against (subject)")
    args = parser.parse_args()

    scoreThresholds(args.genomic_hits, args.benchmark_hits, args.m, args.n)

    """
    scoreThresholds("../results/alignments/test_alignments/dfamseq/DF0000001/",
            "../results/alignments/test_alignments/benchmark/DF0000001/")
    scoreThresholds("../results/alignments/test_alignments/dfamseq/DF0000002/",
            "../results/alignments/test_alignments//benchmark/DF0000002/",
            query_size=consensusSize("../data/consensus/ex_hg38_cons.fa_/DF0000002.fa"))
    scoreThresholds("../results/alignments/test_alignments/dfamseq/DF0000004/",
            "../results/alignments/test_alignments/benchmark/DF0000004/",
            query_size=consensusSize("../data/consensus/ex_hg38_cons.fa_/DF0000004.fa"))
    scoreThresholds("../results/alignments/test_alignments/dfamseq/DF0000244/",
            "../results/alignments/test_alignments/benchmark/DF0000244/",
            query_size=consensusSize("../data/consensus/ex_hg38_cons.fa_/DF0000244.fa"))
    """
    # generateScoreThreshold("../results/alignments/dfamseq/DF0000001_25p35g.sc",
    #                         "../results/alignments/benchmark/DF0000001_25p35g.sc")

