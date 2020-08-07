#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
cache_hits.py - Creates a JSON file from the hits produced by an
alignment file and the score thresholds. The JSON will have the
following structure:

{
    genomic: {
        xxpxxg: [ ... ],
        xxpxxg: [ ... ],
        ...
    },
    benchmark: { ... },
}

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import argparse
import os
import json

from score_thresholds import readScoresFromFile

def writeToFile(hits, fname):
    """
    writeToFile(hits, f) - Write the given json to the given output
    file.

    Args:
        hits - JSON object of hits and thresholds to write
        fname - path to file to write json to.
    """
    with open(fname, "w") as f:
        f.write(json.dumps(hits))

def writeJson(genomic_hits, benchmark_hits, cache_dir):
    """
    writeJson(genomic_hits, benchmark_hits, cache_dir) -

    Reads the scores from genomic hits and benchmark hits of a
    consensus sequence for each GC-tuned matrix of a consensus
    sequence, storing these scores in a sorted array in decreasing
    order.

    Creates a JSON of all the genomic/benchmark hits for a given
    consenus sequence for each GC background and writes the JSON to
    cache_dir.

    Args:
        genomic_hits - path to directory containing all alignments
            of a consensus sequence against the genome.
        benchmark_hits - path to directory containing all alignments
            of a consensus sequnece against the benchmark genome.
        cache_dir - output directory where json will be written to.
    """
    hits = {"genomic": {}, "benchmark": {}}
    print(hits)
    consensus = os.path.dirname(genomic_hits).split("/")[-1]
    print(os.listdir(genomic_hits))
    for sc in os.listdir(genomic_hits):
        matrix = sc.split("_")[1][:-3]
        hits["genomic"][matrix] = readScoresFromFile(
                    os.path.join(genomic_hits, sc))
        hits["benchmark"][matrix] = readScoresFromFile(
                    os.path.join(benchmark_hits, sc))
    writeToFile(hits, cache_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("genomic_hits",
            help="path to directory containing genomic hits")
    parser.add_argument("benchmark_hits",
            help="path to directory containing benchmark hits")
    parser.add_argument("cache_dir",
            help="path to directory to write output jsons")
    args = parser.parse_args()

    print(os.listdir(args.benchmark_hits))
    for consensus in os.listdir(args.benchmark_hits):
        if consensus[:2] == "DF":
            writeJson(os.path.join(args.genomic_hits, consensus),
                    os.path.join(args.benchmark_hits, consensus),
                    os.path.join(args.cache_dir, consensus + ".json"))
