#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
sequence_util.py -

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import sqlite3

DIV_VALUES = [14, 18, 20, 25, 30]

def nearestDivergence(div):
    """
    nearestDivergence(div) - Given a divergence value, rounds to the
    nearest divergence value that we use for our matrices and returns
    that value.

    Args:
        div - Divergence value for a TE family

    Returns: div rounded to the nearest value in DIV_VALUES
    """
    minDist = 100.0
    divVal = -1

    for v in DIV_VALUES:
        if abs(v - div) < minDist:
            minDist = abs(v - div)
            divVal = v

    return divVal

def consensusSize(fa_file):
    """
    consensusSize(fa_file) - Given the name of a consensus fa file,
    return the size of the consensus sequence and update the database
    with the sequence size.

    The fa file should contain a single consensus sequence. This
    function will store not only the number of nucleotides in the
    sequence, but also the family's divergence (rounded).

    Args:
        fa_file - path to consensus fa file containing a single
            sequence.

    Returns: size of the given consensus
    """
    consensus = open(fa_file, "r")
    line = consensus.readline()
    splitLine = line[1:].split()
    name = splitLine[0]
    div = nearestDivergence(float(splitLine[1].split("=")[1]))
    size = 0
    while line != "":
        if line[0] != ">":
            size += len(line) - 1
        line = consensus.readline()
    return size

def genomeSize(fa_file):
    """
    genomeSize(fa_file) - Given the name of a genome fa file, return
    the size of the sequence and update the database with this
    genome's size.

    The fa file should contain a single genome (which may consist of
    multiple sequences). This function will count the number of
    nucleotides in the sequence.

    Args:
        fa_file - path to genome fa file containing sequence(s)

    Returns: size of the given genome
    """
    genome = open(fa_file, "r")
    line = genome.readline()
    size = 0
    while line != "":
        if line[0] != ">":
            size += len(line) - 1
        line = genome.readline()
    return size

if __name__ == '__main__':
    #genomeSize("../data/genomes/dfamseq.mask")
    consensusSize("../data/consensus/ex_hg38_cons.fa_/DF0000001.fa")
