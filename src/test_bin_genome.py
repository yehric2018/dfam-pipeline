#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
test_bin_genome.py

A quick test suite for bin_genome.py.

AUTHOR(S):
    Eric Yeh
"""

import os
import math

from bin_genome import binGenome

GC_BINS = [35, 37, 39, 41, 43, 45, 47, 49, 51, 53]

def failTest(err):
    print("FAILURE: " + err)

def clearDirectory(path):
    """
    clearDirectory(path) - Quickly remove (permanently) all the files
    in the given directory
    """
    filelist = [ f for f in os.listdir(path) ]
    for f in filelist:
        os.remove(os.path.join(path, f))

def countGCAT(seq):
    """
    countGCAT(seq) - Counts the number of GC bases and AT bases in
    seq.

    Args:
        seq: A string of DNA nucleotides.

    Returns:
        A tuple (gc, at), where 'gc' is the count of G and C bases and
        'at' is the count of A and T bases.
    """
    gc = 0
    at = 0
    for base in seq.lower():
        if base == 'g' or base == 'c':
            gc += 1
        elif base == 'a' or base == 't':
            at += 1
    return (gc, at)

def gcBackground(gc, at):
    """
    gcBackground(gc, at) - Computes the GC background from the given
    numbers of GC bases and AT bases. The GC background is computed
    as (100 * # of gc bases / (# of gc bases + # of at bases))
    rounded to the nearest integer.

    If gc + at = 0, return -1 instead.

    Args:
        gc: number of G and C nucleotides
        at: number of A and T nucleotides

    Returns: an int representing the GC background from the given
    parameters, or -1 if gc + at = 0.
    """
    gcat = gc + at
    if gcat == 0:
        return -1
    gcb = 100.0 * gc / gcat
    minDist = 100.0
    binNum = -1

    # Can probably use a binary search algorithm later
    #    O(1)/O(N) => O(1)/O(logN)
    for b in GC_BINS:
        if abs(b - gcb) < minDist:
            minDist = abs(b - gcb)
            binNum = b
    return binNum

def test_binGenome(fa_file, batch_length = 60000, batch_overlap = 2000):
    clearDirectory("testbins")
    binGenome(fa_file, "testbins/", batch_length, batch_overlap)
    filelist = [ os.path.join("testbins", f) for f in os.listdir("testbins") ]
    seqDict = {}

    # Check that batches are correctly binned, while constructing dict
    # to use for later tests.
    for fname in filelist:
        gcContent = int(fname[12:-3])
        with open(fname, "r") as f:
            lines = f.readlines()
            if len(lines) % 2 != 0:
                failTest(fname + " has incorrect format")
            for i in range(int(len(lines) / 2)):
                batchName = lines[2 * i][1:-1]
                batch = batchName.split(":")
                seqName = batch[0]
                seqStart = int(batch[1].split("-")[0])
                seq = lines[2 * i + 1][:-1]
                if seqName not in seqDict:
                    seqDict[seqName] = {}
                seqDict[seqName][seqStart] = seq
                gcat = countGCAT(seq)
                if gcContent != gcBackground(gcat[0], gcat[1]):
                    failTest(batchName + "in wrong bin, found actual"
                            + " GC content = " + str(gcBackground(gcat[0], gcat[1])))

    # Make sure batch sizes are followed, and overlapping regions match
    # Also construct a new dict for seq_name to sequence, appending seqs together
    for item in seqDict.items():
        seqName = item[0]
        sequence = item[1]
        frags = sorted(sequence.keys())
        for i in range(len(frags) - 1):
            if len(sequence[frags[i]]) != batch_length:
                failTest(seqName + " has incorrect batch length at i = " + str(i) +
                        ", actual length is " + str(len(sequence[frags[i]])))
            if sequence[frags[i]][-batch_overlap:] != sequence[frags[i + 1]][:batch_overlap]:
                failTest(seqName + " has incorrect batch overlap after i = " + str(i))

    # Compare to sequences of original fa file

    print("Tests finished for " + fa_file)
    clearDirectory("testbins")



if __name__ == '__main__':
    test_binGenome("../data/test_data/short-human.fa")
    test_binGenome("../data/test_data/short-human.fa", 8, 3)
    test_binGenome("../data/test_data/split_human-1mb.fa")
    test_binGenome("../data/test_data/human-1mb.fa")
    test_binGenome("../data/test_data/split_human-1mb.fa", 1234, 50)
