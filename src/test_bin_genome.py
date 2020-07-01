"""
test_bin_genome.py

A quick test suite for bin_genome.py.

AUTHOR(S):
    Eric Yeh
"""

import os
import math

from split_genome import gcBackground, countGCAT
from bin_genome import binGenome

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

def test_splitGenome(fa_file, batch_length = 60000, batch_overlap = 2000):
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
    test_splitGenome("../data/test_data/short-human.fa")
    test_splitGenome("../data/test_data/short-human.fa", 8, 3)
    test_splitGenome("../data/test_data/split_human-1mb.fa")
    test_splitGenome("../data/test_data/human-1mb.fa")
    test_splitGenome("../data/test_data/split_human-1mb.fa", 1234, 50)