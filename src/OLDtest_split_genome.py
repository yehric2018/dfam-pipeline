"""
test_split_genome.py

A quick test suite for split_genome.py. Verifies that batches are
properly produced for given batch size and overlap size, and that
each batch is placed into the correct bin based on their GC content.

AUTHOR(S):
    Eric Yeh
"""

import os
import math

from split_genome import splitGenome, gcBackground, countGCAT

def clearDirectory(path):
    """
    clearDirectory(path) - Quickly remove (permanently) all the files
    in the given directory
    """
    filelist = [ f for f in os.listdir(path) ]
    for f in filelist:
        os.remove(os.path.join(path, f))

class GenomeTester:
    """
    Stores all DNA sequences listed in the given genome file
    """
    def __init__(self, fa_file, batch_length, batch_overlap, row_length):
        self.fa_file = fa_file
        self.batch_length = batch_length
        self.batch_overlap = batch_overlap
        self.row_length = row_length
        self.seqs = []
        genome = open(fa_file, "r")
        lines = genome.readlines()
        seq_name = ""
        seq = ""
        for line in lines:
            if line[0] == ">":
                if seq_name != "":
                    self.seqs.append(DNASequence(seq_name, seq))
                    seq = ""
                seq_name = line[1:-1]
            else:
                seq += line[:-1]
        self.seqs.append(DNASequence(seq_name, seq))

    def test_produceBins(self):
        """
        test_produceBins(self): Creates bin files for this genome fa
        file so other tests can be run.

        IMPORTANT: Run this test before other tests, and only run this
        test once.
        """
        try:
            splitGenome(self.fa_file, self.batch_length,
                    self.batch_overlap, self.row_length)
            print("SUCCESS: test_produceBins")
        except IOError:
            print("FAILURE: test_produceBins")

    def test_countBatches(self):
        """
        test_countBatches(self): Makes sure the correct number of
        batches are created for the fa file.

        The number of batches created should depend on the batch
        length, batch overlap, and length of the sequence. Due to the
        overlap, only batch_length - batch_overlap new bases will be
        covered for the next batch. Thus, the number of batches
        produced from a given sequence can be computed as follow:

        ceil(seq_length / (batch_length - batch_overlap))

        We round up to the nearest int since the remainder sequence
        at the end will likely be its own batch.

        TODO: Should also account for the edge case where an entire
        batch is made up of N's or non-gcat bases. In this case, a
        batch would not get made for that group of bases.
        """

    def printSeqs(self):
        for seq in self.seqs:
            print(seq)

class DNASequence:
    def __init__(self, seq_name, seq):
        self.seq_name = seq_name
        self.seq = seq

    def getName(self):
        return self.seq_name

    def getSequence(self):
        return self.seq

    def __str__(self):
        return ">" + self.seq_name + "\n" + self.seq + "\n"

def test_numBatches(tester, batch_length, batch_overlap):
    return False

clearDirectory("bins")
tester = GenomeTester("../data/test_data/short-human.fa", 100, 5, 50)
tester.test_produceBins()
tester.printSeqs()
clearDirectory("bins")
