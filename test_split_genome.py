"""
test_split_genome.py

A quick test suite for split_genome.py

AUTHOR(S):
    Eric Yeh
"""

from split_genome import splitGenome, gcBackground, countGCAT

splitGenome("split_human-1mb.fa", batch_length=1000, batch_overlap=50,
            row_length=50)
