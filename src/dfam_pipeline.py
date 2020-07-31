#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
dfam_pipeline.py - Given a set of consensus sequences, compute 10
score thresholds for each sequence for the given genome, one for each
GC tuned matrix. Also generates alignments for the consensus sequence
against the given genome.

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import argparse
import os
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("consensus",
        help="a fa file containing one or more consensus sequences")
    args = parser.parse_args()

    consensus_fa = args.consensus

    dir_names = [consensus_fa + "_1", consensus_fa + "_2", consensus_fa + "_3",
                consensus_fa + "_4", consensus_fa + "_5"]
    for dir_name in dir_names:
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        else:
            filelist = os.listdir(dir_name)
            for f in filelist:
                os.remove(os.path.join(dir_name, f))

    f = open(consensus_fa, "r")
    g = None
    i = 0
    line = f.readline()
    while line != "":
        if line[0] == ">":
            if g != None:
                g.close()
            g = open(os.path.join(dir_names[i % 5], line[1:].split()[0] + ".fa"), "a")
            i += 1
        g.write(line)
        line = f.readline()

    f.close()
    g.close()

if __name__ == '__main__':
    main()
