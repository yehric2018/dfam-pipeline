#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
generate_alignments.py: Code used for creating alignments for
mammalian genomes, which will be used for developing consensus
sequence score thresholds and for the Hood-Price adjudication
project.

Contains variety of util functions for:
 - separating consensus sequences into their own fa file
 - generating alignments between consensus sequences and genome GC
   bins by running rmblastn
 - Extracting scores from alignment output

You can run this script directly to split a given consensus file:
$ python3 generate_alignments.py fa_file
where fa_file is the path to the fa_file.

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import argparse
import sys
import os
import subprocess

from sequence_util import nearestDivergence

DIV_VALUES = [14, 18, 20, 25, 30]
GAP_PARAMS = {
        14: {"open": -35, "ext": -6},
        18: {"open": -33, "ext": -5},
        20: {"open": -30, "ext": -5},
        25: {"open": -27, "ext": -5}
    }

def splitConsensus(fa_file):
    """
    splitConsensus(fa_file): Opens the given fa_file and separates
    the consensus sequences in the file into separate files, placed
    into a new directory. This new directory will be in the same
    directory as the given fa file.

    If the new directory already exists, its contents will be
    deleted. Otherwise, the directory will be created.

    Args:
        fa_file - name of fa file containing consensus sequences to
            be aligned.

    Returns - Name of directory where consensus fa files are saved.
    """
    dir_name = fa_file + "_"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    else:
        filelist = [ f for f in os.listdir(dir_name) ]
        for f in filelist:
            os.remove(os.path.join(dir_name, f))
    f = open(fa_file, "r")
    g = None
    line = f.readline()
    while line != "":
        if line[0] == ">":
            if g != None:
                g.close()
            g = open(os.path.join(dir_name,
                    line[1:].split()[0] + ".fa"), "a")
        g.write(line)
        line = f.readline()
    f.close()
    return dir_name

class ConsensusSequence:
    """
    Stores information on a given consensus sequence, such as the
    sequence name, average kimura divergence, and the sequence
    itself.

    Fields:
        fname - Name of .fa file this consensus sequence is derived
            from.
        name - Name of the consensus sequence as given in original
            .fa file.
        seq - Sequence of bases making up this consensus sequence.
        divergence - Rounded average kimura divergence as given in
            original .fa file.
        gi - gap init parameter
        ge - gap extension parameter
    """
    def __init__(self, fa_file):
        """
        Initializes this ConsensusSequence with the given parameters.
        Retrieves the name and divergence from the given fa file.

        The first line should be in the following format:
        >DF######## avg_kimura=#### [name]
        The first token in this line will be the name of this
        consensus and the avg_kimura will be the divergence,
        rounded to the nearest divergence value in DIV_VALUES.

        Args:
            fa_file - name of the consensus fa file following the fa
            format.
        """
        f = open(fa_file, "r")
        lines = f.readlines()
        description = lines[0][1:].split()
        self.fname = fa_file
        self.name = description[0]
        self.seq = "".join([x[:-1] for x in lines[1:]])
        self.__setDivergence__(float(description[1].split("=")[1]))
        f.close()

    def __setDivergence__(self, div):
        """
        Helper function called by __init__ to set self.divergence,
        along with other parameters needed for running RMBlast.

        Args:
            div - avg kimura divergence from fa file.
        """
        self.divergence = nearestDivergence(div)
        self.gi = GAP_PARAMS[self.divergence]["open"]
        self.ge = GAP_PARAMS[self.divergence]["ext"]

def runRMBlast(consensus, bin_file, output_dir):
    """
    runRMBlast(consensus, bin_file, output_dir) - Run RMBlast against
    bin_file, generating all the alignments of the given
    ConsensusSequence against that bin and placing the results in
    output_dir.

    The format of the output file produced will be:
        [consensus_name]_[##]p[##]g.sc

        Ex: DF00001_25p49g.sc

    The first column of each entry in the output file is the score
    for that particular alignment, which will be extracted in later
    steps to calculate E-values and false discovery rate.

    Args:
        consensus - A ConsensusSequence generated from a fa file for
            a single consensus sequence.
        bin_file - File containing batches to align against, produced
            from bin_genome.py.
        output_dir - Directory to place output alignment files.
    """
    bin_num = bin_file[-5:-3]
    fname = (consensus.name + "/" + consensus.name + "_" +
            str(consensus.divergence) + "p" + bin_num + "g.sc")
    matrix_file = ("../data/matrices/" + str(consensus.divergence) +
                "p" + bin_num + "g.matrix")

    params = [ "/home/rhubley/scripts/rbn",
                bin_file, consensus.fname,
                "-matrix", matrix_file,
                "-gi", str(consensus.gi),
                "-ge", str(consensus.ge),
                "-minmatch", "7",
                "-masklevel", "101",
                "-minscore", "50",
                "-a", "-r" ]
    print(" ".join(params))

    if not os.path.exists(os.path.join(output_dir, consensus.name)):
        os.mkdir(os.path.join(output_dir, consensus.name))
    fStdout = open(os.path.join(output_dir, fname), "w")
    fStderr = open(os.path.join(output_dir, "stderr"), "w")
    proc = None
    try:
        proc = subprocess.check_call(params, stdout=fStdout,
                                        stderr=fStderr)
    except:
        fStdout.write("rmblast exception: " + str(sys.exc_info()[0]) )
        pass
    fStdout.close()
    fStderr.close()

def generateAlignments(consensus_file, bins_dir, output_dir):
    """
    generateAlignments(consensus_file, bins_dir, output_dir) -
    Wrapper for runRMBlast that generates alignments for every bin in
    bins_dir, printing them all to files in output_dir.

    Produces a ConsensusSequence from the given path to a consensus
    file, then for each bin file in bins_dir, runs them in RMBlast.

    Args:
        consensus_file - Path to a .fa file containing a single
            consensus sequence.
        bins_dir - Path to directory containing bins from a genome,
            produced from bin_genome.py.
        output_dir - Directory to place output alignments files.
    """
    cs = ConsensusSequence(consensus_file)
    binlist = [ b for b in os.listdir(bins_dir) ]
    for b in binlist:
        if b[:9] != "ncResults":
            runRMBlast(cs, os.path.join(bins_dir, b), output_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("fa_file",
            help="path to the fa_file containing consensus " +
                "sequence(s), use -m flag to indicate multiple " +
                "consensuses in the given file")
    parser.add_argument("-m", action="store_true",
            help="find alignments for all consensus sequences in a " +
                    "single .fa file")
    parser.add_argument("bins_dir",
            help="path to directory containing genome bins")
    parser.add_argument("output_dir",
            help="path to directory to put alignment output files")
    args = parser.parse_args()

    if args.m:
        dir_name = splitConsensus(args.fa_file)

        filelist = [ f for f in os.listdir(dir_name) ]
        for f in filelist:
            generateAlignments(os.path.join(dir_name, f),
                            args.bins_dir, args.output_dir)
    else:
        generateAlignments(args.fa_file, args.bins_dir, args.output_dir)
