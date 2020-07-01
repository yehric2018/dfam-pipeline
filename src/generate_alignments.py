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
import os

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
            g = open(os.path.join(dir_name,
                    line[1:].split()[0] + ".fa"), "a")
        g.write(line)
        line = f.readline()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("fa_file",
            help="path to the fa_file containing consensus sequences")
    args = parser.parse_args()

    splitConsensus(args.fa_file)
