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

DIV_VALUES = [14, 18, 20, 25, 30]

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
        self.name = description[0]
        self.seq = "".join([x[:-1] for x in lines[1:]])
        self.__setDivergence__(float(description[1].split("=")[1]))
        f.close()

    def __setDivergence__(self, div):
        """
        Helper function called by __init__ to set self.divergence.

        Args:
            div - avg kimura divergence from fa file.
        """
        minDist = 100.0
        divVal = -1

        for v in DIV_VALUES:
            if abs(v - div) < minDist:
                minDist = abs(v - div)
                divVal = v
        self.divergence = divVal

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("fa_file",
            help="path to the fa_file containing consensus sequences")
    args = parser.parse_args()

    dir_name = splitConsensus(args.fa_file)

    filelist = [ f for f in os.listdir(dir_name) ]
    for f in filelist:
        cs = ConsensusSequence(os.path.join(dir_name, f))
        print(cs.name + ": " + str(cs.divergence))
        print(cs.seq + "\n")
