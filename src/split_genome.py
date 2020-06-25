"""
split_genome.py: Reads in a genome fa file and splits sequences into
bins based on GC-background.

This Python3 script reads in a fa file containing a genome and
creates various bin files for each percentage of GC background.
60 kilobase (default) batches will be read and placed into these
bin files based on their GC background. There will be a 2 kilobase
(default) overlap between batches to capture sequences that might
extend past the batch boundaries.

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import argparse
import io

BATCH_LENGTH = 60000
BATCH_OVERLAP = 2000
MIN_BATCH_LENGTH = 40000

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

def batchName(seq_name, start, end):
    """
    batchName(seq_name, start, end) - Constructs a batch name from
    the given seq_name and start and end indices.

    The arguments passed in should use 0-based indexing (start is
    inclusive and end is exclusive). However, the resulting batch
    name converts this to indexing starting at 1, causing start to be
    incremented by 1 in the batch name. The end index in the batch
    name also becomes inclusive, but since we are using 1-based
    indexing, end will be identical.

    Args:
        seq_name: string of the original sequence name
        start: int starting index within original sequence
        end: int ending index within original sequence

    Returns: a string describing the batch.
        ex. if start = 0 and end = 60000, the returned batch-name
        will be "chr1: 1-60000"
    """
    return seq_name + ": " + str(start + 1) + "-" + str(end)

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
    if gc + at == 0:
        return -1
    return int(round(100.0 * gc / (gc + at)))

def writeBatchToBin(batch_name, batch, bin_num, row_length=-1):
    """
    writeBatchToBin(batch_name, batch, bin_num)

    Appends the given batch to a fa file of the given bin_num if
    bin_num >= 0.
    
    Args:
        batch_name: a string describing this batch. It should include
            the title of the original sequence and the start/end
            indices of the batch. (ex. chr1: 1-60000)
        batch: the string of bases making of this batch.
        bin_num: an integer representing the percentage of GC content
            within the batch. If bin_num < 0, don't write this batch.
        row_length: indicates the number of bases per output line. If
            row_length <= 0, the whole batch is printed to one row.
    """
    if bin_num < 0:
        return
    bin_file = open("bins/bin" + str(bin_num) + ".fa", "a")
    bin_file.write(">" + batch_name + "\n")
    if row_length > 0:
        for i in range(int((len(batch) - 1) / row_length) + 1):
            bin_file.write(batch[row_length * i:row_length * (i + 1)]
                            + "\n")
    else:
        bin_file.write(batch + "\n")
    bin_file.close()

def splitGenome(fa_file, batch_length=BATCH_LENGTH,
                batch_overlap=BATCH_OVERLAP, row_length=-1):
    """
    splitGenome(fa_file) - Reads given genome FA file, splits it into
    batches, and sorts the batches into bins based on the batch's GC
    background.

    For each present GC background in the given batches, an fa file
    will be produced for that bin titled "bin[GC content].fa", placed
    in a folder called "bins". Each of these bins will contain the
    batches from the input genome that match that corresponding GC
    background.

    Args:
        fa_file: the name of the fa file containing the genome to be
            split into bins.
        batch_length: the maximum number of bases per batch.
        batch_overlap: number of bases in overlapping regions.
        row_length: number of bases per line in output bin files.
    """
    genome = open(fa_file, "r")
    seq_name = genome.readline()[1:-1]
    batch = ""
    line = ""
    start = end = gc = at = 0
    while True:
        if line == "":
            line = genome.readline()[:-1]
        if line == "":
            writeBatchToBin(batchName(seq_name, start, end),
                            batch, gcBackground(gc, at), row_length)
            break
        if line[0] == ">":
            writeBatchToBin(batchName(seq_name, start, end),
                            batch, gcBackground(gc, at), row_length)
            seq_name = line[1:] # truncate '#' and '\n'
            start = 0
            gc = at = end = 0
            batch = ""
            line = ""
            continue
        if len(batch) + len(line) >= batch_length:
            read_bases = batch_length - len(batch)
            end += read_bases
            batch += line[:read_bases]
            extra = line[read_bases:]
            gcat = countGCAT(line)
            gc += gcat[0]
            at += gcat[1]
            writeBatchToBin(batchName(seq_name, start, end),
                            batch, gcBackground(gc, at), row_length)
            start = end - batch_overlap
            end = start
            line = batch[-batch_overlap:] + extra
            batch = ""
            gc = at = 0
        else:
            batch += line # truncate '\n'
            end += len(line)
            gcat = countGCAT(line)
            gc += gcat[0]
            at += gcat[1]
            line = ""

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_fa_file",
            help="fa file containing genome to parse")
    parser.add_argument("-c", type=int, default=-1,
            help="number of bases per row in output")
    args = parser.parse_args()

    splitGenome(args.genome_fa_file, row_length=args.c)
