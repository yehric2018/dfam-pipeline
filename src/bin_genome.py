"""
bin_genome.py: Reads in a genome fa file and splits sequences into
bins based on GC-background.

This Python3 script reads in a fa file containing a genome and
creates various bin files for each percentage of GC background.
60 kilobase (default) batches will be read and placed into these
bin files based on their GC background. There will be a 2 kilobase
(default) overlap between batches to capture sequences that might
extend past the batch boundaries.

You can run this script directly to generate bins for a given genome
fa file:
$ python3 bin_genome.py fa_file output_dir [-c C]
where fa_file is the fa file containing the genome, output_dir is the
directory where the batches should be written, and C is the number of
bases per line in the output files (default unlimited bases per line)

AUTHOR(S):
    Eric Yeh
"""

#
# Module imports
#
import argparse
import io
import os

BATCH_LENGTH = 60000
BATCH_OVERLAP = 2000
MIN_BATCH_LENGTH = 40000
ROW_LENGTH = 50
GC_BINS = [35, 37, 39, 41, 43, 45, 47, 49, 51, 53]

class GenomeBatch:
    """
    GenomeBatch stores a batch of bases from a genome, which can be
    written to a corresponding bin based on its GC background.
    """
    def __init__(self, seq_name, start, output_dir, seq="",
            bl=BATCH_LENGTH, bo=BATCH_OVERLAP):
        """
        Initializes this batch with the given name and start index,
        as well as the given sequence if provided.

        Args:
            seq_name: The name of the sequence the batch is being
                recorded from (usually follows the ">" in a fa file.
            start: The index where this batch starts, using zero-based
                indexing.
            seq: initialize this batch with some bases, usually
                overlap from the previous sequence.
        """
        self.seq_name = seq_name.split()[0]
        self.start = start
        self.output_dir = output_dir
        self.seq = seq

        self.batch_length = bl
        self.batch_overlap = bo

        self.gc = self.gcat = 0
        self.__countGCAT__(seq)

    def batchName(self):
        """
        batchName(seq_name, start, end) - Constructs a batch name from
        this batch's sequence.

        The resulting name outputs the name of the original sequence
        within the genome along with the start and end indexes of the
        sequence, using one-based indexing and with both start and end
        being inclusive.

        Returns: a string describing the batch.
            ex. if start = 0 and end = 60000, the returned batch-name
            will be "chr1:1-60000"
        """
        return (self.seq_name + ":" + str(self.start + 1) + "-"
                    + str((self.start + len(self.seq))))

    def gcBackground(self):
        """
        gcBackground(self) - Returns the bin number for this batch,
        which is determined by this batch's gc content. The GC
        background is computed by the following:
            100 * [count(G) + count(C)] /
            [count(G) + count(C) + count(A) + count(T)]

        Additionally, the GC background will be rounded to the
        nearest GC bin number. See GC_BINS for the possible bin
        numbers.

        Returns - The GC bin that this batch will be assigned to when
        written.
        """
        if self.gcat == 0:
            return -1
        gcb = 100.0 * self.gc / self.gcat
        minDist = 100.0
        binNum = -1

        # Can probably use a binary search algorithm later
        #    O(1)/O(N) => O(1)/O(logN)
        for b in GC_BINS:
            if abs(b - gcb) < minDist:
                minDist = abs(b - gcb)
                binNum = b
        return binNum

    def addToBatch(self, seq):
        """
        addToBatch(self, seq) - Appends the given seq to the end of
        this batch. Sequences can be added to the batch until the
        batch holds BATCH_LENGTH bases.
        
        If adding to the sequence reaches or exceeds the BATCH_LENGTH
        limit, this function will return the overlap for the next
        sequence plus the remainder sequence if it did not get added
        to this sequence. Otherwise, this function returns "".

        Args:
            seq - string that will be appended to this batch.

        Returns:
            overlap (last BATCH_OVERLAP bases) + remainder of seq if
                this batch reached BATCH_LENGTH bases.
            "" otherwise
        """
        append = seq[:self.batch_length - len(self.seq)]
        self.__countGCAT__(append)
        self.seq += append
        if len(self.seq) == self.batch_length:
            return (self.seq[-self.batch_overlap:] +
                    seq[len(append):])
        return ""

    def writeToBin(self, row_length=ROW_LENGTH):
        """
        writeToBin(self) - Write this batch to its corresponding bin.
        
        Write to bin corresponding to the batch's estimated GC
        background. The batch's name will first be printed, followed
        by the sequence. For example:
            >seq1:1-60000
            ATCG...

        This method should not be called until the next batch is
        finished being read in case the next batch has to be combined
        with this batch.

        Args:
            row_length - number of bases per row in output file.
        """
        bin_num = self.gcBackground()
        if bin_num == -1:
            return
        bin_file = open(os.path.join(self.output_dir,
                    "bin" + str(bin_num) + ".fa"), "a")
        bin_file.write(">" + self.batchName() + "\n")
        if row_length > 0:
            for i in range(int((len(self.seq) - 1) / row_length) + 1):
                bin_file.write(self.seq[row_length * i:row_length * (i + 1)]
                                + "\n")
        else:
            bin_file.write(self.seq + "\n")
        bin_file.close()

    def combineBatch(self, other):
        """
        combineBatch(self, other) - Adds subsequent batch's sequence
        to this batch's sequence. Using this method allows the
        sequence length of this batch to exceed BATCH_LENGTH.

        The GenomeBatch passed in MUST occur directly after this
        GenomeBatch in the original genome sequence so the overlap
        can be merged.

        Args:
            other - another GenomeBatch that will be appended to
                self.seq
        """
        append = other.seq[self.batch_overlap:]
        self.__countGCAT__(append)
        self.seq += append

    def __countGCAT__(self, seq):
        """
        __countGCAT__(self, seq) - Counts the number of GCAT in the
        given sequence and modifies the GC background accordingly.
        
        Here, the GC background is computed by the following:
            100 * [count(G) + count(C)] /
            [count(G) + count(C) + count(A) + count(T)]

        Each sequence within this batch should only use this method
        once, immediately when it is added to the batch, to avoid
        overcounting the sequence's GC content.

        Args:
            seq: A string of DNA nucleotides.

        Modifies: GC-background of this batch.
        """
        for i in seq.lower():
            if i in "gcat":
                if i in "gc":
                    self.gc += 1
                self.gcat += 1

def binGenome(fa_file, output_dir, batch_length=BATCH_LENGTH,
                batch_overlap=BATCH_OVERLAP, row_length=-1):
    """
    binGenome(fa_file) - Reads given genome FA file, splits it into
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
        output_dir: name of directory to place bin files in.
        batch_length: the maximum number of bases per batch.
        batch_overlap: number of bases in overlapping regions.
        row_length: number of bases per line in output bin files.
    """
    print("Binning " + fa_file + " into " + output_dir + " with " +
        "batch_length=" + str(batch_length) + " and batch_overlap="
        + str(batch_overlap))

    genome = open(fa_file, "r")
    seq_name = genome.readline()[1:-1]
    prev_batch = None
    batch = GenomeBatch(seq_name, 0, output_dir, bl=batch_length,
                        bo=batch_overlap)
    line = genome.readline()
    while line != "":
        line = line[:-1]
        if line != "" and line[0] == ">":
            seq_name = line[1:]
            if prev_batch == None:
                batch.writeToBin(row_length)
                batch = GenomeBatch(seq_name, 0, output_dir,
                        bl=batch_length, bo=batch_overlap)
            elif len(batch.seq) < 40000:
                prev_batch.combineBatch(batch)
                prev_batch.writeToBin(row_length)
                prev_batch = None
                batch = GenomeBatch(seq_name, 0, output_dir,
                        bl=batch_length, bo=batch_overlap)
            else:
                prev_batch.writeToBin(row_length)
                prev_batch = None
                batch.writeToBin(row_length)
                batch = GenomeBatch(seq_name, 0, output_dir,
                        bl=batch_length, bo=batch_overlap)
            line = genome.readline()
            continue
        rem = batch.addToBatch(line)
        if rem != "":
            if prev_batch != None:
                prev_batch.writeToBin(row_length)
            prev_batch = batch
            batch = GenomeBatch(seq_name, batch.start +
                    batch_length - batch_overlap, output_dir, rem,
                    bl=batch_length, bo=batch_overlap)
        line = genome.readline()
    if len(batch.seq) < 40000:
        if prev_batch == None:
            batch.writeToBin(row_length)
        else:
            prev_batch.combineBatch(batch)
            prev_batch.writeToBin(row_length)
            prev_batch = None
    else:
        if prev_batch != None:
            prev_batch.writeToBin(row_length)
            prev_batch = None
        batch.writeToBin(row_length)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_fa_file",
            help="fa file containing genome to parse")
    parser.add_argument("output_dir",
            help="dir where bin files are placed")
    parser.add_argument("-c", type=int, default=-1,
            help="number of bases per row in output")
    args = parser.parse_args()

    binGenome(args.genome_fa_file, args.output_dir, row_length=args.c)
