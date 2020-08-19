# Python scripts

The scripts found in this directory are used to generate consensus sequence score thresholds. Check /results for some examples of the output produced.

## Splitting genome into GC bins
Since we use a different set of parameters for each GC background, we want to split portions of a genome into their respective GC categories. To do this, we split a genome into 60k batches, compute the GC content of each batch, and place it into the proper bin. The human genome GC bins can be found in /data/hg38bins (genomic and benchmark).

To produce 10 GC bins containing all of a genome from a .fa file, use the /src/bin\_genome.py script.

`$ python3 bin_genome.py fa_file output_dir [-c C]`
- fa\_file: Path to the genome's fa file.
- output\_dir: Path to the directory to place outputted GC bins.
- C: Optional parameter denoting number of bases per row (Shouldn't matter too much).

If changes need to be made to this script, test\_bin\_genome.py can be used to ensure the correctness of the script.

## Generating alignments and score thresholds for a single consensus sequence
Once the GC bins are created, you can start producing alignments and thresholds for a consensus sequence quickly using run\_job.py.

Note that this pipeline is currently only compatible with human genome consensus sequences. If you want to run this for a different genomes, lines 39 and 42 contains the paths to the genomic and benchmark bins. Simply change these to the path to the genome bins of interest.

`$ python3 run_job.py consensus`
- consensus: Path to the consensus sequence's fa file. This fa file should only contain a single sequence.

The alignments will be outputted to /results/genomic\_hits and /results/benchmark\_hits. The score thresholds will be outputted to /results/thresholds.
