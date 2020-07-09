# Poster presentation
## Background
- Transposable elements (TEs) are sequences of DNA that encode only the function to replicate and move themselves around the genome. Due to their lack of contribution to function, TEs have been formerly been dubbed as "selfish" or "junk" DNA, only having deleterious effects like inserting themselves into "useful" DNA segments. However, scientists are becoming increasingly interested in TEs due to their potential involvement in evolution. The mutations cause by TE insertion can spawn mutations leading to greater genomic diversity, and studying nuances between different TE families can help us deduce the structure of phylogenetic trees.
- Many tools have recently been developed to study these interesting DNA sequences, which make up a significant portion of the human genome. One of these if Dfam, a database which stores TE families, each with a consensus sequence and profile hidden Markov model, in an attempt to classify and annotate these families within genomes.
## Motivation
- Since DNA sequences mutate over time, there are often slight variations between TEs of the same family. To account for these differences, we classify TEs with the Smith-Waterman/BLAST algorithm, which aligns the consensus sequence of a TE family with a sequence of DNA and computes a similarity score. As of now, we have a fixed score threshold, where TE alignments that surpass this threshold are considered part of the TE family. By keeping this threshold high, we can stay safe from potential false matches/hits, but we risk missing potential matches as a consequence. Therefore, in order to maintain high sensitivity in our search, we should find a minimal score threshold that can allow for more matches while keeping the false discovery rate reasonably low.
- RMBlast is a version of BLAST for RepeatMasker
## Objective
- Create a workflow that takes in a genome (fa format) and several hundreds of consensus sequences (fa format)
- Run RMBlast alignment for the consensus sequence against the genome file to get consensus alignments
- Use consensus alignments to generate a consensus score threshold for that family for each GC value
- Questions: How sensitive are the model thresholds to the assembly used? Is it sufficient to use on set of threhsolds for all mammals in which the model is applicable? How significant is the threshold differences of other isochores? What is the expected improvement to annotation if we switch to using per-family thresholds?
## Workflow
- Use Robert's picture to demonstrate the workflow
- Step 1: Split the genome into batches and place into bins based on GC background
- Step 2: Split consensus sequence, run each one once in RMBlast against all bins for its corresponding genome to get alignments
- Step 3: Use regex to parse alignment file and grab scores, calculate E-values using a given equation
- Give a list of parameters used (gap params, Gumbel params, etc)
## Future Work
- Have to calculate score thresholds using both empirical and theoretical methods
- Run several consensus sequences through this workflow to develop the score thresholds
- Analyze these score thresholds and see if we can draw any conclusions
- Find ways to optimize code for running several consensus sequences
- Alignments used in adjudication project
## References/Acknowledgements
- Dfam paper
- Transposable elements paper
- RepeatMasker paper?
