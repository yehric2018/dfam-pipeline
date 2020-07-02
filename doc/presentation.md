# Poster presentation
## Background
- Transposable elements (TEs) are sequences of DNA that encode only the function to replicate and move themselves around the genome. Due to their lack of contribution to function, TEs have been formerly been dubbed as "selfish" or "junk" DNA, only having deleterious effects like inserting themselves into "useful" DNA segments. However, scientists are becoming increasingly interested in TEs due to their potential involvement in evolution. The mutations cause by TE insertion can spawn mutations leading to greater genomic diversity, and studying nuances between different TE families can help us deduce the structure of phylogenetic trees.
- Many tools have recently been developed to study these interesting DNA sequences, which make up a significant portion of the human genome. One of these if Dfam, a database which stores TE families, each with a consensus sequence and profile hidden Markov model, in an attempt to classify and annotate these families within genomes.
- Since DNA sequences mutate over time, there are often slight variations between TEs of the same family. To account for these differences, we classify TEs with the Smith-Waterman algorithm, which aligns the consensus sequence of a TE family with a sequence of DNA and computes a similarity score. As of now, we have a fixed score threshold, where TE alignments that surpass this threshold are considered part of the TE family. By keeping this threshold high, we can stay safe from potential false matches/hits, but we risk missing potential matches as a consequence. Therefore, in order to maintain high sensitivity in our search, we should find a minimal score threshold that can allow for more matches while keeping the false discovery rate reasonably low.
## Goals/Questions
text goes here
## Design
text goes here
## Future Work
text goes here
## References/Acknowledgements
text goes here
