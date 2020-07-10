# Poster presentation
## Background
- Transposable elements (TEs) are sequences of DNA that encode only the function to replicate and move themselves around the genome. Due to their lack of contribution to function, TEs have been formerly been dubbed as "selfish" or "junk" DNA, only having deleterious effects like inserting themselves into "useful" DNA segments. However, scientists are becoming increasingly interested in TEs due to their potential involvement in evolution. The mutations cause by TE insertion can spawn mutations leading to greater genomic diversity, and studying nuances between different TE families can help us deduce the structure of phylogenetic trees.
- Many tools have recently been developed to identify these repeat sequences, which make up a significant portion of the human genome, and store them for further analysis. One of these is Dfam, a database which stores TE family entries consisting of a profile hidden Markov model in an attempt to classify and annotate these families within genomes with greater sensitivity. The use of Dfam allowed for a 2.9% increase in the percentage of the human genome annotated. `Should I define consensus sequence here as well?`
- To maintain a low false discovery rate, we annotate a simulated sequence preserving dinucleotide frequencies which should only yield false hits. We estimate that each hit in the simulated sequence equates to one false hit in the genomic sequence.
## Motivation
- Since DNA sequences mutate over time, there are often slight variations between TEs of the same family. To account for these differences, we classify TEs with the Smith-Waterman/BLAST algorithm, which aligns the consensus sequence of a TE family with a sequence of DNA and computes a similarity score. In the past, we have a used fixed score threshold, where TE alignments that surpass this threshold are considered part of the TE family. By keeping this threshold high, we can stay safe from potential false matches/hits, but we risk missing older and more diverged matches as a consequence. To maintain high sensitivity in our annotation, we should find a minimal score threshold that can allow for more matches while keeping the false discovery rate reasonably low, in our case below 0.2%.
- <span style="color:red">Do we need a hypothesis?</span>
## Objective
- Our goal is to create a workflow that can compute sensitive but conservative score thresholds for each TE family consensus sequence per GC background. From the resulting thresholds we can answer the following questions:
  - How sensitive are the model thresholds to the assembly used?
  - Is it sufficient to use one set of threhsolds for all mammals in which the model is applicable?
  - How significant is the threshold differences of other isochores?
  - What is the expected improvement to annotation if we switch to using per-family thresholds?
## Workflow
- Use Roberts picture to demonstrate the workflow
- Step 1: Split the genome into batches and place into bins based on GC background
- Step 2: Split consensus sequence, run each one once in RMBlast against all bins for its corresponding genome to get alignments
- Step 3: Use regex to parse alignment file and grab scores, calculate E-values using a given equation
- Give a list of parameters used (gap params, Gumbel params, etc)
## Future Work
- Run several consensus sequences through this workflow to develop the score thresholds
- Analyze these score thresholds and see if we can draw any conclusions
- Find ways to optimize code for running several consensus sequences in a cluster
- Alignments used in adjudication project `Should I describe this in the background section?`
## References/Acknowledgements
- Dfam paper
- Transposable elements paper
- RepeatMasker paper?
