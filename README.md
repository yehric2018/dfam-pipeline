# dfam-pipeline

## Overview
The goal of this project is to create a pipeline that can compute sensitive but conservative score thresholds for each TE family consensus sequence per GC-tuned scoring matrix. By selecting score thresholds maintaining a false discovery rate below 0.2%, we can develop a robust dataset that can be used for Dfam and other projects.

From the resulting thresholds we can answer the following questions:
- Is it sufficient to use one set of thresholds for all taxa in which the model is applicable?
- How significant are the differences between isochores?
- What is the expected improvement to annotation if we switch to using per-family thresholds?

### Transposable elements (TEs)
TEs are DNA sequences that encode the ability to replicate themselves (autonomous) or rely on the products of other TE families to replicate (non-autonomous). Due to their minimal direct contribution to host function, TEs have formerly been labeled as "selfish" or "junk" DNA, often having deleterious effects like inserting themselves into functional DNA segments. However, the success of TEs in many organisms (over 50% of the human genome is derived from TE copies) is an indication of the major role they have played in genome evolution, and their individual copies provide a rich dataset for many areas of research.

### Dfam
Dfam is a database developed at the Hood-Price lab that stores TE family entries consisting of consensus sequences and profile hidden Markov models (HMMs) to classify and annotate families within genomes. The use of Dfam with profile HMMs identified an additional 2.9% of the human genome as being derived from TEs.

### False Discovery Rate (FDR)
The proportion of sequence alignments (annotations) that are expected to be false positives (Type I errors). To account for mutations over time, we identify TE copies with a sequence similarity search algorithm, which aligns the consensus sequence of a TE family with a sequence of DNA and computes a similarity score. In the past, we have used a small set of conservative fixed score thresholds to maintain a low overall FDR, but applying the same set of thresholds to large groups of TE families risks being overly conservative with some families and misses true matches. To maintain high sensitivity in our annotation, we should find a minimal score threshold(s) that can find more matches while keeping the FDR reasonably low, in our case below 0.2%. To estimate the alignment threshold for a fixed FDR, we align a family consensus sequence against a benchmark genome (a realistic sequence devoid of TE content, generated using GARLIC), and a reference genome. The threshold would then be the lowest reference genome alignment score which satisfies the target FDR.

## Scripts
**/src/dfam\_pipeline.py:** Split fa file of consensus sequences into 5 directories, placing each consensus sequence in its own individual fa file.<br />
**/src/jobs\_batch.py:** Runs the pipeline for a directory of consensus sequences. Open 5 new screen sessions and run this script for each of the directories produced by dfam\_pipeline.py so thresholds can be generated concurrently.<br />
<br />
**/src/bin\_genome.py:** pass in genome .fa file, splits genome into batches and categorizes them based on GC-background. <br />
**/src/generate\_alignments.py:** pass in a directory of GC bins and a consensus sequence fa file, produce alignments of the consensus sequence against the bins by running RMBlast. Alignments are stored in files based on GC background of the bins aligned against and gives each alignment a score. <br />
**/src/score\_thresholds.py:** pass in alignment files for a consensus sequence against both the genomic and benchmark bins. Uses the scores from these alignment files to compute a score threshold that maintains a false discovery rate below 0.2%. <br />
**/src/sequence\_utils.py:** a variety of util functions that can be used to manage sequences throughout the pipeline. <br />
**/src/test\_bin\_genome.py:** maintainability suite to ensure correctness of bin\_genome.py works properly.
