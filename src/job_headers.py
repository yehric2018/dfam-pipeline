#! /usr/bin/env python
import os

f = open("../data/consensus/hg38list.txt", "r")
seqs = f.read().splitlines()
f.close()

for seq in [s[:-3] for s in seqs]:
    job = '''#!/bin/bash
#SBATCH --partition wheeler_lab_large_cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --job-name={name}
#SBATCH --output=%x-slurm-%j.out
#SBATCH --error=%x-slurm-%j.err

python3 ../src/run_job.py {fpath}'''

    g = open(os.path.join("../bin/jobs/", seq + ".sh"), "w")

    g.write(job.format(name = seq, fpath = os.path.join(
            "../data/consensus/Dfam_HG38_Families.fa_", seq + ".fa")))

    g.close()

f = open("../bin/batch_run.sh", "w")
for seq in [s[:-3] for s in seqs]:
    f.write("batch ../bin/jobs/{name}.sh\n".format(name = seq))
f.close()
