#!/bin/bash
#SBATCH --partition wheeler_lab_large_cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --job-name=DF0000230
#SBATCH --output=%x-slurm-%j.out
#SBATCH --error=%x-slurm-%j.err

python3 ../src/run_job.py ../data/consensus/Dfam_HG38_Families.fa_/DF0000230.fa