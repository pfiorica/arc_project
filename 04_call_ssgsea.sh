#!/bin/bash
#SBATCH --job-name=ssgsea
#SBATCH --output=logs/R-%x.o.%j.log
#SBATCH --error=logs/R-%x.e%j.err
#SBATCH --time=23:59:00
#SBATCH --partition=general-compute --qos=general-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64000



module load gcc foss java samtools r-bundle-bioconductor


cd /projects/rpci/songyao/pnfioric/arc_project

Rscript 04_ssGSEA_Script.R
