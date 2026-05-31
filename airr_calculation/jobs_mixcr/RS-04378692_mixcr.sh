#!/bin/bash
#SBATCH --job-name=RS-04378692_mixcr
#SBATCH --output=logs_mixcr/%x.o.%j.log
#SBATCH --error=logs_mixcr/%x.e%j.err
#SBATCH --time=23:59:00
#SBATCH --partition=general-compute --qos=general-compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16000

# module load gcc foss java samtools r-bundle-bioconductor
cd $SLURM_SUBMIT_DIR 

module load gcc foss java



/projects/rpci/songyao/pnfioric/software/mixcr_new/mixcr analyze rna-seq --species hsa --threads 12 -Xmx16g \
    /projects/rpci/songyao/pnfioric/arc_project/airr_calculation/fastq_files/RS-04378692/*RS-04378692*_R1_*.fastq.gz \
    /projects/rpci/songyao/pnfioric/arc_project/airr_calculation/fastq_files/RS-04378692/*RS-04378692*_R2_*.fastq.gz \
    /projects/rpci/songyao/pnfioric/arc_project/airr_calculation/tcr_quant/results/mixcr_RS-04378692