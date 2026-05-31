# Running Snakemake for all 864 samples
# Peter Fiorica
# 14 May 2026

cd /projects/rpci/songyao/pnfioric/software/Pantry/phenotyping

module load gcc foss snakemake bedtools samtools star r-bundle-bioconductor subread
export PATH=/projects/rpci/songyao/pnfioric/software/gffread:$PATH
export PATH=$PATH:/projects/rpci/songyao/pnfioric/tils_r01_RNAseq/splicing/regtools/build
export SNAKEMAKE_TMPDIR=/vscratch/grp-songyao/pnfioric/temp_dir
export LD_LIBRARY_PATH=/cvmfs/soft.ccr.buffalo.edu/versions/2023.01/easybuild/software/Core/gcccore/11.2.0/lib64:$LD_LIBRARY_PATH

snakemake -s /projects/rpci/songyao/pnfioric/arc_project/pantry/Snakefile -j 950 --latency-wait 60 --profile slurm 
