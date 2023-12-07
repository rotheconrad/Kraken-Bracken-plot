#!/bin/bash
#SBATCH -J k2
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=10
#SBATCH --mem=160G                     	
#SBATCH -t 2:00:00           
#SBATCH -q inferno

cd ${SLURM_SUBMIT_DIR}
env_path=path/to/kraken_env
db_path=path/to/kraken_db

module load anaconda3  
conda activate ${env_path}                                                                                         

# PAIRED FASTQ:
kraken2 --db ${db_path} --paired ${input}.[12].fastq --output "-" --report ${input}.k2 --threads 10
# should you want to provide gzipped FASTQ, then change the input call and add the --gzip-compressed flag

# INTERLEAVED FASTA
# kraken2 --db ${db_path} ${input} --output "-" --report ${output}.k2 --threads 10

# USAGE:
# for file in *.1.fastq; do sampleid=$(echo ${file} | cut -f 1 -d '.');	sbatch --export	input=${sampleid} kraken2.sbatch; done

# if files exist in your working directory as:
# sample1.1.fastq
# sample1.2.fastq
# sample2.1.fastq
# sample2.2.fastq
# etc.
