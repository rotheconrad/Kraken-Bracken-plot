#!/bin/bash
#SBATCH -J bracken
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem=16G 
#SBATCH -t 2:30:00 
#SBATCH -q inferno

cd ${SLURM_SUBMIT_DIR}
db_path=path/to/kraken_db
bracken=path/to/Bracken/bracken

# below generates bracken results for species, genus, class, and family levels. 
# comment out whichever you don't want:
${bracken} -d ${db_path} -i ${input} -o ${input}_S.bracken -r 150 -l 'S'
${bracken} -d ${db_path} -i ${input} -o ${input}_G.bracken -r 150 -l 'G'
${bracken} -d ${db_path} -i ${input} -o ${input}_C.bracken -r 150 -l 'C'
${bracken} -d ${db_path} -i ${input} -o ${input}_F.bracken -r 150 -l 'F'

# USAGE:
# for report in *.k2; do sbatch --export input=${report} bracken.sbatch; done

# feel free to change the variables passed from sbatch call if you don't like
# a growing tail of file extensions on your final product.
