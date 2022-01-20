#!/bin/bash

#SBATCH --job-name=micropipe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=s%A.micropipe_guppy3.6.1_cpu_12samples_72h.out
#SBATCH --error=s%A.micropipe_guppy3.6.1_cpu_12samples_72h.err
#SBATCH --time=96:00:00
#SBATCH --partition='longq'

module load nextflow/20.07.1-multi
#source activate nextflow
module load singularity/3.6.4 

#Cloud9: It is recommended to run the nextflow command in the background inside a tmux/screen session to avoid potential issues when submitting the pipeline in a batch script

#directory containing the nextflow.config file and the main.nf script
dir=/scratch/director2172/vmurigneux/micropipe
cd ${dir}
#datadir=${dir}/test_data
datadir=${dir}/Illumina
#out_dir=${dir}/results_3.6.1_gpu
out_dir=${dir}/results_3.6.1_cpu
#Run A, B or C depending on whether you are starting with ONT fast5 (A or B) or fastq files (C or D) 

#A) Workflow including basecalling, demultiplexing and assembly
fast5_dir=${dir}/fast5_pass
#fast5_dir=${dir}/fast5_pass/test
csv=${dir}/test_data/samples_all_basecalling.csv
#csv=${dir}/test_data/samples_1_basecalling.csv
#nextflow main.nf --gpu true --basecalling  -profile zeus --slurm_account='director2172' --demultiplexing --samplesheet ${csv} --outdir ${out_dir} --fast5 ${fast5_dir} --datadir ${datadir}
nextflow main.nf --gpu false --basecalling --guppy_num_callers 16 -profile zeus --slurm_account='director2172' --demultiplexing --samplesheet ${csv} --outdir ${out_dir} --fast5 ${fast5_dir} --datadir ${datadir}

#B) Workflow including basecalling and assembly (skip demultiplexing step)
#fast5_dir=${dir}/fast5_pass
#csv=${dir}/test_data/samples_1_basecalling_single_isolate.csv
#nextflow main.nf --basecalling --samplesheet ${csv} --outdir ${out_dir} --fast5 ${fast5_dir} --datadir ${datadir}

#C) Workflow including demultiplexing and assembly
#fastq_dir=${dir}/fastq
#csv=${dir}/test_data/samples_1_basecalling.csv
#nextflow main.nf --demultiplexing --samplesheet ${csv} --outdir ${out_dir} --fastq ${fastq_dir} --datadir ${datadir}

#D) Assembly workflow (skip basecalling and demultiplexing step)
#csv=${dir}/test_data/samples_1.csv
#nextflow main.nf --samplesheet ${csv} --outdir ${out_dir} --datadir ${datadir}

#to restart the pipeline if something failed, use the -resume flag after correcting the issue
#nextflow main.nf -resume --samplesheet ${csv} --outdir ${out_dir} --datadir ${datadir}
