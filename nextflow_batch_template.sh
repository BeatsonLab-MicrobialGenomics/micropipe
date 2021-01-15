#!/bin/bash

#SBATCH --job-name=pipeline
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=s%A.pipeline_assembly.out
#SBATCH --error=s%A.pipeline_assembly.err

source activate nextflow

#Cloud9: It is recommended to run the nextflow command in the background inside a tmux/screen session to avoid potential issues when submitting the pipeline in a batch script

#directory containing the nextflow.config file and the main.nf script
dir=/scratch/uqvmurig/ST131_03c/pipeline_v15
cd ${dir}
out_dir=${dir}/results

#Run A, B or C depending on whether you are starting with ONT fast5 (A or B) or fastq files (C or D) 

#A) Workflow including basecalling, demultiplexing and assembly
#fast5_dir=${dir}/fast5_pass
#csv=${dir}/samplesheet/samples_1_basecalling.csv
#nextflow main.nf --basecalling --demultiplexing --samplesheet ${csv} --outdir ${out_dir} --fast5 ${fast5_dir}

#B) Workflow including basecalling and assembly (skip demultiplexing step)
#fast5_dir=${dir}/fast5_pass
#csv=${dir}/samplesheet/samples_1_basecalling_single_isolate.csv
nextflow main.nf --basecalling --samplesheet ${csv} --outdir ${out_dir} --fast5 ${fast5_dir}

#C) Workflow including demultiplexing and assembly
#fastq_dir=${dir}/fastq
#csv=${dir}/samplesheet/samples_1_basecalling.csv
#nextflow main.nf --demultiplexing --samplesheet ${csv} --outdir ${out_dir} --fastq ${fastq_dir}

#D) Assembly workflow (skip basecalling and demultiplexing step)
csv=${dir}/samples_1.csv
nextflow main.nf --samplesheet ${csv} --outdir ${out_dir}

#to restart the pipeline if something failed, use the -resume flag after correcting the issue
#nextflow main.nf -resume --samplesheet ${csv} --outdir ${out_dir}
