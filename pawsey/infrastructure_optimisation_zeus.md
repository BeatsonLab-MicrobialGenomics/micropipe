microPIPE on Zeus/Topaz @ Pawsey
===========

---

# Accessing tool/workflow

The workflow can be downloaded from the GitHub page https://github.com/BeatsonLab-MicrobialGenomics/micropipe using the command: 
```
git clone https://github.com/BeatsonLab-MicrobialGenomics/micropipe.git
```
---

# Installation

* **[Nextflow](https://www.nextflow.io/)**  
A modified version of Nextflow, capable of submitting jobs to Zeus, Topaz and Magnus, has been installed as a system module and can be accessed with the command:
```
module load nextflow/20.07.1-multi
```
* **[Singularity](https://singularity.lbl.gov/install-linux)**  
Singularity has been installed as a system module and can be accessed with the command:
```
module load singularity/3.6.4 
```
* **Guppy** (3.6.1 was the latest working version)   
Due to the Oxford Nanopore Technologies terms and conditions, we are not allowed to redistribute the Guppy software either in its binary form or packaged form e.g. Docker or Singularity images. Therefore users will have to either install Guppy, provide a container image or start the pipeline from the basecalled fastq files.  See [Usage](https://github.com/BeatsonLab-MicrobialGenomics/micropipe#usage) section for instructions. 
---

# Quickstart tutorial

A tutorial is available on the GitHub page: https://github.com/BeatsonLab-MicrobialGenomics/micropipe#usage. The steps are summarised below including the specific instructions required to run the pipeline at Pawsey Zeus.  

**1. Prepare the Nextflow configuration file (nextflow.config)**  
Use the configuration file to run microPIPE at Pawsey Zeus [here](./nextflow.config).

**2. Prepare the samplesheet file (csv)**  
See instructions at the microPIPE [GitHub page](https://github.com/BeatsonLab-MicrobialGenomics/micropipe#usage), section 2. Prepare the samplesheet file. 

**3. Prepare the slurm script (e.g. nextflow_batch_template.sh)**  
The pipeline will be launched using a Slurm script submitted to Zeus. This script will load the required modules, define the input/output directories and files, and include the nextflow command line with optional parameters. Note that the configuration profile definition for the Zeus cluster should be specified when launching the pipeline execution by using the "-profile zeus" command line option, as well as the slurm account allocation by using the "--slurm_account='director2172'" command line option (replace 'director2172' by your account identifier). 
```
#!/bin/bash

#SBATCH --job-name=micropipe
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=s%A.micropipe_guppy3.6.1_gpu_12samples.out
#SBATCH --error=s%A.micropipe_guppy3.6.1_gpu_12samples.err
#SBATCH --time=24:00:00

module load nextflow/20.07.1-multi
module load singularity/3.6.4 

#directory containing the nextflow.config file and the main.nf script
dir=/scratch/director2172/vmurigneux/micropipe
cd ${dir}
datadir=${dir}/Illumina
out_dir=${dir}/results_3.6.1_gpu

#Run A, B or C depending on whether you are starting with ONT fast5 (A or B) or fastq files (C or D) 

#A) Workflow including basecalling, demultiplexing and assembly
fast5_dir=${dir}/fast5_pass
csv=${dir}/test_data/samples_all_basecalling.csv
nextflow main.nf --gpu true --basecalling  -profile zeus --slurm_account='director2172' --demultiplexing --samplesheet ${csv} --outdir ${out_dir} --fast5 ${fast5_dir} --datadir ${datadir}
#nextflow main.nf --gpu false --basecalling --guppy_num_callers 16 -profile zeus --slurm_account='director2172' --demultiplexing --samplesheet ${csv} --outdir ${out_dir} --fast5 ${fast5_dir} --datadir ${datadir}

#B) Workflow including basecalling and assembly (skip demultiplexing step)
#fast5_dir=${dir}/fast5_pass
#csv=${dir}/test_data/samples_1_basecalling_single_isolate.csv
#nextflow main.nf --basecalling --samplesheet ${csv} --outdir ${out_dir} --fast5 ${fast5_dir} --datadir ${datadir} -profile zeus --slurm_account='director2172'

#C) Workflow including demultiplexing and assembly
#fastq_dir=${dir}/fastq
#csv=${dir}/test_data/samples_1_basecalling.csv
#nextflow main.nf --demultiplexing --samplesheet ${csv} --outdir ${out_dir} --fastq ${fastq_dir} --datadir ${datadir} -profile zeus --slurm_account='director2172'

#D) Assembly workflow (skip basecalling and demultiplexing step)
#csv=${dir}/test_data/samples_1.csv
#nextflow main.nf --samplesheet ${csv} --outdir ${out_dir} --datadir ${datadir} -profile zeus --slurm_account='director2172'

#to restart the pipeline if something failed, use the -resume flag after correcting the issue
#nextflow main.nf -resume --samplesheet ${csv} --outdir ${out_dir} --datadir ${datadir} -profile zeus --slurm_account='director2172'
```

**4. Run the pipeline by submitting a job at Pawsey Zeus**
```
sbatch nextflow_batch_template.sh
```
---

# Optimisation required

MicroPIPE was originally developed on a cluster for which the jobs could be submitted to both CPU nodes and a GPU node.  At Pawsey, the CPU and GPU nodes are accessed from different clusters ie Zeus (CPU) and Topaz (GPU). Therefore, a modified version of Nextflow, capable of submitting jobs to Zeus, Topaz and Magnus, has been installed as a system module.  
* The modified Nextflow module should be loaded prior to running the main nextflow command by using ```module load nextflow/20.07.1-multi```.   
* The MicroPIPE pipeline will be launched using a Slurm script submitted to Zeus.  
* As a result, Nextflow will automatically submit the GPU tasks to Topaz and the CPU tasks to Zeus.  

Here is a template script to hack Nextflow for multiple clusters (thanks to [@marcodelapierre](https://github.com/marcodelapierre)):   
https://github.com/marcodelapierre/toy-gpu-nf/blob/master/extra/install-nextflow-hack-slurm-multi-cluster.sh


---

# Infrastructure usage and benchmarking

---

## Summary

## Exemplar 1: Assembly of 12 *E.coli* ST131 samples using GPU and CPU resources
You can collect usage metrics from your Canu run using the NCI Gadi optimised workflow using scripts available on the Sydney Informatics Hub, University of Sydney GitHub repository.
* We used the *E.coli* data from the [microPIPE publication](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07767-z) available from the NCBI SRA [BioProject PRJNA679678](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA679678/) (Oxford Nanopore) and the [BioProject PRJEB2968](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB2968) (Illumina). 

* See Nextflow configuration file used [here](./nextflow.config) and slurm submission script [here](./nextflow_batch_template.sh). 
* See Nextflow [HTML execution report](./micropipe_ecoli_ST131_pawsey_guppy3.6.1_gpu.report.html), [trace report](./micropipe_ecoli_ST131_pawsey_guppy3.6.1_gpu.trace.txt) and [HTML processes execution timeline](./micropipe_ecoli_ST131_pawsey_guppy3.6.1_gpu.timeline.html). 

* The table below summarised the assembly results for each strain. 

|Strain|Chromosome/plasmid|Size (bps)|Circularised?|
|-------|:-----:|:-----:|:-----:|
|S24EC| Chromosome <br> Plasmid A | 5078304 <br> 114708 | Yes <br> Yes |    
|S34EC| Chromosome <br> Plasmid A <br> Plasmid B | 5050427 <br> 153321 <br> 108135 | Yes <br> Yes <br> Yes |    
|S37EC| Chromosome <br> Plasmid A <br> Plasmid B | 4981928 <br> 157642 <br> 61072 | Yes <br> Yes <br> Yes |    
|S39EC| Chromosome <br> Plasmid A <br> Plasmid B <br> Plasmid C <br> Plasmid D <br> Plasmid E <br> Plasmid F | 5054402 <br> 141007 <br> 94979 <br> 68049 <br> 62085 <br> 2070 <br> 1846 | Yes <br> Yes <br> Yes <br> Yes <br> Yes <br> Yes <br> Yes |    
|S65EC| Chromosome <br> Plasmid A | 5205011 <br> 147412 | Yes <br> Yes |    
|S96EC| Chromosome <br> Plasmid A <br> Plasmid B <br> Plasmid C <br> Plasmid D | 5069496 <br> 164355 <br> 115965 <br> 14479 <br> 4184 | Yes <br> Yes <br> Yes <br> Yes <br> Yes |  
|S97EC| Chromosome <br> Plasmid A <br> Plasmid B <br> Plasmid C <br> Plasmid D | 5178868 <br> 166099 <br> 96788 <br> 4092 <br> 3209 | Yes <br> Yes <br> Yes <br> Yes <br> Yes |  
|S112EC| Chromosome <br> Plasmid A <br> Plasmid B <br> Plasmid C <br> Plasmid D | 5020013 <br> 161028 <br> 68847 <br> 5338 <br> 4136 | Yes <br> Yes <br> Yes <br> Yes <br> Yes |  
|S116EC| Chromosome <br> Plasmid A <br> Plasmid B <br> Plasmid C <br> Plasmid D | 4989207 <br> 66792 <br> 5263 <br> 4257 <br> 4104 | Yes <br> Yes <br> Yes <br> Yes <br> Yes |  
|S129EC| Chromosome <br> Plasmid A <br> Plasmid B <br> Plasmid C <br> Plasmid D <br> Plasmid E <br> Plasmid F <br> Plasmid G | 5193964 <br> 163681 <br> 93505 <br> 33344 <br> 4087 <br> 2401 <br> 2121 <br> 1571 | Yes <br> Yes <br> Yes <br> Yes <br> Yes <br> Yes <br> Yes <br> Yes |  
|EC958| Chromosome <br> Plasmid A <br> Plasmid B <br> Plasmid C | 5126816 <br> 136157 <br> 4145 <br> 1830 | Yes <br> Yes <br> Yes <br> Yes |    
|HVM2044| Chromosome <br> Plasmid A <br> Plasmid B <br> Plasmid C | 5003288 <br> 142959 <br> 18716 <br> 18345 | Yes <br> Yes <br> Yes <br> Yes |    


## Exemplar 2: Assembly of 12 *E.coli* ST131 samples using CPU resources 

* See Nextflow configuration file used [here](./nextflow.config) and slurm submission script [here](./nextflow_batch_template.sh). 

* See Nextflow [HTML execution report](./micropipe_ecoli_ST131_pawsey_guppy3.6.1_cpu.report.html), [trace report](./micropipe_ecoli_ST131_pawsey_guppy3.6.1_cpu.trace.txt) and [HTML processes execution timeline](./micropipe_ecoli_ST131_pawsey_guppy3.6.1_cpu.timeline.html). 

---

# Acknowledgements / citations / credits

- The deployment of the workflow at the Pawsey Supercomputing Centre was supported by the Australian BioCommons via funding from Bioplatforms Australia, the Australian Research Data Commons (https://doi.org/10.47486/PL105) and the Queensland Government RICF programme. Bioplatforms Australia and the Australian Research Data Commons are funded by the National Collaborative Research Infrastructure Strategy (NCRIS).
- Marco de la Pierre (Pawsey Supercomputing Centre) [@marcodelapierre](https://github.com/marcodelapierre)
- Johan Gustafsson (Australian BioCommons) [@supernord](https://github.com/supernord)
```
Any attribution information that is relevant to the workflow being documented, or the infrastructure being used.
```

---
