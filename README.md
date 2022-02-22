<p align="center"><img src="docs/micropipe_logo.png" alt="Logo" width="45%"></p>

microPIPE: a pipeline for high-quality bacterial genome construction using ONT and Illumina sequencing
================

  - [Description](#description)
  - [Diagram](#diagram)
  - [User guide](#user-guide)
      - [Quick start guide](#quick-start-guide)
      - [Step by step user guide](#step-by-step-user-guide)
      - [Optional parameters](#optional-parameters)
      - [Structure of the output folders](#structure-of-the-output-folders)
      - [Test Data](#example-data)
      - [Infrastructure usage and recommendations](#infrastructure-usage-and-recommendations)
      - [Compute resource usage across tested infrastructures](#compute-resource-usage-across-tested-infrastructures)
  - [Benchmarking](#benchmarking)
  - [Workflow summaries](#workflow-summaries)
      - [Metadata](#metadata)
      - [Component tools](#component-tools)
      - [Third party tools /dependencies](#third-party-tools--dependencies)
  - [Additional notes](#additional-notes)
  - [Help/FAQ/Troubleshooting](#help--faq--troubleshooting)
  - [Licence(s)](#licences)
  - [Acknowledgements/citations/credits](#acknowledgements--citations--credits)

# Description

microPIPE was developed to automate high-quality complete bacterial genome assembly using Oxford Nanopore Sequencing in combination with Illumina sequencing. 

To build microPIPE we evaluated the performance of several tools at each step of bacterial genome assembly, including basecalling, assembly, and polishing. Results at each step were validated using the high-quality ST131 *Escherichia coli* strain EC958 (GenBank: HG941718.1). After appraisal of each step, we selected the best combination of tools to achieve the most consistent and best quality bacterial genome assemblies.

Please note that this pipeline does not perform extensive quality assessment of the input sequencing data. Contamination and sequencing read quality should be assessed independently to avoid problems with assembly.

Micropipe has been written in Nextflow and uses Singularity containers. It can use both GPU and CPU resources. 

For more information please see our publication here: https://doi.org/10.1186/s12864-021-07767-z.

-----

# Diagram

The diagram below summarises the different steps of the pipeline (with each selected tool) and the approximate run time (using GPU basecalling, averaged over 12 E. coli isolates sequenced on a R9.4 MinION flow cell). Dashed boxes correspond to optional steps in the pipeline.
 <p align="center">
 <img src="https://github.com/BeatsonLab-MicrobialGenomics/micropipe/blob/main/docs/Fig_workflow.png" alt="Workflow" width="400"/>
 </p>

-----

# User guide

## Quick start guide

1. Basecalling, demultiplexing and assembly workflow

`nextflow main.nf --basecalling --demultiplexing --samplesheet /path/to/samples.csv --fast5 /path/to/fast5/directory/ --datadir /path/to/datadir/ --outdir /path/to/outdir/`

2. Demultiplexing and assembly workflow (basecalling already complete)

`nextflow main.nf --demultiplexing --samplesheet /path/to/samples.csv --fastq /path/to/fastq/directory/ --datadir /path/to/datadir/ --outdir /path/to/outdir/`

3. Assembly only workflow (basecalling and demultiplexing already complete)

`nextflow main.nf --samplesheet /path/to/samples.csv --fastq /path/to/fastq/directory/ --datadir /path/to/datadir/ --outdir /path/to/outdir/`

An infrastructure specific guide for Zeus @ [Pawsey Supercomputing Centre](https://pawsey.org.au/) (Perth, Western Australia) is provided [here](https://github.com/vmurigneu/micropipe_pawsey/blob/master/docs/infrastructure_optimisation_zeus.md#quickstart-tutorial).

## Step by step user guide

**0. Requirements**

microPIPE has been built using Nextflow and Singularity to enable ease of use and installation across different platforms.

* [Nextflow](https://www.nextflow.io/) >= 20.10.0

Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and Java 8 (or later, up to 15) to be installed.

To install Nextflow, run the command:

`wget -qO- https://get.nextflow.io | bash` or `curl -s https://get.nextflow.io | bash`

It will create the nextflow main executable file in the current directory. Optionally, move the nextflow file to a directory accessible by your $PATH variable. 

* [Singularity](https://singularity.lbl.gov/install-linux) >= 2.3 (microPIPE has been tested with version 3.4.1, 3.5.0 and 3.6.3)

* Guppy (4.4.1 was the latest working version)
 
Due to the Oxford Nanopore Technologies terms and conditions, we are not allowed to redistribute the Guppy software either in its binary form or packaged form e.g. Docker or Singularity images. Therefore users will have to either install Guppy, provide a container image or start the pipeline from the basecalled fastq files. 

**1. Installing microPIPE**

Download the microPIPE repository using the command:
``` 
git clone https://github.com/BeatsonLab-MicrobialGenomics/micropipe.git
```
microPIPE requires the files `main.nf`, `nexflow.config` and a samplesheet file to run.

**2. Prepare the Nextflow configuration file**

When a Nexflow pipeline script is launched, Nextflow looks for a file named **nextflow.config** in the current directory. The configuration file defines default parameters values for the pipeline and cluster settings such as the executor (e.g. "slurm", "local") and queues to be used (https://www.nextflow.io/docs/latest/config.html). 

The pipeline uses separated [Singularity containers](https://www.nextflow.io/docs/latest/singularity.html) for all processes. Nextflow will automatically pull the singularity images required to run the pipeline and cache those images in the singularity directory in the pipeline work directory by default or in the singularity.cacheDir specified in the [nextflow.config](./nextflow.config) file: 

```
singularity {
  enabled = true
  singularity.cacheDir = '/path/to/cachedir'
}
```

The **nextflow.config** file should be modified to specify the location of Guppy using one of the following options:

* Download and unpack the Guppy .tar.gz archive file. Provide the path to the Guppy binary folder in the params section and comment the following lines in the process section: 
   ```
   params {
   //Path to the Guppy GPU and CPU binary folder. Change this as appropriate when providing Guppy as a binary folder and do not forget the "/" at the end of the path
   guppy_gpu_folder= "/scratch/ont-guppy/bin/"
   guppy_cpu_folder = "/scratch/ont-guppy-cpu/bin/"
   }
   ```
   ```
   process {
   //Path to the Guppy GPU and CPU container images. Uncomment and change this as appropriate if providing Guppy as a container image.
   //withLabel: guppy_gpu { container = '' }
   //withLabel: guppy_cpu { container = '' }
   }
   ```
 
* Provide the link to a Guppy container in the process section and uncomment the two following lines in the params section:
   ```  
   params {
   //Uncomment the two following lines when providing Guppy container images (and comment the two previous lines)
   guppy_gpu_folder = ""
   guppy_cpu_folder = ""
   }
   ```
   ```
   process { 
   //Path to the Guppy GPU and CPU container images. Uncomment and change this as appropriate if providing Guppy as a container image. 
   withLabel: guppy_gpu { container = '' }
   withLabel: guppy_cpu { container = '' }
   }
   ```

An example configuration file can be found in this [repository](https://github.com/BeatsonLab-MicrobialGenomics/micropipe/blob/main/nextflow.config). 

Two versions of the configuration file are available and correspond to microPIPE v0.8 (utilizing Guppy v3.4.3) and v0.9 (utilizing Guppy v3.6.1), as referenced in the [paper](https://doi.org/10.1186/s12864-021-07767-z). 

**3. Prepare the samplesheet file (csv)**

The samplesheet file (comma-separated values) defines the input fastq files (Illumina [short] and Nanopore [long], and their directory path), barcode number, sample ID, and the estimated genome size (for Flye assembly). The header line should match the header line in the examples below:

1. If using demultiplexing:

```
barcode_id,sample_id,short_fastq_1,short_fastq_2,genome_size
barcode01,S24,S24EC.filtered_1P.fastq.gz,S24EC.filtered_2P.fastq.gz,5.5m
barcode02,S34,S34EC.filtered_1P.fastq.gz,S34EC.filtered_2P.fastq.gz,5.5m

```
2. If not using demultiplexing (single isolate):

```
sample_id,short_fastq_1,short_fastq_2,genome_size
S24,S24EC.filtered_1P.fastq.gz,S24EC.filtered_2P.fastq.gz,5.5m
```

3. If using assembly only:

```
barcode_id,sample_id,long_fastq,short_fastq_1,short_fastq_2,genome_size
barcode01,S24,barcode01.fastq.gz,S24EC.filtered_1P.fastq.gz,S24EC.filtered_2P.fastq.gz,5.5m
barcode02,S34,barcode02.fastq.gz,S34EC.filtered_1P.fastq.gz,S34EC.filtered_2P.fastq.gz,5.5m
```

4. If Illumina reads are not available (--skip_illumina), do not include the two columns with the Illumina files:

```
barcode_id,sample_id,long_fastq,genome_size
barcode01,S24,barcode01.fastq.gz,5.5m
barcode02,S34,barcode02.fastq.gz,5.5m
```

**4. Run the pipeline**

The pipeline can be used to run:

* **Basecalling, demultiplexing and assembly workflow**

The entire workflow from basecalling to polishing will be run. The input files will be the ONT fast5 files and the Illumina fastq files. 

`nextflow main.nf --basecalling --demultiplexing --samplesheet /path/to/samples.csv --fast5 /path/to/fast5/directory/ --datadir /path/to/datadir/ --outdir /path/to/outdir/`
```
--samplesheet: samplesheet file
--basecalling: flag to run the basecalling step 
--demultiplexing: flag to run the demultiplexing step 
--fast5: directory containing the ONT fast5 files
--outdir: path to the output directory to be created
--datadir: path to the directory containing the Illumina fastq files
--guppy_config_gpu: Guppy configuration file name for basecalling using GPU resources (default=dna_r9.4.1_450bps_hac.cfg suitable if the Flow Cell Type = FLO-MIN106 and Kit = SQK-RBK004)
--guppy_config_cpu: Guppy configuration file name for basecalling using CPU resources (default=dna_r9.4.1_450bps_fast.cfg)
--medaka_model: Medaka model (default=r941_min_high, Available models: r941_min_fast, r941_min_high, r941_prom_fast, r941_prom_high, r10_min_high, r941_min_diploid_snp), see [details](https://github.com/nanoporetech/medaka#models)
```
**NOTE:** to use **GPU** resources for basecalling and demultiplexing, use the `--gpu` flag.

Example of samplesheet file: 
```
barcode_id,sample_id,short_fastq_1,short_fastq_2,genome_size
barcode01,S24,S24EC.filtered_1P.fastq.gz,S24EC.filtered_2P.fastq.gz,5.5m
barcode02,S34,S34EC.filtered_1P.fastq.gz,S34EC.filtered_2P.fastq.gz,5.5m
```

* **Basecalling and assembly workflow (single isolate)**

The entire workflow from basecalling to polishing will be run (excluding demultiplexing). The input files will be the ONT fast5 files and the Illumina fastq files. 

`nextflow main.nf --basecalling --samplesheet /path/to/samples.csv --fast5 /path/to/fast5/directory/ --datadir /path/to/datadir/ --outdir /path/to/outdir/`
```
--samplesheet: path to the samplesheet file
--basecalling: flag to run the basecalling step
--fast5: path to the directory containing the ONT fast5 files
--outdir: path to the output directory to be created
--datadir: path to the directory containing the Illumina fastq files
--guppy_config_gpu: Guppy configuration file name for basecalling using GPU resources (default=dna_r9.4.1_450bps_hac.cfg suitable if the Flow Cell Type = FLO-MIN106 and Kit = SQK-LSK109)
--guppy_config_cpu: Guppy configuration file name for basecalling using CPU resources (default=dna_r9.4.1_450bps_fast.cfg)
--medaka_model: name of the Medaka model (default=r941_min_high, Available models: r941_min_fast, r941_min_high, r941_prom_fast, r941_prom_high, r10_min_high, r941_min_diploid_snp), see [details](https://github.com/nanoporetech/medaka#models)
```
Example of samplesheet file (containing only one sample): 
```
sample_id,short_fastq_1,short_fastq_2,genome_size
S24,S24EC.filtered_1P.fastq.gz,S24EC.filtered_2P.fastq.gz,5.5m
```

* **Demultiplexing and assembly workflow**

The entire workflow from demultiplexing to polishing will be run. The input files will be the ONT fastq files and the Illumina fastq files. 

`nextflow main.nf --demultiplexing --samplesheet /path/to/samples.csv --fastq /path/to/fastq/directory/ --datadir /path/to/datadir/ --outdir /path/to/outdir/`
```
--samplesheet: path to the samplesheet file
--demultiplexing: flag to run the demultiplexing step
--fastq: path to the directory containing the ONT fastq files (gzip compressed)
--outdir: path to the output directory to be created
--datadir: path to the directory containing the Illumina fastq files
--guppy_config_gpu: Guppy configuration file name for basecalling using GPU resources (default=dna_r9.4.1_450bps_hac.cfg suitable if the Flow Cell Type = FLO-MIN106 and Kit = SQK-LSK109)
--guppy_config_cpu: Guppy configuration file name for basecalling using CPU resources (default=dna_r9.4.1_450bps_fast.cfg)
--medaka_model: name of the Medaka model (default=r941_min_high, available models: r941_min_fast, r941_min_high, r941_prom_fast, r941_prom_high, r10_min_high, r941_min_diploid_snp), see [details](https://github.com/nanoporetech/medaka#models)
```
Example of samplesheet file: 
```
barcode_id,sample_id,short_fastq_1,short_fastq_2,genome_size
barcode01,S24,S24EC.filtered_1P.fastq.gz,S24EC.filtered_2P.fastq.gz,5.5m
barcode02,S34,S34EC.filtered_1P.fastq.gz,S34EC.filtered_2P.fastq.gz,5.5m
```

* **Assembly only workflow**

The assembly workflow from adapter trimming to polishing will be run. The input files will be the ONT fastq files and the Illumina fastq files. 

`nextflow main.nf --samplesheet /path/to/samples.csv --fastq /path/to/fastq/directory/ --datadir /path/to/datadir/ --outdir /path/to/outdir/`
```
--samplesheet: path to the samplesheet file
--fastq: path to the directory containing the ONT fastq files (gzip compressed)
--outdir: path to the output directory to be created
--datadir: path to the directory containing the Illumina fastq files
--medaka_model: name of the Medaka model (default=r941_min_high, Available models: r941_min_fast, r941_min_high, r941_prom_fast, r941_prom_high, r10_min_high, r941_min_diploid_snp), see [details](https://github.com/nanoporetech/medaka#models)
```
Example of samplesheet file: 
```
barcode_id,sample_id,long_fastq,short_fastq_1,short_fastq_2,genome_size
barcode01,S24,barcode01.fastq.gz,S24EC.filtered_1P.fastq.gz,S24EC.filtered_2P.fastq.gz,5.5m
barcode02,S34,barcode02.fastq.gz,S34EC.filtered_1P.fastq.gz,S34EC.filtered_2P.fastq.gz,5.5m
```

## Optional parameters

Some parameters can be added to the command line in order to include or skip some steps and modify some parameters:

Basecalling 
* `--gpu`: use the GPU node to run the Guppy basecalling and/or demultiplexing step (default=true) 
* `--guppy_basecaller_args`: Guppy basecaller parameters (default="--recursive --trim_barcodes -q 0")
* `--guppy_num_callers`: number of parallel basecallers to create when running guppy basecalling (default=8)
* `--guppy_cpu_threads_per_caller`: number of CPU worker threads per basecaller (default=1). The number of CPU threads (num_callers * cpu_threads_per_caller ) used should generally not exceed the number of logical CPU cores your machine has.
*	`--guppy_gpu_device`:	Basecalling device for Guppy: "auto" or "cuda:<device_id>" (default="auto")
*	`--flowcell`:	Name of the ONT flow cell used for sequencing (default=false). Ignored if '--guppy_config_gpu' or '--guppy_congif_cpu' is specified
*	`--kit`: Name of the ONT kit used for sequencing (default=false). Ignored if '--guppy_config_gpu' or '--guppy_congif_cpu' is specified

Quality control:
* `--skip_pycoqc`: skip the pycoQC step to generate a quality control html report, when --basecalling (default=false)

Demultiplexing:
* `--demultiplexer`: demultiplexing tool: "qcat" or "guppy" (default=`--demultiplexer "guppy"`)
* `--qcat_args`: qcat optional parameters (default=""), see [details](https://github.com/nanoporetech/qcat#full-usage)
* `--guppy_barcoder_args`: Guppy barcoder parameters (default="--recursive --trim_barcodes -q 0")
* `--guppy_barcode_kits`: Space separated list of barcoding kit(s) to detect against (default="SQK-RBK004")
* `--guppy_barcoder_threads`: number of worker threads to spawn for Guppy barcoder to use. Increasing this number will allow Guppy barcoder to make better use of multi-core CPU systems, but may impact overall system performance (default=2)

Adapter trimming:
* `--skip_porechop`: skip the Porechop trimming step (default=false)
* `--porechop_threads`: number of threads for Porechop (default=4)
* `--porechop_args`: Porechop optional parameters (default=""), see [details](https://github.com/rrwick/Porechop#full-usage)

Filtering:
* `--skip_filtering`: skip the filtering step (default=false)
* `--filtering`: filtering tool: "japsa" or "filtlong" (default=`--filtering "japsa"`)
* `--japsa_args`: Japsa optional parameters (default="--lenMin 1000 --qualMin 10"), see [details](https://japsa.readthedocs.io/en/latest/tools/jsa.np.filter.html)
* `--filtlong_args`: Filtlong optional parameters (default="--min_length 1000 --keep_percent 90"), see [details](https://github.com/rrwick/Filtlong#full-usage)
* `--skip_rasusa`: Skip the sub-sampling Rasusa step (default=true)
* `--rasusa_coverage`: The desired coverage to sub-sample the reads to (default=100), see [details](https://github.com/mbhall88/rasusa#-c---coverage)

Assembly:
* `--flye_args`: Flye optional parameters (default=`--flye_args "--plasmids"`), see [details](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md)
* `--flye_threads`: number of threads for Flye (default=4)

Polishing:
* `--polisher`: Long-read polishing tool: "medaka" (racon followed by medaka) or "nextpolish" (default="medaka")
* `--racon_nb`: number of Racon long-read polishing iterations (default=4)
* `--racon_args`: Racon optional parameters (default="-m 8 -x -6 -g -8 -w 500"), see [details](https://github.com/isovic/racon#usage)
* `--racon_threads`: number of threads for Racon (default=4)
* `--medaka_threads`: number of threads for Medaka (default=4)
* `--skip_illumina`: skip the short-read polishing step if Illumina reads are not available (not recommended, default=false)
* `--nextpolish_threads`: number of threads for nextPolish (default=4)
* `--nextpolish_task_SR`: task to run for nextPolish short-read polishing ("12" or "1212", default="1212"), see [details](https://nextpolish.readthedocs.io/en/latest/OPTION.html#cmdoption-arg-task)
* `--nextpolish_task_LR`: task to run for nextPolish long-read polishing ("5" or "55",  default="55"), see [details](https://nextpolish.readthedocs.io/en/latest/OPTION.html#cmdoption-arg-task)
* `--skip_fixstart`: skip the Circlator fixstart step (default=false), see [details](https://github.com/sanger-pathogens/circlator/wiki/Task:-fixstart)
* `--fixstart_args`: Circlator fixstart optional parameters (default=""). Example `--fixstart_args "--genes_fa /path/to/datadir/fasta"` (the file should be located in the nextflow launch directory or in the datadir).

Assembly evaluation:
* `--skip_quast`: skip the QUAST assembly assessment step (default=false)
* `--quast_args`: QUAST optional parameters (default=""), see [details](http://quast.sourceforge.net/docs/manual.html#sec2.3). Example: `--quast_args "-r /path/to/datadir/fasta"` (the file should be located in the nextflow launch directory or in the datadir).
* `--quast_threads`: number of threads for QUAST (default=1)

## Structure of the output folders

The pipeline will create several folders corresponding to the different steps of the pipeline. 
The main output folder (`--outdir`) will contain the following folders:
* **0_basecalling:** Fastq files containing the basecalled reads, Guppy sequencing_summary.txt file
* **0_pycoQC:** Quality control report (pycoQC.html, see for [example]( https://a-slide.github.io/pycoQC/pycoQC/results/Guppy-2.1.3_basecall-1D_RNA.html))
* **a folder per sample:** see content below (the folder is named as in the column sample_id in the samplesheet file)

Each sample folder will contain the following folders:
* **1_filtering:** Fastq files containing filtered reads (sample_id_filtered.fastq.gz) 
* **2_assembly:** Flye assembly output files (.fasta, .gfa, .gv, .info.txt), see [details](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-flye-output)
* **3_polishing_long_reads:** Long-read polished assembly fasta file (sample_id_flye_polishedLR.fasta)
* **4_polishing_short_reads:** Final polished assembly fasta file (sample_id_flye_polishedLR_SR_fixstart.fasta)
* **5_quast:** QUAST quality assessment report, see [details](http://quast.sourceforge.net/docs/manual.html)

## Example data

To test the pipeline, we have provided some [test data](https://github.com/BeatsonLab-MicrobialGenomics/micropipe/tree/main/test_data). In this directory you will find:

File | Description
---|---
S24EC_1P_test.fastq.gz | Illumina reads 1st pair
S24EC_2P_test.fastq.gz | Illumina reads 2nd pair
barcode01.fastq.gz | ONT fastq reads 
samples_1.csv | sample sheet for running assembly-only pipeline

To test the assembly-only pipeline, edit the `sample_1.csv` samplesheet to point to the correct test files. Then run:

`nextflow main.nf --samplesheet /path/to/samples_1.csv --outdir /path/to/test_outdir/`

## Infrastructure usage and recommendations

### General recommendations for using microPIPE

When using microPIPE to run the Oxford Nanopore data basecalling and demultiplexing, it is recommended to use the GPU resources. As a result, the basecalling step will be performed using the high accuracy model (instead of the fast model) and the workflow will complete faster than with only the CPU resources. 
 
To use GPU resources for basecalling and demultiplexing, use the `--gpu` flag in the main nextflow command:  
```
nextflow main.nf --gpu true --basecalling --demultiplexing --samplesheet /path/to/samples.csv --fast5 /path/to/fast5/directory/ --datadir /path/to/datadir/ --outdir /path/to/outdir/
```

-----

## Compute resource usage across tested infrastructures

The table below summarised the basecalling run time depending on the resources used at the Pawsey Supercomputing Centre. 

|Resources (Cluster)|Basecalling model|Guppy Configuration file|Run time|
|-------|:-----:|:-----:|:-----:|
|GPU (Pawsey Topaz)| high-accuracy | dna_r9.4.1_450bps_hac.cfg | 10h 17m 17s |    
|CPU (Pawsey Zeus)| fast | dna_r9.4.1_450bps_fast.cfg | 3d 19h 21m 31s |    

-----

# Benchmarking

## Summary

### Exemplar 1: Assembly of 12 *E.coli* ST131 samples using GPU and CPU resources @ Pawsey

* We used the *E.coli* data from the [microPIPE publication](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07767-z) available from the NCBI SRA [BioProject PRJNA679678](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA679678/) (Oxford Nanopore) and the [BioProject PRJEB2968](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB2968) (Illumina). 

* See Nextflow configuration file used [here](./pawsey/nextflow.config) and slurm submission script [here](./pawsey/nextflow_batch_template.sh). 
* See Nextflow [HTML execution report](./pawsey/micropipe_ecoli_ST131_pawsey_guppy3.6.1_gpu.report.html), [trace report](./pawsey/micropipe_ecoli_ST131_pawsey_guppy3.6.1_gpu.trace.txt) and [HTML processes execution timeline](./pawsey/micropipe_ecoli_ST131_pawsey_guppy3.6.1_gpu.timeline.html). 

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


### Exemplar 2: Assembly of 12 *E.coli* ST131 samples using CPU resources @ Pawsey

* See Nextflow configuration file used [here](./pawsey/nextflow.config) and slurm submission script [here](./nextflow_batch_template.sh). 

* See Nextflow [HTML execution report](./pawsey/micropipe_ecoli_ST131_pawsey_guppy3.6.1_cpu.report.html), [trace report](./pawsey/micropipe_ecoli_ST131_pawsey_guppy3.6.1_cpu.trace.txt) and [HTML processes execution timeline](./pawsey/micropipe_ecoli_ST131_pawsey_guppy3.6.1_cpu.timeline.html). 

-----

# Workflow summaries

## Metadata

| metadata field   | workflow\_name / workflow\_version |
| ---------------- | :--------------------------------: |
| Version          |         v0.9          |
| Maturity         |               stable               |
| Creators         |        Valentine Murigneux, Leah W Roberts,  Brian M Forde, Minh-Duy Phan, Nguyen Thi Khanh Nhu, Adam D Irwin, Patrick N A Harris, David L Paterson, Mark A Schembri, David M Whiley, Scott A Beatson          |
| Source           |                 https://github.com/BeatsonLab-MicrobialGenomics/micropipe                 |
| License          |                 https://github.com/BeatsonLab-MicrobialGenomics/micropipe/blob/main/LICENSE                 |
| Workflow manager |              NextFlow              |
| Container        |                Singularity                |
| Install method   |               Manual               |
| GitHub           |                 https://github.com/BeatsonLab-MicrobialGenomics/micropipe                 |
| bio.tools        |                 NA                 |
| BioContainers    |                 NA                 |
| bioconda         |                 NA                 |

-----

## Component tools

| Workflow element | Workflow element version | Workflow title |
| ---------------- | :----------------------: | :------------: |
| Guppy           |  v3.6.1   | microPIPE v0.9 |
| qcat           |   v1.0.1   | microPIPE v0.9 |
| rasusa           |   v0.3.0   | microPIPE v0.9 |
| pycoQC           |   v2.5.0.23   | microPIPE v0.9 |
| Porechop           |  v0.2.3  | microPIPE v0.9 |
| Filtlong           |  v0.2.0  | microPIPE v0.9 |
| Japsa           |  v1.9-01a   | microPIPE v0.9 |
| Flye           |   v2.5   | microPIPE v0.9 |
| Racon           |   v1.4.9   | microPIPE v0.9 |
| Medaka           |   v0.10.0 | microPIPE v0.9 |
| NextPolish           |   v1.1.0   | microPIPE v0.9 |
| Circlator           |   v1.5.5   | microPIPE v0.9 |
| QUAST           |   v5.0.2   | microPIPE v0.9 |

-----

## Third party tools / dependencies

* [Nextflow](https://www.nextflow.io/) >= 20.10.0

Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and Java 8 (or later, up to 15) to be installed.

To install Nextflow, run the command:

`wget -qO- https://get.nextflow.io | bash` or `curl -s https://get.nextflow.io | bash`

It will create the nextflow main executable file in the current directory. Optionally, move the nextflow file to a directory accessible by your $PATH variable. 

* [Singularity](https://singularity.lbl.gov/install-linux) >= 2.3 (microPIPE has been tested with version 3.4.1, 3.5.0 and 3.6.3)

* Guppy (4.4.1 was the latest working version)
 
Due to the Oxford Nanopore Technologies terms and conditions, we are not allowed to redistribute the Guppy software either in its binary form or packaged form e.g. Docker or Singularity images. Therefore users will have to either install Guppy, provide a container image or start the pipeline from the basecalled fastq files. 

-----

# Additional notes

* The pipeline has been tested using the following grid based executors: SLURM, PBS Pro and LSF.

* Do not forget to delete the /work directory created by Nextflow once the pipeline has completed.

* Planned upgrades:

    Enabling GPU resource for Racon and Medaka processes.

-----

# Help / FAQ / Troubleshooting
* In versions greater than Guppy v4.5.2, the default Guppy parameters have changed. If you wish to use Guppy > v4.5.2, please modify the `nexflow.config` to run Guppy with the "--disable_qscore_filtering" flag: 
``` 
params {
        guppy_basecaller_args = "--recursive --trim_barcodes -q 0 --disable_qscore_filtering"
}
```

-----

# Licence(s)
 https://github.com/BeatsonLab-MicrobialGenomics/micropipe/blob/main/LICENSE  
 
-----

# Acknowledgements / citations / credits

## Citations

- Murigneux, V., Roberts, L.W., Forde, B.M. et al. MicroPIPE: validating an end-to-end workflow for high-quality complete bacterial genome construction. BMC Genomics 22, 474 (2021). [https://doi.org/10.1186/s12864-021-07767-z](https://doi.org/10.1186/s12864-021-07767-z)
- Murigneux, V. (2021). microPIPE: a pipeline for high-quality bacterial genome construction using ONT and Illumina sequencing. WorkflowHub. [https://doi.org/10.48546/WORKFLOWHUB.WORKFLOW.140.1](https://doi.org/10.48546/WORKFLOWHUB.WORKFLOW.140.1)

## Acknowledgements

This work was supported by funding from the Queensland Genomics Health Alliance (now Queensland Genomics), Queensland Health, the Queensland Government.

The deployment of the workflow at the Pawsey Supercomputing Centre was supported by the Australian BioCommons via funding from Bioplatforms Australia, the Australian Research Data Commons (https://doi.org/10.47486/PL105) and the Queensland Government RICF programme. Bioplatforms Australia and the Australian Research Data Commons are funded by the National Collaborative Research Infrastructure Strategy (NCRIS).

-----

