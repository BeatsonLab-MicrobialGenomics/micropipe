<p align="center"><img src="docs/micropipe_logo.png" alt="Logo" width="50%"></p>

**microPIPE: a pipeline for high-quality bacterial genome construction using ONT and Illumina sequencing**
======

microPIPE was developed to automate high-quality complete bacterial genome assembly using Oxford Nanopore Sequencing in combination with Illumina sequencing. 

To build microPIPE we evaluated the performance of several tools at each step of bacterial genome assembly, including basecalling, assembly, and polishing. Results at each step were validated using the high-quality ST131 *Escherichia coli* strain EC958 (GenBank: HG941718.1). After appraisal of each step, we selected the best combination of tools to achieve the most consistent and best quality bacterial genome assemblies.

The workflow below summarises the different steps of the pipeline (with each selected tool) and the approximate run time (using GPU basecalling, averaged over 12 *E. coli* isolates sequenced on a R9.4 MinION flow cell). Dashed boxes correspond to optional steps in the pipeline. 

Micropipe has been written in Nextflow and uses Singularity containers. It can use both GPU and CPU resources. 

For more information please see our preprint here: **URL for BioRxiv**

<p align="center">
  <img src="docs/Fig_workflow.png" alt="Workflow" width="400"/>
</p>

Please note that this pipeline does not perform extensive quality assessment of the input sequencing data. Contamination and sequencing read quality should be assessed independently to avoid problems with assembly. 

# Contents

* [Quickstart](#quickstart)
* [Installation](#installation)
* [Usage](#usage)
* [Test data](#example-data)
* [Optional parameters](#optional-parameters)
* [Structure of the output folders](#structure-of-the-output-folders)
* [Comments](#comments)


# Quickstart

1. Basecalling, demultiplexing and assembly workflow

`nextflow main.nf --basecalling --demultiplexing --samplesheet /path/to/samples.csv --fast5 /path/to/fast5/directory/ --outdir /path/to/outdir/`

2. Demultiplexing and assembly workflow (basecalling already complete)

`nextflow main.nf --demultiplexing --samplesheet /path/to/samples.csv --fastq /path/to/fastq/directory/ --outdir /path/to/outdir/`

3. Assembly only workflow (basecalling and demultiplexing already complete)

`nextflow main.nf --samplesheet /path/to/samples.csv --outdir /path/to/outdir/`


# Installation

microPIPE has been built using Nextflow and Singularity to enable ease of use and installation across different platforms. 

**0. Requirements**

* [Nextflow](https://www.nextflow.io/) >= 20.10.0

Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and Java 8 (or later, up to 15) to be installed.

To install Nextflow, run the command:

`wget -qO- https://get.nextflow.io | bash` or `curl -s https://get.nextflow.io | bash`

It will create the nextflow main executable file in the current directory. Optionally, move the nextflow file to a directory accessible by your $PATH variable. 

* [Singularity](https://singularity.lbl.gov/install-linux) >= 2.3 (microPIPE has been tested with version 3.4.1, 3.5.0 and 3.6.3)


**1. Installing microPIPE**

microPIPE only requires the `main.nf` and `nexflow.config` files to run. You will also need to provide a samplesheet (explained below). 

# Usage

**1. Prepare the Nextflow configuration file**

When a Nexflow pipeline script is launched, Nextflow looks for a file named **nextflow.config** in the current directory. The configuration file defines default parameters values for the pipeline and cluster settings such as the executor (e.g. "slurm", "local") and queues to be used (https://www.nextflow.io/docs/latest/config.html). 

The pipeline uses separated Singularity containers for all processes. Nextflow will automatically pull the singularity images required to run the pipeline and cache those images in the singularity directory in the pipeline work directory by default or in the singularity.cacheDir specified in the [nextflow.config](https://www.nextflow.io/docs/latest/singularity.html) file: 

```
singularity {
  enabled = true
  autoMounts = true
  singularity.cacheDir = '/path/to/cachedir'
}
```

An example configuration file can be found in this [repository](https://github.com/BeatsonLab-MicrobialGenomics/micropipe/blob/main/nextflow.config). 

Two versions of the configuration file are available and corresponds to MicroPIPE v0.8 (utilizing Guppy v3.4.3) and v0.9 (utilizing Guppy v3.6.1) as referenced in the paper.  

**NOTE:** to use **GPU** resources for basecalling and demultiplexing, use the `--gpu` flag.


**2. Prepare the samplesheet file (csv)**

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

**3. Run the pipeline**

The pipeline can be used to run:

* **Basecalling, demultiplexing and assembly workflow**

The entire workflow from basecalling to polishing will be run. The input files will be the ONT fast5 files and the Illumina fastq files. 

`nextflow main.nf --basecalling --demultiplexing --samplesheet /path/to/samples.csv --fast5 /path/to/fast5/directory/ --outdir /path/to/outdir/`
```
--samplesheet: samplesheet file
--outdir: name of the output directory to be created
--basecalling: flag to run the basecalling step 
--demultiplexing: flag to run the demultiplexing step 
--fast5: directory containing the ONT fast5 files
--guppy_config_gpu: Guppy configuration file name for basecalling using GPU resources (default=dna_r9.4.1_450bps_hac.cfg suitable if the Flow Cell Type = FLO-MIN106 and Kit = SQK-RBK004)
--guppy_config_cpu: Guppy configuration file name for basecalling using CPU resources (default=dna_r9.4.1_450bps_fast.cfg)
--medaka_model: Medaka model (default=r941_min_high, Available models: r941_min_fast, r941_min_high, r941_prom_fast, r941_prom_high, r10_min_high, r941_min_diploid_snp), see [details](https://github.com/nanoporetech/medaka#models)
```
Example of samplesheet file: 
```
barcode_id,sample_id,short_fastq_1,short_fastq_2,genome_size
barcode01,S24,S24EC.filtered_1P.fastq.gz,S24EC.filtered_2P.fastq.gz,5.5m
barcode02,S34,S34EC.filtered_1P.fastq.gz,S34EC.filtered_2P.fastq.gz,5.5m
```

* **Basecalling and assembly workflow (single isolate)**

The entire workflow from basecalling to polishing will be run (excluding demultiplexing). The input files will be the ONT fast5 files and the Illumina fastq files. 

`nextflow main.nf --basecalling --samplesheet /path/to/samples.csv --fast5 /path/to/fast5/directory/ --outdir /path/to/outdir/`
```
--samplesheet: path to the samplesheet file
--basecalling: flag to run the basecalling step
--fast5: path to the directory containing the ONT fast5 files
--outdir: path to the output directory to be created
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

`nextflow main.nf --demultiplexing --samplesheet /path/to/samples.csv --fastq /path/to/fastq/directory/ --outdir /path/to/outdir/`
```
--samplesheet: path to the samplesheet file
--demultiplexing: flag to run the demultiplexing step
--fastq: path to the directory containing the ONT fastq files (gzip compressed)
--outdir: path to the output directory to be created
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

`nextflow main.nf --samplesheet /path/to/samples.csv --outdir /path/to/outdir/`
```
--samplesheet: path to the samplesheet file
--outdir: path to the output directory to be created
--medaka_model: name of the Medaka model (default=r941_min_high, Available models: r941_min_fast, r941_min_high, r941_prom_fast, r941_prom_high, r10_min_high, r941_min_diploid_snp), see [details](https://github.com/nanoporetech/medaka#models)
```
Example of samplesheet file: 
```
barcode_id,sample_id,long_fastq,short_fastq_1,short_fastq_2,genome_size
barcode01,S24,barcode01.fastq.gz,S24EC.filtered_1P.fastq.gz,S24EC.filtered_2P.fastq.gz,5.5m
barcode02,S34,barcode02.fastq.gz,S34EC.filtered_1P.fastq.gz,S34EC.filtered_2P.fastq.gz,5.5m
```

# Example data

To test the pipeline, we have provided some [test data](https://github.com/BeatsonLab-MicrobialGenomics/micropipe/tree/main/test_data). In this directory you will find:

File | Description
---|---
S24EC_1P_test.fastq.gz | Illumina reads 1st pair
S24EC_2P_test.fastq.gz | Illumina reads 2nd pair
barcode01.fastq.gz | ONT fastq reads 
samples_1.csv | sample sheet for running assembly-only pipeline
samples_1_basecalling.csv | sample sheet for full pipeline
samples_1_basecalling_single_isolate.csv | sample sheet for a single isolate

To test the assembly-only pipeline, edit the `sample_1.csv` samplesheet to point to the correct test files. Then run:

`nextflow main.nf --samplesheet /path/to/samples_1.csv --outdir /path/to/test_outdir/`


# Optional parameters

Some parameters can be added to the command line in order to include or skip some steps and modify some parameters:

Basecalling 
* `--gpu`: use the GPU node to run the Guppy basecalling and/or demultiplexing step (default=false) 
* `--guppy_basecaller_args`: Guppy basecaller parameters (default="--recursive --trim_barcodes -q 0")
* `--guppy_num_callers`: number of parallel basecallers to create when running guppy basecalling (default=8)
* `--guppy_cpu_threads_per_caller`: number of CPU worker threads per basecaller (default=1). The number of CPU threads (num_callers * cpu_threads_per_caller ) used should generally not exceed the number of logical CPU cores your machine has.
*	`--guppy_gpu_device`:	Basecalling device for Guppy: "auto" or "cuda:<device_id>" (default="auto")
*	`--flowcell`:	Name of the ONT flow cell used for sequencing (default=false). Ignored if '--guppy_config_gpu' or '--guppy_congif_cpu' is specified
*	`--kit`: Name of the ONT kit used for sequencing (default=false). Ignored if '--guppy_config_gpu' or '--guppy_congif_cpu' is specified

Quality control:
* `--skip_pycoqc`: skip the pycoQC step to generate a quality control html report (when --basecalling)

Demultiplexing:
* `--demultiplexer`: demultiplexing tool: "qcat" or "guppy" (default="qcat")
* `--qcat_args`: qcat optional parameters (default="")
* `--guppy_barcoder_args`: Guppy barcoder parameters (default="--recursive --trim_barcodes -q 0")
* `--guppy_barcode_kits`: Space separated list of barcoding kit(s) to detect against (default="SQK-RBK004")
* `--guppy_barcoder_threads`: number of worker threads to spawn for Guppy barcoder to use. Increasing this number will allow Guppy barcoder to make better use of multi-core CPU systems, but may impact overall system performance (default=2)

Adapter trimming:
* `--skip_porechop` : skip the Porechop trimming step 
* `--porechop_threads`: number of threads for Porechop (default=4)
* `--porechop_args`: Porechop optional parameters (default=""), see [details](https://github.com/rrwick/Porechop#full-usage)

Filtering:
* `--skip_filtering` : skip the filtering step 
* `--filtering`: filtering tool: "japsa" or "filtlong" (default="japsa")
* `--japsa_args`: Japsa optional parameters (default="--lenMin 1000 --qualMin 10"), see [details](https://japsa.readthedocs.io/en/latest/tools/jsa.np.filter.html)
* `--filtlong_args`: Filtlong optional parameters (default="--min_length 1000 --keep_percent 90"), see [details](https://github.com/rrwick/Filtlong#full-usage)
* `--skip_rasusa`: Skip the sub-sampling Rasusa step
* `--rasusa_coverage`: The desired coverage to sub-sample the reads to (default=100), see [details](https://github.com/mbhall88/rasusa#-c---coverage)

Assembly:
* `--flye_args`: Flye optional parameters (default="--plasmids")
* `--flye_threads`: number of threads for Flye (default=4)

Polishing:
* `--polisher`: Long-read polishing tool: "medaka" (racon followed by medaka) or "nextpolish" (default="medaka")
* `--racon_nb`: number of Racon long-read polishing iterations (default=4)
* `--racon_args`: Racon optional parameters (default="-m 8 -x -6 -g -8 -w 500")
* `--racon_threads`: number of threads for Racon (default=4)
* `--medaka_threads`: number of threads for Medaka (default=4)
* `--skip_illumina`: skip the short-read polishing step if Illumina reads are not available (not recommended)
* `--nextpolish_threads`: number of threads for Nextpolish (default=4)
* `--skip_fixstart`: skip the Circlator fixstart step, see [details](https://github.com/sanger-pathogens/circlator/wiki/Task:-fixstart)
* `--fixstart_args`: Circlator fixstart optional parameters (default=""). Example "--genes_fa /path/to/fasta".

Assembly evaluation:
* `--skip_quast`: skip the QUAST assembly assessment step
* `--quast_args`: QUAST optional parameters (default=""), see [details](http://quast.sourceforge.net/docs/manual.html#sec2.3)
* `--quast_threads`: number of threads for QUAST (default=1)

# Structure of the output folders

The pipeline will create several folders corresponding to the different steps of the pipeline. 
The main output folder (`--outdir`) will contain the following folders:
* **0_basecalling:** Fastq files containing the basecalled reads, Guppy sequencing_summary.txt file
* **0_pycoQC:** Quality control report (pycoQC.html, see for [example]( https://a-slide.github.io/pycoQC/pycoQC/results/Guppy-2.1.3_basecall-1D_RNA.html))
* **a folder per sample:** see content below (the folder is named as in the column sample_id in the samplesheet file)

Each sample folder will contain the following folders:
* **1_filtering:** Fastq files containing trimmed reads (sample_id_trimmed.fastq.gz) and filtered reads (sample_id_filtered.fastq.gz) 
* **2_assembly:** Flye assembly output files (.fasta, .gfa, .gv, .info.txt), see [details](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-flye-output)
* **3_polishing_long_reads:** Long-read polished assembly fasta file (sample_id_flye_polishedLR.fasta)
* **4_polishing_short_reads:** Final polished assembly fasta file (sample_id_flye_polishedLR_SR.fasta)
* **5_quast:** QUAST quality assessment report, see [details](http://quast.sourceforge.net/docs/manual.html)

# Comments

The pipeline has been tested using the following grid based executors: SLURM, PBS Pro and LSF.  

Planned upgrades:
- Enabling GPU resource for Racon and Medaka processes.

