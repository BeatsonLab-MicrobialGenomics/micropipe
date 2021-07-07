#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/*
========================================================================================
        microPIPE - Bacterial genome construction pipeline using ONT sequencing
========================================================================================
 #### Documentation
https://github.com/BeatsonLab-MicrobialGenomics/micropipe
 #### Authors
 Valentine Murigneux <v.murigneux@uq.edu.au>
========================================================================================
*/

def helpMessage() {
	log.info"""
	=========================================
	microPIPE v${workflow.manifest.version}
	=========================================
	Usage:
	Basecalling, demultiplexing and assembly workflow:
	nextflow main.nf --basecalling --demultiplexing --samplesheet /path/to/samples.csv --fast5 /path/to/fast5/directory/ --outdir /path/to/outdir/ --datadir /path/to/datadir/

	Basecalling and assembly workflow (single isolate):
	nextflow main.nf --basecalling --samplesheet /path/to/samples.csv --fast5 /path/to/fast5/directory/ --outdir /path/to/outdir/ --datadir /path/to/datadir/

	Demultiplexing and assembly workflow (basecalling already complete):
	nextflow main.nf --demultiplexing --samplesheet /path/to/samples.csv --fastq /path/to/fastq/directory/ --outdir /path/to/outdir/ --datadir /path/to/datadir/

	Assembly only workflow (basecalling and demultiplexing already complete):
	nextflow main.nf --samplesheet /path/to/samples.csv --fastq /path/to/fastq/directory/ --outdir /path/to/outdir/ --datadir /path/to/datadir/

	Required arguments:
		--samplesheet				Path to the samplesheet file
		--outdir				Path to the output directory
		--datadir				Path to the directory containing the Illumina fastq files
    
	Basecalling:
		--basecalling				Flag to run the basecalling step (default=false)
		--fast5					Path to the directory containing the ONT fast5 files
		--gpu					Use the GPU node to run the Guppy basecalling and/or demultiplexing step (default=false) 
		--guppy_basecaller_args			Guppy basecaller parameters (default="--recursive --trim_barcodes -q 0")
		--guppy_num_callers 			Number of parallel basecallers to create when running guppy basecalling (default=8)
		--guppy_cpu_threads_per_caller 		Number of CPU worker threads per basecaller (default=1). The number of CPU threads (num_callers * cpu_threads_per_caller ) used should generally not exceed the number of logical CPU cores your machine has.
		--guppy_gpu_device			Basecalling device for Guppy: "auto" or "cuda:<device_id>" (default="auto")
		--guppy_gpu_folder			Path to the Gupppy GPU binary folder (default="/scratch/ont-guppy/bin")
		--guppy_cpu_folder			Path to the Gupppy CPU binary folder (default="/scratch/ont-guppy-cpu/bin")
		--guppy_config_gpu			Guppy configuration file name for basecalling using GPU resources (default=dna_r9.4.1_450bps_hac.cfg suitable if the Flow Cell Type = FLO-MIN106 and Kit = SQK-RBK004)
		--guppy_config_cpu 			Guppy configuration file name for basecalling using CPU resources (default=dna_r9.4.1_450bps_fast.cfg)
		--flowcell            			Name of the ONT flow cell used for sequencing (default=false). Ignored if '--guppy_config_gpu' or '--guppy_congif_cpu' is specified
		--kit                 			Name of the ONT kit used for sequencing (default=false). Ignored if '--guppy_config_gpu' or '--guppy_congif_cpu' is specified

	Quality control:
		--skip_pycoqc				Skip the pycoQC step to generate a quality control html report (when --basecalling)

	Demultiplexing:
		--demultiplexing			Flag to run the demultiplexing (default=false)
		--fastq					Path to the directory containing the ONT fastq files (gzip compressed)
		--demultiplexer 			Demultiplexing tool: "qcat" or "guppy" (default="qcat")
		--qcat_args 				Qcat optional parameters (default="")
		--guppy_barcoder_args 			Guppy barcoder parameters (default="--recursive --trim_barcodes -q 0")
		--guppy_barcode_kits 			Space separated list of barcoding kit(s) to detect against (default="SQK-RBK004")
		--guppy_barcoder_threads 		Number of worker threads to spawn for Guppy barcoder to use. Increasing this number will allow Guppy barcoder to make better use of multi-core CPU systems, but may impact overall system performance (default=2)

	Adapter trimming:
		--skip_porechop 			Skip the Porechop trimming step 
		--porechop_threads			Number of threads for Porechop (default=4)
		--porechop_args 			Porechop optional parameters (default=""), see https://github.com/rrwick/Porechop#full-usage

	Filtering:
		--skip_filtering 			Skip the filtering step 
		--filtering				Filtering tool: "japsa" or "filtlong" (default="japsa")
		--japsa_args				Japsa optional parameters (default="--lenMin 1000 --qualMin 10"), see https://japsa.readthedocs.io/en/latest/tools/jsa.np.filter.html
		--filtlong_args 			Filtlong optional parameters (default="--min_length 1000 --keep_percent 90"), see https://github.com/rrwick/Filtlong#full-usage
		--skip_rasusa				Skip the sub-sampling Rasusa step
		--rasusa_coverage			The desired coverage to sub-sample the reads to (default=100)

	Assembly:
		--flye_args				Flye optional parameters (default="--plasmids")
		--flye_threads 				Number of threads for Flye (default=4)

	Polishing:
		--polisher				Long-read polishing tool: "medaka" (racon followed by medaka) or "nextpolish" (default="medaka")
		--racon_nb				Number of Racon long-read polishing iterations (default=4)
		--racon_args 				Racon optional parameters (default="-m 8 -x -6 -g -8 -w 500")
		--racon_threads 			Number of threads for Racon (default=4)
		--medaka_threads 			Number of threads for Medaka (default=4)
		--medaka_model				Medaka model (default=r941_min_high, Available models: r941_min_fast, r941_min_high, r941_prom_fast, r941_prom_high, r10_min_high, r941_min_diploid_snp), see https://github.com/nanoporetech/medaka#models
		--skip_illumina 			Skip the short-read polishing step if Illumina reads are not available (not recommended)
		--nextpolish_threads 			Number of threads for Nextpolish (default=4)
		--skip_fixstart 			Skip the Circlator fixstart step, see https://github.com/sanger-pathogens/circlator/wiki/Task:-fixstart
		--fixstart_args 			Circlator fixstart optional parameters (default=""). Example "--genes_fa /path/to/fasta"

	Assembly evaluation:
		--skip_quast 				SKip the QUAST assembly assessment step
		--quast_args 				QUAST optional parameters (default=""), see http://quast.sourceforge.net/docs/manual.html#sec2.3
		--quast_threads				Number of threads for QUAST (default=1)	
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

process basecalling {
	cpus "${params.guppy_num_callers}"
	label "gpu"
	label "guppy_gpu"
	containerOptions '--nv'
	publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*.txt'
	publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*.log'
	input:
		path(fast5_dir)
	output:
		path "sequencing_summary.txt", emit: sequencing_summary
		path "fastq_runid*fastq", emit: basecalled_fastq
		path("*.log")
		path("guppy_basecaller_version.txt")
	when:
	params.basecalling & params.gpu & params.demultiplexer == 'qcat'
	script:
	"""
	set +eu
	if [[ "${params.guppy_config_gpu}" != "false" ]] ; then
		${params.guppy_gpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --device ${params.guppy_gpu_device} --config ${params.guppy_config_gpu} ${params.guppy_basecaller_args}
	elif [[ "${params.flowcell}" != "false" ]] && [[ "${params.kit}" != "false" ]]; then
		${params.guppy_gpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --device ${params.guppy_gpu_device} --flowcell ${params.flowcell} --kit ${params.kit} --num_callers ${params.guppy_num_callers} ${params.guppy_basecaller_args}
	fi
	cp .command.log guppy_basecaller.log
	${params.guppy_gpu_folder}guppy_basecaller --version > guppy_basecaller_version.txt
	"""
}

process basecalling_single_isolate {
	cpus "${params.guppy_num_callers}"
	label "gpu"
	label "guppy_gpu"
	containerOptions '--nv'
	publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*.txt'
	publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*.log'
	publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*fastq.gz'
	input:
		tuple path(fast5_dir), val(sample)
	output:
		path "sequencing_summary.txt", emit: sequencing_summary
		path "*fastq.gz", emit: basecalled_fastq
		path("*.log")
		path("guppy_basecaller_version.txt")
	when:
	params.basecalling & params.gpu & !params.demultiplexing
	script:
	"""
	set +eu
	if [[ "${params.guppy_config_gpu}" != "false" ]] ; then
		${params.guppy_gpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --device ${params.guppy_gpu_device} --config ${params.guppy_config_gpu} ${params.guppy_basecaller_args}
	elif [[ "${params.flowcell}" != "false" ]] && [[ "${params.kit}" != "false" ]]; then
		${params.guppy_gpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --device ${params.guppy_gpu_device} --flowcell ${params.flowcell} --kit ${params.kit} --num_callers ${params.guppy_num_callers} ${params.guppy_basecaller_args}
	fi
	cp .command.log guppy_basecaller.log
	cat *.fastq > ${sample}.fastq
	gzip ${sample}.fastq
	${params.guppy_gpu_folder}guppy_basecaller --version > guppy_basecaller_version.txt
	"""
}

process basecalling_cpu {
    cpus "${params.guppy_num_callers}"
    label "cpu"
    label "guppy_cpu"
    publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*.txt'
    publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*.log'
    input:
        path(fast5_dir)
    output:
    	path "sequencing_summary.txt", emit: sequencing_summary
    	path "fastq_runid*fastq", emit: basecalled_fastq
		path("*.log")
		path("guppy_basecaller_version.txt")
    when:
	params.basecalling & !params.gpu & params.demultiplexer == 'qcat'
	script:
	"""
	set +eu
	if [[ "${params.guppy_config_cpu}" != "false" ]] ; then
		${params.guppy_cpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --config ${params.guppy_config_cpu} --num_callers ${params.guppy_num_callers} --cpu_threads_per_caller ${params.guppy_cpu_threads_per_caller} ${params.guppy_basecaller_args}
	elif [[ "${params.flowcell}" != "false" ]] && [[ "${params.kit}" != "false" ]]; then
		${params.guppy_cpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --flowcell ${params.flowcell} --kit ${params.kit} --num_callers ${params.guppy_num_callers} --cpu_threads_per_caller ${params.guppy_cpu_threads_per_caller} ${params.guppy_basecaller_args}
	fi
	cp .command.log guppy_basecaller.log
	${params.guppy_cpu_folder}guppy_basecaller --version > guppy_basecaller_version.txt
	"""
}

process basecalling_cpu_single_isolate {
	cpus "${params.guppy_num_callers}"
	label "cpu"
	label "guppy_cpu"
	publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*.txt'
	publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*.log'
	publishDir "$params.outdir/0_basecalling",  mode: 'copy', pattern: '*fastq.gz'
	input:
		tuple path(fast5_dir), val(sample)
	output:
		path "sequencing_summary.txt", emit: sequencing_summary
		path "*fastq.gz", emit: basecalled_fastq
		path("*.log")
		path("guppy_basecaller_version.txt")
	when:
	params.basecalling & !params.gpu & !params.demultiplexing
	script:
	"""
	set +eu
	if [[ "${params.guppy_config_cpu}" != "false" ]] ; then
		${params.guppy_cpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --config ${params.guppy_config_cpu} --num_callers ${params.guppy_num_callers} --cpu_threads_per_caller ${params.guppy_cpu_threads_per_caller} ${params.guppy_basecaller_args}
	elif [[ "${params.flowcell}" != "false" ]] && [[ "${params.kit}" != "false" ]]; then
		${params.guppy_cpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --flowcell ${params.flowcell} --kit ${params.kit} --num_callers ${params.guppy_num_callers} --cpu_threads_per_caller ${params.guppy_cpu_threads_per_caller} ${params.guppy_basecaller_args}
	fi
	cp .command.log guppy_basecaller.log
	cat *.fastq > ${sample}.fastq
	gzip ${sample}.fastq
	${params.guppy_cpu_folder}guppy_basecaller --version > guppy_basecaller_version.txt
	"""
}

process demultiplexing_qcat {
	cpus 1
	label "cpu"
	publishDir "$params.outdir/0_basecalling",  mode: 'copy'
	input:
		path(fastq)
	output:
		path "*fastq.gz", emit: demultiplexed_fastq
		path("qcat.log")
	path("qcat_version.txt")
	when:
	params.demultiplexer == 'qcat'
	script:
	"""
	set +eu
	if [[ -d "${fastq}" ]] ; then
		zcat ${fastq}/*q.gz | qcat -b \$PWD ${params.qcat_args}
	else
		cat ${fastq} | qcat -b \$PWD ${params.qcat_args}
	fi
	cp .command.log qcat.log
	gzip barcode*fastq 
	if [[ -f "none.fastq" ]] ; then
		gzip none.fastq
	fi
	qcat --version > qcat_version.txt
	"""
}

process basecalling_demultiplexing_guppy {
	cpus "${params.guppy_num_callers}"
	label "gpu"
	label "guppy_gpu"
	containerOptions '--nv'
	publishDir "$params.outdir/0_demultiplexing", mode: 'copy'
	input:
		path(fast5_dir)
	output:
		path "sequencing_summary.txt", emit: sequencing_summary
		path "*fastq.gz", emit: demultiplexed_fastq
		path("*log")
		path("guppy_basecaller_version.txt")
	when:
	params.basecalling & params.gpu & params.demultiplexer == 'guppy'
	script:
	"""
	set +eu
	if [[ "${params.guppy_config_gpu}" != "false" ]]; then
		${params.guppy_gpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --device ${params.guppy_gpu_device} --config "${params.guppy_config_gpu}" --compress_fastq --num_callers ${params.guppy_num_callers} ${params.guppy_basecaller_args} --barcode_kits ${params.guppy_barcode_kits}
	elif [[ "${params.flowcell}" != "false" ]] && [[ "${params.kit}" != "false" ]]; then
		${params.guppy_gpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --device ${params.guppy_gpu_device} --flowcell ${params.flowcell} --kit ${params.kit} --compress_fastq --num_callers ${params.guppy_num_callers} ${params.guppy_basecaller_args} --barcode_kits ${params.guppy_barcode_kits}	
	fi
	cp .command.log guppy_basecaller.log
	for dir in barc*/ uncl*/; do
		barcode_id=\${dir%*/}
		cat \${dir}/*.fastq.gz > \${barcode_id}.fastq.gz
	done
	${params.guppy_gpu_folder}guppy_basecaller --version > guppy_basecaller_version.txt
	"""
}

process basecalling_demultiplexing_guppy_cpu {
	cpus "${params.guppy_num_callers}"
	label "cpu"
	label "guppy_cpu"
	publishDir "$params.outdir/0_demultiplexing", mode: 'copy'
	input:
		path(fast5_dir)
	output:
		path "sequencing_summary.txt", emit: sequencing_summary
		path "*fastq.gz", emit: demultiplexed_fastq
		path("*log")
		path("guppy_basecaller_version.txt")
	when:
	params.basecalling & !params.gpu & params.demultiplexer == 'guppy'
	script:
	"""
	set +eu
	if [[ "${params.guppy_config_gpu}" != "false" ]] ; then
		${params.guppy_cpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD --config "${params.guppy_config_cpu}" --compress_fastq --num_callers ${params.guppy_num_callers} --cpu_threads_per_caller ${params.guppy_cpu_threads_per_caller} ${params.guppy_basecaller_args} --barcode_kits ${params.guppy_barcode_kits}
	elif [[ "${params.flowcell}" != "false" ]] && [[ "${params.kit}" != "false" ]]; then
		${params.guppy_cpu_folder}guppy_basecaller -i ${fast5_dir} -s \$PWD  --flowcell ${params.flowcell} --kit ${params.kit} --compress_fastq --num_callers ${params.guppy_num_callers} --cpu_threads_per_caller ${params.guppy_cpu_threads_per_caller} ${params.guppy_basecaller_args} --barcode_kits ${params.guppy_barcode_kits}
	fi
	cp .command.log guppy_basecaller.log
	for dir in barc*/ uncl*/; do
		barcode_id=\${dir%*/}
		cat \${dir}/*.fastq.gz > \${barcode_id}.fastq.gz
	done
	${params.guppy_cpu_folder}guppy_basecaller --version > guppy_basecaller_version.txt
	"""
}

process demultiplexing_guppy {
	cpus "${params.guppy_barcoder_threads}"
	label "gpu"
	label "guppy_gpu"
	containerOptions '--nv'
	publishDir "$params.outdir/0_demultiplexing", mode: 'copy'
	input:
		path(fastq_dir)
	output:
		path "*fastq.gz", emit: demultiplexed_fastq
		path("*log")
		path "barcoding_summary.txt"
		path("guppy_barcoder_version.txt")
	when:
	params.demultiplexer == 'guppy' & params.demultiplexing & params.gpu
	script:
	"""
	set +eu
	${params.guppy_gpu_folder}guppy_barcoder -i ${fastq_dir} -s \$PWD --device ${params.guppy_gpu_device} --compress_fastq ${params.guppy_barcoder_args} --barcode_kits ${params.guppy_barcode_kits} --worker_threads ${params.guppy_barcoder_threads}
	cp .command.log guppy_barcoder.log
	for dir in barc*/ uncl*/; do
		barcode_id=\${dir%*/}
		cat \${dir}/*.fastq.gz > \${barcode_id}.fastq.gz
	done
	${params.guppy_gpu_folder}guppy_barcoder --version > guppy_barcoder_version.txt
	"""
}

process demultiplexing_guppy_cpu {
    cpus "${params.guppy_barcoder_threads}"
    label "cpu"
    label "guppy_cpu"
    publishDir "$params.outdir/0_demultiplexing", mode: 'copy'
    input:
		path(fastq_dir)
    output:
    	path "*fastq.gz", emit: demultiplexed_fastq
		path("*log")
		path "barcoding_summary.txt"
		path("guppy_barcoder_version.txt")
	when:
	params.demultiplexer == 'guppy' & params.demultiplexing & !params.gpu
	script:
	"""
	set +eu
	${params.guppy_cpu_folder}guppy_barcoder -i ${fastq_dir} -s \$PWD --compress_fastq ${params.guppy_barcoder_args} --barcode_kits ${params.guppy_barcode_kits} --worker_threads ${params.guppy_barcoder_threads}
	cp .command.log guppy_barcoder.log
	for dir in barc*/ uncl*/; do
		barcode_id=\${dir%*/}
		cat \${dir}/*.fastq.gz > \${barcode_id}.fastq.gz
	done
	${params.guppy_cpu_folder}guppy_barcoder --version > guppy_barcoder_version.txt
	"""
}

process pycoqc {
	cpus 1
	label "cpu"
	label "pycoqc"
	publishDir "$params.outdir/0_pycoQC",  mode: 'copy'
	input:
		path(sequencing_summary)
	output:
		path("pycoQC.html")
		path("pycoqc_version.txt")
	when:
	params.basecalling & !params.skip_pycoqc
	script:
	"""
	set +eu
	pycoQC -f ${sequencing_summary} -o \$PWD/pycoQC.html
	pycoQC --version > pycoqc_version.txt
	"""
}

process rasusa {
	cpus 1
	tag "${sample}"
	label "cpu"
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: "*_version.txt"
	input:
		tuple val(barcode), file(long_reads), val(sample), val(genome_size)
	output:
		tuple val(barcode), file("subsampled.fastq.gz"), val(sample), val(genome_size), emit: subsampled_fastq
		path("rasusa.log")
		path("rasusa_version.txt")
	when:
	!params.skip_rasusa
	script:
	"""
	set +eu
	rasusa --coverage ${params.rasusa_coverage} --genome-size ${genome_size} --input ${long_reads} --output subsampled.fastq.gz
	cp .command.log rasusa.log
	rasusa --version > rasusa_version.txt
	"""
}

process porechop {
	cpus "${params.porechop_threads}"
	tag "${sample}"
	label "cpu"
	label "big_mem"
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: "*_version.txt"
	input:
		tuple val(barcode), file(long_reads), val(sample), val(genome_size)
	output:
		tuple val(barcode), file("trimmed.fastq.gz"), val(sample), val(genome_size), emit: trimmed_fastq
		path("porechop.log")
		path("porechop_version.txt")
	when:
	!params.skip_porechop
	script:
	"""
	set +eu
	porechop -i ${long_reads} -t ${params.porechop_threads} -o trimmed.fastq.gz ${params.porechop_args}
	cp .command.log porechop.log
	porechop --version > porechop_version.txt
	"""
}

process japsa {
	cpus 1
	tag "${sample}"
	label "cpu"
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: '*filtered.fastq.gz', saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(barcode), path(trimmed), val(sample), val(genome_size)
	output:
		tuple val(barcode), path("filtered.fastq.gz"), val(sample), val(genome_size), emit: filtered_fastq
	when:
	!params.skip_filtering & params.filtering == 'japsa'
	script:
	"""
	set +eu
	jsa.np.filter --input ${trimmed} ${params.japsa_args} --output filtered.fastq.gz
	"""
}

process filtlong {
	cpus 1
	tag "${sample}"
	label "cpu"
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: '*filtered.fastq.gz', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: 'filtlong.log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/1_filtering",  mode: 'copy', pattern: "*_version.txt"
	input:
		tuple val(barcode), path(trimmed), val(sample), val(genome_size)
	output:
		tuple val(barcode), path("filtered.fastq.gz"), val(sample), val(genome_size), emit: filtered_fastq
		path("*.log")
		path("filtlong_version.txt")
	when:
	!params.skip_filtering & params.filtering == 'filtlong'
	script:
	"""
	set +eu
	filtlong ${params.filtlong_args} ${trimmed} | gzip > filtered.fastq.gz
	cp .command.log filtlong.log
	filtlong --version > filtlong_version.txt
	"""
}

process flye {
	cpus "${params.flye_threads}"
	tag "${sample}"
	label "cpu"
	label "big_mem"
	publishDir "$params.outdir/$sample/2_assembly",  mode: 'copy', pattern: 'assembly*', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/2_assembly",  mode: 'copy', pattern: 'flye.log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/2_assembly",  mode: 'copy', pattern: "*_version.txt"
	input:
	tuple val(barcode), path(filtered), val(sample), val(genome_size)
	output:
	tuple val(barcode), path(filtered), val(sample), path("assembly.fasta"), path("assembly_info.txt"), path("assembly_graph.gfa"), path("assembly_graph.gv"), emit: assembly_out
	path("flye.log")
	path("flye_version.txt")
	script:
	"""
	set +eu
	flye --nano-raw ${filtered} --genome-size ${genome_size} --threads ${params.flye_threads} --out-dir \$PWD ${params.flye_args}
	flye -v 2> flye_version.txt
	"""
}

prefix="flye"
prefix_lr="flye_polishedLR"
prefix_lr_sr="flye_polishedLR_SR"
raconv="racon"
medakav="medaka"

process racon_cpu {
	cpus "${params.racon_threads}"
	tag "${sample}"
	label "cpu"
	label "racon"
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: '*fasta', saveAs: { filename -> "${sample}_${prefix}_${raconv}_${params.racon_nb}.fasta"}
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: "*_version.txt"
	input:
		tuple val(barcode), path(filtered), val(sample), path(assembly), path(info), path(gfa), path(gv)
	output:
		tuple val(barcode), path(filtered), val(sample), path("${prefix}_${raconv}_${params.racon_nb}.fasta"), emit: polished_racon
		path("racon.log")
		path("racon_version.txt")
	when:
	params.polisher == 'medaka'
	script:
	"""
	set +eu
	ln -s ${assembly} ${prefix}_${raconv}_0.fasta
	for i in `seq 1 ${params.racon_nb}`; do
		ii=\$((\$i-1))
		minimap2 -t ${params.racon_threads} -ax map-ont ${prefix}_${raconv}_\$ii.fasta ${filtered} > ${prefix}.gfa\$i.sam
		racon ${params.racon_args} -t ${params.racon_threads} ${filtered} ${prefix}.gfa\$i.sam ${prefix}_${raconv}_\$ii.fasta --include-unpolished > ${prefix}_${raconv}_\$i.fasta
		rm ${prefix}.gfa\$i.sam
		python3 $projectDir/bin/rotate_circular_fasta.py ${prefix}_${raconv}_\$i.fasta ${info} ${prefix}_${raconv}_\$i.tmp.fasta
		cp ${prefix}_${raconv}_\$i.tmp.fasta ${prefix}_${raconv}_\$i.fasta
	done
	cp .command.log racon.log
	racon --version > racon_version.txt
	"""
}

process medaka_cpu {
	cpus "${params.medaka_threads}"
	tag "${sample}"
	label "cpu"
	label "medaka"
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: '*fasta', saveAs: { filename -> "${sample}_${prefix_lr}.fasta"}
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: "*_version.txt"	
    input:
	tuple val(barcode), path(filtered), val(sample), path(draft)
	output:
	tuple val(barcode), path(filtered), val(sample), path ("consensus.fasta"), emit: polished_medaka
	path("medaka.log")
	path("medaka_version.txt")
	when:
	params.polisher == 'medaka'
	script:
	"""
	set +eu
	medaka_consensus -i ${filtered} -d ${draft} -o \$PWD -t ${params.medaka_threads} -m ${params.medaka_model}
	rm consensus_probs.hdf calls_to_draft.bam calls_to_draft.bam.bai
	cp .command.log medaka.log
	medaka --version > medaka_version.txt
	"""
}

process nextpolish_LR {
	cpus "${params.nextpolish_threads}"
	tag "${sample}"
	label "cpu"
	label "nextpolish"
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: '*fasta', saveAs: { filename -> "${sample}_${prefix_lr}.fasta"}
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: "*_version.txt"
	input:
		tuple val(barcode), path(filtered), val(sample), path(assembly), path(info), path(gfa), path(gv)
	output:
		tuple val(barcode), path(filtered), val(sample), path ("${sample}_${prefix_lr}.fasta"), emit: polished_LR
		path("nextpolish_LR.log")
		path("nextpolish_version.txt")
	when:
	params.polisher == 'nextpolish'
	script:
	"""
	set +eu
	ls ${filtered} > lgs.fofn
	echo -e "task = ${params.nextpolish_task_LR}\ngenome = ${assembly}\nmultithread_jobs = ${task.cpus}\nlgs_fofn = lgs.fofn\nlgs_minimap2_options = -x map-ont -t ${params.nextpolish_threads}" > nextpolish.cfg
	nextPolish nextpolish.cfg
	if [[ "${params.nextpolish_task_LR}" == "55" ]] || [[ "${params.nextpolish_task_LR}" == "best" ]] ; then
		cat 01.lgs_polish/*polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > ${sample}_${prefix_lr}.fasta
		rm -r 00.lgs_polish 01.lgs_polish
	elif [[ "${params.nextpolish_task_LR}" == "5" ]]; then
		cat 00.lgs_polish/*polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > ${sample}_${prefix_lr}.fasta
		rm -r 00.lgs_polish
	fi
	rm input.lgspart.*.gz
	cp .command.log nextpolish_LR.log
	nextPolish --version 2> nextpolish_version.txt
	"""
}

process nextpolish {
	cpus "${params.nextpolish_threads}"
	tag "${sample}"
	label "cpu"
	label "nextpolish"
	publishDir "$params.outdir/$sample/4_polishing_short_reads",  mode: 'copy', pattern: '*fasta', saveAs: { filename -> "${sample}_${prefix_lr_sr}.fasta"}
	publishDir "$params.outdir/$sample/4_polishing_short_reads",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/4_polishing_short_reads",  mode: 'copy', pattern: "*_version.txt"
	input:
		tuple val(barcode), path(filtered), val(sample), path(draft), path(reads_1), path(reads_2)
	output:
		tuple val(barcode), path(filtered), val(sample), path ("${sample}_${prefix_lr_sr}_2.fasta"), emit: polished_SR
		path("nextpolish.log")
		path("nextpolish_version.txt")
	when:
	!params.skip_illumina
	script:
	"""
	set +eu
	ls ${reads_1} ${reads_2} > sgs.fofn
	echo -e "task = ${params.nextpolish_task_SR}\ngenome = ${draft}\nsgs_fofn = sgs.fofn\nmultithread_jobs = ${params.nextpolish_threads}" > nextpolish.cfg
	nextPolish nextpolish.cfg
	if [[ "${params.nextpolish_task_SR}" == "1212" ]] || [[ "${params.nextpolish_task_SR}" == "best" ]] ; then
		cat 01.kmer_count/*polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > ${sample}_${prefix_lr_sr}_1.fasta
		cat 03.kmer_count/*polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > ${sample}_${prefix_lr_sr}_2.fasta
		rm -r 00.score_chain 01.kmer_count 02.score_chain 03.kmer_count
	elif [[ "${params.nextpolish_task_SR}" == "12" ]]; then	
		cat 01.kmer_count/*polish.ref.sh.work/polish_genome*/genome.nextpolish.part*.fasta > ${sample}_${prefix_lr_sr}_2.fasta
		rm -r 00.score_chain 01.kmer_count
	fi
	rm input.sgspart*.fastq.gz
	cp .command.log nextpolish.log
	nextPolish --version 2> nextpolish_version.txt
	"""
}

process fixstart {
	cpus 1
	label "circlator"
	tag "${sample}"
	publishDir "$params.outdir/$sample/4_polishing_short_reads",  mode: 'copy', pattern: '*fixstart.fasta', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/4_polishing_short_reads",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(barcode), path(filtered), val(sample), path(polished)
	output:
		tuple val(barcode), path(filtered), val(sample), path ("${prefix_lr_sr}_fixstart.fasta"), emit: polished_fixstart 
		path("*log")
	script:
	"""
	set +eu
	circlator fixstart ${params.fixstart_args} ${polished} ${prefix_lr_sr}_fixstart
	"""
}

process fixstart_LR {
	cpus 1
	label "circlator"
	tag "${sample}"
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: '*fixstart.fasta', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/3_polishing_long_reads",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
	input:
		tuple val(barcode), path(filtered), val(sample), path(polished)
	output:
		tuple val(barcode), path(filtered), val(sample), path ("${prefix_lr}_fixstart.fasta"), emit: polished_fixstart 
		path("*log")
	script:
	"""
	set +eu
	circlator fixstart ${params.fixstart_args} ${polished} ${prefix_lr}_fixstart
	"""
}

process quast {
	cpus "${params.quast_threads}"
	tag "${sample}"
	label "cpu"
	publishDir "$params.outdir/$sample/5_quast",  mode: 'copy', pattern: 'report*', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_quast",  mode: 'copy', pattern: 'quast.log', saveAs: { filename -> "${sample}_$filename" }
	publishDir "$params.outdir/$sample/5_quast",  mode: 'copy', pattern: "*_version.txt"
	input:
		tuple val(barcode), path(filtered), val(sample), path(polished)
	output:
		tuple path("report.txt"), path("report.html"), path("report.tsv"), path("report.pdf"), path("quast.log"), emit: quast_out
		path("quast_version.txt")
	script:
	"""
	set +eu
	quast.py -o \$PWD -t ${params.quast_threads} -l ${sample} ${polished} ${params.quast_args}
	quast --version > quast_version.txt
	"""
}

workflow assembly {
	take: 
	ch_samplesheet
	ch_samplesheet_illumina
	main:
	if (!params.skip_porechop & !params.skip_filtering) {
		if (!params.skip_rasusa) {
			rasusa(ch_samplesheet)
			porechop(rasusa.out.subsampled_fastq)
		} else if (params.skip_rasusa) {
			porechop(ch_samplesheet)
		}
		if (params.filtering == "japsa") {
			japsa(porechop.out.trimmed_fastq)
			flye(japsa.out.filtered_fastq)
		} else if (params.filtering == "filtlong") {
			filtlong(porechop.out.trimmed_fastq)
			flye(filtlong.out.filtered_fastq)
		}
	} else if (!params.skip_porechop & params.skip_filtering) {
		if (!params.skip_rasusa) {
			rasusa(ch_samplesheet)
			porechop(rasusa.out.subsampled_fastq)
		} else if (params.skip_rasusa) {
			porechop(ch_samplesheet)
		}
		flye(porechop.out.trimmed_fastq)
	} else if (params.skip_porechop & !params.skip_filtering) {
		if (params.filtering == "japsa") {
			if (!params.skip_rasusa) {
				rasusa(ch_samplesheet)
				japsa(rasusa.out.subsampled_fastq)
				flye(japsa.out.filtered_fastq)
			} else if (params.skip_rasusa) {
				japsa(ch_samplesheet)
				flye(japsa.out.filtered_fastq)
			}
		} else if (params.filtering == "filtlong") {
			if (!params.skip_rasusa) {
				rasusa(ch_samplesheet)
				filtlong(rasusa.out.subsampled_fastq)
				flye(filtlong.out.filtered_fastq)
			} else if (params.skip_rasusa) {
				filtlong(ch_samplesheet)
				flye(filtlong.out.filtered_fastq)
            }
		}
	} else {
		if (!params.skip_rasusa) {
			rasusa(ch_samplesheet)
			flye(rasusa.out.subsampled_fastq)
		} else if (params.skip_rasusa) {
			flye(ch_samplesheet)
		}
	}
	if (params.polisher == 'medaka') {
		racon_cpu(flye.out.assembly_out)
		medaka_cpu(racon_cpu.out.polished_racon)
		if (!params.skip_illumina) {
			nextpolish(medaka_cpu.out.polished_medaka.combine (ch_samplesheet_illumina, by: 0))
			if (params.skip_fixstart) {
			quast(nextpolish.out.polished_SR)
			}
			else if (!params.skip_fixstart) {
				fixstart(nextpolish.out.polished_SR)
				quast(fixstart.out.polished_fixstart)
			}
		}
		else if (params.skip_illumina) {
			if (params.skip_fixstart) {
				quast(medaka_cpu.out.polished_medaka)
			}
			else if (!params.skip_fixstart) {
				fixstart_LR(medaka_cpu.out.polished_medaka)
				quast(fixstart_LR.out.polished_fixstart)
			}
		}
	}
	else if (params.polisher == 'nextpolish') {
		nextpolish_LR(flye.out.assembly_out)
		if (!params.skip_illumina) {
			nextpolish(nextpolish_LR.out.polished_LR.combine (ch_samplesheet_illumina, by: 0))
			if (params.skip_fixstart) {
				quast(nextpolish.out.polished_SR)
			}
			else if (!params.skip_fixstart) {
				fixstart(nextpolish.out.polished_SR)
				quast(fixstart.out.polished_fixstart)
			}
		}
		else if (params.skip_illumina) {
			if (params.skip_fixstart) {
				quast(nextpolish_LR.out.polished_LR)
			}
			else if (!params.skip_fixstart) {
				fixstart_LR(nextpolish_LR.out.polished_LR)
				quast(fixstart_LR.out.polished_fixstart)
			}
		}
	}
}

workflow {
	//basecalling, demultiplexing and assembly workflow
	if( params.basecalling && params.demultiplexing) {
		Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
		.splitCsv(header:true, sep:',')
		.map { row -> tuple(row.barcode_id, row.sample_id, row.genome_size) }
		.set { ch_samplesheet_basecalling }
		ch_samplesheet_basecalling.view()
		if ( !params.skip_illumina ) {
			Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
			.splitCsv(header:true, sep:',')          
			.map { row -> tuple(row.barcode_id, file(row.short_fastq_1, checkIfExists: true), file(row.short_fastq_2, checkIfExists: true)) }
			.set { ch_samplesheet_illumina }
			ch_samplesheet_illumina.view()
		}
		fast5 = Channel.fromPath("${params.fast5}", checkIfExists: true )
		if( params.demultiplexer == "qcat") {
			if( params.gpu ) {
				basecalling(fast5)
				demultiplexing_qcat(basecalling.out.basecalled_fastq)
				pycoqc(basecalling.out.sequencing_summary)
			} else {
				basecalling_cpu(fast5)
				demultiplexing_qcat(basecalling_cpu.out.basecalled_fastq)
				pycoqc(basecalling_cpu.out.sequencing_summary)
			}
			ch_fastq=demultiplexing_qcat.out.demultiplexed_fastq.map { file -> tuple(file.simpleName, file) }.transpose()
			ch_fastq.view()
			ch_data=ch_fastq.combine(ch_samplesheet_basecalling, by: 0)
			ch_data.view()
		} else if (params.demultiplexer == "guppy") {
			if( params.gpu ) {
				basecalling_demultiplexing_guppy(fast5)
				pycoqc(basecalling_demultiplexing_guppy.out.sequencing_summary)
				ch_fastq=basecalling_demultiplexing_guppy.out.demultiplexed_fastq.map { file -> tuple(file.simpleName, file) }.transpose()
			} else {
				basecalling_demultiplexing_guppy_cpu(fast5)
				pycoqc(basecalling_demultiplexing_guppy_cpu.out.sequencing_summary)
				ch_fastq=basecalling_demultiplexing_guppy_cpu.out.demultiplexed_fastq.map { file -> tuple(file.simpleName, file) }.transpose()
			}
			ch_fastq.view()
			ch_data=ch_fastq.combine(ch_samplesheet_basecalling, by: 0)
		}
		if ( !params.skip_illumina ) {
			assembly( ch_data, ch_samplesheet_illumina)
		} else {
			assembly( ch_data, Channel.empty() )
		}
	//basecalling and assembly workflow (single isolate)
	} else if( params.basecalling && !params.demultiplexing) {
		Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
		.splitCsv(header:true, sep:',')
		.map { row -> tuple(row.sample_id, row.genome_size) }
		.set { ch_samplesheet_basecalling }
		ch_samplesheet_basecalling.view()
		if ( !params.skip_illumina ) {
			Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
			.splitCsv(header:true, sep:',')          
			.map { row -> tuple(row.sample_id, file(row.short_fastq_1, checkIfExists: true), file(row.short_fastq_2, checkIfExists: true)) }
			.set { ch_samplesheet_illumina }
			ch_samplesheet_illumina.view()
		}
		fast5 = Channel.fromPath("${params.fast5}", checkIfExists: true )
		ch_sample = ch_samplesheet_basecalling.first().map { it[0] }
		ch_fast5 = fast5.concat( ch_sample ).collect()
		ch_fast5.view()
		if( params.gpu ) {
			basecalling_single_isolate(ch_fast5)
			pycoqc(basecalling_single_isolate.out.sequencing_summary)
			ch_fastq=basecalling_single_isolate.out.basecalled_fastq.map { file -> tuple(file.simpleName, file) }.transpose()
		} else {
			basecalling_cpu_single_isolate(ch_fast5)
			pycoqc(basecalling_cpu_single_isolate.out.sequencing_summary)
			ch_fastq=basecalling_cpu_single_isolate.out.basecalled_fastq.map { file -> tuple(file.simpleName, file) }.transpose()
		}
		ch_fastq.view()
		if ( !params.skip_illumina ) {
			ch_data = ch_fastq.concat( ch_samplesheet_basecalling ).collect()
			ch_data.view()
			assembly( ch_data, ch_samplesheet_illumina )
		} else {
			ch_data = ch_fastq.concat( ch_samplesheet_basecalling ).collect()
			ch_data.view()
			assembly( ch_data, Channel.empty() )
		}
	//demultiplexing and assembly workflow
	} else if ( !params.basecalling && params.demultiplexing ){
		Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
		.splitCsv(header:true, sep:',')
		.map { row -> tuple(row.barcode_id, row.sample_id, row.genome_size) }
		.set { ch_samplesheet_basecalling }
		ch_samplesheet_basecalling.view()
		if ( !params.skip_illumina ) {
			Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
			.splitCsv(header:true, sep:',')          
			.map { row -> tuple(row.barcode_id, file(row.short_fastq_1, checkIfExists: true), file(row.short_fastq_2, checkIfExists: true)) }
			.set { ch_samplesheet_illumina }
			ch_samplesheet_illumina.view()
		}
		fastq = Channel.fromPath("${params.fastq}", checkIfExists: true )
		if( params.demultiplexer == "qcat") {
			demultiplexing_qcat(fastq)
			ch_fastq=demultiplexing_qcat.out.demultiplexed_fastq.map { file -> tuple(file.simpleName, file) }.transpose()
		} else if (params.demultiplexer == "guppy") {
			if( params.gpu ) {
				demultiplexing_guppy(fastq)
				ch_fastq=demultiplexing_guppy.out.demultiplexed_fastq.map { file -> tuple(file.simpleName, file) }.transpose()
			} else {
				demultiplexing_guppy_cpu(fastq)
				ch_fastq=demultiplexing_guppy_cpu.out.demultiplexed_fastq.map { file -> tuple(file.simpleName, file) }.transpose()
			}
		}
		ch_fastq.view()
		ch_data=ch_fastq.combine(ch_samplesheet_basecalling, by: 0)
		ch_data.view()
		if ( !params.skip_illumina ) {
			assembly( ch_data, ch_samplesheet_illumina)
		} else {
			assembly( ch_data, Channel.empty() )
		}
	//assembly only workflow
	} else if ( !params.basecalling && !params.demultiplexing ) {
		Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
		.splitCsv(header:true, sep:',')
		.map { row -> tuple(row.barcode_id, file(row.long_fastq, checkIfExists: true), row.sample_id, row.genome_size) }
		.set { ch_samplesheet }
		ch_samplesheet.view()
		if ( !params.skip_illumina ) {
			Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
			.splitCsv(header:true, sep:',')          
			.map { row -> tuple(row.barcode_id, file(row.short_fastq_1, checkIfExists: true), file(row.short_fastq_2, checkIfExists: true)) }
			.set { ch_samplesheet_illumina }
			ch_samplesheet_illumina.view()
			assembly( ch_samplesheet, ch_samplesheet_illumina )
		} else {
			assembly( ch_samplesheet, Channel.empty() )
		}
	}
}
