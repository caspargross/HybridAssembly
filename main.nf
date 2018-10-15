#!/usr/bin/env/ nextflow

/* 
===============================================================================
   M I C R O B I A L   H Y B R I D   A S S E M B L Y   P I P E L I N E 
===============================================================================
Nextflow pipeline for hybrid assembly, quality check and plasmid finding
of bacterial genomes.
-------------------------------------------------------------------------------
@ Author
Caspar Groß <mail@caspar.one>
-------------------------------------------------------------------------------
@ Documentation
https://github.com/caspargross/hybridassembly/README.md
------------------------------------------------------------------------------
Processes overview:
... to be completed
------------------------------------------------------------------------------
*/


/* 
------------------------------------------------------------------------------
                       C O N F I G U R A T I O N 
------------------------------------------------------------------------------
*/
//Define valid run modes:
validModes = ['spades_simple', 'spades', 'spades_plasmid', 'canu', 'unicycler', 'flye', 'miniasm', 'all']


if (params.help) exit 0, helpMessage()
if (!params.mode) exit 0, helpMessage()

// check if mode input is valid
modes = params.mode.tokenize(',')
if (!modes.every{validModes.contains(it)}) {
    exit 1,  log.info "Wrong execution mode, should be one of " + validModes
}

// Target coverage for long reads before assembly
target_sr_number = (params.target_shortread_cov * params.genome_size) / params.shortread_length
target_lr_length = params.target_longread_cov * params.genome_size

// inputFiles
files = Channel.fromPath(params.pathFile)
    .ifEmpty {error "Cannot find file with path locations in ${params.files}"}
    .splitCsv(header: true)
    .view()


// Shorthands for conda environment activations
PY27 = "source activate ha_py27"
PY36 = "source activate ha_py36"

/* 
------------------------------------------------------------------------------
                           P R O C E S S E S 
------------------------------------------------------------------------------
*/

process seqpurge {
// Trim adapters on short read files
    tag{id}
    
    input:
    set id, sr1, sr2, lr from files

    output:
    set id, file('sr1.fastq.gz'), file('sr2.fastq.gz'), lr into files_purged
    file("${id}_readQC.qcml")
    
    script:
    """
    $PY27  
    SeqPurge -t ${params.cpu} -gc ${id}_readQC.qcml -in1 ${sr1} -in2 ${sr2} -threads ${params.cpu} -out1 sr1.fastq.gz -out2 sr2.fastq.gz
    """
}

process sample_shortreads {
// Subset short reads
    tag{id}

    input:
    set id, sr1, sr2, lr from files_purged

    output:
    set id, file('sr1_filt.fastq'), file('sr2_filt.fastq'), lr into files_filtered
    
    script:
    """
    $PY27
    seqtk sample -s100 ${sr1} ${target_sr_number} > sr1_filt.fastq 
    seqtk sample -s100 ${sr2} ${target_sr_number} > sr2_filt.fastq 
    """
}

process porechop {
// Trim adapter sequences on long read nanopore files
    tag{id}
        
    input:
    set id, sr1, sr2, lr from files_filtered
    
    output:
    set id, sr1, sr2, file('lr_porechop.fastq') into files_porechop
    set id, lr, val("raw") into files_nanoplot_raw
    
    script:
    // Join multiple longread files if possible 
    """
    $PY36
    cat ${lr} > nanoreads.fastq
    porechop -i nanoreads.fastq -t ${params.cpu} -o lr_porechop.fastq
    """
}

process filtlong {
// Quality filter long reads
    tag{id}

    input: 
    set id, sr1, sr2, lr from files_porechop
    
    output:
    set id, sr1, sr2, file("lr_filtlong.fastq") into files_pre_unicycler, files_pre_spades, files_pre_canu, files_pre_miniasm, files_pre_flye,  files_fastqc
    set id, file("lr_filtlong.fastq"), val('filtered') into files_nanoplot_filtered
    
    script:
    """
    $PY36
    filtlong -1 ${sr1} -2 ${sr2} \
    --min_length 1000 \
    --keep_percent 90 \
    --target_bases  ${target_lr_length} \
    ${lr} > lr_filtlong.fastq
    """
}

process nanoplot {
// Quality check for nanopore reads and Quality/Length Plots
    tag{id}
    publishDir "${params.outDir}/${id}_${params.assembly}/nanoplot/", mode: 'copy'
    
    input:
    set id, lr, type from files_nanoplot_raw .mix(files_nanoplot_filtered)

    output:
    file '*'
    
    script:
    """
    $PY36
    NanoPlot -t ${params.cpu} -p ${type}_  --title ${id}_${type} -c darkblue --fastq ${lr}
    """
}

/*process fastqc{
// Create FASTQC quality check on short reads
    tag{id}
    
    publishDir "${params.outDir}/${id}_${params.assembly}/fastQC/", mode: 'copy'

    input: 
    set id, sr1, sr2, lr from files_fastqc

    output:
    file "fastqc/\*"

    script: 
    """
    mkdir -p fastqc
    ${FASTQC} ${sr1} ${sr2} -o fastqc
    """
} */

process unicycler{
// complete bacterial hybrid assembly pipeline
    tag{id}
    publishDir "${params.outDir}/${id}_${params.assembly}/unicycler", mode: 'copy'   
   
    input:
    set id, sr1, sr2, lr from files_pre_unicycler

    output:
    set id, sr1, sr2, lr, file("${id}/assembly.fasta"), val('unicycler') into assembly_unicycler
    file("${id}/*")

    when:
    params.assembly in ['unicycler', 'all']

    script:
    """ 
    $PY36
    unicycler -1 ${sr1} -2 ${sr2} -l ${lr} -o ${id} -t ${params.cpu}
    """
}

/*
* SPAades assembler
*
*
*/
process spades{
    tag{id}
    publishDir "${params.outDir}/${id}_${params.assembly}", mode: 'copy'   

    input:
    set id, forward, reverse, longread from files_pre_spades  

    output:
    set id, forward, reverse, longread, file("spades/scaffolds.fasta"), val('spades_plasmid') into files_spades_sspace, files_spades_links, files_spades_plasmid 
    file("spades/scaffolds.fasta")
    file("spades/${id}_spades_graph.gfa")


    when:
    params.assembly in ['spades_sspace','spades_links','all']
     
    script:
    if ( params.plasmid == true) 
        """
        ${SPADES} -t ${params.cpu} -m ${params.mem} \
        --phred-offset 33 --careful \
        --pe1-1 ${forward} \
        --pe1-2 ${reverse} \
        --nanopore ${longread} \
        --plasmid \
        -o spades
        cp spades/assembly_graph_with_scaffolds.gfa spades/${id}_spades_graph.gfa
        """
    else  
        """
        ${SPADES} -t ${params.cpu} -m ${params.mem} \
        --phred-offset 33 --careful \
        --pe1-1 ${forward} \
        --pe1-2 ${reverse} \
        --nanopore ${longread} \
        -o spades
        cp spades/assembly_graph_with_scaffolds.gfa spades/${id}_spades_graph.gfa
        """
}


/*
*  SSPACE scaffolder 
*
*
*/
process sspace_scaffolding{
    tag{data_id}

    input:
    set data_id, forward, reverse, longread, scaffolds, plasmid from files_spades_sspace  

    output:
    set data_id, forward, reverse, longread, file("sspace/scaffolds.fasta"), val('spades_sspace') into files_sspace 
    
    when:
    params.assembly in ['spades_sspace','all'] && !params.plasmid

    script:
    """
    perl ${SSPACE} -c ${scaffolds} -p ${longread} -b sspace -t ${params.cpu}
    """
}
    
/*
* Links scaffolder
*
*
*/
process links_scaffolding{
    tag{data_id}
    
    input:
    set data_id, forward, reverse, longread, scaffolds, plasmid from files_spades_links
    
    output:
    set data_id, forward, reverse, longread, file("${data_id}_links.fasta"), val('spades_links') into files_links

    when:
    params.assembly in ['spades_links', 'all'] && !params.plasmid
    
    script:
    """
    echo ${longread} > longreads.txt
    perl ${LINKS} -f ${scaffolds} -s longreads.txt -b links
    mv links.scaffolds.fa ${data_id}_links.fasta
    """
}

process gapfiller{
   tag{data_id}
   
   input:
   set data_id, forward, reverse, longread, scaffolds, type from files_sspace .mix(files_links)
          
   output:
   set data_id, forward, reverse, longread, file("${data_id}_gapfiller.fasta"), type into assembly_gapfiller

   script:
   """
   echo 'Lib1GF bowtie '${forward} ${reverse} '500 0.5 FR' > gapfill.lib
   perl ${GAPFILLER} -l gapfill.lib -s ${scaffolds} -m 32 -t 10 -o 2 -r 0.7 -d 200 -n 10 -i 15 -g 0 -T 5 -b out
   mv out/out.gapfilled.final.fa ${data_id}_gapfiller.fasta
   """
}

/*
*  Canu assembler
*
*
*/
process canu_parameters {

    output: 
    file('canu_settings.txt') into canu_settings

    """
    echo \
    'genomeSize=$params.genome_size 
    minReadLength=1000
    maxMemory=$params.mem 
    maxThreads=$params.cpu' > canu_settings.txt
    """
}

process canu{
    tag{id}
    publishDir "${params.outDir}/${id}_${params.assembly}/canu", mode: 'copy'

    input:
    set id, sr1, sr2, lr from files_pre_canu
    file canu_settings
    
    output: 
    set id, sr1, sr2, lr, file("${id}.contigs.fasta"), val('canu') into files_unpolished_canu
    file("${id}.report")
    file("${id}_canu.gfa")

    when:
    params.assembly in ['canu','all']

    script:
    """
    ${CANU} -s ${canu_settings} -p ${id} -nanopore-raw ${lr}
    cp ${id}.unitigs.gfa ${id}_canu.gfa
    """
}

/*
*  miniasm assembler
*
*
*/
process miniasm{
    tag{id}
    publishDir "${params.outDir}/${id}_${params.assembly}/miniasm", mode: 'copy'

    input:
    set id, sr1, sr2, lr from files_pre_miniasm
    
    output:
    set id, sr1, sr2, lr, file("miniasm_assembly.fasta") into files_noconsensus
    file("${id}_miniasm_graph.gfa")

    when:
    params.assembly in ['miniasm', 'all']

    script:
    """
    ${MINIMAP2} -x ava-ont -t ${params.cpu} ${lr} ${lr} > ovlp.paf
    ${MINIASM} -f ${lr} ovlp.paf > ${id}_miniasm_graph.gfa
    awk '/^S/{print ">"\$2"\\n"\$3}' ${id}_miniasm_graph.gfa | fold > miniasm_assembly.fasta
    """
}

process racon {
    tag{id}
    publishDir "${params.outDir}/${id}_${params.assembly}/racon", mode: 'copy'
    
    input:
    set id, sr1, sr2, lr, assembly from files_noconsensus

    output:
    set id, sr1, sr2, lr, file("assembly_consensus.fasta"), val("miniasm") into files_unpolished_racon
    file("assembly_consensus.fasta")

    script:
    """
    ${MINIMAP2} -x map-ont -t ${params.cpu} ${assembly} ${lr} > assembly_map.paf
    ${RACON} -t ${params.cpu} ${lr} assembly_map.paf ${assembly} assembly_consensus.fasta
    """
}

process flye {
// Assembly step using Flye assembler
    errorStrategy 'ignore'
    tag{id}
    publishDir "${params.outDir}/${id}_${params.assembly}", mode: 'copy'

    input:
    set id, sr1, sr2, lr from files_pre_flye

    output:
    set id, sr1, sr2, lr, file("flye/scaffolds.fasta"), val('flye') into files_unpolished_flye
    file("flye/assembly_info.txt")
    file("flye/${id}_flye_graph.gfa")

    
    when:
    params.assembly in ['flye', 'all']

    script:
    """
    ${FLYE} --nano-raw ${lr} --out-dir flye \
    --genome-size ${params.genome_size} --threads ${params.cpu} -i 0
    cp flye/2-repeat/graph_final.gfa flye/${id}_flye_graph.gfa
    """
}

// Create channel for all unpolished files to be cleaned with Pilon
files_unpolished = Channel.create()
files_pilon = files_unpolished.mix(files_unpolished_canu, files_unpolished_racon, files_unpolished_flye)

/*
* Pilon polisher
*
*
*/
process pilon{
    tag{id}

    input:
    set id, sr1, sr2, lr, contigs, type from files_pilon

    output:
    set id, sr1, sr2, lr, file("after_polish.fasta"), type into assembly_pilon

    script:
    """
    ${BOWTIE2_BUILD} ${contigs} contigs_index.bt2 

    ${BOWTIE2} --local --very-sensitive-local -I 0 -X 2000 -x contigs_index.bt2 \
    -1 ${sr1} -2 ${sr2} | ${SAMTOOLS} sort -o alignments.bam -T reads.tmp 
    
    ${SAMTOOLS} index alignments.bam

    java -jar $PILON --genome ${contigs} --frags alignments.bam --changes \
    --output after_polish --fix all
    """
}

// Merge channel output from different assembly paths
assembly=Channel.create()
if (params.plasmid) {
    assembly_merged = assembly.mix(files_spades_plasmid, assembly_unicycler, assembly_pilon)
} else {
    assembly_merged = assembly.mix(assembly_gapfiller, assembly_unicycler, assembly_pilon)
}

/*
* Length filter trimming of contigs < 2000bp from the final assembly
* Creates a plot of contig lenghts in the assembly
*/
process length_filter {
    publishDir "${params.outDir}/${id}_${params.assembly}/", mode: 'copy'

    input:
    set id, sr1, sr2, lr, contigs, type from assembly_merged

    output:
    set id, type into complete_status
    file("${id}_${type}_final.fasta")
    
    // Uses python2 
    script:
    """
    #!/usr/bin/env python

    import sys
    import os
    import numpy as np
    from Bio import SeqIO
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    long_contigs = []
    input_handle=open('${contigs}', 'rU')
    output_handle=open('${id}_${type}_final.fasta', 'w')
    
    for index, record in enumerate(SeqIO.parse(input_handle, 'fasta')):
        if len(record.seq) >= ${min_contig_length}:
            record.id = "${id}." + str(index+1)
            record.description = "assembler=${type} length=" + str(len(record.seq))
            long_contigs.append(record)
    
    SeqIO.write(long_contigs, output_handle, "fasta")
    
    input_handle.close()
    output_handle.close()

    """

}


/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def helpMessage() {
  // Display help message
  // this.pipelineMessage()
  log.info "    Usage:"
  log.info "       nextflow run caspargross/hybridAssembly --samples <file.csv> --mode [mode] [options] "
  log.info "    --input <file.tsv>"
  log.info "       Specify a TSV file containing paths to sample files."
  log.info "    --mode ${validModes}"
  log.info "       Default: none, choose one or multiple modes to run the pipeline "
  log.info:   
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Cont Engine : " + workflow.containerEngine
  log.info "Out Dir     : " + params.outDir
  log.info "TSV file    : ${tsvFile}"
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Step        : " + step
  log.info "Containers"
  if (params.repository != "") log.info "  Repository   : " + params.repository
  if (params.containerPath != "") log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
  log.info "Reference files used:"
  log.info "  dbsnp       :\n\t" + referenceMap.dbsnp
  log.info "\t" + referenceMap.dbsnpIndex
  log.info "  genome      :\n\t" + referenceMap.genomeFile
  log.info "\t" + referenceMap.genomeDict
  log.info "\t" + referenceMap.genomeIndex
  log.info "  bwa indexes :\n\t" + referenceMap.bwaIndex.join(',\n\t')
  log.info "  intervals   :\n\t" + referenceMap.intervals
  log.info "  knownIndels :\n\t" + referenceMap.knownIndels.join(',\n\t')
  log.info "\t" + referenceMap.knownIndelsIndex.join(',\n\t')
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def pipelineMessage() {
  // Display hybridAssembly info  message
  log.info "hybridAssembly Pipeline ~  version ${workflow.manifest.version} - revision " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def startMessage() {
  // Display start message
  this.nextflowMessage()
  this.minimalInformationMessage()
}

workflow.onComplete {
  // Display complete message
  this.nextflowMessage()
  this.minimalInformationMessage()
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError {
  // Display error message
  this.nextflowMessage()
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}

def isMode = {
  // returns whether a given list of arguments contains at least one valid mode
it.any {modes.contains(it)}
}
