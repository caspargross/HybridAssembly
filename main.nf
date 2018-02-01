#!/usr/bin/env/ nextflow

params.cpu = 10
params.mem = 100

//params.shortreadName = '*_*_L001_R{1,2}*.fastq.gz'
params.shortreadName = 'ESBL2048_S2_L001_R{1,2}*.fastq.gz'
params.shortreadFolder = '/mnt/projects/external/Microbiome/Citrobacter/samples/Illumina/2016-07-20-IMP'

//params.longreadName = '*/*.fastq'
params.longreadName = 'ESBL2048/*.fastq'
params.longreadFolder = '/mnt/projects/external/Microbiome/Citrobacter/samples/Nanopore'

params.outFolder = 'mnt/projects/external/Microbiome/Citrobacter/analysis'

// Input channel for short read (Illumina) files
println("Short read input files:")
Channel
    .fromFilePairs("$params.shortreadFolder/$params.shortreadName", flat: true)
    .set{readPairs}

// Input channel for nanopore reads
println("Long read input files:")
Channel
    .fromPath("${params.longreadFolder}/${params.longreadName}")
    .set{longreads}


process assembly{
    tag{data_id}

    // Write spades output to folder
    publishDir "${params.outFolder}/${data_id}/spades", mode: 'copy'

    input:
    set data_id, file(forward), file(reverse) from readPairs
    file(longread) from longreads

    output:
    set data_id, file("${data_id}SpadesScaffolds.fasta") into spadesScaffolds

    script:
    """
    ${SPADES} -t ${params.cpu} -m ${params.mem} \
    --phred-offset 33 --careful \
    --pe1-1 ${forward} \
    --pe1-2 ${reverse} \
    --nanopore ${longread} \
    -o spades
    mv spades/scaffolds.fasta ${data_id}SpadesScaffolds.fasta
    """
}


process sspace_scaffolding{
    tag{data_id}

    input:
    set data_id, file(forward), file(reverse) from readPairs
    file(scaffolds) from spadesScaffolds

    output:
    file("${data_id}sspace.final.scaffolds.fasta") into sspaceScaffolds

    script:
    """
    echo 'Lib1 bowtie '${forward}, ${reverse} '500 0.5 FR' > sspace.lib
    perl $SSPACE -l sspace.lib -s scaffolds -g 0 -x 0 -T ${params.cpu} -k 3 -a 0.7 -n 20 -z 500 -b ${data_id} -p 1
    """
}


process gapfiller{
   tag{data_id}
   
   publishDir "${params.outFolder}/${data_id}/gapfiller", mode: 'copy' 
   
   input:
   set data_id, file(forward), file(reverse) from readPairs
   file(scaffolds) from sspaceScaffolds
   
   output:
   file("${data_id}.gapfilled.final.fa") into finalScaffolds

   script:
   """
   echo 'Lib1GF bowtie '${forward} ${reverse} '500 0.5 FR' > gapfill.lib
   perl $GAPFILLER -l gapfill.lib -s scaffolds -m 32 -t 10 -o 2 -r 0.7 -d 200 -n 10 -i 15 -g 0 -T 5 -b ${data_id}
   """
}


