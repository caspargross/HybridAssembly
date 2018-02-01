#!/usr/bin/env/ nextflow

params.cpu = 10
params.mem = 100

//params.shortreadName = '*_*_L001_R{1,2}*.fastq.gz'
params.shortreadName = 'ESBL2048_S2_L001_R{1,2}*.fastq.gz'
params.shortreadFolder = '/mnt/projects/external/Microbiome/Citrobacter/samples/Illumina/2016-07-20-IMP'

//params.longreadName = '*/*.fastq'
params.longreadName = 'ESBL2048/*.fastq'
params.longreadFolder = '/mnt/projects/external/Microbiome/Citrobacter/samples/Nanopore'

params.outFolder = '/mnt/projects/external/Microbiome/Citrobacter/analysis'
params.reference = 

// Input channel for short read (Illumina) files
Channel
    .fromFilePairs("$params.shortreadFolder/$params.shortreadName", flat: true)
    .set{readPairs}

// Multiply read pairs for all subprocesses
readPairs.into {
    readPairs1
    readPairs2
    readPairs3
}

// Input channel for nanopore reads
Channel
    .fromPath("${params.longreadFolder}/${params.longreadName}")
    .set{longreads}


process assembly{
    tag{data_id}

    // Write spades output to folder
    publishDir "${params.outFolder}/${data_id}/spades", mode: 'copy'

    input:
    set data_id, file(forward), file(reverse) from readPairs1
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
    set data_id, file(forward), file(reverse) from readPairs2
    file(scaffolds) from spadesScaffolds

    output:
    file("${data_id}sspace.final.scaffolds.fasta") into sspaceScaffolds

    script:
    """
    echo 'Lib1 bowtie '${forward}, ${reverse} '500 0.5 FR' > sspace.lib
    perl ${SSPACE} -l sspace.lib -s ${scaffolds} -g 0 -x 0 -T ${params.cpu} -k 3 -a 0.7 -n 20 -z 500 -b ${data_id} -p 1
    """
}


process gapfiller{
   tag{data_id}
   
   publishDir "${params.outFolder}/${data_id}/gapfiller", mode: 'copy' 
   
   input:
   set data_id, file(forward), file(reverse) from readPairs3
   file(scaffolds) from sspaceScaffolds
   
   output:
   file("${data_id}.gapfilled.final.fa") into finalScaffolds

   script:
   """
   echo 'Lib1GF bowtie '${forward} ${reverse} '500 0.5 FR' > gapfill.lib
   perl ${GAPFILLER} -l gapfill.lib -s ${scaffolds} -m 32 -t 10 -o 2 -r 0.7 -d 200 -n 10 -i 15 -g 0 -T 5 -b ${data_id}
   """
}

process reference_alignment{
    tag{data_id}
    
    publishDir "${params.outFolder}/${data_id}/mummer", mode: 'copy'

    input:
    file(gapfilled) from finalScaffolds

    script:
    """
    ${MUMMER}/nucmer --mum -l 100 -c 150 -p={$data_id}mummer ${params.reference} ${gapfilled}
    ${MUMMER}/delta-filter -m {$data_id}mummer.delta > {$data_id}mummer.fdelta
    ${MUMMER}/delta-filter -q {$data_id}mummer.delta > {$data_id}mummer.qdelta
    ${MUMMER}/delta-filter -1 {$data_id}mummer.delta > {$data_id}mummer.1delta
    ${MUMMER}/show-coords -lrcT {$data_id}mummer.fdelta | sort -k13 -k1n -k2n > {$data_id}mummer.coords
    ${MUMMER}/show-tiling -c -l 1 -i 0 -V 0 {$data_id}mummer.fdelta > {$data_id}mummer.tiling
    ${MUMMER}/show-snps -ClrTH {$data_id}mummer.1delta > {$data_id}mummer.snps
    ${MUMMER}/mummerplot {$data_id}mummer.qdelta -R ${params.reference} -Q ${gapfilled} -p gapfilled {$data_id}mummer/out --filter --layout -postscript
    """ 
}
