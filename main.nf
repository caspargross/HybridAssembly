#!/usr/bin/env/ nextflow
params.assembly = ''
params.return_all = false 


/* Pipeline paths:

spades_sspace
spades_links
canu
unicycler

*/


//inputFiles
files = Channel.fromPath(params.pathFile)
    .ifEmpty {error "Cannot find file with path locations in ${params.files}"}
    .splitCsv(header: true)
    .view()


// Validate assembly protocol choice:
if (!(params.assembly in ['spades_sspace', 'spades_links', 'canu', 'unicycler' ])){
    exit 1, "Invalid assembly protocol: (${params.assembly}) \nMust be one of the following:\n    'canu'\n    'spades_links'\n    'spades_sspace'\n    'unicycler'"
}



// Trim adapter sequences on long read nanopore files
process porechop {
    tag{id}
        
    input:
    set id, sr1, sr2, lr from files
    
    output:
    set id, sr1, sr2, file('lr_porechop.fastq') into files_porechop

    // Join multiple longread files if possible 
    script:
    """
    cat ${lr} > nanoreads.fastq
    $PORECHOP -i nanoreads.fastq -t ${params.cpu} -o lr_porechop.fastq
    """
}


// Quality filter long reads
process filtlong {
    tag{id}

    input: 
    set id, sr1, sr2, lr from files_porechop
    
    output:
    set id, sr1, sr2, file("lr_filtlong.fastq") into files_filtlong, files_fastqc
    
    script:
    """
    $FILTLONG -1 ${sr1} -2 ${sr2} \
    --min_length 1000 \
    --keep_percent 90 \
    --target_bases  100000000 \
    ${lr} > lr_filtlong.fastq
    """
    // Expected genome size: 5.3Mbp --> Limit to 100Mbp for approx 20x coverage
}


// Create FASTQC quality check on short reads
process fastqc{
    tag{id}
    
    publishDir "${params.outFolder}/${id}_${params.assembly}/fastQC/", mode: 'copy'

    input: 
    set id, sr1, sr2, lr from files_fastqc

    output:
    files("fastqc/*")

    script: 
    """
    mkdir -p fastqc
    ${FASTQC} sr1 sr2 -o fastqc
    """
}

if (params.assembly == 'unicycler') {

    process unicycler{
    tag{id}
    publishDir "${params.outFolder}/${id}_${params.assembly}/unicycler", mode: 'copy'   
   
    input:
    set id, sr1, sr2, lr from files_filtlong

    output:
    set id, sr1, sr2, lr, file("${id}/assembly.fasta") into files_assembled
    file("${id}/*")

    script:
    """ 
    python3 /mnt/users/ahgrosc1/tools/Unicycler/unicycler-runner.py \
    -1 ${sr1} -2 ${sr2} -l ${lr}\
    -o ${id} -t ${params.cpu}\
    --spades_path ${SPADES}\
    --racon_path ${RACON}
    --pilon_path ${PILON}\
    --bowtie2_build_path ${BOWTIE2_BUILD}\
    --bowtie2_path ${BOWTIE2}\
    --samtools_path ${SAMTOOLS}
    """
    }
}


if (params.assembly == "spades_sspace" || params.assembly == "spades_links") {
    
    // Run SPADes assembly
    process spades{
        tag{data_id}

        input:
        set data_id, forward, reverse, longread from files_filtlong  

        output:
        set data_id, forward, reverse, longread, file("spades/scaffolds.fasta") into files_spades
        file("spades/contigs.fasta")

        script:
        """
        ${SPADES} -t ${params.cpu} -m ${params.mem} \
        --phred-offset 33 --careful \
        --pe1-1 ${forward} \
        --pe1-2 ${reverse} \
        --nanopore ${longread} \
        -o spades
        """
    }
}



// Scaffold using SSPACE
if(params.assembly == 'spades_sspace'){

    process sspace_scaffolding{
        tag{data_id}

        input:
        set data_id, forward, reverse, longread, scaffolds from files_spades  

        output:
        set data_id, forward, reverse, longread, file("sspace/scaffolds.fasta") into files_sspace 

        script:
        """
        perl ${SSPACE} -c ${scaffolds} -p ${longread} -b sspace -t ${params.cpu}
        """
    }
    
    process gapfiller{
       tag{data_id}
       
       input:
       set data_id, forward, reverse, longread, scaffolds from files_sspace
              
       output:
       set data_id, forward, reverse, longread, file("${data_id}_gapfiller.fasta") into files_assembled

       script:
       """
       echo 'Lib1GF bowtie '${forward} ${reverse} '500 0.5 FR' > gapfill.lib
       perl ${GAPFILLER} -l gapfill.lib -s ${scaffolds} -m 32 -t 10 -o 2 -r 0.7 -d 200 -n 10 -i 15 -g 0 -T 5 -b out
       mv out/out.gapfilled.final.fa ${data_id}_gapfiller.fasta
       """
    }
}


if(params.assembly == 'spades_links'){
    process links_scaffolding{
        tag{data_id}
        
        input:
        set data_id, forward, reverse, longread, scaffolds from files_spades
        
        output:
        set data_id, forward, reverse, longread, file("${data_id}_links.fasta") into files_assembled

        script:
        """
        echo ${longread} > longreads.txt
        perl ${LINKS} -f ${scaffolds} -s longreads.txt -b links
        mv links.scaffolds.fa ${data_id}_links.fasta
        """
    }

}

if (params.assembly == 'canu') {
    
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

        input:
        set id, sr1, sr2, lr from files_filtlong
        file canu_settings
        
        output: 
        set id, sr1, sr2, lr, file("${id}.contigs.fasta") into files_canu
        file("${id}.report")

        script:
        """
        ${CANU} -s ${canu_settings} -p ${id} -nanopore-raw ${lr}
        """
    }
}


if (params.assembly == 'canu'){
    
    process pilon{
        tag{id}

        input:
        set id, sr1, sr2, lr, contigs from files_canu

        output:
        set id, sr1, sr2, lr, file("after_polish.fasta") into files_assembled

        script:
        """
        ${BOWTIE2} --local --very-sensitive-local -I 0 -X 2000 -x ${contigs} \
        -1 ${sr1} -2 ${sr2} | samtools sort -o alignments.bam -T reads.tmp 
        
        samtools index alignments.bam

        java -jar $PILON --genome ${contigs} --frags alignments.bam --changes \
        --output after_polish --fix all
        """

    }
}



process contig_length {
    // This script filters contigs by length (standard 200bp)
    // and creates a plot of the read length distribution
    // and writes the final fasta file to the disk
    publishDir "${params.outFolder}/${id}_${params.assembly}/", mode: 'copy'

    input:
    set id, sr1, sr2, lr, contigs from files_assembled

    output:
    set id, file("${id}_final.fasta") into analysis_mummer, analysis_quast
    file("${id}_contig_lengthDist.pdf")
    
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
    output_handle=open('${id}_final.fasta', 'w')
    
    for record in SeqIO.parse(input_handle, 'fasta'):
        if len(record.seq) >= 2000 :
            long_contigs.append(record)
    
    
    SeqIO.write(long_contigs, output_handle, "fasta")
    
    input_handle.close()
    output_handle.close()

    lengths = map(len, long_contigs)
    
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(lengths, np.zeros_like(lengths)+1, 'bs')
    
    title = 'Contig lengths for ${id}'
    ax.set_title(title)
    ax.set_xscale('symlog')
    ax.set_xlabel("Nucleotides")
    ax.tick_params(
        axis = 'y',
        which = 'both',
        left = 'off',
        right = 'off',
        labelleft = 'off'
    )
    ax.xaxis.grid(False)
    fig.savefig("${id}_contig_lengthDist.pdf", format='pdf')
    """

}


process align_reference{
    
    publishDir "${params.outFolder}/${data_id}_${params.assembly}/mummer/", mode: 'copy'
    
    input: 
    set data_id, contigs from analysis_mummer

    output:
    file("${data_id}_mummerplot.ps")
    file("${data_id}.dnadiff.report")
    

    script:
    """
    ${MUMMER}/dnadiff -p ${data_id}.dnadiff ${params.reference} ${contigs}
    
    ${MUMMER}/nucmer --mum -l 100 -c 150 -p ${data_id} ${params.reference} ${contigs}
    ${MUMMER}/delta-filter -m ${data_id}.delta > ${data_id}.fdelta
    ${MUMMER}/delta-filter -q ${data_id}.delta > ${data_id}.qdelta
    ${MUMMER}/delta-filter -1 ${data_id}.delta > ${data_id}.1delta
    ${MUMMER}/show-coords -lrcT ${data_id}.fdelta | sort -k13 -k1n -k2n > ${data_id}.coords
    ${MUMMER}/show-tiling -c -l 1 -i 0 -V 0 ${data_id}.fdelta > ${data_id}.tiling
    ${MUMMER}/show-snps -ClrTH ${data_id}.1delta > ${data_id}.snps
    ${MUMMER}/mummerplot ${data_id}.qdelta -R ${params.reference} -Q ${contigs} -p ${data_id}_mummerplot --filter --layout -postscript

    """
}

process quast{
    tag{id}
    
    publishDir "${params.outFolder}/${data_id}_${params.assembly}/quast/", mode: 'copy'

    input: 
    set id, assembly from analysis_quast

    output:
    file("quast/*")

    script:
    """
    ${QUAST} ${assembly} -t ${params.cpu} -o quast -R ${params.reference}
    """
}
