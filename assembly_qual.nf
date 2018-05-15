/**
*
*   HYBRID ASSEMBLY QUALITY CHECK  
*
*   This pipeline runs several quality checks on the 
*   assembled genomes. 
*   Caspar Gross 2018
* 
**/

params.reference = false


// Arguments:
samples = Channel.fromPath(params.pathFile)
    .ifEmpty {error "Cannot find file with path locations in ${params.pathFile}"}
    .splitCsv(header: true)
    

// Creates channels for all possible assemblers per sample
process buildChannel {
   tag{id}

   input:
   set id, sr1, sr2, lr, ref from samples 

   output:
   set id, file("*.fasta"), ref into samples_listed
   set id, file("*.fasta"), sr1, sr2 into samples_listed_bwa
    
   script:
   """ 
   cp ${params.outFolder}${id}*/*.fasta .
   """

}

// Split Channels for different assemblers
genomes = samples_listed
    .transpose()
    .map{it ->  [it[0], (it[1] =~ /_(.+)_/)[0][1], it[1], it[2]]}
//    .view()
    .into{genomes_mummer; genomes_bwa; genomes_dnadiff; genomes_quast; genomes_bwa}

genomes_bwa = samples_listed_bwa
    .transpose()
    .map{it ->  [it[0], (it[1] =~ /_(.+)_/)[0][1], it[1], it[2], it[3] ]}
//    .view()

// Mummerplot against reference
process mummer_alignment {
    tag{id}

    input:
    set id, type, genome, ref from genomes_mummer
    
    output:
    set id, type, file("${id}.coords"), ref into mummer_coords 

	when:
	params.reference
    
	script:
    """
    ${MUMMER}/nucmer --mum -l 100 -c 150 -p ${id} ${ref} ${genome}
    ${MUMMER}/delta-filter -m ${id}.delta > ${id}.delta.m
    ${MUMMER}/show-coords -c ${id}.delta.m > ${id}.coords

    """

}

// Calculate DNADiff
process mummer_dnadiff {
    tag{id}
    publishDir "${params.analysisFolder}/${id}/MummerPlot/", mode: 'copy'

    input:
    set id, type, genome, ref from genomes_dnadiff
    
    output:
    file("${id}_${type}_dnadiff.report")

	when:
    params.reference

    script:
    """
    ${MUMMER}/dnadiff -p ${id}_${type}_dnadiff ${ref} ${genome}

    """

}


process mummer_plot {
    tag{id}
    publishDir "${params.analysisFolder}/${id}/MummerPlot/", mode: 'copy'

    input:
    set id, type, coords, ref from mummer_coords

    output:
    file "*" optional true


    script:
    """
    Rscript ${baseDir}/scripts/mummerCoordsDotPlotly.R -i ${coords} -o ${id}_${type}_mplot -p 9 -s -t -m 100 -q 100 -l
    """
}


// Move quast analysis from hybridAssembly to over here
process quast{
    tag{id}
    publishDir "${params.analysisFolder}/${id}/quast/", mode: 'copy'

    input:
    set id, type, genome, ref from genomes_quast

    output:
    file "quast_${type}/*"

    script:
    if (params.reference) 
        """
        ${QUAST} ${genome} -t ${params.cpu} -o quast_${type} -R ${ref} --labels ${id}_${type} --min-identity 85
        """
    else
        """
        ${QUAST} ${genome} -t ${params.cpu} -o quast_${type} --labels ${id}_${type} --min-identity 85
        """

}

process bwa_aln {
    publishDir "${params.analysisFolder}/${id}/aligned_reads/", mode: 'copy'

    input:
    set id, type, genome, sr1, sr2 from genomes_bwa

    output:
    set id, type, genome, file('*.bam') into alignment_bwa
    file('*.bai')

    script:
    """
    ${BWA} index ${genome}
    ${BWA} aln ${genome} ${sr1} > R1.sai
    ${BWA} aln ${genome} ${sr2} > R2.sai
    ${BWA} sampe ${genome} R1.sai R2.sai ${sr1} ${sr2} \
    | ${SAMTOOLS} sort -o ${genome.baseName}_${type}.bam
    ${SAMTOOLS} index ${genome.baseName}_${type}.bam ${genome.baseName}_${type}.bai
    """

}


process bwa_call {
    tag{id}
    publishDir "${params.analysisFolder}/${id}/variants/", mode: 'copy'

    input:
    set id, type, genome, aln from alignment_bwa

    output:
    file("${id}_${type}_snp.bcf")
    set id, type, genome, aln, file("${id}_${type}_nsnp.txt") into snp_number

    script:
    """
    ${BCFTOOLS} mpileup -Ou -f ${genome} ${aln} | \
    ${BCFTOOLS} call -mv -Ob -o ${id}_${type}_snp.bcf

    ${BCFTOOLS} view ${id}_${type}_snp.bcf | wc -l > ${id}_${type}_nsnp.txt
    """
}

