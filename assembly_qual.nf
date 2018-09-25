/**

    HYBRID ASSEMBLY QUALITY CHECK.

    This pipeline runs several quality checks on the 
    finished hybrid assembly and can be run with or
    without reference genome. It uses three tools to 
    create reference statistics, check the assembly
    using mapped paired end reads and functional 
    validation using core genes with checkM.

    written by
    Caspar Gross

    ------------------------------
    Process overview:
    - Quast - Calculate sumamry statistics with quast
    - checkM - Run full lineage analysis with checkM
    - checkM plot - Plot additional checkM figures
    - REAPR

*/

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



