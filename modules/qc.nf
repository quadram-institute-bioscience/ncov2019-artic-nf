process makeQCCSV {
    tag { sampleName }
    
    label 'largecpu'
    
    publishDir "${params.outdir}/qc_plots", pattern: "${sampleName}.depth.png", mode: 'copy'
    
    
    input:
    tuple val(sampleName), path(bam), path(fasta), path(ref)

    output:
    path "${params.prefix}.${sampleName}.qc.csv", emit: csv
    path "${sampleName}.depth.png"

    script:
    if ( params.illumina ) {
       qcSetting = "--illumina"
    } else {
       qcSetting = "--nanopore"
    }

    """
    samtools index ${bam}
    qc.py ${qcSetting} --outfile ${params.prefix}.${sampleName}.qc.csv --sample ${sampleName} --ref ${ref} --bam ${bam} --fasta ${fasta}
    """
}


process writeQCSummaryCSV {
    tag { params.prefix }

    input:
    val lines

    exec:
    file("${params.outdir}/${params.prefix}.qc.csv").withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
}

process remove_human_dna {
    publishDir "${params.outdir}/cleanedFastq", mode: 'copy'
    
    tag "$sample_name"
    
    label 'largemem'
    
    input:
    path(params.kraken2_human_db)
    tuple val(sample_name), path(reads)

    output:
    if (params.se) {
        tuple val(sample_name), path("*.fastq")    
    } else {
        tuple val(sample_name), path("*_1.fastq"), path("*_2.fastq")
    }
    path("*.tsv")
    
    script:
    def kraken_input = params.se ? "--unclassified-out ${sample_name}.fastq --classified-out ${sample_name}_human.fastq" : "--paired-end --unclassified-out ${sample_name}_#.fastq --classified-out {}_human_#.fastq"
    """
    kraken2 \
    --threads $task.cpus \
    --db ${params.kraken2_human_db} \
    --gzip-compress \
    ${kraken_input} \
    --output read_reports_${sample_name}.tsv \
    --memory-mapping \
    $reads
    """ 
}