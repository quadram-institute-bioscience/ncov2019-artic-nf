process bam2fastq {
    publishDir "${params.outdir}/ENA/${params.prefix}", pattern: "${sampleName}", mode: 'copy'
    tag "$sample_id"
    
    cpus 4

    input:
    // set sample_id, file(bam) from ch_mapped_bam
    path(sample_folder)
    
    output:
    // set sample_id, file("*.fastq.gz") into ch_fastq
    path("*.fastq.gz")

    script:
    def output_reads = params.pe ? "-1 ${sample_id}_R1.fq.gz -2 ${sample_id}_R2.fq.gz -n" : "-0 ${sample_id}.fastq.gz"
    def repair = params.repair ? "repair.sh -Xmx8g in1=${sample_id}_R1.fq.gz in2=${sample_id}_R2.fq.gz out1=${sample_id}_R1.fastq.gz out2=${sample_id}_R2.fastq.gz" : ""
    // samtools view -@ $task.cpus -b -f 0x2 $bam > proper.bam 
    // $repair
    """
    samtools fastq -@ $task.cpus $output_reads $bam
    """
}


process collateSamples {
    tag { sampleName }

    publishDir "${params.outdir}/qc_climb_upload/${params.lib_name}", pattern: "${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(fasta)

    output:
    path("${sampleName}")

    script:
    """
    mkdir ${sampleName}
    mv ${bam} ${fasta} ${sampleName}
    """
}

process prepareUploadDirectory {
    tag { params.prefix }

    input:
    path("${params.prefix}/*")

    output:
    path("${params.prefix}")

    script:
    """
    echo "dummy" > dummyfile
    """
}


process uploadToCLIMB {
    tag { params.prefix }

    input:
    tuple(path(sshkey), path(uploadDir))

    output:

    script:
    """
    rsync -Lav -e "ssh -i ${sshkey} -l ${params.CLIMBUser}" ${uploadDir} ${params.CLIMBHostname}:upload/
    """
}

