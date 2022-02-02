// ARTIC processes

process articDownloadScheme{
    tag params.schemeRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scheme", mode: "copy"

    output:
    path "scheme/**/${params.schemeVersion}/*.reference.fasta" , emit: reffasta
    path "scheme/**/${params.schemeVersion}/${params.scheme}.bed" , emit: bed
    path "scheme" , emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} scheme
    """
}

process articGuppyPlex {
    tag { samplePrefix + "-" + fastqDir }

    label 'largemem'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.prefix}*.fastq*", mode: "copy"

    input:
    path(fastqDir)

    output:
    path "${params.prefix}*.fastq*", emit: fastq

    script:
    """
    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --prefix ${params.prefix} \
    --directory ${fastqDir}
    """
}

process barcodeToCOG {
    tag { fastq.baseName }

    errorStrategy 'ignore'

    input:
    path fastq
    
    output:
    path "*.fastq", emit: fastq
    
    script:
   
    """
    hostname > hostname
    rename.py ${params.barcode_lookup} ${fastq}
    """
}


process articMinIONMedaka {
    tag { sampleName }
    
    validExitStatus 0,20
    
    errorStrategy 'ignore'

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy", 
    saveAs: {filename -> "${sampleName}/${filename}"}

    input:
    tuple file(fastq), file(schemeRepo)

    output:
    file("${sampleName}.*")
    
    tuple val(sampleName), file("${sampleName}.primertrimmed.rg.sorted.bam"), optional: true, emit: ptrim
    tuple val(sampleName), file("${sampleName}.sorted.bam"), optional: true, emit: mapped
    tuple val(sampleName), file("${sampleName}.consensus.fasta"), optional: true, emit: consensus_fasta

    script:
    // Make an identifier from the fastq filename
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')
    coronahit = params.coronahit ? "--coronahit" : ""
    // Configure artic minion pipeline
    minionRunConfigBuilder = []

    if ( params.normalise ) {
    minionRunConfigBuilder.add("--normalise ${params.normalise}")
    }
    
    if ( params.bwa ) {
    minionRunConfigBuilder.add("--bwa")
    } else {
    minionRunConfigBuilder.add("--minimap2")
    }

    minionFinalConfig = minionRunConfigBuilder.join(" ")

    """
    artic minion --medaka \
    --medaka-model r941_min_high_g360 \
    ${coronahit} \
    ${minionFinalConfig} \
    --threads ${task.cpus} \
    --scheme-directory ${schemeRepo}/${params.schemeDir} \
    --read-file ${fastq} \
    ${params.scheme}/${params.schemeVersion} \
    ${sampleName} 2>&1 > ${sampleName}_minion.log
    """
}
// samtools bam2fq ${sampleName}.sorted.bam | seqkit fq2fa | seqkit fx2tab -l -n > ${sampleName}.tab

process articMinIONNanopolish {
    tag { sampleName }
    errorStrategy 'ignore'
    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy",
    saveAs: {filename -> "${sampleName}/${filename}"}

    input:
    tuple file(fastq), file(schemeRepo), file(fast5Pass), file(seqSummary)

    output:
    file("${sampleName}.*")
    
    tuple val(sampleName), file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple val(sampleName), file("${sampleName}.sorted.bam"), emit: mapped
    tuple val(sampleName), file("${sampleName}.consensus.fasta"), emit: consensus_fasta

    script:
    // Make an identifier from the fastq filename
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')
    coronahit = params.coronahit ? "--coronahit" : ""
    // Configure artic minion pipeline
    minionRunConfigBuilder = []

    if ( params.normalise ) {
    minionRunConfigBuilder.add("--normalise ${params.normalise}")
    }
    
    if ( params.bwa ) {
    minionRunConfigBuilder.add("--bwa")
    } else {
    minionRunConfigBuilder.add("--minimap2")
    }

    minionFinalConfig = minionRunConfigBuilder.join(" ")

    """
    artic minion ${minionFinalConfig} \
    ${coronahit} \
    --threads ${task.cpus} \
    --scheme-directory ${schemeRepo}/${params.schemeDir} \
    --read-file ${fastq} \
    --fast5-directory ${fast5Pass} \
    --sequencing-summary ${seqSummary} \
    ${params.scheme}/${params.schemeVersion} \
    ${sampleName}
    """
}

process articRemoveUnmappedReads {
    tag { sampleName }

    cpus 1

    input:
    tuple val(sampleName), path(bamfile)

    output:
    tuple val(sampleName), path("${sampleName}.mapped.sorted.bam")

    script:
    """
    samtools view -F4 -o ${sampleName}.mapped.sorted.bam ${bamfile} 
    """
}

