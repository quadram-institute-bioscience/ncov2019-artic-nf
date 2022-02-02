
process typeVariants {

    tag { sampleName }
    label 'genotyping'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/variants", pattern: "${sampleName}.variants.csv", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/vcf", pattern: "${sampleName}.csq.vcf", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/typing", pattern: "${sampleName}.typing.csv", mode: 'copy'

    input:
    tuple val(sampleName), path(variants), path(gff), path(ref), path(yaml)

    output:
    path "${sampleName}.variants.csv", optional: true, emit: variants_csv
    path "${sampleName}.typing.csv", optional: true, emit: typing_csv
    path "${sampleName}.csq.vcf", emit: csq_vcf

    script:
    if( params.illumina )
        """
        type_vcf.py \
        -i ${sampleName} \
        -y ${yaml} \
        -ov ${sampleName}.csq.vcf \
        -ot ${sampleName}.typing.csv \
        -os ${sampleName}.variants.csv \
        -dp ${params.csqDpThreshold} \
        -af ${params.csqAfThreshold} \
        -t ${variants} \
        ${gff} ${ref}
        """
    else
        """
        type_vcf.py \
        -i ${sampleName} \
        -y ${yaml} \
        -ov ${sampleName}.csq.vcf \
        -ot ${sampleName}.typing.csv \
        -os ${sampleName}.variants.csv \
        -dp ${params.csqDpThreshold} \
        -af ${params.csqAfThreshold} \
        -v ${variants} \
        ${gff} ${ref}
        """
}

process mergeTypingCSVs {

    tag { params.prefix }

    publishDir "${params.outdir}", pattern: "${params.prefix}.typing_summary.csv", mode: 'copy'
    publishDir "${params.outdir}", pattern: "${params.prefix}.variant_summary.csv", mode: 'copy'
    label 'genotyping'
    input:
    tuple path('typing/*'), path('variant/*')

    output:
    path "${params.prefix}.typing_summary.csv", emit: typing_summary_csv
    path "${params.prefix}.variant_summary.csv", emit: variant_summary_csv

    script:
    """
    #!/usr/bin/env python3
    import glob
    import csv

    dirs = ['typing', 'variant']

    for dir in dirs:
        globstring = dir + '/*.csv'
        files = glob.glob(globstring)

        header_written = False
        out_fn = "${params.prefix}." +dir+ '_summary.csv'
        with open(out_fn, 'w') as outfile:
            for fl in files:
                with open(fl, 'r' ) as csvfile:
                    csvreader = csv.DictReader(csvfile)
                    for row in csvreader:
                        if not header_written:
                            writer = csv.DictWriter(outfile, fieldnames=list(row.keys()))
                            writer.writeheader()
                            header_written = True

                        writer.writerow(row)

    """
}

process send_discord {
    // publishDir "${params.outdir}/pangolin", mode: "copy"
    
    conda '/home/ubuntu/miniconda3/envs/monitor'

    input: 
    path(typing_summary_csv)
    output:
    
    script:
    """
    parse_variant_typing.py --name ${params.prefix} ${typing_summary_csv}
    """
 }


 process collect_fasta {
    // publishDir "${params.output}", mode: "copy"
    tag {"💦"}
    executor 'slurm'
    
    // clusterOptions '-w slurm-node-112cpus-161tb'
    
    conda '/home/ubuntu/miniconda3'
    
    cpus 4
       
    input:
    
    file(fasta)
    
    output:
    file("combined.tmp")

    script:
    """
    cat *.fa* | seqkit replace -p '(Consensus_|_S[0-9]{1,5}.primertrimmed.consensus_threshold_0.75_quality_20|/ARTIC/medaka MN908947.3)' - | seqkit grep -r -s -p "A|C|T|G|N" - > combined.tmp
    """
}


process nextclade {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"

    tag {"🦠 Running 🥕"}

    executor 'slurm'
       
    // container '/beegfs/singularity/nextclade/nextclade.sif'
    
    cpus 64
    
    memory '128.GB'
    
    input:
    
    path("combined.fasta")
    
    output:
    path("${task.process.replaceAll(":","_")}.csv"), emit: csv
    path("${params.name}.*")

    script:
    """
    nextclade run -j ${task.cpus} \
    -i combined.fasta \
    --input-root-seq=${params.reference} \
    --input-gene-map=${params.genemap} \
    --genes=E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
    --input-qc-config=${params.qc} \
    --input-tree=${params.tree} \
    --input-pcr-primers=${params.primers} \
    --output-dir=nextclade \
    --output-basename=${params.name} \
    -c ${task.process.replaceAll(":","_")}.csv
    mv nextclade/*.* .
    """
}

//  nextclade --input-fasta=sequences.fasta --input-root-seq=reference.fasta --genes=E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --input-gene-map=genemap.gff --input-tree=tree.json --input-qc-config=qc.json --input-pcr-primers=primers.csv --output-json=output/nextclade.json --output-csv=output/nextclade.csv --output-tsv=output/nextclade.tsv --output-tree=output/nextclade.auspice.json --output-dir=output/ --output-basename=nextclade

process pangolin {
    publishDir "${params.outdir}", mode: "copy"

    tag {"Running 😀"}

    executor 'slurm'
    
    // clusterOptions '-w slurm-node-112cpus-161tb'
    
    container '/beegfs/singularity/pangolin/pangolin.sif'
    
    cpus 32
    
    // memory '800G'
    
    input:
    
    file("combined.fasta")
    
    output:
    file("pangolin")

    script:
    """
    pangolin combined.fasta --no-temp --max-ambig 0.7 -t ${task.cpus} -o pangolin
    """
}
