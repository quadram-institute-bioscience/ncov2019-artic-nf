process fastqMergeLanes {
  publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: 'copy'

  tag { sampleName }

  input:
  tuple val(sampleName), file(forward), file(reverse)
  
  output:
  tuple val(sampleName), file("${sampleName}_R1.fastq.gz"), file("${sampleName}_R2.fastq.gz")
  
  script:
  """
  zcat $forward | gzip > ${sampleName}_R1.fastq.gz
  zcat $reverse | gzip > ${sampleName}_R2.fastq.gz
  """
}


process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    // publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 2

    input:
    tuple val(sampleName), path(forward), path(reverse)

    output:
    tuple val(sampleName), path("*_val_1.fq.gz"), path("*_val_2.fq.gz") optional true

    script:
    """
    if [[ \$(gunzip -c ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      exit 0
    else
      trim_galore --paired $forward $reverse
    fi
    """
}

process indexReference {
    /**
    * Indexes reference fasta file in the scheme repo using bwa.
    */

    tag { ref }

    input:
        path(ref)

    output:
        tuple path('ref.fa'), path('ref.fa.*')

    script:
        """
        ln -s ${ref} ref.fa
        bwa index ref.fa
        """
}

process readMapping {
    /**
    * Maps trimmed paired fastq using BWA (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */
    tag { sampleName }

    label 'largecpu'

    // publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted.bam", mode: 'copy'

    input:
        tuple val(sampleName), path(forward), path(reverse), path(ref), path("*")

    output:
        tuple val(sampleName), path("${sampleName}.sorted.bam")

    script:
        """
        bwa mem -t ${task.cpus} ${ref} ${forward} ${reverse} | \
        samtools sort -o ${sampleName}.sorted.bam
        """
}

process trimPrimerSequences {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.bam", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(bedfile)

    output:
    tuple val(sampleName), path("${sampleName}.mapped.bam"), emit: mapped
    tuple val(sampleName), path("${sampleName}.mapped.primertrimmed.sorted.bam" ), emit: ptrim

    script:
    if (params.allowNoprimer){
        ivarCmd = "ivar trim -e"
    } else {
        ivarCmd = "ivar trim"
    }
    if (params.schemeVersion != 'V4') {
        """
        samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
        samtools index ${sampleName}.mapped.bam
        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out
        samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        """
    } else {
        """
        samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
        samtools index ${sampleName}.mapped.bam
        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out
        #samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        mask_dgy_bases.py -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        """
    }
        
}

process callVariants {

    tag { sampleName }
    
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.tsv", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(ref)

    output:
    tuple val(sampleName), path("${sampleName}.variants.tsv"), emit: variants

    script:
        """
        samtools mpileup -A -d 0 --reference ${ref} -B -Q 0 ${bam} |\
        ivar variants -r ${ref} -m ${params.ivarMinDepth} -p ${sampleName}.variants -q ${params.ivarMinVariantQuality} -t ${params.ivarMinFreqThreshold}
        """
}

process makeConsensus {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.fa", mode: 'copy'

    input:
        tuple val(sampleName), path(bam)

    output:
        tuple val(sampleName), path("${sampleName}.primertrimmed.consensus.fa")

    script:
        """
        samtools mpileup -aa -A -B -d ${params.mpileupDepth} -Q0 ${bam} | \
        ivar consensus -t ${params.ivarFreqThreshold} -m ${params.ivarMinDepth} \
        -n N -p ${sampleName}.primertrimmed.consensus
        """
}

//Copied from https://github.com/jts/ncov2019-artic-nf/blob/master/modules/illumina.nf#L174
process callConsensusFreebayes {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}_consensus", pattern: "${sampleName}.consensus.fasta", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}_variant", pattern: "${sampleName}.variants.norm.vcf", mode: 'copy'
    
    cpus 8
    memory '32.GB'

    input:
    tuple val(sampleName), path(bam), path(ref)

    output:
    tuple val(sampleName), path("${sampleName}.consensus.fasta"), emit: fasta
    tuple val(sampleName), path("${sampleName}.variants.norm.vcf"), emit: vcf

    script:
        """
        # the sed is to fix the header until a release is made with this fix
        # https://github.com/freebayes/freebayes/pull/549
        freebayes -p 1 \
                  -f ${ref} \
                  -F 0.2 \
                  -C 1 \
                  --pooled-continuous \
                  --min-coverage ${params.ivarMinDepth} \
                  --gvcf --gvcf-dont-use-chunk true ${bam} | sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > ${sampleName}.gvcf
        # make depth mask, split variants into ambiguous/consensus
        # NB: this has to happen before bcftools norm or else the depth mask misses any bases exposed during normalization
        process_gvcf.py -d ${params.ivarMinDepth} \
                        -l ${params.ivarMinFreqThreshold} \
                        -u ${params.ivarFreqThreshold} \
                        -m ${sampleName}.mask.txt \
                        -v ${sampleName}.variants.vcf \
                        -c ${sampleName}.consensus.vcf ${sampleName}.gvcf
        # normalize variant records into canonical VCF representation
        for v in "variants" "consensus"; do
            bcftools norm -f ${ref} ${sampleName}.\$v.vcf > ${sampleName}.\$v.norm.vcf
        done
        # split the consensus sites file into a set that should be IUPAC codes and all other bases, using the ConsensusTag in the VCF
        for vt in "ambiguous" "fixed"; do
            cat ${sampleName}.consensus.norm.vcf | awk -v vartag=ConsensusTag=\$vt '\$0 ~ /^#/ || \$0 ~ vartag' > ${sampleName}.\$vt.norm.vcf
            bgzip -f ${sampleName}.\$vt.norm.vcf
            tabix -f -p vcf ${sampleName}.\$vt.norm.vcf.gz
        done
        # apply ambiguous variants first using IUPAC codes. this variant set cannot contain indels or the subsequent step will break
        bcftools consensus -f ${ref} -I ${sampleName}.ambiguous.norm.vcf.gz > ${sampleName}.ambiguous.fasta
        # apply remaninng variants, including indels
        bcftools consensus -f ${sampleName}.ambiguous.fasta -m ${sampleName}.mask.txt ${sampleName}.fixed.norm.vcf.gz | sed s/MN908947.3/${sampleName}/ > ${sampleName}.consensus.fasta
        """
}


process compareConsensus {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*", mode: 'copy'

    tag {sample_id}

    cpus 4

    errorStrategy 'ignore'
    
    shell = ['/bin/bash', '-uo', 'pipefail']

    input:
    tuple val(sample_id), path(ivar_fasta), path(freebayes_fasta), path(ref)

    output:
    path("${sample_id}_ivar_vs_freebayes.vcf"), optional: true

    script:
    """
    seqkit replace -p '(.*)' -r '${sample_id}_ivar' ${ivar_fasta}> ivar
    seqkit replace -p '(.*)' -r '${sample_id}_freebayes' ${freebayes_fasta} > freebayes
    
    cat ivar freebayes > all.fasta
    nextalign -r ${ref} -i all.fasta -o ${sample_id}.all.aln
    snp-sites -v -o ${sample_id}_ivar_vs_freebayes.vcf ${sample_id}.all.aln >> /dev/null 2>&1
    
    if [[ ! -f ${sample_id}_ivar_vs_freebayes.vcf ]];then
        touch ${sample_id}_ivar_vs_freebayes.vcf
    fi
    """
}


process cramToFastq {
    /**
    * Converts CRAM to fastq (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to CRAM, to FastQ (http://www.htslib.org/doc/samtools.html)
    * @input
    * @output
    */

    input:
        tuple val(sampleName), file(cram)

    output:
        tuple val(sampleName), path("${sampleName}_1.fastq.gz"), path("${sampleName}_2.fastq.gz")

    script:
        """
        samtools collate -u ${cram} -o tmp.bam
        samtools fastq -1 ${sampleName}_1.fastq.gz -2 ${sampleName}_2.fastq.gz tmp.bam
        rm tmp.bam
        """
}

