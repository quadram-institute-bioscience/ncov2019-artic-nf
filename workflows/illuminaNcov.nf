#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// import modules
include {articDownloadScheme } from '../modules/artic.nf' 
include {readTrimming} from '../modules/illumina.nf' 
include {indexReference} from '../modules/illumina.nf'
include {readMapping} from '../modules/illumina.nf' 
include {trimPrimerSequences} from '../modules/illumina.nf' 
include {callVariants} from '../modules/illumina.nf'
include {callVariants as callVariants_1X} from '../modules/illumina.nf' params(ivarMinDepth: '1', ivarMinVariantQuality: 20, ivarMinFreqThreshold: '0.25', outdir: params.outdir)
include {makeConsensus} from '../modules/illumina.nf'
include {makeConsensus as makeConsensus_1X} from '../modules/illumina.nf' params(ivarMinDepth: '1', mpileupDepth: '1000', ivarFreqThreshold: '0.75', outdir: params.outdir)
include {cramToFastq} from '../modules/illumina.nf'
include {fastqMergeLanes} from '../modules/illumina.nf'
include {callConsensusFreebayes} from '../modules/illumina.nf'
include {compareConsensus} from '../modules/illumina.nf'
include {makeQCCSV} from '../modules/qc.nf'
include {writeQCSummaryCSV} from '../modules/qc.nf'

include {collateSamples} from '../modules/upload.nf'

// import subworkflows
include {CLIMBrsync} from './upload.nf'
include {Genotyping} from './typing.nf'
include {Genotyping as Genotyping_1X} from './typing.nf'

workflow fastqMergeFourLanes {
  take:
    ch_filePairs
  
  main:
    // Group 4 lane elements into 1 element in the channel
    ch_filePairs
            .map {
            it -> [it[0].replaceAll(~/\_L00[1,2,3,4]/,""), it[1], it[2]]
            }
            .groupTuple(by:0)
            .set { ch_reads_four_lanes }

    fastqMergeLanes(ch_reads_four_lanes)
  
  emit: fastqMergeLanes.out
}

workflow prepareReferenceFiles {
    // Get reference fasta
    if (params.ref) {
      Channel.fromPath(params.ref)
              .set{ ch_refFasta }
    } else {
      articDownloadScheme()
      articDownloadScheme.out.reffasta
                          .set{ ch_refFasta }
    }


    /* Either get BWA aux files from reference 
       location or make them fresh */
    
    if (params.ref) {
      // Check if all BWA aux files exist, if not, make them
      bwaAuxFiles = []
      refPath = new File(params.ref).getCanonicalPath()
      new File(refPath).getParentFile().eachFileMatch( ~/.*.bwt|.*.pac|.*.ann|.*.amb|.*.sa/) { bwaAuxFiles << it }
     
      if ( bwaAuxFiles.size() == 5 ) {
        Channel.fromPath( bwaAuxFiles )
               .set{ ch_bwaAuxFiles }

        ch_refFasta.combine(ch_bwaAuxFiles.collect().toList())
                   .set{ ch_preparedRef }
      } else {
        indexReference(ch_refFasta)
        indexReference.out
                      .set{ ch_preparedRef }
      }
    } else {
      indexReference(ch_refFasta)
      indexReference.out
                    .set{ ch_preparedRef }
    }
  
    /* If bedfile is supplied, use that,
       if not, get it from ARTIC github repo */ 
 
    if (params.bed ) {
      Channel.fromPath(params.bed)
             .set{ ch_bedFile }

    } else {
      articDownloadScheme.out.bed
                         .set{ ch_bedFile }
    }

    emit:
      bwaindex = ch_preparedRef
      bedfile = ch_bedFile
      reffasta = ch_refFasta
}


workflow sequenceAnalysis {
    take:
      ch_filePairs
      ch_preparedRef
      ch_bedFile

    main:
      if (params.readTrimming) {
        readTrimming(ch_filePairs)
        readTrimming.out.set{ ch_preprocess }
      } else {
        ch_filePairs.set {ch_preprocess }
      }
      

      readMapping(ch_preprocess.combine(ch_preparedRef))

      trimPrimerSequences(readMapping.out.combine(ch_bedFile))
 
      callVariants(trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] }))     
      
      callConsensusFreebayes(trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] }))

      makeConsensus(trimPrimerSequences.out.ptrim)

      //Compare consensus
      compareConsensus(makeConsensus.out.join(callConsensusFreebayes.out.fasta).combine(ch_preparedRef.map{ it[0] }))

      //1X call
      callVariants_1X(trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] }))
      makeConsensus_1X(trimPrimerSequences.out.ptrim)

      makeQCCSV(trimPrimerSequences.out.ptrim.join(makeConsensus.out, by: 0)
                                   .combine(ch_preparedRef.map{ it[0] }))

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
    		       }
                       .set { qc }

      writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

      if (params.get_all) {
          makeQCCSV.out.csv.splitCsv()
                           .unique()
                           .set {collate_ch}
      } else {
        qc.pass.set {collate_ch}
      }
      // Pass
      // makeQCCSV.out.csv.splitCsv()
      //                  .unique().subscribe{ println "$it" } 
      // collate_ch.subscribe{ println "$it" }

      collateSamples(collate_ch.map{ it[0] }
                           .join(makeConsensus.out, by: 0)
                           .join(trimPrimerSequences.out.mapped))

    emit:
      qc_pass = collateSamples.out
      variants = callVariants.out.variants
      variants_1X = callVariants_1X.out.variants
}

workflow ncovIllumina {
    take:
      ch_filePairs

    main:
      // Build or download fasta, index and bedfile as required
      prepareReferenceFiles()
      //If fastqs are separated in four lane files, need to merge first
      if (params.fourLanes) {
        fastqMergeFourLanes(ch_filePairs)
        fastqMergeFourLanes.out.set { ch_filePairs_new }
      } else {
        ch_filePairs.set { ch_filePairs_new }
      }
      // Actually do analysis
      sequenceAnalysis(ch_filePairs_new, prepareReferenceFiles.out.bwaindex, prepareReferenceFiles.out.bedfile)
 
      // Do some typing if we have the correct files
    //   if ( params.gff ) {
    //       Channel.fromPath("${params.gff}")
    //              .set{ ch_refGff }

    //       Channel.fromPath("${params.yaml}")
    //              .set{ ch_typingYaml }

    //       Genotyping(sequenceAnalysis.out.variants, ch_refGff, prepareReferenceFiles.out.reffasta, ch_typingYaml) 
    //       // Genotyping_1X(sequenceAnalysis.out.variants_1X, ch_refGff, prepareReferenceFiles.out.reffasta, ch_typingYaml) 
    //   }

      // Upload files to CLIMB
      if ( params.upload ) {
        
        Channel.fromPath("${params.CLIMBkey}")
               .set{ ch_CLIMBkey }
      
        CLIMBrsync(sequenceAnalysis.out.qc_pass, ch_CLIMBkey )
      }

}

workflow ncovIlluminaCram {
    take:
      ch_cramFiles
    main:
      // Convert CRAM to fastq
      cramToFastq(ch_cramFiles)

      // Run standard pipeline
      ncovIllumina(cramToFastq.out)
}

