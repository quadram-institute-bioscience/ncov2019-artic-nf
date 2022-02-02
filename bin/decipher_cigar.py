import pysam
import pathlib

def decipher(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam:
        cigar = [x[0] for x in read.cigar]
        print(cigar)

if __name__ == "__main__":
    bam = "/beegfs/Covid-19_Seq/result.illumina.20200701/qc_climb_upload/NORW-20200701/BLANK20-1_S171/BLANK20-1_S171.mapped.bam"
    decipher(bam)