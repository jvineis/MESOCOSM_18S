#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import argparse

parser = argparse.ArgumentParser(description='''Read in a pair of fastq files and pull out sequence\
 pairs if they contain defline ids provided in a given file''')
parser.add_argument('--R1', help = 'The R1 (read1) of paired reads in fastq format')
parser.add_argument('--R2', help = 'The R2 (read2) of paired reads in fastq format')
parser.add_argument('--outfile1', help = 'The name of the R1 output file containing your target seqs')
parser.add_argument('--outfile2', help = 'The name of the R2 output file containing your target seqs')
parser.add_argument('--infile', help = 'A txt file with a single column of defline ids')
args = parser.parse_args()

headers = []
infile = open(args.infile, 'rU')
for line in infile:
    x =line.strip()
    headers.append(x)

read1 = open(args.R1, 'rU')
read2 = open(args.R2, 'rU')
print("Hi microbiologist, Im searching for your %d sequecnces and writing them to read 1 file [ %s] and read2 file [ %s]. Its gonna be awesome :)" %(len(headers), str(args.outfile1), str(args.outfile2)))

def select_fastq(fastqfile, output):
    outfile = open(output, 'w')
    for seq_record in SeqIO.parse(fastqfile, "fastq"):
        print(seq_record.id)
        if seq_record.id in headers:
            outfile.write(seq_record.format("fasta"))
    return outfile

select_fastq(read1,args.outfile1)
select_fastq(read2,args.outfile2)
