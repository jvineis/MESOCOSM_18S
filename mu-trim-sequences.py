#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import argparse

parser = argparse.ArgumentParser(description='''Read in a fastq file and trim sequnces from whatever end you request.  Its the best''')
parser.add_argument('--f', help = 'the fastq file')
parser.add_argument('--side', help = 'either F or L: eg. trim the First x number of sequences or the Last x number of sequences.')
parser.add_argument('--n', help = 'The number of sequences you wish to remove')
parser.add_argument('--out', help = 'The outfile to write your trimmed sequences')
args = parser.parse_args()

outfile = open(args.out, 'w')

if args.side == "F":
    for seq in SeqIO.parse(open(args.f, 'rU'),"fastq"):

        outfile.write("@"+seq.id+'\n'+str(seq.seq[0:int(args.n)])+'\n'+"+"+'\n'+str(seq.letter_annotations['phred_quality'])+'\n')
elif args.side == "R":
    for seq in SeqIO.parse(open(args.f, 'rU'), "fastq"):
        print(int(len(seq.seq)))
        outfile.write("@"+seq.id+'\n'+str(seq.seq[int(args.n):int(len(seq.seq))])+'\n'+"+"+'\n'+str(seq.letter_annotations['phred_quality'])+'\n')


