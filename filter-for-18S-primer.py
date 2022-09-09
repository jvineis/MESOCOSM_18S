#!/usr/bin/env python

from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description='''find the primer sequences for 18S used in the DeVargas paper in each sequence and create a new fasta containing only sequences that have both the fwd and reverse primer''')
parser.add_argument('--i', help='a merged fasta file')
parser.add_argument('--o', help='the name you want to give to your beautiful fasta that contains only the best sequences that contain the primer')
args=parser.parse_args()

outfile = open(args.o, 'w')

def filter_left_primers(records, Fprimer,Fprimer1,Fprimer2,Fprimer3):
    left_side_found_ids = []
    for record in records: # look for the fwd and rev primers in the sequence
        if record.seq.lower().startswith(Fprimer) or record.seq.lower().startswith(Fprimer1) or record.seq.startswith(Fprimer2) or record.seq.startswith(Fprimer3):
            left_side_found_ids.append(record)
    print("found this many good left", len(left_side_found_ids))
    return(left_side_found_ids)


def filter_right_primers(seq_list, Fprimer, Fprimer1, Fprimer2,Fprimer3):
    both_primers_found_ids = []
    for record in seq_list:
        if record.seq.lower().endswith(Fprimer) or record.seq.lower().endswith(Fprimer1) or record.seq.lower().endswith(Fprimer2) or record.seq.endswith(Fprimer3):
            both_primers_found_ids.append(record)
    print("found this many right", len(both_primers_found_ids))
    return(both_primers_found_ids)
                                          
                               
original_reads = SeqIO.parse(args.i, "fasta")
left_trim = filter_left_primers(original_reads, "gtgccagcagccgcggtaa","gtgccagccgccgcggtaa","gtgtcagcagccgcggtaa","gtgtcagccgccgcggtaa")
both_found = filter_right_primers(left_trim, "gcgaaagcatttgccaa","gcgaaagcattcgccaa", "gcgaaagcatttaccaa", "gcgaaagcattcaccaa")
SeqIO.write(both_found, args.o, "fasta")
