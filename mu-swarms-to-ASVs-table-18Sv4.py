#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''This script takes the output returned from swarm-v3.1 -w flag, which is a list of lists '\n' containing the sequence ids for each swarm, and converts it into a '\n' ASV count matrix. Details of chimeric, taxonomic, swarm statistics for each ASV are also provided in the resulting table''')
parser.add_argument('-s' , help = "The swarm output of the -o flag: a list of ids for each swarm (therefore a list of lists)  -  commonly run like this: swarm -d 1 -f -t 10 -z pooled-samples-derep.fa -s pooled-samples-derep-stats.txt -w pooled-samples-node-representatives.fa -o pooled-samples-node-table.txt")
parser.add_argument('-o' , help = "The prefix for the two output files. One for the anvio output with ALL the information(PREFIX-swarm-matrix-for-anvio.txt), the othet a bare matrix of ASV ids and counts for tree building and flashweave analysis (PREFIX-swarm-matrix-for-R.txt)")
parser.add_argument('-l' , help = "A single column list of the fasta files from each of the dereplicated fastas that were pooled and then went into the swarm analysis. Each line in the file should contain the nmame of a fasta e.g. B3D5T1A_18S-primer-derep-for-swarm.fa")
parser.add_argument('-n' , help = "The swarm representative nodes which is the -w flag in swarm")
parser.add_argument('-c', help = "The output of running vseach for chimeric detection")
parser.add_argument('-t', help = "The output of running the taxonomic assingments using vsearch followed by combine-node-hits-with-tax-strings.py")
parser.add_argument('-st', help = "The stats file produced buy the -s flag in SWARM. usually calld  pooled-samples-derep-stats.txt")
#parser.add_argument('-min', help = "The threshold value to keep a swarm,  anything observed fewer than this many times will be removed")
args = parser.parse_args()

outfile_for_anvio = open(args.o+"-swarm-matrix-for-anvio.txt", 'w')
outfile_for_R = open(args.o+"-swarm-matrix-for-R.txt", 'w')
swarms = open(args.s, 'r')
samples = open(args.l, 'r')
node_seqs = open(args.n, 'r')
chimeras = open(args.c, 'r')
taxa = open(args.t, 'r')
stats = open(args.st, 'r')

####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#############################

## We want to store more information in our output, including the taxononmy and details for each swarm. To do this, we will create several other dictionaries keeping the swarm as the key for all, so that in the end, we can loop through all of the ids to write the data contained in each unique dict.

## A dictionary to hold the node ids, the size of the node, and the sequence
node_dict = {}
for seq in SeqIO.parse(node_seqs, "fasta"):
    x = seq.id.split(';')
    node_dict[x[0]] = x

## A dictionary for the uchime chimeric seqs    
chimera_dict = {}
for c in chimeras:
    x = c.strip().split('\t')
    chimera_dict[x[1].split(';')[0]] = x[0:len(x)]

## A dictionary for the taxonomy
taxa_dict = {}
for t in taxa:
    x = t.strip().split('\t')
    taxa_dict[x[0]] = x[0:len(x)] ## This is what the contents of the stats file looks like "31af4e133d2c982e7197d55e6d3e546cf012a2f1	49097	57.9	Eukaryota|Opisthokonta|Metazoa|Arthropoda|Hexapoda|Insecta|Ceratocombus|Ceratocombus_australiensis    AY252300.1.1835_U". The swarm id, abundance of the swarm, percent identity of the hit, the taxonomy of te hit, and the ref id. 

## A dictionary for the stats - these include the 
stats_dict = {}
for s in stats:
    x = s.strip().split('\t')
    stats_dict[x[2]] = x[0:len(x)] ## The stats file contains 0: number of unique amplicons in the cluster, 1:total abundance of amplicons in the cluster, 2: label of the intial seed, 3: abundance of the initial seed, 4:number of amplicons with 1 as an abundance in the cluster, 5: max number of iterations before the cluster reached its limit, 6: steps in the path. 

swarm_dict = {} ## A dictionary to hold the swarm data
for line in swarms: ## for each swarm
    x = line.strip().split(' ') ## separate all the sequences and their size information
    swarm_list = []
    if len(x) > 5 and int(x[0].split("=")[1]) > 9: ## separate all the sequences and their size information for the SWARMs that have more than 5 representative sequences with the most abundant sequence found at least ten times. 
        for element in x:
            swarm_list.append(element.split(";")[0])
        swarm_dict[x[0].split(";")[0]] = swarm_list ## make the first sequence in the swarm the key and all other sequences/sizes as the values.. values look like this bb3f24804f300e8478f59542ce4465f7f7b36ddc;size=1
 
## Now we are going to bild a dictionary with each of the samples as the keys and the deflines with abundance information as the values. That way we can look through this dictionary on a sample basis and record the abundance for each of the swarms.
samples_dict = {} ## here is the holder for our new dictionary
for fa in samples: ## loop through each of the fasta files
    rec_id_list = [] ## create a list to hold all of the seq ids
    x = fa.strip() ## collect the fasta file name
    for rec in SeqIO.parse(x,"fasta"): ## open the fasta file
        fa_rec = rec.id ## collect the defline of the individual sequence
        rec_id_list.append(fa_rec) ## append this to the list of individual sequence ids. This will also contain the abundance information d0ffa3d4f07d4702cff1b39ca02089a8607554a6;size=1;
    samples_dict[x.split(".")[0]] = rec_id_list ## create the dictionary entry for the sample as the key B1D2T1_18S-primer-derep-for-swarm and each sequence id as values. 

####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#############################

def synthesize_swarm_information (swarm, sample): ## This function will return the number of observations for a particular sample and swarm combination. We will need to run it for each swarm across each sample
    count = 0 # The base count for each of the swarm sample combinations.
    for seq in samples_dict[sample]: # For each sequence in a sample
        x = seq.split(";")[0] ## get the swarm name
        y = seq.split(";")[1].split("=")[1] ## get the count information
        if x in swarm_dict[swarm]: ## if the swarm name is found in that set of sequences in the swarm, add the "y" value (abundance of the sequence found in the sample) to the 
            count += int(y)
    return count

####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#############################

synthesized_swarm_count_dict = {} ## Create an empty dictionary to store the swarm as a key and the counts for each sample as the value
for swarm in swarm_dict.keys(): ## loop through each of the swarms
    sample_counts = [] ## set up an empty list to store the counts of the swarm for each sample
    for sample in samples_dict.keys(): ## startlooping through each sample contained in the dictionary that we made containing the sample name and list of amplicon ids.
        x = synthesize_swarm_information(swarm, sample) ## use the sample and the swarm ID
        sample_counts.append(x) ## append the sum of all tabulated amplicon counts for each swarm in each sample
    synthesized_swarm_count_dict[swarm] = sample_counts ## store this information in the dictionary.

####################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#############################

## Now we want to write all the SWARM information to a single output file. Lets open it, then write the headers and finally complie all of the information for each SWARM. The "outfile" variable was defined above on line 17. "This will need to be altered depending on your samples and the order that they are listed in your "-l" file. They should be in the order that apears if you list them in your terminal.. eg ls *derep-for-swarm.fa

outfile_for_anvio.write("uid"+'\t'+"SWARM"+'\t'+"unique_amplicons"+'\t'+"total_abundance"+'\t'+"singletons_in_cluster"+'\t'+"B1D2T1"+'\t'+"B1D3T1"+'\t'+"B1D4T1A"+'\t'+"B1D5T1A"+'\t'+"B1D5T2A"+'\t'+"B1D6T1B"+'\t'+"B2D2T1"+'\t'+"B2D3T1"+'\t'+"B2D4T1A"+'\t'+"B2D5T1A"+'\t'+"B3D2T1"+'\t'+"B3D3T1"+'\t'+"B3D4T1A"+'\t'+"B3D5T1A"+'\t'+"B3D5T2C"+'\t'+"Kingdom"+'\t'+"division"+'\t'+"phylum"+'\t'+"class"+'\t'+"order"+'\t'+"family"+'\t'+"genus"+'\t'+"species"+'\t'+"uid"+'\t'+"percent_id"+'\t'+"chmeric"+'\n')
count = 0
for key in swarm_dict.keys():
    outfile_for_anvio.write("s_"+str(count)+'\t'+key+'\t'+## The swarm id
                             stats_dict[key][0]+'\t'+ ## The number of unique amplicons in the swarm
                             stats_dict[key][1]+'\t'+ ## The total abundance of all amplicon counts in the swarm
                             stats_dict[key][3]+'\t'+ ## The unmber of singletons in the swarm cluster
                             '\t'.join(str(x) for x in synthesized_swarm_count_dict[key])+'\t'+ ## All the swarm count data
                             '\t'.join(taxa_dict[key][3:11])+'\t'+ ## All the taxonomic classification
                             taxa_dict[key][1]+'\t'+ ## The uid for the PR2 database
                             taxa_dict[key][2]+'\t'+ ## The percent identity of the hit
                             chimera_dict[key][17]+'\n') ## yes or no on the chimeric status
    count += 1

outfile_for_R.write("uid"+'\t'+"B1D2T1"+'\t'+"B1D3T1"+'\t'+"B1D4T1A"+'\t'+"B1D5T1A"+'\t'+"B1D5T2A"+'\t'+"B1D6T1B"+'\t'+"B2D2T1"+'\t'+"B2D3T1"+'\t'+"B2D4T1A"+'\t'+"B2D5T1A"+'\t'+"B3D2T1"+'\t'+"B3D3T1"+'\t'+"B3D4T1A"+'\t'+"B3D5T1A"+'\t'+"B3D5T2C"+'\n')
count = 0
for key in swarm_dict.keys():
    outfile_for_R.write("s_"+str(count)+'\t'+
                         '\t'.join(str(x) for x in synthesized_swarm_count_dict[key])+'\n') ## All the swarm count data
                
    count += 1


        # The stats file contains 0: number of unique amplicons in the cluster, 1:total abundance of amplicons in the cluster, 2: label of the intial seed, 3: abundance of the initial seed, 4:number of amplicons with 1 as an abundance in the cluster, 5: max number of iterations before the cluster reached its limit, 6: steps in the path.
