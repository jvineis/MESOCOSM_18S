#!/usr/bin/env python

import sys

fw_dict = {}
tax_dict = {}

## open the flashweave edges file with the first two lines removed, should be called "x_SWARM-counts-for-flashweave-with-MAGs-output.edgelist"
## and the lines should look something like this
## s_12	s_13	0.7492479085922241
count = 0
for line in open(sys.argv[1],'r'):
    x = line.strip().split()
    if "ASV" not in x[0]:
        fw_dict[count] = [x[0],x[1],x[2]]
    count += 1
    
## open the "x_SWARM-contingency-table.txt-swarm-matrix-for-anvio.txt" file
for line in open(sys.argv[2], 'r'):
    x = line.strip().split()   
    tax_dict[x[0]] = x[0:len(x)]


outfile = open(sys.argv[3],'w')
outfile.write("pair"+'\t'+"iud_1"+'\t'+"iud_2"+'\t'+"interaction"+'\t'+"tax_1"+'\t'+"tax_2"+'\t'+"per_id1"+'\t'+"per_id2"+'\n')

for key in fw_dict.keys():
    taxa1 = ";".join(tax_dict[fw_dict[key][0]][20:28])
    taxa2 = ";".join(tax_dict[fw_dict[key][1]][20:28])
    per_id1 = tax_dict[fw_dict[key][0]][29]
    per_id2 = tax_dict[fw_dict[key][1]][29]
    assoc = fw_dict[key][2]
    print(taxa1, taxa2, per_id1, per_id2)
    outfile.write(str(key) +'\t'+fw_dict[key][0]+'\t'+fw_dict[key][1]+'\t'+assoc+'\t'+taxa1+'\t'+taxa2+'\t'+per_id1+'\t'+per_id2+'\n')

