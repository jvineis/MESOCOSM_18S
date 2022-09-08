#!/usr/bin/env python

import sys
outfile =open(sys.argv[2],'w')
for line in open(sys.argv[1]):
    x = line.strip().split('\t')
    outfile.write(x[0]+'\t'+x[1][0:8]+'\n')

