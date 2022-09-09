#!/usr/bin/env R

dat = read.table("x_SWARM-counts-for-flashweave.tsv", header = TRUE, row.names = 1, sep  = "\t")
dat1 = dat[rowSums(dat[])>50,]
write.table(t(dat1),"x_SWARM-counts-for-flashweave-transposed-min50.csv", sep = ",")


