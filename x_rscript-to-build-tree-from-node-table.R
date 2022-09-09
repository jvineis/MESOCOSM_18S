#!/usr/bin R

library("ggplot2")
library("ape")
library("vegan")

# Read in the counts table produced by convert-node-hits-to-tax-node-table.py - you may need to edit the file name depending on what you named your file. Most everything else will work with no changes at all.
dat = read.table("x_SWARMS-and-tax-for-anvio-counts.txt", header = TRUE, row.names = 1, sep = '\t')

# This is the number of quality filtered reads, we may want to filter based on that someday.
#reads = read.table("~/Dropbox/PARASITES-BBW-JHV/JENNA-OCEANS/x_Quality-and-primer-filtered-reads-per-sample.txt", header = TRUE, row.names = 1, sep = '\t')

# This creates a table of relative abundance for each ASV within each sample using the sum of the columns 
relative_abund = scale(dat, center = FALSE, scale = colSums(dat))
# Fix the header of the relative abundance table so that column 1 has a header "OTU"
ra = data.frame("OTU"=rownames(relative_abund), relative_abund)

# Write the relative abundance to a table with the rows as ASVs and the columns as samples
write.table(ra,"x_SWARMS-and-tax-for-anvio-relative-abundance.txt", sep = '\t', row.names = FALSE)

# Calculate the distance among ASVs based on the relative abundance and then write the tree
asv_dist = vegdist(relative_abund, "bray")
asv_clust = hclust(asv_dist)
write.tree(as.phylo(asv_clust),"x_SWARMS-and-tax-for-anvio-relative-abundance.tre")

#Calculate the distance among samples based on the relative abundance and then write the tree
sample_dist = vegdist(t(relative_abund), "bray")
sample_clust = hclust(sample_dist)
write.tree(as.phylo(sample_clust),"x_SWARMS-and-tax-for-anvio-relative-abundance-samples.tre")