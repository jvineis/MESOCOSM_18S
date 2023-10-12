#!/usr/bin/env python

import sys

microfile = open(sys.argv[1], 'r') # these are your microscopy cell count estimates "cell-counts-for-biovolume-scripts.txt"
biomassfile = open(sys.argv[2], 'r') # this is the biomass database "diatom-biovolume-for-scripts.txt"
outfile = open(sys.argv[3], 'w') # the outfile "
outfile1 = open(sys.argv[4],'w') # outfile of the average per species biovolume "diatom-database-speceies-averages.txt"
# a dictionary of the microscopy data
micro_dict = {} 

# a dictionary for the biomass database
biomass_dict = {}


for i in microfile:
    x = i.strip().split('\t')
    if x[0] == "Index":
        outfile.write('\t'.join(x[0:len(x)])+'\t'+"biovolume"+'\n')
    else:
        micro_dict[x[0]] = x[0:len(x)]

for i in biomassfile:
    x = i.strip().split('\t')
    biomass_dict[x[0]] = x[0:len(x)]


# We want to get an estimate of the biovolume for each observation in the dataset and write the table
# to a file that we can import into R etc.

# To do this we will use the species classification if available, otherwise we will need to use the mean
# biovolume estimate for genus.

# Lets get a dictionary of the average values for genus and species

# an empty list for the genera in the biomass database
genera = []
# populate the genera list
for key in biomass_dict.keys():
    print(key, biomass_dict[key])
    genera.append(biomass_dict[key][1])

# an empty list for the species in the biomass database
species = []
for key in biomass_dict.keys():
    species.append('_'.join([biomass_dict[key][1],biomass_dict[key][2]]))

# make a unique list of the genera

unique_genera = set(genera)
unique_gen_dict = {}
print("there are %d genera in your biomass database" %(len(unique_genera)))

# make a unique list of the species

unique_species = set(species)
unique_species_dict = {}
print("there are %d species in your biomass database" %(len(unique_species)))


# A dictionary to hold the values for the average biomass of the unique GENERA in the biomass database. 
for gen in unique_genera:
    vol = []
    for key in biomass_dict.keys():
        if biomass_dict[key][1] == gen:
            vol.append(float(biomass_dict[key][3]))
    unique_gen_dict[gen] = sum(vol)/len(vol)

# print the genus and the mean biovolume
#for key in unique_gen_dict.keys():
#    print(key,unique_gen_dict[key])

    
# A dictionary to hold the values for the average biomass of the unique SPECIES in the biomass database
for spe in unique_species:
    svol = []
    for key in biomass_dict.keys():
        if biomass_dict[key][1] == spe.split('_')[0] and biomass_dict[key][2] == spe.split('_')[1]:
#            print(biomass_dict[key][3])
            svol.append(float(biomass_dict[key][3]))
    unique_species_dict[spe] = sum(svol)/len(svol)


# print the species and the mean biovolume
for key in unique_species_dict.keys():
    outfile1.write(key+'\t'+str(unique_species_dict[key])+'\n')
#    print(key, unique_species_dict[key])

# lets loop throught the microscopy data to find the biomass estimate for each of the species and if we don't find the species
# we will use the average biomass for the genus.
micro_genus_not_found = [] # Keep track of the genera that we don't find in eithe the unique genus or species lists
micro_found = [] # to hold the keys found so I don't look for the genus in both the unique genus and species lists 
for key in micro_dict.keys():
    mgenus = micro_dict[key][2] # microscopy genus
    mspecies = micro_dict[key][3] # microscopy species
    for skey in unique_species_dict.keys():
        bgenus = skey.split('_')[0]
        bspecies = skey.split('_')[1]
        if mgenus == bgenus and mspecies == bspecies:
#            print("found species",mgenus, mspecies, unique_species_dict[skey])  
            outfile.write('\t'.join(micro_dict[key])+'\t'+ str(unique_species_dict[skey])+'\n')
            micro_found.append(key)

            
    for gkey in unique_gen_dict.keys():
        if mgenus == gkey and key not in micro_found:
 #           print("mgenus", gkey, unique_gen_dict[gkey])
            outfile.write('\t'.join(micro_dict[key])+'\t'+ str(unique_gen_dict[mgenus])+'\n')
            micro_found.append(key)
    if mgenus not in unique_gen_dict.keys():
        micro_genus_not_found.append(key)

for item in micro_genus_not_found:
#    print(micro_dict[item])
    outfile.write('\t'.join(micro_dict[item])+'\t'+"unknown"+'\n')
