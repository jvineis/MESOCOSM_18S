# MESOCOSM_18S
This git contains the code and steps to process the v4 region of the 18S rRNA gene and other supplemental data associated with the manuscript Co-occurrence and successional patterns among diatoms, dinoflagellates, and potential parasites in a coastal upwelling experiment. The steps for amplicon library construction are outlined in the publication links to this git and detailed steps are contained in the "16S v4 18s Amplicon Protocol for Illumina Sequencing.doc", located in the code tab of this git repository. Most of the code below contains information for the SLURM scheduler as well. Note to self, you are working here /scratch/gpfs/WARD/JOE/MOSS_BLOOM/18S_20220502_Data.

### Code required for this tutorial and links to the conda install or more specific install instructions. 
    illumina-utilities: conda install -c bioconda illumina-utils
    anvio-7.1: https://anvio.org/install/
    vsearch: conda install -c bioconda vsearch
    oligotyping:install the latest stable release (artistic mode) https://merenlab.org/2014/08/16/installing-the-oligotyping-pipeline/

### Quality filtering: 

#### Following the sequencing of samples, the libraries require demultiplexing, unless this was handled by the core facility. The demultiplexed raw data for the MLML v4 sequences can be found here (link to SRA accession numbers) and the unknown fastq files used for the demultiplexing steps are available upon request. The sample to barcode file used below can be found in the code tab in this repository. 

#### First I had to trim the barcodes in the barcode_to_sample.txt file because the length of the barcodes in the sequencing run were eitht bases, but the barcodes in the illumina samplesheet are 10 bases. Here is a simple script to do this which is also found in the code tab of this repository.

    python x_trim-barcodes.py barcode_to_sample.txt barcode_to_sample_8index.txt
    
#### Second, I ran the demultiplexing step. 

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=00:30:00
    #SBATCH --mem=80Gb

    iu-demultiplex -s barcode_to_sample_8index.txt --r1 4491__Ward-MOSS-18s-amplicons-and-metagenomics_pool_A-for-516-cycles-000000000-DG7K7_1_Read_1_passed_filter.fastq --r2 4491__Ward-MOSS-18s-amplicons-and-metagenomics_pool_A-for-516-cycles-000000000-DG7K7_1_Read_4_passed_filter.fastq --index 4491__Ward-MOSS-18s-amplicons-and-metagenomics_pool_A-for-516-cycles-000000000-DG7K7_1_Read_2_Index_Read_passed_filter.fastq -o DEMULTIPLEX
    
    
#### Third, merge the demultiplexed sequences

##### before you can merge, you need .ini files which can be made like this, but will depend on your sample names.. Be careful here to make sure you create the "x_18S-samples.txt" name right so it will work properly with your array command. 

    iu-gen-configs 00_DEMULTIPLEXING_18S --r1-prefix ^GTG[C,T]CAGC[A,C]GCCGCGGTAA --r2-prefix ^TTGG[C,T][A,G]AATGCTTTCGC
    ls *.ini | sed 's/\.ini//g' | grep 18S> x_18S-samples.txt
    
##### now you can merge the sequnces for each sample and use vsearch to dereplicate each sample separately. 

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=00:30:00
    #SBATCH --mem=80Gb
    #SBATCH --array=1-15

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_18S-samples.txt)
    iu-merge-pairs ${SAMPLE}.ini
    vsearch --quiet --derep_fulllength ${SAMPLE}_MERGED --sizeout --fasta_width 0 --relabel_sha1 --output ${SAMPLE}-primer-derep.fa
    
   

#### THE SWARM CLUSTERING WITH PRIMERS TRIMMED FROM AMPLICONS - change the time and memory depending on which step you are at. SWARM will likely take the longest, but likely not more than 1hr.
    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=20Gb
    #SBATCH --time=00:30:00

    ### These steps replicate the work here
    ### https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline#global-dereplication-clustering-and-chimera-detection
    ### The cat and vsearch steps should be run with the conda vsearch envionment 
    ### The sarm should be run with this environment "conda activate /home/jv2474/.conda/envs/swarm-v3.1"
    ### The python script should be run with the bioconda environment.

    ## 1. Concatenate the merged and filtered sequences
    cat *primer-derep-for-swarm.fa > pooled-samples.fa
    ## 2. Dereplicaete the concatenated sequences
    vsearch --derep_fulllength pooled-samples.fa --fasta_width 0 --sizeout --sizein --output pooled-samples-derep.fa
    ## 3. Cluster the sequences
    #swarm -d 1 -f -t 40 -z pooled-samples-derep.fa -i pooled-samples-derep-struct.txt -s pooled-samples-derep-stats.txt -w pooled-samples-node-representatives.fa -o pooled-samples-node-table.txt
    ## 4. Sort the clustered node representatives
    #vsearch --fasta_width 0 --sortbysize pooled-samples-node-representatives.fa --output pooled-samples-node-representatives-sorted.fa
    ## 5. Chimera cheking
    #vsearch --uchime_denovo pooled-samples-node-representatives-sorted.fa --uchimeout pooled-samples-node-representatives-sorted.uchime

### The database for the PR2 taxonomy that you will want to use for the taxonomic assignment of your SWARM ASVs can be found here
https://github.com/pr2database/pr2database/releases.  You will need to have the fasta "pr2_version_4.14.0_SSU_mothur.fasta" and the tax id file "SILVA_138.1_SSURef_NR99_tax_silva-fixed.tax" found in this git for the taxonomic assignment. Below is the general slurm script to run the taxonomy. You will need to make sure that vsearch is active. 

    #!/bin/bash

    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=00:30:00
    #SBATCH --mem=8Gb

    vsearch --usearch_global pooled-samples-node-representatives-sorted.fa --db /scratch/gpfs/WARD/JOE/DBs/pr2_version_4.14.0_SSU_mothur.fasta --blast6out NODE-HITS-PR2.txt --id 0.4
    
 #### The vsearch command above provides taxonomy ids, but not taxonomy strings for each of the unique swarm IDs in your file. The script below will create a new file that contains the swarm ID and the taxonoimic string.. because the string is what we are really after. Position one in the script is the output from the command above, position 2 is the location of the tax file you downloaded for PR2 and position 3 is the name of the output file name of your choosing.. I would keep it consistent with what I have named it in order to simplify life downstream. 
 
    #!/bin/bash

    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=00:30:00
    #SBATCH --mem=80Gb
    python combine-node-hits-with-tax-strings.py NODE-HITS-PR2.txt /scratch/gpfs/WARD/JOE/DBs/pr2_version_4.14.0_SSU_mothur.tax NODE-HITS-PR2-tax_strings.txt
    
### Look for chimeric sequences in in the pooled-samples-node-representatives-sorted.fa file using vsearch

 #!/bin/bash

    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=00:30:00
    #SBATCH --mem=8Gb
    
    vsearch --uchime_denovo pooled-samples-node-representatives-sorted.fa --uchimeout pooled-samples-node-representatives-sorted.uchime
    
### Now we have a lot of files. Including, the stats, structure, swarms, and representative sequences produced by SWARM, taxonomic hits, and chimera tables produced by VSEARCH. We want to combine this information into a single table to enable visualization and understanding of the community according to our 18S sequences. There is a single script to do this, but it won't work unless you have done everything above exactly the same as I have. This is where keeping the file names the same will be really helpful. 


    python ~/scripts/mu-swarms-to-ASVs-table-18Sv4.py -s pooled-samples-node-table.txt -o x_SWARM-contingency-table.txt -l x_fasta-names-for-swarm-otu-table-construction.txt -n pooled-samples-node-representatives.fa -t NODE-HITS-PR2-tax_strings.txt -c pooled-samples-node-representatives-sorted.uchime -st pooled-samples-derep-stats.txt
    
#### You need to transpose the matrix for flashweave analysis below. To do this, run the script "x_transpose-and-filter-counts.R" through the "x_transpose-and-filter-counts.shx" sbatch script followed by the "x_filter-flashweave-input-table.R" using the same sbatch script. Its important to also remove the sample names from the resulting ".csv file. I have been doing this step manually. 


### Now we can look at a network of the v4 sequences..I think that we shoud include all of the ASVs for this.. even the low abundance SWARMS. The x_SWARM-contingency-table.txt produced above will be most helpful to making this happen. 
### Run the flashweave analysis (I edited this section on 01042023). First create the julia script to run Flashweave. It should look like this. and the input file "x_SWARM-contingency-table.txt-swarm-matrix-for-flashweave.csv" was produced by the script above. 

    #!/usr/bin/env julia
                    
    using FlashWeave
    println("Hello climber, I'm running Flashweave for you")
                                                                                    
    v4_count_data = "x_SWARM-contingency-table.txt-swarm-matrix-for-flashweave.csv"
    v4_network_count = learn_network(v4_count_data,max_k=2, n_obs_min=3)
    save_network("x_SWARM-counts-for-flashweave-with-MAGs-output-test.edgelist", v4_network_count)                                            
    save_network("x_SWARM-counts-for-flashweave-with-MAGs-output-test.gml", v4_network_count)

### The output matrix can be merged with taxonomy etc. so that the interacting pairs contain the taxonomic identity and details that you would like in order to understand who is interacting. I copy the "x_SWARM-counts-for-flashweave-with-MAGs-output-test.edgelist" to "x_SWARM-counts-for-flashweave-with-MAGs-output.edgelist" then remove the first two lines from the "x_SWARM-counts-for-flashweave-with-MAGs-output-test.edgelist" file. Then I run the "select-associations-from-flashweave.py" script on the head node because it is a very small job. Here ia an example of how to run it.

    python ~/scripts/select-associations-from-flashweave.py x_SWARM-counts-for-includes-chimera-flashweave-output_test-23.edgelist x_SWARM-contingency-table-includes-chimera-swarm-matrix-for-anvio.txt x_SWARM-counts-for-includes-chimera-flashweave-output_test-23-edgelist-with-taxa.txt
    
    
    
 
    

    
