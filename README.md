# MESOCOSM_18S
This git contains the code and steps to process the v4 region of the 18S rRNA gene and other supplemental data associated with the manuscript - "Co-occurrence and successional patterns among diatoms, dinoflagellates, and potential parasites in a coastal upwelling experiment". 

### Code required for this tutorial and links to the conda install or more specific install instructions. 
    illumina-utilities: conda install -c bioconda illumina-utils
    anvio-7.1: https://anvio.org/install/
    vsearch: conda install -c bioconda vsearch
    oligotyping:install the latest stable release (artistic mode) https://merenlab.org/2014/08/16/installing-the-oligotyping-pipeline/

### First download the sequences from the Short Read Archive https://dataview.ncbi.nlm.nih.gov/object/PRJNA1060340?reviewer=ciehuil2a33fpn566l3kq3clsq. 

##### before you can merge, you need .ini files which can be made according to the command below and the 00_DEMULTIPLEXING_18S can be found in the Code section of this git. 

    iu-gen-configs 00_DEMULTIPLEXING_18S --r1-prefix ^GTG[C,T]CAGC[A,C]GCCGCGGTAA --r2-prefix ^TTGG[C,T][A,G]AATGCTTTCGC
    ls *.ini | sed 's/\.ini//g' | grep 18S > x_18S-samples.txt
    
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


    
    
    
 
    

    
