# MESOCOSM_18S
This git contains the code and steps to process the v4 region of the 18S rRNA gene. The steps for amplicon library construction are outlined in publication links to this git and detailed steps are contained in the "16S v4 18s Amplicon Protocol for Illumina Sequencing.doc", located in the code tab of this git repository. Most of the code below contains information for the SLURM scheduler as well. Note to self, you are working here /scratch/gpfs/WARD/JOE/MOSS_BLOOM/18S_20220502_Data.

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

    iu-gen-configs 00_DEMULTIPLEXING_REPORT
    ls *.ini | sed 's/\.ini//g' | grep 18S> x_18S-samples.txt
    
##### now you can merge the sequnces for each sample

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=00:30:00
    #SBATCH --mem=80Gb
    #SBATCH --array=1-15

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_18S-samples.txt)
    iu-merge-pairs ${SAMPLE}.ini

#### IF YOU CHOOSE TO TRIM SEQUENCES PRIOR TO MERGING:

##### First run prinseq on your set of samples in order to trim both read1 and read2 to the desired length. The standard prinseq command to do this looks like this

    prinseq-lite.pl -fastq B1D2T1_18S-R1.fastq -trim_to_len 240 -out_good B1D2T1_18S-R1-prinseq
    prinseq-lite.pl -fastq B1D2T1_18S-R2.fastq -trim_to_len 240 -out_good B1D2T1_18S-R2-prinseq
    
##### however, you will want to run it in a loop using the array of the SLURM scheduler which will look like this.
    
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=00:05:00
    #SBATCH --mem=80Gb
    #SBATCH --array=1-15

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_18S-samples.txt)
    prinseq-lite.pl -fastq ${SAMPLE}-R1.fastq -trim_to_len 240 -out_good ${SAMPLE}-R1-prinseq
    prinseq-lite.pl -fastq ${SAMPLE}-R2.fastq -trim_to_len 240 -out_good ${SAMPLE}-R2-prinseq
    
##### Then you can run the merging script on the trimmed sequences. But first you will need to create the .ini files for this. You can do this by editing your 00_DEMULTIPLEXING_REPORT and then generating the ini files. you will just have to be careful that you activate the proper conda environments when running these commands. 

    sed 's/R1/R1-prinseq/g' 00_DEMULTIPLEXING_REPORT | sed 's/R2/R2-prinseq/g' > 00_DEMULTIPLEXING_REPORT-prinseq
    
###### Don't forget to add these lines of text to your 00_DEMULTIPLEXING_REPORT-prinseq

    sample	num_indexes_found	num_reads_stored	r1	r2

###### The 00_DEMULTIPLEXING_REPORT-prinseq should look something like this

    sample	num_indexes_found	num_reads_stored	r1	r2
    B3D2T1_18S-prinseq  91686	91686	B3D2T1_18S-R1-prinseq.fastq	    B3D2T1_18S-R2-prinseq.fastq
    B1D5T2A_18S-prinseq 86714	86714	B1D5T2A_18S-R1-prinseq.fastq	B1D5T2A_18S-R2-prinseq.fastq
    B3D5T2C_18S-prinseq	80791	80791	B3D5T2C_18S-R1-prinseq.fastq	B3D5T2C_18S-R2-prinseq.fastq
    B1D5T1A_18S-prinseq	59856	59856	B1D5T1A_18S-R1-prinseq.fastq	B1D5T1A_18S-R2-prinseq.fastq
    B1D4T1A_18S-prinseq	56279   56279	B1D4T1A_18S-R1-prinseq.fastq	B1D4T1A_18S-R2-prinseq.fastq
    
##### Now you can generate the ini files and run the merging. 
    
    iu-gen-configs 00_DEMULTIPLEXING_REPORT-prinseq
    
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=00:30:00
    #SBATCH --mem=80Gb
    #SBATCH --array=1-15

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_samples.txt)
    iu-merge-pairs ${SAMPLE}-prinseq.ini
    
#### Next we filter reads that don't have the primers. Here is the slurm for that

#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --time=00:30:00
#SBATCH --mem=80Gb
#SBATCH --array=1-30

SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_merged-samples.txt)
python ~/scripts/filter-for-18S-primer.py --i ${SAMPLE} --o ${SAMPLE}-primer-filtered.fa
    
#### Now you will want to run SWARM using the quality and primer filtered sequences. I use this general script below to run swarm, but I comment out the parts of the slurm and activate different parts of the command one at a time using the "#" character. 

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=200Gb
    #SBATCH --time=04:00:00

    ### These steps replicate the work here
    ### https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline#global-dereplication-clustering-and-chimera-detection
    ### The cat and vsearch steps should be run with the conda vsearch envionment
    ### The sarm should be run with this environment "conda activate /home/jv2474/.conda/envs/swarm-v3.1"
    ### The python script should be run with the bioconda environment.

    ## 1. Concatenate the merged and filtered sequences
    cat *-primer-derep.fa > pooled-samples.fa
    ## 2. Dereplicaete the concatenated sequences
    vsearch --derep_fulllength pooled-samples.fa --fasta_width 0 --sizeout --sizein --output pooled-samples-derep.fa
    ## 3. Cluster the sequences
    swarm -d 1 -f -t 40 -z pooled-samples-derep.fa -i pooled-samples-derep-struct.txt -s pooled-samples-derep-stats.txt -w pooled-samples-node-representatives.fa -o pooled-samples-node-table.txt
    ## 4. Sort the clustered node representatives
    vsearch --fasta_width 0 --sortbysize pooled-samples-node-representatives.fa --output pooled-samples-node-representatives-sorted.fa


### The database for the PR2 taxonomy that you will want to use for the taxonomic assignment of your SWARM ASVs can be found here
https://github.com/pr2database/pr2database/releases.  You will need to have the fasta "pr2_version_4.14.0_SSU_mothur.fasta" and the tax id file "SILVA_138.1_SSURef_NR99_tax_silva-fixed.tax" found in this git for the taxonomic assignment. Below is the general slurm script to run the taxonomy. You will need to make sure that vsearch is active.

    #!/bin/bash

    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=00:30:00
    #SBATCH --mem=80Gb

    vsearch --usearch_global pooled-samples-node-representatives-sorted.fa --db /scratch/gpfs/WARD/JOE/DBs/pr2_version_4.14.0_SSU_mothur.fasta --blast6out NODE-HITS-PR2.txt --id 0.4
 
 #### Now you will want to create a file that contains meaningful taxonomy for your taxonomic hits contained in the "NODE-HITS.txt" file. Here is the command that will get you there. This script will produce two files.. 1. contains the metadata includin the ASV and the taxonomic string and the other is the matrix of counts for each sample
 
    #!/bin/bash

    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=00:30:00
    #SBATCH --mem=80Gb
 
    python ~/scripts/convert-node-hits-to-tax-node-table.py -n NODE-HITS-PR2.txt -o x_SWARMS-and-tax-for-anvio-pr2 -r /scratch/gpfs/WARD/JOE/DBs/pr2_version_4.14.0_SSU_mothur.tax -s x_SWARM-contingency-table.txt -min 50
 
 
 
    
    
