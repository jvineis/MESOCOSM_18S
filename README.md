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

##### before you can merge, you need .ini files which can be made like this

    iu-gen-configs 00_DEMULTIPLEXING_REPORT
    ls *.ini | sed 's/\.ini//g' > x_samples.txt
    
##### now you can merge the sequnces for each sample

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --time=00:30:00
    #SBATCH --mem=80Gb
    #SBATCH --array=1-30

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
    #SBATCH --array=1-30

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p x_18S-samples.txt)
    prinseq-lite.pl -fastq ${SAMPLE}-R1.fastq -trim_to_len 240 -out_good ${SAMPLE}-R1-prinseq
    prinseq-lite.pl -fastq ${SAMPLE}-R2.fastq -trim_to_len 240 -out_good ${SAMPLE}-R2-prinseq


    
