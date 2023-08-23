#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cm
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=best_gg
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rj931@nyu.edu
#SBATCH --output=slurm_%j.out



# DIRECTIONS TO ALTER THE INPUT OF THE PAIR END READS AND ALSO THE BACKGROUND GENOMIC DATA (BG)

# this is what the pipeline has by defualt , only change bg regions if you have your own. which is required

# MOCK READS
# params.reads = "/scratch/rj931/tf_sle_project/all_sle_data/45*-Mock*cut*_{R1,R2}*.fastq.gz"
# params.filts = "filt_files/45*-Mock*cut*_{R1,R2}*.filt*"
# params.bg_regions = "/scratch/rj931/tf_sle_project/all_sle_data/461-IgG*cut*_{R1,R2}*.fastq.gz"

# IAV READS
# reads = "/scratch/rj931/tf_sle_project/all_sle_data/45*-IAV*cut*_{R1,R2}*.fastq.gz"
# filts = "filt_files/45*-IAV*cut*_{R1,R2}*.filt*"
# bg_regions = "/scratch/rj931/tf_sle_project/all_sle_data/461-IgG*cut*_{R1,R2}*.fastq.gz"

# BLEO READS
# reads = "/scratch/rj931/tf_sle_project/all_sle_data/45*-Bleo*cut*_{R1,R2}*.fastq.gz"
# filts = "filt_files/45*-Bleo*cut*_{R1,R2}*.filt*"
# bg_regions = "/scratch/rj931/tf_sle_project/all_sle_data/461-IgG*cut*_{R1,R2}*.fastq.gz"

#####################
# IRF1 input parameters

# IRF1 MOCK READS
#  params.reads = "/scratch/rj931/tf_sle_project/all_sle_data/46*-Mock*cut*_{R1,R2}*.fastq.gz"
# params.filts = "filt_files/46*-Mock*cut*_{R1,R2}*.filt*"
# params.bg_regions = "/scratch/rj931/tf_sle_project/all_sle_data/461-IgG*cut*_{R1,R2}*.fastq.gz"

# IRF1 IAV READS
# reads = "/scratch/rj931/tf_sle_project/all_sle_data/46*-IAV*cut*_{R1,R2}*.fastq.gz"
# filts = "filt_files/46*-IAV*cut*_{R1,R2}*.filt*"
# bg_regions = "/scratch/rj931/tf_sle_project/all_sle_data/461-IgG*cut*_{R1,R2}*.fastq.gz"

# IRF1 BLEO READS
# reads = "/scratch/rj931/tf_sle_project/all_sle_data/46*-Bleo*cut*_{R1,R2}*.fastq.gz"
# filts = "filt_files/46*-Bleo*cut*_{R1,R2}*.filt*"
# bg_regions = "/scratch/rj931/tf_sle_project/all_sle_data/461-IgG*cut*_{R1,R2}*.fastq.gz"

module purge

module load nextflow/23.04.1

# i have to run this twice since the files are being created and sent to the dir
# so now that they exist the pipeline can move on, how to fix this and make it run
# or check the directories again instead

nextflow run -resume pe_sle_pipeline.nf  --reads "/scratch/rj931/tf_sle_project/all_sle_data/45*-Bleo*cut*_{R1,R2}*.fastq.gz" --filts "filt_files/45*-Bleo*cut*_{R1,R2}*.filt*"
nextflow run -resume pe_sle_pipeline.nf  --reads "/scratch/rj931/tf_sle_project/all_sle_data/45*-Bleo*cut*_{R1,R2}*.fastq.gz" --filts "filt_files/45*-Bleo*cut*_{R1,R2}*.filt*"





find . -name *fastqc.zip > fastqc_files.txt

module load multiqc/1.9

multiqc -force --file-list fastqc_files.txt --filename 'multiqc_report.html'

#nextflow run -resume pe_sle_pipeline.nf

# here i want to out put the genomic regions
# specify the path and the glob pattern for the --reads
# specify the the out dir for this set of filtered outputs
# specify the location of the bg filtered output fq files and the same glob pattern in the reads

#nextflow run -resume pe_sle_pipeline.nf \
#--reads "/scratch/rj931/tf_sle_project/all_sle_data/461-IgG*cut*_{R1,R2}*.fastq.gz" \
#--outdir "bg_filt" \
#--filts "bg_filt/461-IgG*cut*_{R1,R2}*.fastq.gz" \



# this section is to run the make_homer.sh script part of the pipeline

#sbatch make_homer.sh

# here I want to combine the tag directories for knockdown and control to make one peak file that contains 
# differentially expressed peaks from each replicate
#sbatch homer_anno_combined.sh

# this allows for the second script to only run after the first is finished
JOBID1=$(sbatch --parsable --array=1-6 make_homer.sh)
sbatch --dependency=afterok:$JOBID1 homer_anno_combined.sh