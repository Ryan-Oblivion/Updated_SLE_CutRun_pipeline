#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cm
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=test_nf
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=rj931@nyu.edu
#SBATCH --output=slurm_%j.out
#SBATCH --array=1-6


# THIS IS A TEST TO USE THE IgG AS THE CONTROL INSTEAD OF THE CONTROL. I WILL DO THIS TWICE ALSO
# ONCE FOR KNOCKDOWN WITH IGG AND THEN WITH CONTROL AND IGG

cd store_normal_bam_files

#rm kd_bam.txt ctr_bam.txt bg_region.txt kd_ctr_bg_bam.txt


ls 454*.bam > kd_bam.txt
ls 455*.bam > ctr_bam.txt
ls ../bg_sort_bam_files/*.bam > bg_region.txt
paste kd_bam.txt ctr_bam.txt bg_region.txt > kd_ctr_bg_bam.txt

path_to_bg='../bg_sort_bam_files/'

pair_bams='kd_ctr_bg_bam.txt'

line="$(less "$pair_bams" | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)"
bam_kd="$(printf "%s" "${line}"| cut -f1)"
bam_ctr="$(printf "%s" "${line}"| cut -f2)"
bam_bg="$(printf "%s" "${line}"| cut -f3)"

module load homer/4.11
module load samtools/intel/1.14

makeTagDirectory $bam_kd'_tag_dir/' $bam_kd -tbp 1
makeTagDirectory $bam_ctr'_tag_dir/' $bam_ctr -tbp 1
makeTagDirectory $bam_bg'_tag_dir/' $path_to_bg$bam_bg -tbp 1


# here are the changes. targets are kd and ctr tag directories, while in both the input is the IgG tag directories
findPeaks $bam_kd'_tag_dir/' -style factor -i $bam_bg'_tag_dir/' -o $bam_kd'_peaks.txt'

findPeaks $bam_ctr'_tag_dir/' -style factor -i $bam_bg'_tag_dir/' -o $bam_ctr'_peaks.txt'



findPeaks $bam_bg'_tag_dir/' -style factor -o $bam_bg'_bg_peaks.txt'

#now i want to convert peaks.txt files into bed files, so i can use them in R CHIC package

pos2bed.pl $bam_kd'_peaks.txt' > $bam_kd'_peakfile.bed'

# added this one line for the conversion of the new control peak to bed
pos2bed.pl $bam_ctr'_peaks.txt' > $bam_ctr'_peakfile.bed'

pos2bed.pl $bam_bg'_bg_peaks.txt' > $bam_bg'_bg_peakfile.bed'


echo $bam_bg
echo $path_to_bg$bam_bg

ref="/scratch/work/courses/BI7653/hw3.2023/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa"

findMotifsGenome.pl $bam_kd'_peaks.txt' $ref $bam_kd'_motifOutput/' -size 200 -bg $bam_bg'_bg_peaks.txt' 

# added this new line to see the motifs in the control also. maybe if the above works, this would be needed
findMotifsGenome.pl $bam_ctr'_peaks.txt' $ref $bam_ctr'_motifOutput/' -size 200 -bg $bam_bg'_bg_peaks.txt'
