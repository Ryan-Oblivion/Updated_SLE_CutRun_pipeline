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


findPeaks $bam_kd'_tag_dir/' -style factor -i $bam_ctr'_tag_dir/' -o $bam_kd'_peaks.txt'

findPeaks $bam_bg'_tag_dir/' -style factor -o $bam_bg'_bg_peaks.txt'

echo $bam_bg
echo $path_to_bg$bam_bg

ref="/scratch/work/courses/BI7653/hw3.2023/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa"

findMotifsGenome.pl $bam_kd'_peaks.txt' $ref $bam_kd'_motifOutput/' -size 200 -bg $bam_bg'_bg_peaks.txt' 
