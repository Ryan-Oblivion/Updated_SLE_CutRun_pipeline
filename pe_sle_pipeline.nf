// this works to get the read pairs in a new tuple form


// run this to see how it looks again
params.reads = 
"/scratch/rj931/tf_sle_project/all_sle_data/454-Mock*cut*_{R1,R2}*.fastq.gz"


FASTP='fastp/intel/0.20.1'
FASTQC='fastqc/0.11.9'
BWA='bwa/intel/0.7.17'
SAMTOOLS='samtools/intel/1.14'



params.outdir = 'filt_files'

params.filts = "filt_files/454-Mock*cut*_{R1,R2}*.filt*"

nextflow.enable.dsl=2

process fastp {

//tag
publishDir params.outdir, mode: 'copy'

input:
tuple val(pair_id), path(PE_reads)
// input in fastp would be ${reads[0]} and ${reads[1]}
// something i can do in fastp is to now concatenate the $pair_id'_R1.filt.fastq.gz' 
// same for output 2 in fastp use $pair_id'_R2.filt.fastq.gz'

output:
path "${pair_id}_R1.filt.fq.gz", emit: fastp_out_f
path "${pair_id}_R2.filt.fq.gz", emit: fastp_out_r

"""
#!/bin/env bash

module load $FASTP

fastp \
-i ${PE_reads[0]} \
-I ${PE_reads[1]} \
-o $pair_id'_R1.filt.fq.gz' \
-O $pair_id'_R2.filt.fq.gz' \
--detect_adapter_for_pe \
--trim_front1 7 \
--trim_front2 7 \
"""

}

process bwa {

input:
val ref
tuple val(pair_id), path(filt_pe)


output:
path "${pair_id}.bam", emit: bam_file

"""
#!/bin/env bash

module load $BWA
module load $SAMTOOLS

bwa index -a bwtsw $ref

bwa aln -t 8 $ref ${filt_pe[0]} > ${filt_pe[0]}'reads_1.sai'
bwa aln -t 8 $ref ${filt_pe[1]} > ${filt_pe[1]}'reads_2.sai'

bwa sampe $ref ${filt_pe[0]}'reads_1.sai' ${filt_pe[1]}'reads_2.sai' ${filt_pe[0]} ${filt_pe[1]} \
> ${pair_id}'aligned_reads.sam'

samtools view -b -h -q 20 ${pair_id}'aligned_reads.sam' -o $pair_id'.bam'
"""
}


PE_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
filt_pe = Channel.fromFilePairs(params.filts, checkIfExists: true)

workflow{

// look what channel.value does. it lets me use the ref multiple times
ref = Channel.value(params.ref)

//PE_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)

//PE_reads.view()
//filt_pe.view()
main:

fastp(PE_reads)
bwa(ref, filt_pe)
//fastp.out.fastp_out_f.view()
//fastp.out.fastp_out_r.view()
bwa.out.bam_file.view()
}
