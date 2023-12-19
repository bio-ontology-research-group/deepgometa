#!/bin/bash

HOSTNAME=""
HOSTGENOME=""
phiX=""

#Trimming (fastp)
fastp -i ${1}_1.fastq.gz \
	-I ${1}_2.fastq.gz \
	-o ${1}_trimmed_1.fastq.gz \
	-O ${1}_trimmed_2.fastq.gz \
	-q 30 \
	--thread 6

#Host removal (Bowtie2)
bowtie2-build $phiX phiX
bowtie2-build $HOSTGENOME $HOSTNAME

bowtie2 --threads 6 --quiet -x $HOSTNAME \
-1 ${1}_trimmed_1.fastq.gz \
-2 ${1}_trimmed_2.fastq.gz \
--un-conc-gz ${1}_trimmed.nohost.%.fq.gz

#Assembly
megahit \
-1 ${1}_trimmed.nohost.1.fq.gz \
-2 ${1}_trimmed.nohost.2.fq.gz \
-o ${1}_megahitout \
--mem-flag 128 \
-t 18

#Prodigal
prodigal -i ${1}_megahitout/final.contigs.fa \
-a ${1}_megahitout/${1}.prodigal.fa \
-p meta
