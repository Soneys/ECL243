#!/usr/bin/env bash 

Trinity --seqType fq \
	--max_memory 80G \
	--trimmomatic \
	--SS_lib_type RF \
	--CPU 8 \
	--full_cleanup \
	--left SRR6466772_1.fastq.gz,SRR6466760_1.fastq.gz,SRR6466774_1.fastq.gz,SRR6466761_1.fastq.gz \
	--right SRR6466772_2.fastq.gz,SRR6466760_2.fastq.gz,SRR6466774_2.fastq.gz,SRR6466761_2.fastq.gz
