#!/usr/bin/env bash

for i in SRR6466778_2.fastq.gz SRR6466778_1.fastq.gz \
	 SRR6466777_2.fastq.gz SRR6466777_1.fastq.gz \
	 SRR6466776_2.fastq.gz SRR6466776_1.fastq.gz \
	 SRR6466775_2.fastq.gz SRR6466775_1.fastq.gz \
	 SRR6466774_2.fastq.gz SRR6466774_1.fastq.gz \
         SRR6466773_2.fastq.gz; do
	echo "$i"
	cp "$i" "$i~" &&
	gunzip -cd "$i~" | sed 's/ /:/g' | gzip > "$i"
done
