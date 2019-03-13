#!/usr/bin/env bash

for i in *.fq.gz
do
	echo $i
	kallisto quant --single -l 250 -s 50 -i dev_transcriptome.idx -o kallisto_outs/${i%%_z} $i
done
