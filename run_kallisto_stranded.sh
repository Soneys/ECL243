#!/usr/bin/env bash

for f1 in *_1_val_1.fq.gz
do
	echo ${f1%%_1_val_1.fq.gz}
	f2=${f1%%_1_val_1.fq.gz}"_2_val_2.fq.gz"
	kallisto quant --rf-stranded -b 100 -i dev_transcriptome.idx -o kallisto_outs/${f1%%_1_val_1.fq.gz} $f1 $f2
done
