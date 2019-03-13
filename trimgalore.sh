#!/bin/bash 
for f1 in *_1.fastq.gz 
do
        f2=${f1%%_1.fastq.gz}"_2.fastq.gz"       
	echo $f1 $f2 \n
	trim_galore --paired -o trim_ena/ $f1 $f2 
done
