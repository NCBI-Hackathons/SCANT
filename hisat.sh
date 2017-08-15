#!/bin/bash

hisat2 -p $1 -x $2 --sra-acc $3 --novel-splicesite-outfile $4 | samtools sort -@ $1 -o $5
