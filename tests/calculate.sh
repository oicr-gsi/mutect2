#!/bin/bash
cd $1
for v in *.vcf.gz;do zcat $v | grep -v ^# | md5sum; done | sort -V
