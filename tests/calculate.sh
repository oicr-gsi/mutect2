#!/bin/bash
cd $1
find *vcf* -type f -size +0 | sed 's/.*\.//' | sort | uniq -c