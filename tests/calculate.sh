#!/bin/bash
cd $1
find *vcf* -xtype f -size +0 | sed 's/.*\.//' | sort | uniq -c
