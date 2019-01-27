#!/bin/bash

# name : Leila Esmaeili
# date : December 5, 2018
# contact info:  lesmaeili@dons.usfca.edu
#  This is a bash script to do a fastq dump to access data on the NCBI website
# for the final project article.

# set the PATH for fastq-dump and other NCBI tools
export PATH=$PATH:/Users/leilaesmaeili/sratoolkit.2.9.2-mac64/bin

# exclude first line
for SRA_number in $(cut -f 10 data/metadata/SraRunTable.txt | tail -n +2)
do
   echo Downloading "$SRA_number"
   fastq-dump -v "SRA_number"-O data/raw_data
done