#!/bin/bash

folder=$0
folder=${folder%/*}

wget -O ${folder}/SRR3192412_1.fastq.gz https://www.encodeproject.org/files/ENCFF001RFF/@@download/ENCFF001RFF.fastq.gz
wget -O ${folder}/SRR3192412_2.fastq.gz https://www.encodeproject.org/files/ENCFF001RFE/@@download/ENCFF001RFE.fastq.gz

wget -O ${folder}/SRR3192396_1.fastq.gz https://www.encodeproject.org/files/ENCFF001RFH/@@download/ENCFF001RFH.fastq.gz
wget -O ${folder}/SRR3192396_2.fastq.gz https://www.encodeproject.org/files/ENCFF001RFG/@@download/ENCFF001RFG.fastq.gz

wget -O ${folder}/ERR188297_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188297/ERR188297_1.fastq.gz
wget -O ${folder}/ERR188297_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188297/ERR188297_2.fastq.gz
