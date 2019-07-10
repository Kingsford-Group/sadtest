#!/bin/bash

folder=$0
folder=${folder%/*}

while read -r line; do

	IFS=$'\t' read -a x <<< ${line}
	sra=${x[16]}
	id=${x[6]}
	echo -e ${id}"\t"${sra}
	if [[ ! -e ${sra}.fastq ]] && [[ ! -e ${sra}.fastq.gz ]]; then
		echo ${sra}
		fastq-dump -O ${folder}/ ${sra}
	fi
done < ${folder}/lr_Metadata.txt

linecount=0
while read -r line; do
	((linecount++))
	
	if ((linecount==1)); then
		continue
	fi

	IFS=$'\t' read -ra x <<< ${line}
	id=${x[5]}
	url=${x[0]}
	if [[ ${url} == *"_1.fastq.gz" ]]; then
		wget ${url} -O ${folder}/${id}_1.fastq.gz
	elif [[ ${url} == *"_2.fastq.gz" ]]; then
		wget ${url} -O ${folder}/${id}_2.fastq.gz
	else
		echo "Error"
	fi
done < ${folder}/sr_RunTable.txt
