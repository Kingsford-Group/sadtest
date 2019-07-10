#!/bin/bash

DataFolder=/mnt/disk33/user/congm1/RawData/LongRead1000G
OutFolder=/home/congm1/savanna/savannacong33/SADrealdata/LongRead1000G/

SalmonIndex=/mnt/disk33/user/congm1/NCBI/gencode.v26.full
GTFfile=/mnt/disk33/user/congm1/NCBI/gencode.v26.annotation.gtf
GenomeFasta=/mnt/disk33/user/congm1/NCBI/GRCh38.primary_assembly.genome.fa
TransFasta=/mnt/disk33/user/congm1/NCBI/gencode.v26.full.transcripts.fa

count=0
while read -r line; do

	IFS=$'\t' read -a x <<< ${line}
	sra=${x[16]}
	id=${x[6]}
	echo ${id}

	((count++))
	if ((count <= 0)) || ((count > 9)); then
		continue
	fi

	if [[ ! -e ${OutFolder}/${id}/sr_salmon/quant.sf ]]; then
		echo -e "\tSALMON QUANTIFYING..."
		mkdir -p ${OutFolder}/${id}/sr_salmon
		salmon quant -p 4 -l A -i ${SalmonIndex} -1 ${DataFolder}/${id}_1.fastq.gz -2 ${DataFolder}/${id}_2.fastq.gz --gcBias --seqBias --posBias --dumpEqWeights -o ${OutFolder}/${id}/sr_salmon --writeMappings=${OutFolder}/${id}/sr_salmon/mapping.sam
	fi

	if [[ ! -e ${OutFolder}/${id}/lr_mummer/nucmer_${sra}.delta ]]; then
		echo -e "\tMUMMER ALIGNING LONG READ..."
		mkdir -p ${OutFolder}/${id}/lr_mummer
		nucmer  ${GenomeFasta} ${DataFolder}/${sra}.fastq -p ${OutFolder}/${id}/lr_mummer/nucmer_${sra}
	fi

	if [[ ! -e ${OutFolder}/${id}/lr_mummer/nucmer_${sra}.txt ]]; then
		echo -e "\tSHOW-COORDS..."
		show-coords ${OutFolder}/${id}/lr_mummer/nucmer_${sra}.delta > ${OutFolder}/${id}/lr_mummer/nucmer_${sra}.txt
	fi

	if [[ ! -e ${OutFolder}/${id}/sr_salmon/sad_unadjustable_pvalue.tsv ]]; then
		echo python3 /home/congm1/savanna/savannacong33/Code/sad/scripts/SADpipe.py -t ${TransFasta} -a ${GTFfile} -s ${OutFolder}/${id}/sr_salmon/ -o ${OutFolder}/${id}/sr_salmon/sad
	fi

	# analysis: compare SAD unadjustable anomalies and long reads
	if [[ ! -e ${OutFolder}/${id}/sr_salmon/sad_unadjustable_validation_${sra}.txt ]] && [[ -e ${OutFolder}/${id}/sr_salmon/sad_unadjustable_pvalue.tsv ]] && [[ -e ${OutFolder}/${id}/lr_mummer/nucmer_${sra}.txt ]]; then
		echo -e "\tCOMPARING SAD UNADJUSTABLE ANOMALIES AND LONG READ ALIGNMENTS..."
		python /home/congm1/savanna/savannacong33/Code/SAD_ori/src/LongreadsNucmer.py ${GTFfile} ${OutFolder}/${id}/sr_salmon/sad_unadjustable_pvalue.tsv ${OutFolder}/${id}/lr_mummer/nucmer_${sra}.txt ${OutFolder}/${id}/sr_salmon/sad_unadjustable_validation_${sra}.txt
	fi

	# merge the hits in long reads
	if [[ ! -e ${OutFolder}/${id}/sr_salmon/sad_unadjustable_validation_merge.txt ]]; then
		python /home/congm1/savanna/savannacong33/Code/SAD_ori/src/MergeLongreadsNucmer.py ${OutFolder}/${id}/sr_salmon/sad_unadjustable_validation_ERR ${OutFolder}/${id}/sr_salmon/sad_unadjustable_validation_merge.txt
	fi

done < lr_Metadata.txt

