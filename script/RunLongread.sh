#!/bin/bash

args=("$@")
if [[ ${#args[@]} == 0 ]]; then
	echo "Usage: RunLongread.sh <sad directory> <prepare data output directory>"
	echo -e "\tsad code directory: the top directory of SAD that contains scripts and bin folder"
	echo -e "\tprepare data output directory: the output directory specified when running prepare_data.sh, which will also be used as storing the output directory of this script."
	exit
fi

scriptdir=$0
scriptdir=${scriptdir%/*}
codedir=$1
prepdir=$2

MetaFile="${prepdir}/data/LongRead1000G/lr_Metadata.txt"

GTFfile="${prepdir}/gencode.v26.annotation.gtf"
TransFasta="${prepdir}/gencode.v26.transcripts.fa"
SalmonIndex="${prepdir}/gencode.v26.full"
GenomeFasta="${prepdir}/GRCh38.primary_assembly.genome.fa"
ReadFolder="${prepdir}/data/LongRead1000G"
OutFolder="${prepdir}/LongRead1000G"

count=0
while read -r line; do

	IFS=$'\t' read -a x <<< ${line}
	sra=${x[16]}
	id=${x[6]}
	echo ${id}

	((count++))
	#if ((count <= 0)) || ((count > 9)); then
	#	continue
	#fi

	if [[ ! -e ${OutFolder}/${id}/sr_salmon/quant.sf ]]; then
		echo -e "\tSALMON QUANTIFYING..."
		mkdir -p ${OutFolder}/${id}/sr_salmon
		salmon quant -p 4 -l A -i ${SalmonIndex} -1 ${ReadFolder}/${id}_1.fastq.gz -2 ${ReadFolder}/${id}_2.fastq.gz --gcBias --seqBias --posBias --dumpEqWeights -o ${OutFolder}/${id}/sr_salmon --writeMappings=${OutFolder}/${id}/sr_salmon/mapping.sam
	fi

	if [[ ! -e ${OutFolder}/${id}/lr_mummer/nucmer_${sra}.delta ]]; then
		echo -e "\tMUMMER ALIGNING LONG READ..."
		mkdir -p ${OutFolder}/${id}/lr_mummer
		nucmer  ${GenomeFasta} ${ReadFolder}/${sra}.fastq -p ${OutFolder}/${id}/lr_mummer/nucmer_${sra}
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
		python3 ${scriptdir}/src/LongreadsNucmer.py ${GTFfile} ${OutFolder}/${id}/sr_salmon/sad_unadjustable_pvalue.tsv ${OutFolder}/${id}/lr_mummer/nucmer_${sra}.txt ${OutFolder}/${id}/sr_salmon/sad_unadjustable_validation_${sra}.txt
	fi

	# merge the hits in long reads
	if [[ ! -e ${OutFolder}/${id}/sr_salmon/sad_unadjustable_validation_merge.txt ]]; then
		python3 ${scriptdir}/src/MergeLongreadsNucmer.py ${OutFolder}/${id}/sr_salmon/sad_unadjustable_validation_ERR ${OutFolder}/${id}/sr_salmon/sad_unadjustable_validation_merge.txt
	fi

done < ${MetaFile}

