#!/bin/bash

args=("$@")
if [[ ${#args[@]} == 0 ]]; then
	echo "Usage: RunRealGEUDAVIS.sh <sad directory> <prepare data output directory>"
	echo -e "\tsad code directory: the top directory of SAD that contains scripts and bin folder"
	echo -e "\tprepare data output directory: the output directory specified when running prepare_data.sh, which will also be used as storing the output directory of this script."
	exit
fi

scriptdir=$0
scriptdir=${scriptdir%/*}
codedir=$1
prepdir=$2

MetaFile="${prepdir}/data/GEUVADIS/Metadata.txt"

Type=("Full")
GTFfiles=("${prepdir}/gencode.v26.annotation.gtf")
TransFastas=("${prepdir}/gencode.v26.transcripts.fa")
SalmonIndex=("${prepdir}/gencode.v26.full")
StarIndex=("${prepdir}/StarIndex")
StarTransIndex=("${prepdir}/gencode.v26.StarIndex")
RSEMIndex="${prepdir}/gencode.v26.RSEM/ref"
ReadFolder="${prepdir}/data/GEUVADIS"
OutDirectory="${prepdir}/GEUVADIS"
SADFolder="sad"

count=0

i=0
#for ((i=0; i<${#Type[@]}; i++)); do
	t=${Type[${i}]}
	gtffile=${GTFfiles[${i}]}
	transfasta=${TransFastas[${i}]}
	salmonindex=${SalmonIndex[${i}]}
	starindex=${StarIndex[${i}]}
	startransindex=${StarTransIndex[${i}]}
	while read -r line; do

		((count++))
		#if ((count<=0)) || ((count>6)); then
		#	continue
		#fi

		read -ra x <<< ${line}
		ID=${x[${#x[@]}-1]}
		read1=${ReadFolder}/${ID}/${ID}"_1.fastq.gz"
		read2=${ReadFolder}/${ID}/${ID}"_2.fastq.gz"

		echo "${OutDirectory}/salmon_${t}_${ID}/"

		# running salmon quantification
		if [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/quant.sf ]]; then
			mkdir -p ${OutDirectory}/salmon_${t}_${ID}/
			echo -e "\tsalmon quantification"

			salmon quant -p 4 -l A -i ${salmonindex} -1 ${read1} -2 ${read2} --gcBias --seqBias --posBias --dumpEqWeights --numBootstraps 100 -o ${OutDirectory}/salmon_${t}_${ID}/ --writeMappings=${OutDirectory}/salmon_${t}_${ID}/mapping.sam
			samtools view -Shb ${OutDirectory}/salmon_${t}_${ID}/mapping.sam -o ${OutDirectory}/salmon_${t}_${ID}/mapping.bam
		fi

		# run SADpipe
		if [[ ! -e ${prepdir}/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue.tsv ]]; then
			python3 ${codedir}/scripts/SADpipe.py -m 0 -t ${transfasta} -a ${gtffile} -s ${OutDirectory}/salmon_${t}_${ID}/ -o ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad
		fi

		# post-process SAD
		if [[ -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue.tsv ]] && [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv ]]; then
			echo -e "\tpost-process SAD p values"
			python3 ${codedir}/scripts/AdjustPValue.py 1 ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue.tsv ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv
			python3 ${codedir}/scripts/Pred_Trans2Gene.py ${gtffile} ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv
		fi

		# STAR alignment for transcriptome assembly
		if [[ ! -e ${OutDirectory}/star_${t}_${ID}/Aligned.sortedByCoord.out.bam ]]; then
			mkdir -p ${OutDirectory}/star_${t}_${ID}/
			echo -e "\tSTAR aligning"
			STAR --runThreadN 4 --genomeDir ${starindex} --readFilesIn ${read1} ${read2} --readFilesCommand gunzip -c --outFileNamePrefix ${OutDirectory}/star_${t}_${ID}/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --sjdbGTFfile ${gtffile}
		fi

		# transcriptome assembly: Scallop
		if [[ ! -e ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf ]]; then
			echo -e "\tScallop"
			scallop -i ${OutDirectory}/star_${t}_${ID}/Aligned.sortedByCoord.out.bam --verbose 0 -o ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf
			gffcompare -o ${OutDirectory}/star_${t}_${ID}/gffcompscallop -r ${gtffile} ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf
		fi

		# transcriptome assembly StringTie
		if [[ ! -e ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf ]]; then
			echo -e "\tStringTie"
			stringtie ${OutDirectory}/star_${t}_${ID}/Aligned.sortedByCoord.out.bam -G ${gtffile} -o ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf
			gffcompare -o ${OutDirectory}/star_${t}_${ID}/gffcompstringtie -r ${gtffile} ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf
		fi

		# analysis: compare SAD and scallop
		if [[ -e ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf ]] && [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_scallopcomp.txt ]]; then
			echo -e "\tcomparing SAD and Scallop prediction"
			python3 ${codedir}/FindSADpredinAssembly.py 1 ${gtffile} ${OutDirectory}/star_${t}_${ID}/scallopgenes.gtf ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_scallopcomp.txt
		fi

		# analysis: compare SAD and stringtie
		if [[ -e ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf ]] && [[ ! -e ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_stringtiecomp.txt ]]; then
			echo -e "\tcomparing SAD and StringTie prediction"
			python3 ${codedir}/FindSADpredinAssembly.py 1 ${gtffile} ${OutDirectory}/star_${t}_${ID}/stringtiegenes.gtf ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv ${OutDirectory}/salmon_${t}_${ID}/${SADFolder}/sad_unadjustable_stringtiecomp.txt
		fi


		# About RSEM
		# Run RSEM
		if [[ ! -e ${OutDirectory}/rsem_${t}_${ID}/rsem.isoforms.results ]]; then
			echo -e "\tquantifying using RSEM"
			mkdir -p ${OutDirectory}/rsem_${t}_${ID}/
			rsem-calculate-expression -p 4 --star --star-gzipped-read-file --estimate-rspd --paired-end ${read1} ${read2} ${RSEMIndex} ${OutDirectory}/rsem_${t}_${ID}/rsem
		fi

		# SAD for RSEM
		if [[ ! -e ${OutDirectory}/rsem_${t}_${ID}/${SADFolder}/sad_unadjustable_pvalue.tsv ]]; then
			python3 ${codedir}/scripts/SADpipe.py -m 1 -t ${transfasta} -a ${gtffile} -s ${OutDirectory}/rsem_${t}_${ID}/rsem -o ${OutDirectory}/rsem_${t}_${ID}/${SADFolder}/sad
		fi


		# Identifiability using eXpress
		# star alignment onto the transcript index
		if [[ ! -e ${OutDirectory}/star_transcript_${t}_${ID}/Aligned.sortedByCoord.out.bam ]]; then
			echo "STAR transcriptome aligning..."
			mkdir -p ${OutDirectory}/star_transcript_${t}_${ID}
			STAR --runThreadN 4 --genomeDir ${startransindex} --readFilesIn ${read1} ${read2} --readFilesCommand gunzip -c --outFileNamePrefix ${OutDirectory}/star_transcript_${t}_${ID}/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 100
		fi

		# sort by read name
		if [[ ! -e ${OutDirectory}/star_transcript_${t}_${ID}/Aligned.sortedByName.out.bam ]]; then
			echo "SORTING by read name..."
			samtools sort -n ${OutDirectory}/star_transcript_${t}_${ID}/Aligned.sortedByCoord.out.bam -o ${OutDirectory}/star_transcript_${t}_${ID}/Aligned.sortedByName.out.bam
		fi

		# run eXpress
		if [[ ! -e ${OutDirectory}/star_transcript_${t}_${ID}/results.xprs ]]; then
			echo "RUNNING eXpress..."
			express -o ${OutDirectory}/star_transcript_${t}_${ID}/ ${transfasta} ${OutDirectory}/star_transcript_${t}_${ID}/Aligned.sortedByName.out.bam
		fi
			

	done < ${MetaFile}
#done

