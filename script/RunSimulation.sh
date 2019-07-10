#!/bin/bash

args=("$@")
if [[ ${#args[@]} == 0 ]]; then
        echo "Usage: preparing_data.sh <sad directory> <prepare data output directory>"
        exit
fi

scriptdir=$0
scriptdir=${scriptdir%/*}
codedir=$1
prepdir=$2

echo "SADTEST SCRIPT DIRECTORY: ${scriptdir}"
echo "SAD DIRECTORY: ${codedir}"
echo "PREPARED DATA DIRECTORY: ${preddir}"

Type=("PC" "Full")
GTFfiles=("${prepdir}/gencode.v26.annotation.pc.gtf" "${prepdir}/gencode.v26.annotation.gtf")
TransFastafiles=("${prepdir}/gencode.v26.pc.transcripts.fa" "${prepdir}/gencode.v26.transcripts.fa")
NumEvents=(200 500 1000 1500)
RefQuantPrefix=("${prepdir}/data/BaseExpression/ERR188297" "${prepdir}/data/BaseExpression/SRR3192396" "${prepdir}/data/BaseExpression/SRR3192412")
RefQuantID=(GEU GM12878 K562)
SADFolder="sad"

#i=0
for ((i=0; i<${#Type[@]}; i++)); do
	t=${Type[${i}]}
	gtffile=${GTFfiles[${i}]}
	for n in ${NumEvents[@]}; do
		for ((j=0; j<${#RefQuantPrefix[@]}; j++)); do
			ID=${RefQuantID[${j}]}

			echo ${prepdir}/simu_${t}_${ID}_${n}

			# Aligning reads to reference genome
			if [[ ! -e ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/Aligned.sortedByCoord.out.bam ]]; then
				echo -e "\tAligning reads to reference genome"
				mkdir -p ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/
				STAR --runThreadN 4 --genomeDir ${prepdir}/StarIndex/ --readFilesIn ${prepdir}/data/simulation/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz ${prepdir}/data/simulation/polyester_${t}_${ID}_${n}/sample_01_2.fasta.gz --readFilesCommand gunzip -c --outFileNamePrefix ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --sjdbGTFfile ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_annotation.gtf --limitBAMsortRAM 32416692217
			fi

			# transcriptome assembly: stringtie
			if [[ ! -e ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtiegenes.gtf ]]; then
				echo -e "\tstringtie assembly"
				stringtie ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/Aligned.sortedByCoord.out.bam -G ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_annotation.gtf -o ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtiegenes.gtf
			fi

			# transcriptome assembly: scallop
			if [[ ! -e ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallopgenes.gtf ]]; then
				echo -e "\tscallop assembly"
				scallop --verbose 0 -i ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/Aligned.sortedByCoord.out.bam -o ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallopgenes.gtf
			fi

			# summarize assembly result
			if [[ ! -e ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtie_multiexon_novelisoforms.txt ]] || [[ ! -e ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallop_multiexon_novelisoforms.txt ]]; then
				echo -e "\tsummarizing transcript assembly result for stringtie"
				gffcompare -o ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtiegffcomp -r ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_annotation.gtf ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtiegenes.gtf
				${scriptdir}/../bin/assemblypost ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_annotation.gtf ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtiegenes.gtf ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtiegffcomp.stringtiegenes.gtf.tmap ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtie
				echo -e "\tsummarizing transcript assembly result for scallop"
				gffcompare -o ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallopgffcomp -r ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_annotation.gtf ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallopgenes.gtf
				${codedir}/bin/assemblypost ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_annotation.gtf ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallopgenes.gtf ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallopgffcomp.scallopgenes.gtf.tmap ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallop
			fi

			# evaluation assembly result
			if [[ ! -e ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/evaluation_scallop_existing_novelisoforms_existjunc ]] || [[ ! -e ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/evaluation_stringtie_existing_novelisoforms_existjunc ]]; then
				echo -e "\tevaluating transcript assembly result for stringtie"
				python3 ${codedir}/scripts/Pred_Trans2Gene.py ${gtffile} ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtie_existing_novelisoforms.txt ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtie_existing_novelisoforms_uniqgene.txt
				python3 ${scriptdir}/Simulation.py evatrans ${prepdir}/data/simulation/groundtruth_${t}_${ID}_${n}/removelist_existjunc.txt ${prepdir}/data/simulation/groundtruth_${t}_${ID}_${n}/fusion_existjunc.txt ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/stringtie_existing_novelisoforms_uniqgene.txt ${gtffile} ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/evaluation_stringtie_existing_novelisoforms_existjunc

				echo -e "\tevaluating transcript assembly result for scallop"
				python3 ${codedir}/scripts/Pred_Trans2Gene.py ${gtffile} ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallop_existing_novelisoforms.txt ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallop_existing_novelisoforms_uniqgene.txt
				python3 ${scriptdir}/Simulation.py evatrans ${prepdir}/data/simulation/groundtruth_${t}_${ID}_${n}/removelist_existjunc.txt ${prepdir}/data/simulation/groundtruth_${t}_${ID}_${n}/fusion_existjunc.txt ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/scallop_existing_novelisoforms_uniqgene.txt ${gtffile} ${prepdir}/data/simulation/Star/${t}_${ID}_${n}/evaluation_scallop_existing_novelisoforms_existjunc
			fi

			# salmon indexing
			if [[ ! -e ${prepdir}/data/simulation/IndexSalmon/${t}_${ID}_${n}/hash.bin ]]; then
				echo -e "\tsalmon indexing"
				mkdir -p ${prepdir}/data/simulation/IndexSalmon/${t}_${ID}_${n}/
				salmon index -t ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_reference.fa -i ${prepdir}/data/simulation/IndexSalmon/${t}_${ID}_${n}/
			fi

			# salmon quantification
			if [[ ! -e ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/quant.sf ]]; then
				echo -e "\tsalmon quantification"
				mkdir -p ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/
				salmon quant -p 4 -l A -i ${prepdir}/data/simulation/IndexSalmon/${t}_${ID}_${n}/ -1 ${prepdir}/data/simulation/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz -2 ${prepdir}/data/simulation/polyester_${t}_${ID}_${n}/sample_01_2.fasta.gz --gcBias --seqBias --posBias --dumpEqWeights --numBootstraps 100 -o ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/ --writeMappings=${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/mapping.sam
				samtools view -Shb ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/mapping.sam -o ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/mapping.bam
			fi

			# expression difference: comparing simulated reads info and salmon quant
			if [[ ! -e ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/expression_difference.txt ]]; then
				echo -e "\tcomparing simulated ground truth with salmon quant"
				python3 ${scriptdir}/GroupTruthExp.py ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_target.fa ${prepdir}/data/simulation/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/quant.sf ${gtffile} ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/
			fi

			# run SADpipe
			if [[ ! -e ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue.tsv ]]; then
				python3 ${codedir}/scripts/SADpipe.py -t ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_reference.fa -a ${prepdir}/data/simulation/simu_${t}_${ID}_${n}_annotation.gtf -s ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/ -o ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad
			fi

			# sorting and adjust pvalues
			if [[ -e ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue.tsv ]] && [[ ! -e ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv ]]; then
				echo -e "\tsorting SAD unadjustable anomalies"
				python3 ${codedir}/scripts/AdjustPValue.py 1 ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue.tsv ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv
			fi

			# convert to gene-level predictions
			if [[ -e ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv ]] && [[ ! -e ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv ]]; then
				echo -e "\tconverting to gene-level predictions"
				python3 ${codedir}/scripts/Pred_Trans2Gene.py ${gtffile} ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue_sorted.tsv ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv
			fi

			# evaluate SAD new isoform prediction
			if [[ -e ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv ]] && [[ ! -e ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/evaluation_genelevel_overall_exist ]]; then
				echo -e "\tevaluating SAD new isoform prediction"
				python3 ${codedir}/scripts/Pred_Trans2Gene.py ${gtffile} ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_pvalue_overall_sorted ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_pvalue_overall_sorted_uniqgene
				python3 ${scriptdir}/Simulation.py evatrans ${prepdir}/data/simulation/groundtruth_${t}_${ID}_${n}/removelist_existjunc.txt ${prepdir}/data/simulation/groundtruth_${t}_${ID}_${n}/fusion_existjunc.txt ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/sad_unadjustable_pvalue_sorted_uniqgene.tsv ${gtffile} ${prepdir}/data/simulation/salmon/${t}_${ID}_${n}/${SADFolder}/evaluation_genelevel_overall_exist
			fi

		done
	done
done
