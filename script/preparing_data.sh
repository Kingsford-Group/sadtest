#!/bin/bash

args=("$@")
if [[ ${#args[@]} == 1 ]]; then
	echo "Usage: preparing_data.sh <output directory>"
	exit
fi

indir=$0
indir=${indir%/*}
outdir=$1

# download reference and build index
if [[ ! -e ${outdir}/gencode.v26.annotation.gtf ]]; then
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz | gunzip > ${outdir}/gencode.v26.annotation.gtf
fi
if [[ ! -e ${outdir}/GRCh38.p10.genome.fa ]]; then
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.p10.genome.fa.gz | gunzip > ${outdir}/GRCh38.p10.genome.fa
fi
if [[ ! -e ${outdir}/gencode.v26.transcripts.fa ]]; then
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz | gunzip > ${outdir}/gencode.v26.transcripts.fa
fi
if [[ ! -e ${outdir}/gencode.v26.pc.transcripts.fa ]]; then
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.pc_transcripts.fa.gz | gunzip > ${outdir}/gencode.v26.pc.transcripts.fa
fi
if [[ ! -e ${outdir}/StarIndex/Genome ]]; then
	mkdir -p ${outdir}/StarIndex
	STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ${outdir}/StarIndex/ --genomeFastaFiles GRCh38.p10.genome.fa && mv Log.out ${outdir}/StarIndex/
fi
if [[ ! -e ${outdir}/gencode.v26.full/hash.bin ]]; then
	salmon index --gencode -t gencode.v26.transcripts.fa -i ${outdir}/gencode.v26.full
fi
if [[ ! -e ${outdir}/gencode.v26.pc/hash.bin ]]; then
	salmon index --gencode -t gencode.v26.pc.transcripts.fa -i ${outdir}/gencode.v26.pc
fi
if [[ ! -e ${outdir}/gencode.v26.Gene_Trans_Map.txt ]]; then
	awk '{if($3=="transcript") print substr($10,2,length($10)-3)"\t"substr($12,2,length($12)-3)}' ${outdir}/gencode.v26.annotation.gtf > ${outdir}/gencode.v26.Gene_Trans_Map.txt
fi

# extract protein-coding only transcript
python ${indir}/ExtractPCAnnotation.py ${outdir}/gencode.v26.pc.transcripts.fa ${outdir}/gencode.v26.annotation.gtf ${outdir}/gencode.v26.annotation.pc.gtf

# download data
cp -r ${indir}/../data ${outdir}/
echo "DOWNLOADING BASE EXPRESSION SAMPLES..."
${outdir}/data/BaseExpression/getdata.sh
echo "DOWNLOADING GEUVADIS SAMPLES..."
${outdir}/data/GEUVADIS/getdata.sh
echo "DOWNLOADING HUMAN BODY MAP SAMPLES..."
${outdir}/data/HumanBodyMap/getdata.sh
echo "DOWNLOADING 1000GP LONG READ AND SHORT READ SAMPLES..."
${outdir}/data/LongRead1000G/getdata.sh

# quantify baseline expression
if [[ ! -e ${outdir}/data/BaseExpression/ERR188297_salmon_full_quant.sf ]]; then
	salmon quant -p 4 -l A -i ${outdir}/gencode.v26.full -1 ${outdir}/data/BaseExpression/ERR188297_1.fastq.gz -2 ${outdir}/data/BaseExpression/ERR188297_2.fasta.gz --gcBias --seqBias --posBias -o ${outdir}/data/BaseExpression/ERR188297_full
	mv ${outdir}/data/BaseExpression/ERR188297_full/quant.sf ${outdir}/data/BaseExpression/ERR188297_salmon_full_quant.sf
fi
if [[ ! -e ${outdir}/data/BaseExpression/ERR188297_salmon_pc_quant.sf ]]; then
	salmon quant -p 4 -l A -i ${outdir}/gencode.v26.pc -1 ${outdir}/data/BaseExpression/ERR188297_1.fastq.gz -2 ${outdir}/data/BaseExpression/ERR188297_2.fasta.gz --gcBias --seqBias --posBias -o ${outdir}/data/BaseExpression/ERR188297_pc
	mv ${outdir}/data/BaseExpression/ERR188297_pc/quant.sf ${outdir}/data/BaseExpression/ERR188297_salmon_pc_quant.sf
fi
if [[ ! -e ${outdir}/data/BaseExpression/SRR3192396_salmon_full_quant.sf ]]; then
        salmon quant -p 4 -l A -i ${outdir}/gencode.v26.full -1 ${outdir}/data/BaseExpression/SRR3192396_1.fastq.gz -2 ${outdir}/data/BaseExpression/SRR3192396_2.fasta.gz --gcBias --seqBias --posBias -o ${outdir}/data/BaseExpression/SRR3192396_full
        mv ${outdir}/data/BaseExpression/SRR3192396_full/quant.sf ${outdir}/data/BaseExpression/SRR3192396_salmon_full_quant.sf
fi
if [[ ! -e ${outdir}/data/BaseExpression/SRR3192396_salmon_pc_quant.sf ]]; then
	salmon quant -p 4 -l A -i ${outdir}/gencode.v26.pc -1 ${outdir}/data/BaseExpression/SRR3192396_1.fastq.gz -2 ${outdir}/data/BaseExpression/SRR3192396_2.fasta.gz --gcBias --seqBias --posBias -o ${outdir}/data/BaseExpression/SRR3192396_pc
	mv ${outdir}/data/BaseExpression/SRR3192396_pc/quant.sf ${outdir}/data/BaseExpression/SRR3192396_salmon_pc_quant.sf
fi
if [[ ! -e ${outdir}/data/BaseExpression/SRR3192412_salmon_full_quant.sf ]]; then
	salmon quant -p 4 -l A -i ${outdir}/gencode.v26.full -1 ${outdir}/data/BaseExpression/SRR3192412_1.fastq.gz -2 ${outdir}/data/BaseExpression/SRR3192412_2.fasta.gz --gcBias --seqBias --posBias -o ${outdir}/data/BaseExpression/SRR3192412_full
	mv ${outdir}/data/BaseExpression/SRR3192412_full/quant.sf ${outdir}/data/BaseExpression/SRR3192412_salmon_full_quant.sf
fi
if [[ ! -e ${outdir}/data/BaseExpression/SRR3192412_salmon_pc_quant.sf ]]; then
	salmon quant -p 4 -l A -i ${outdir}/gencode.v26.pc -1 ${outdir}/data/BaseExpression/SRR3192412_1.fastq.gz -2 ${outdir}/data/BaseExpression/SRR3192412_2.fasta.gz --gcBias --seqBias --posBias -o ${outdir}/data/BaseExpression/SRR3192412_pc
	mv ${outdir}/data/BaseExpression/SRR3192412_pc/quant.sf ${outdir}/data/BaseExpression/SRR3192412_salmon_pc_quant.sf
fi

# simulating deletion, fusion events of transcriptome; simulating reads with polyester
Type=("PC" "Full")
GTFfiles=("${outdir}/gencode.v26.annotation.pc.gtf" "${outdir}/gencode.v26.annotation.gtf")
TransFastafiles=("${outdir}/gencode.v26.pc.transcripts.fa" "${outdir}/gencode.v26.transcripts.fa")
NumEvents=(200 500 1000 1500)
RefQuantPrefix=("${outdir}/data/BaseExpression/ERR188297" "${outdir}/data/BaseExpression/SRR3192396" "${outdir}/data/BaseExpression/SRR3192412")
RefQuantID=(GEU GM12878 K562)
for ((i=0; i<${#Type[@]}; i++)); do
	t=${Type[${i}]}
	gtffile=${GTFfiles[${i}]}
	transfasta=${TransFastafiles[${i}]}
	for n in ${NumEvents[@]}; do
		for ((j=0; j<${#RefQuantPrefix[@]}; j++)); do
			ID=${RefQuantID[${j}]}
			if ((i==0)); then
				refquantfile=${RefQuantPrefix[${j}]}"_salmon_pc_quant.sf"
			else
				refquantfile=${RefQuantPrefix[${j}]}"_salmon_full_quant.sf"
			fi

			echo simu_${t}_${ID}_${n}

			# simulating deletion and fusion
			if [[ ! -e ${outdir}/data/simulation/simu_${t}_${ID}_${n}_theoexp.txt ]]; then
				echo -e "\tSimulating deletion and fusion."
				mkdir -p ${outdir}/data/simulation/simu_${t}_${ID}_${n}
				python3 ${indir}/Simulation.py ${transfasta} ${gtffile} ${refquantfile} ${outdir}/data/simulation/simu_${t}_${ID}_${n} ${n} ${n}
				python3 ${indir}/ExtractPCAnnotation.py ${outdir}/data/simulation/simu_${t}_${ID}_${n}_reference.fa ${outdir}/gencode.v26.annotation.gtf ${outdir}/data/simulation/simu_${t}_${ID}_${n}_annotation.gtf
			fi

			# categorize the simulated events into with novel junctions and without novel junctions
			if [[ ! -e ${outdir}/data/simulation/groundtruth_${t}_${ID}_${n}/removelist_noveljunc.txt ]]; then
				echo -e "\tCategorizing simulated events."
				mkdir -p ${outdir}/data/simulation/groundtruth_${t}_${ID}_${n}
				${indir}/../bin/categorizesimulation ${gtffile} ${outdir}/data/simulation/simu_${t}_${ID}_${n}_removelist.txt ${outdir}/data/simulation/simu_${t}_${ID}_${n}_target.fa ${outdir}/data/simulation/groundtruth_${t}_${ID}_${n}/
			fi

			# simulating reads with polyester
			if [[ ! -e ${outdir}/data/simulation/polyester_${t}_${ID}_${n}/sample_01_1.fasta.gz ]] && [[ ! -e  ${outdir}/data/simulation/polyester_${t}_${ID}_${n}/sample_01_1.fasta ]]; then
				echo -e "\tSimulating reads with polyester."
				mkdir -p ${outdir}/data/simulation/polyester_${t}_${ID}_${n}
				#/opt/local/stow/R-3.3.2/bin/Rscript /mnt/disk33/user/congm1/Code/SAD/src/SimulationReads.R simu_${t}_${ID}_${n}_target.fa simu_${t}_${ID}_${n}_theoexp.txt polyester_${t}_${ID}_${n}/ 1
				Rscript ${indir}/SimulationReads.R ${outdir}/data/simulation/simu_${t}_${ID}_${n}_target.fa ${outdir}/data/simulation/simu_${t}_${ID}_${n}_theoexp.txt ${outdir}/data/simulation/polyester_${t}_${ID}_${n}/ 1
			fi

			# gzipping simulated reads
			if [[ $(find ${outdir}/data/simulation/polyester_${t}_${ID}_${n}/ -name *fasta) != "" ]]; then
				echo -e "\tGzipping simulated reads."
				gzip ${outdir}/data/simulation/polyester_${t}_${ID}_${n}/*fasta
			fi

		done
	done
done
