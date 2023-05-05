#!/bin/bash

#supply genome and transcriptome files from Ensembl, place in the directory where this script is located
#EXAMPLE genome "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
#EXAMPLE transcriptome "Homo_sapiens.GRCh38.cdna.all.fa.gz"

# run this script with the following command `bash salmon_script.sh`

## EDIT ME if the `00_fastq` directory is not in the current directory
read_files="./00_fastq"
## EDIT ME

##LEAVE ALONE
genome="dna.primary_assembly.fa.gz"
transcriptome="cdna.all.fa.gz"
threads="$(grep -c ^processor /proc/cpuinfo)"

grep "^>" <(gunzip -c ./*"$genome") | cut -d " " -f 1 >decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat ./*"$transcriptome" ./*"$genome" >gentrome.fa.gz

salmon index -t gentrome.fa.gz -d decoys.txt -p "$threads" -i ./salmon_index

#quantify transcripts against index
destination_directory="./salmon_output"
index="./salmon_index"
for sample in "$read_files"/*_R1_001.fastq.gz; do
	base_name="${sample##*/}"
	echo "Processing sample ${base_name%%_R1_001*}"
	salmon quant \
		-i "$index" \
		--gcBias \
		-l A \
		-1 "${sample}" \
		-2 "${sample%%_R1*}_R2_001.fastq.gz" \
		-o "$destination_directory"/quants/"${base_name%%_R1*}"
	mapped=$(grep "percent" "$destination_directory"/quants/"${base_name%%_R1*}"/aux_info/meta_info.json | cut -d : -f 2 | cut -d , -f 1)
	echo "${base_name%%_R1*}" "$mapped" >>"$destination_directory"/mapped_percent.txt
done
##LEAVE ALONE
