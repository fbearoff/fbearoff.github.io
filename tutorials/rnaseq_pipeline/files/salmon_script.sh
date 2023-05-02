#!/bin/bash

##EDIT ME
read_files="location/of/your/rnaseq/reads/YOUR_PROJECT_NAME/YOUR_PROJECT_ID/00_fastq"
#EXAMPLE read_files="$HOME/data/rnaseq/YOUR_PROJECT_NAME/YOUR_PROJECT_ID/00_fastq"
##EDIT ME

#supply these files from ensembl, place in the directory where this script is located
genome="*.primary_assembly.fa.gz"
transcriptome="*.cdna.all.fa.gz"

##LEAVE ALONE
threads="$(grep -c ^processor /proc/cpuinfo)"
parentdir="$(dirname "$genome")"

grep "^>" <(gunzip -c "$genome") | cut -d " " -f 1 >"$parentdir"/decoys.txt
sed -i.bak -e 's/>//g' "$parentdir"/decoys.txt

cat "$transcriptome" "$genome" >"$parentdir"/gentrome.fa.gz

salmon index -t "$parentdir"/gentrome.fa.gz -d "$parentdir"/decoys.txt -p "$threads" -i "$parentdir"/salmon_index

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
