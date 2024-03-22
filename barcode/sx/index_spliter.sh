#!/bin/bash

# usage
# index_spliter.sh fasta_genome csvfiles

scriptPath="$(dirname $(realpath $0))"
fasta_genome=$1
shift
for csvfile in "$@"
do
    paste -d "" <("$scriptPath/get_primer.sh" <"$csvfile") <("$scriptPath/get_barcode.sh" <"$csvfile") | paste -d "\n" <("$scriptPath/get_faHead.sh" <"$csvfile") - >"$csvfile.primer+barcode.fa"
    bowtie2-build -q "$csvfile.primer+barcode.fa" "$csvfile.primer+barcode"

    paste -d "" <("$scriptPath/get_adapter.sh" <"$csvfile") <("$scriptPath/get_sgRNA.sh" <"$csvfile") <("$scriptPath/get_scaffold.sh" <"$csvfile") | paste -d "\n" <("$scriptPath/get_faHead.sh" <"$csvfile") - >"$csvfile.adapter+sgRNA+scaffold.fa"
    bowtie2-build -q "$csvfile.adapter+sgRNA+scaffold.fa" "$csvfile.adapter+sgRNA+scaffold"

    "${scriptPath}/get_sgRNA.sh" <"${csvfile}" >"${csvfile}.sgRNA"

    "${scriptPath}/get_ref.sh" "${csvfile}" ${fasta_genome} >"${csvfile}.ref12"
done