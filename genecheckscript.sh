#!/bin/sh

input=$1
genelist=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/knownSCAgeneslist.txt

for gene in $(cat $genelist); do
	grep 'Validated' $input | grep $gene > $gene.tmp
done
wc -l *.tmp > geneinfo.txt
rm *.tmp
