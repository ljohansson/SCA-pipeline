#!/bin/sh

#This code Validates the binding sites found in the Jaspar dataset by comparing it with Chipseq data from the Encode project
#It will add either Validated or a . after the line to show if the binding side was found for that gene.

#The database containing encode data. 
encodedata=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/wgEncodeRegTfbsClusteredV3.bed
input=$1
header=$2
numberofcolumns=$3
tmpdir=$4

#add header to input file.
cat $header $input > ${tmpdir}filewithheader.vcf

#This compares with the encode dataset

bedtools intersect -a ${tmpdir}filewithheader.vcf -wb -b $encodedata > ${tmpdir}encodedatafile.tmp

rm ${tmpdir}filewithheader.vcf


columnfoundsite=$((numberofcolumns+10))
columnvalidatedsite=$((numberofcolumns+16))
#Goes over the output file by line comparing the binding site with a found binding site, if they are the same it will add Validated to the line, if they are not the same it will add a "." 
while read -r line || [[ -n "$line" ]]; do
	Foundsite=$(echo $line | cut -d$' ' -f${columnfoundsite})
	Validatedsite=$(echo $line | cut -d$' ' -f${columnvalidatedsite})
	if [ $Foundsite == $Validatedsite ]
	then
		line="${line}	Validated" 
		echo $line >> ${tmpdir}validated_output.tmp
	else 
		line="${line}	."
                echo $line >> ${tmpdir}nonvalidated_output.tmp
	fi
done < ${tmpdir}encodedatafile.tmp

#This code checks if the binding site already has a validated site found, if there is not it will return the first line found that ends with a dot.

for pos in $(cut -d$'\t' -f2 $input); do
	grep $pos ${tmpdir}nonvalidated_output.tmp > ${tmpdir}nonval.tmp
	val=$(grep $pos ${tmpdir}validated_output.tmp)
	if [ -z "$val" ]
	then		
		head -1 ${tmpdir}nonval.tmp >>  ${tmpdir}nonvalidated.tmp
	else 
		:  
	fi
	
done

#This combines both loops into one output file.
cat ${tmpdir}validated_output.tmp ${tmpdir}nonvalidated.tmp > ${tmpdir}Combined_output.tmp
#This turns the whitespace seperator into tabs and removes any columns that do not add additional information
tr ' ' \\t < ${tmpdir}Combined_output.tmp > ${tmpdir}validatedMain.tmp
#This line removes any lines that are leftover that are entirely the same.
awk '!x[$0]++' ${tmpdir}validatedMain.tmp > ${tmpdir}output.txt


