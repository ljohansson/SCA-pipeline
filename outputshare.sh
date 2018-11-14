#!bin/sh

#This program takes the output file and combines it with the input file to replace the lines it added information to and keep the other lines that were not found in non-coding regions. 
input=$1
wholefile=$2
header=$3
output=$4
filename=$5
tmpdir=$6

#The databases required the chr to be added before the chromosome number but to keep the original file structure it is removed here. 
#The positions of the analysed variants are extracted and all the lines not matching those positions are found. Then the analysed file is combined with the adapted input file.
#After this the header is added back on.
sed 's/chr//' $input | sort -m > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
cut -d$'\t' -f1,2 ${tmpdir}output.txt > ${tmpdir}excludefile.tmp
grep -v -f ${tmpdir}excludefile.tmp $wholefile > ${tmpdir}excludedfile.tmp
cat  ${tmpdir}excludedfile.tmp ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
sed 's/chr//' $input | sort > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
cat $header ${tmpdir}output.txt > ${output}${filename}_output.vcf

