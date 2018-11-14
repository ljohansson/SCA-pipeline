#!/bin/sh

#This script removes duplicate lines by seeking lines with the same exact location in the file and adding the newly added columns together split by a "/" as separator.

input=$1
numberofcolumns=$2
tmpdir=$3

#The total number of columns in the file is determined for the combining of duplicated lines. 
awk '{print NF}' $input | sort -nu | tail -n 1 > ${tmpdir}numberofcolumnsinfile.tmp
numberofcolumnsinfile=$(cat ${tmpdir}numberofcolumnsinfile.tmp)

#This loop goes through the lines based on the position of the found variants where there is more than one line per variant. 
#It compares all the columns and combines the columns that are not the same by adding a "/" between the values. 
for dup in $(cut -d$'\t' -f2 $input | sort -m | uniq -d); do
	grep $dup $input > ${tmpdir}line.tmp
	numberoflines=$(wc -l < ${tmpdir}line.tmp)
	Addedcolumns=""
	for column in $(seq ${numberofcolumnsinfile}); do
		Column=$(cut -d$'\t' -f${column} ${tmpdir}line.tmp)
		
#This checks if the values of both line are the same in the column to prevent adding double information.
		var1=$(echo $Column | cut -d$' ' -f1)
		var2=$(echo $Column | cut -d$' ' -f2)
		if [ "$var1" == "$var2" ]
		then
			Combined=$var1
		else
			Combined="$var1/$var2" 
		fi
		Addedcolumns="${Addedcolumns} ${Combined}"
		
	done
	echo $Addedcolumns >> ${tmpdir}duplicatedline.tmp
done


#This part takes the unique lines in the file and adds them to the output file. 
for nodup in $(cut -d$'\t' -f2 $input | sort -m |  uniq -u); do
	grep $nodup $input >> ${tmpdir}noduplicatedlines.tmp 
done

#This combines the output of both for loops together

cat ${tmpdir}duplicatedline.tmp ${tmpdir}noduplicatedlines.tmp > ${tmpdir}output.txt
tr ' ' \\t < ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt

#These lines contain the locations of the important added columns from the analysis. The locations will differ based on the amount of columns in the input file and this should correct for that. 
column1=$((numberofcolumns+4))
column2=$((numberofcolumns+7))
column3=$((numberofcolumns+12))
column4=$((numberofcolumns+21))
#Here the tabs between columns are turned into | to create one column with only the non-coding information
cut -d$'\t' -f${column1},${column2}-${column3},${column4}- ${tmpdir}output.txt | tr '\t' '|' > ${tmpdir}tempoutput.tmp
cut -d$'\t' -f1-${numberofcolumns} ${tmpdir}output.txt > ${tmpdir}input.tmp
paste ${tmpdir}input.tmp ${tmpdir}tempoutput.tmp > ${tmpdir}output.txt
