#!/bin/bash

#This is a program that takes an input bed or csv file and then compares that with datasets from the JASPAR database to find transcriptionfactor binding sites. 
#When a site is found this program will open the corresponding matrix file and print the full line of scores for every base in this position. 

#input
inputfile=$1
header=$2
numberofcolumns=$3
tmpdir=$4

#variables
Jaspardatafile=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/Jasparbed/
Jasparmatrixfile=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/Jasparmatrix/
Ensembl=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/Ensembler_lexi


bedtools sort -i $inputfile > ${tmpdir}tmp && mv ${tmpdir}tmp $inputfile

cat $header $inputfile > ${tmpdir}bedtoolsfile.vcf

#bedtools sort -i ${tmpdir}bedtoolsfile.vcf > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}bedtoolsfile.vcf

filename=$(basename -- "$file")

filenamefile="${filename%.*}"


for f in ${Jaspardatafile}MA*; do
	
	#Getting the filename without the path of extention to match it with the matrix in a different file with a different extention.

	filename=$(basename -- "$f")	
	filename="${filename%.*}"

	#Perform bedtools over the input file on the given datasets one by one put it in an output file named after the dataset used.
	bedtools intersect -wa -a ${tmpdir}bedtoolsfile.vcf -wb -b $f > ${tmpdir}${filename}_output.tmp
	#Check if file is empty
	if [ -s ${tmpdir}${filename}_output.tmp ]
	then	
		
		#If the file is not empty it goes through the file line by line
		while IFS='' read -r line || [[ -n "$line" ]]; do
			#Find the location of the mutation in the file by substracting the column with the start position of the binding site by the position of the mutation. 
			#The +1 is there because the matrix file starts with the base in column 1, which makes position 1 of the mutation position 2 in the matrix.  
			columnstartpos=$((numberofcolumns + 2))
			columnstartposneg=$((numberofcolumns + 3))
			columnposstring=$(( numberofcolumns + 6 ))
			posornegstring=$( echo $line | awk -v posorneg="$columnposstring" '{ print $posorneg}' )
			if [ "$posornegstring" == "-" ]
			then
				MutationLocation=$(echo "$line" | awk -v val=$columnstartposneg '{ print $val-$2+2}')
			else 
				MutationLocation=$(echo "$line" | awk -v val=$columnstartpos '{ print $2-$val+1}')
			fi
			
			#This check removes any indicated lines that are not actually in the binding site and skips to the next iteration.
			if [ "$MutationLocation" -eq 1 ]
			then
				continue
			fi

			#Adds the name of the transciption factor binding site to the data. 
			line="${line} ${filename}"

			#Matches the right dataset with the right matrix file.
			Matrixfile=${Jasparmatrixfile}${filename}.jaspar


			#This finds the right column from the matrix file.
			Matrixgenename=$(sed '5q;d' $Matrixfile)
			
			#These lines determine the total value of the matrix per line and it pulls the specific mutation from the input file.
			sumofmatrix=$(awk '{s+=$2}END{print s}' $Matrixfile)
			valueofmut=$(echo $line | cut -d' ' -f5)
			valueoforiginal=$(echo $line | cut -d' ' -f4)
			#These lines determine the length of the mutations.
			lengthvalmut=$(echo $valueofmut | wc -m)
			lengthvalorg=$(echo $valueoforiginal | wc -m)
			#Here the length of the mutation is determined. If the mutation is more than a simple substitution no score can currently be given so this data is added to the output file without going over the score code.
			if [ $lengthvalmut -gt 2 ] || [ $lengthvalorg -gt 2 ] 
			then	

				line="${line} . . ${Matrixgenename}"
				echo $line >> ${tmpdir}output.txt

			else
				Columntobeadded=$(awk -v MutLoc="$MutationLocation" '{ print $1 $MutLoc }' $Matrixfile)
				#This substitutes the base for a row number in the matrix file
				#This checks for a negative strand and corrects the bases for a correct score. 
				if [ "$posornegstring" == "-" ]
				then
					valueofmut=$( echo "$valueofmut" | tr ACGT TGCA )
					valueoforiginal=$( echo "$valueoforiginal" | tr ACGT TGCA )
				else
					:
				fi 
				valofmut=$( echo "$valueofmut" | tr ACGT 1234 )
				valoforiginal=$( echo "$valueoforiginal" | tr ACGT 1234 )
				#These 2 lines pull the row and column from the matrix to get the values of both spots in the matrix
				Newval=$(awk -v MutLoc=$MutationLocation -v Val=$valofmut 'FNR == Val {print $MutLoc}' $Matrixfile)
				Oldval=$(awk -v MutLoc=$MutationLocation -v Val=$valoforiginal 'FNR == Val {print $MutLoc}' $Matrixfile)
				if [ -z "$Newval"  ];
				then
					continue
				fi
				 
				#Here the current way of determining the score is done.	
				Newvalue=$(echo "scale=10 ; $Newval / $sumofmatrix + 0.001" | bc)
				Oldvalue=$(echo "scale=10 ; $Oldval / $sumofmatrix + 0.001" | bc)
				totalvalueofmutation=$(echo "scale=10 ; $Oldvalue / $Newvalue" | bc)
				totvalueofmutation=${totalvalueofmutation/.*}
				if [[ "$totvalueofmutation" -ge 1 ]]
				then 
					newmutationvalue=$(echo "scale=2 ; $totalvalueofmutation / 1" | bc)
				else	
					newmutationvalue=$(echo "scale=2 ; -1 / $totalvalueofmutation" | bc)
				fi
				#Removes the whitespaces between the values of the column in the matrix

				Combinedcolumn="$(echo $Columntobeadded | tr ' ' ':' )"
				Combinedcolumnfix="$(echo $Combinedcolumn | cut -d':' -f1-4 )"

				line="${line} ${Combinedcolumnfix} ${newmutationvalue} ${Matrixgenename}"

				#Stores the line in an output file
				echo $line >> ${tmpdir}output.txt
			fi
		done < ${tmpdir}${filename}_output.tmp
	else
		:
	fi
done


#This recreates tabs in the data where it was changes into spaces. 
tr ' ' \\t < ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
#This line clears any empty lines in the file.
sed -i '1{/^[[:space:]]*$/d}' ${tmpdir}output.txt
#Here the file is sorted for the bedtools closest command, the header is added back aswell.
bedtools sort -i ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
cat $header ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp  ${tmpdir}output.vcf
bedtools closest -wb -b $Ensembl -a ${tmpdir}output.vcf -t last > ${tmpdir}Genelist.tmp
#Getting the column of genes from the analysis and adding it to the sorted file.
columnofstartpos=$((numberofcolumns+12))
columnofgenes=$((numberofcolumns+15))
cut -d$'\t' -f${columnofstartpos},${columnofgenes} ${tmpdir}Genelist.tmp > ${tmpdir}genelistfile.tmp
paste ${tmpdir}output.txt ${tmpdir}genelistfile.tmp > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt

