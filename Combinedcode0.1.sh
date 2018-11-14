#!/bin/bash

#################################INPUTFILES AND DATASETS#######################################################
encodedata=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/wgEncodeRegTfbsClusteredV3.bed
genexprss=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct
tissuelist=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/tissuelist.txt
Jaspardatafile=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/Jasparbed/
Jasparmatrixfile=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/Jasparmatrix/
Ensembl=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/Ensembler_lexi

input=$1
tissue=$2

output="IWAN_output/"
tmpdir="IWAN_output/tmp/"

#Creating tmp and output directories to store created files.
if [ ! -d "$output" ]; then mkdir IWAN_output ; fi
if [ ! -d "$tmpdir" ]; then mkdir IWAN_output/tmp ; fi

module load BEDTools

###############################################################################################################

###################################HELP AND INFORMATION########################################################

function showHelp() {
        #
        # Display commandline help on STDOUT.
        #
	cat <<EOH
===============================================================================================================
Script to analyze variants in transcription factor binding sites.
Usage:
	
	



===============================================================================================================
EOH
	trap - EXIT
        exit 0
}


###############################################################################################################

##############################ANALYSISSCRIPT###################################################################


function Analysis_script () {
	inputfile=$1
	header=$2
#	numberofcolumns=$3
#	tmpdir=$4

	bedtools sort -i $inputfile > ${tmpdir}tmp && mv ${tmpdir}tmp $inputfile
	cat $header $inputfile > ${tmpdir}bedtoolsfile.vcf
	filename=$(basename -- "$file")
	filenamefile="${filename%.*}"

	for f in ${Jaspardatafile}MA*; do
	        filename=$(basename -- "$f")
		filename="${filename%.*}"
### Can be done in one line with basename $f .txt(?)

		bedtools intersect -wa -a ${tmpdir}bedtoolsfile.vcf -wb -b $f > ${tmpdir}${filename}_output.tmp
		if [ -s ${tmpdir}${filename}_output.tmp ];
		then
	                while IFS='' read -r line || [[ -n "$line" ]]; do
	                        local columnstartpos=$((numberofcolumns + 2))
	                        local columnstartposneg=$((numberofcolumns + 3))
	                        local columnposstring=$(( numberofcolumns + 6 ))
	                        local posornegstring=$( echo $line | awk -v PosorNeg="$columnposstring" '{ print $PosorNeg}' )
	                        if [ "$posornegstring" == "-" ]
	                        then
	                                MutationLocation=$(echo "$line" | awk -v val=$columnstartposneg '{ print $val-$2+2}')
	                        else
	                                MutationLocation=$(echo "$line" | awk -v val=$columnstartpos '{ print $2-$val+1}')
	                        fi

	                        if [ "$MutationLocation" -eq 1 ] ; then continue ; fi
				line="${line} ${filename}"
	                        Matrixfile=${Jasparmatrixfile}${filename}.jaspar
	                        Matrixgenename=$(sed '5q;d' $Matrixfile)

	                        sumofmatrix=$(awk '{s+=$2}END{print s}' $Matrixfile)
	                        valueofmut=$(echo $line | cut -d' ' -f5)
	                        valueoforiginal=$(echo $line | cut -d' ' -f4)
	                        #These lines determine the length of the mutations.
	                        lengthvalmut=$(echo $valueofmut | wc -m)
	                        lengthvalorg=$(echo $valueoforiginal | wc -m)
	                        if [ $lengthvalmut -gt 2 ] || [ $lengthvalorg -gt 2 ]
	                        then
	
	                                line="${line} . . ${Matrixgenename}"
	                                echo $line >> ${tmpdir}output.txt
				else
	                                Columntobeadded=$(awk -v MutLoc="$MutationLocation" '{ print $1 $MutLoc }' $Matrixfile)
	                                if [ "$posornegstring" == "-" ]
	                                then
	                                        valueofmut=$( echo "$valueofmut" | tr ACGT TGCA )
	                                        valueoforiginal=$( echo "$valueoforiginal" | tr ACGT TGCA )
	                                fi
	                                valofmut=$( echo "$valueofmut" | tr ACGT 1234 )
	                                valoforiginal=$( echo "$valueoforiginal" | tr ACGT 1234 )
	                                #These 2 lines pull the row and column from the matrix to get the values of both spots in the matrix
	                                Newval=$(awk -v MutLoc=$MutationLocation -v Val=$valofmut 'FNR == Val {print $MutLoc}' $Matrixfile)
	                                Oldval=$(awk -v MutLoc=$MutationLocation -v Val=$valoforiginal 'FNR == Val {print $MutLoc}' $Matrixfile)
	                                if [ -z "$Newval"  ]; then continue ; fi

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
					Combinedcolumn="$(echo $Columntobeadded | tr ' ' ':' )"
	                                Combinedcolumnfix="$(echo $Combinedcolumn | cut -d':' -f1-4 )"
	
	                                line="${line} ${Combinedcolumnfix} ${newmutationvalue} ${Matrixgenename}"
	
	                                echo $line >> ${tmpdir}output.txt
	                        fi
	                done < ${tmpdir}${filename}_output.tmp
	        fi
	done


	tr ' ' \\t < ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
	sed -i '1{/^[[:space:]]*$/d}' ${tmpdir}output.txt
	bedtools sort -i ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
	cat $header ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp  ${tmpdir}output.vcf
	bedtools closest -wb -b $Ensembl -a ${tmpdir}output.vcf -t last > ${tmpdir}Genelist.tmp
	columnofstartpos=$((numberofcolumns+12))
	columnofgenes=$((numberofcolumns+15))
	cut -d$'\t' -f${columnofstartpos},${columnofgenes} ${tmpdir}Genelist.tmp > ${tmpdir}genelistfile.tmp
	paste ${tmpdir}output.txt ${tmpdir}genelistfile.tmp > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt




}
###############################################################################################################

##########################################VALIDATIONSCRIPT#####################################################
function Validation_script () {

	input=$1
	header=$2
	numberofcolumns=$3
	tmpdir=$4

	cat $header $input > ${tmpdir}filewithheader.vcf

	bedtools intersect -a ${tmpdir}filewithheader.vcf -wb -b $encodedata > ${tmpdir}encodedatafile.tmp

	rm ${tmpdir}filewithheader.vcf

	columnfoundsite=$((numberofcolumns+10))
	columnvalidatedsite=$((numberofcolumns+16))
	while read -r line || [[ -n "$line" ]]; do
	        Foundsite=$(echo $line | cut -d$' ' -f${columnfoundsite})
	        Validatedsite=$(echo $line | cut -d$' ' -f${columnvalidatedsite})
	        if [ $Foundsite == $Validatedsite ]
	        then
	                line="${line}   Validated"
	                echo $line >> ${tmpdir}validated_output.tmp
	        else
	                line="${line}   ."
	                echo $line >> ${tmpdir}nonvalidated_output.tmp
	        fi
	done < ${tmpdir}encodedatafile.tmp

	for pos in $(cut -d$'\t' -f2 $input); do
	        grep $pos ${tmpdir}nonvalidated_output.tmp > ${tmpdir}nonval.tmp
	        val=$(grep $pos ${tmpdir}validated_output.tmp)
	        if [ -z "$val" ]
	        then
	                head -1 ${tmpdir}nonval.tmp >>  ${tmpdir}nonvalidated.tmp
		fi	
	done

	cat ${tmpdir}validated_output.tmp ${tmpdir}nonvalidated.tmp > ${tmpdir}Combined_output.tmp
	tr ' ' \\t < ${tmpdir}Combined_output.tmp > ${tmpdir}validatedMain.tmp
	awk '!x[$0]++' ${tmpdir}validatedMain.tmp > ${tmpdir}output.txt
 
}

###############################################################################################################

######################################TISSUETYPES##############################################################

function Tissue_types () {

	input=$2
	tissuetypes=$1
	numberofcolumns=$3
	tmpdir=$4

	genexprss=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct
	tissuelist=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/tissuelist.txt

	columnwithgenename=$((numberofcolumns + 12))
	columnwithTFgenename=$((numberofcolumns + 10))

	for number in $( echo $tissuetypes | sed 's/,/ /g' ); do
	        numb=$(($number+2))
	        tissue=$(awk -v var=$numb 'NR==var' $tissuelist)
	        for gene in $(cut -d$'\t' -f${columnwithgenename} $input); do
	                grep  $'\t'${gene}$'\t' $genexprss | cut -d$'\t' -f${numb} > ${tmpdir}geneexpressionvalue.tmp
	                sed "s/^/${tissue}=/" ${tmpdir}geneexpressionvalue.tmp >> ${tmpdir}geneexprssval.tmp
	        done
	        if [ -f ${tmpdir}Geneexprssval.tmp ];
	        then
	                paste ${tmpdir}Geneexprssval.tmp ${tmpdir}geneexprssval.tmp > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}Geneexprssval.tmp
	        else
	                paste ${tmpdir}geneexprssval.tmp > ${tmpdir}Geneexprssval.tmp
	        fi
	        tr '\t' '/' < ${tmpdir}Geneexprssval.tmp > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}Geneexprssval.tmp
	        rm ${tmpdir}geneexprssval.tmp
	        for gene in $(cut -d$'\t' -f${columnwithTFgenename} $input); do
	                grep  $'\t'${gene}$'\t' $genexprss | cut -d$'\t' -f${numb} > ${tmpdir}genebindexpressionvalue.tmp
	                sed "s/^/TF${tissue}=/" ${tmpdir}genebindexpressionvalue.tmp  >> ${tmpdir}genebindexprssval.tmp
	        done
	                if [ -f ${tmpdir}Genebindexprssval.tmp ];
	        then
	                paste ${tmpdir}Genebindexprssval.tmp ${tmpdir}genebindexprssval.tmp > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}Genebindexprssval.tmp
	        else
	                paste ${tmpdir}genebindexprssval.tmp > ${tmpdir}Genebindexprssval.tmp
	        fi
	        tr '\t' '/' < ${tmpdir}Genebindexprssval.tmp > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}Genebindexprssval.tmp
	        rm ${tmpdir}genebindexprssval.tmp
	done

	paste ${tmpdir}Genebindexprssval.tmp ${tmpdir}Geneexprssval.tmp > ${tmpdir}combinedexprssval.tmp
	paste $input ${tmpdir}combinedexprssval.tmp > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt

}

###############################################################################################################

#####################################DUPLICATIONFIX############################################################

function Duplication_fix () {

	input=$1
	numberofcolumns=$2
	tmpdir=$3

	awk '{print NF}' $input | sort -nu | tail -n 1 > ${tmpdir}numberofcolumnsinfile.tmp
	numberofcolumnsinfile=$(cat ${tmpdir}numberofcolumnsinfile.tmp)

	for dup in $(cut -d$'\t' -f2 $input | sort -m | uniq -d); do
	        grep $dup $input > ${tmpdir}line.tmp
	        numberoflines=$(wc -l < ${tmpdir}line.tmp)
	        Addedcolumns=""
	        for column in $(seq ${numberofcolumnsinfile}); do
	                Column=$(cut -d$'\t' -f${column} ${tmpdir}line.tmp)

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
	
	
	for nodup in $(cut -d$'\t' -f2 $input | sort -m |  uniq -u); do
	        grep $nodup $input >> ${tmpdir}noduplicatedlines.tmp
	done
	
	
	cat ${tmpdir}duplicatedline.tmp ${tmpdir}noduplicatedlines.tmp > ${tmpdir}output.txt
	tr ' ' \\t < ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
	
	column1=$((numberofcolumns+4))
	column2=$((numberofcolumns+7))
	column3=$((numberofcolumns+12))
	column4=$((numberofcolumns+21))
	cut -d$'\t' -f${column1},${column2}-${column3},${column4}- ${tmpdir}output.txt | tr '\t' '|' > ${tmpdir}tempoutput.tmp
	cut -d$'\t' -f1-${numberofcolumns} ${tmpdir}output.txt > ${tmpdir}input.tmp
	paste ${tmpdir}input.tmp ${tmpdir}tempoutput.tmp > ${tmpdir}output.txt

}

###############################################################################################################

#####################################OUTPUTSHARE###############################################################

function Output_share () {

	input=$1
	wholefile=$2
	header=$3
	output=$4
	filename=$5
	tmpdir=$6

	sed 's/chr//' $input | sort -m > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
	cut -d$'\t' -f1,2 ${tmpdir}output.txt > ${tmpdir}excludefile.tmp
	grep -v -f ${tmpdir}excludefile.tmp $wholefile > ${tmpdir}excludedfile.tmp
	cat  ${tmpdir}excludedfile.tmp ${tmpdir}output.txt > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
	sed 's/chr//' $input | sort > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt
	cat $header ${tmpdir}output.txt > ${output}${filename}_output.vcf

}

###############################################################################################################

########################################MAIN###################################################################

input=$1
tissue=$2

#This will store the file in the input folder as to leave the original input file intact.

filename=$(basename -- "$input")


#This will check if the file is compressed and it will decompress if needed.
if [[ $input =~ \.t?gz$ ]];
then
        echo "Unzipping data.."
        cp $input ${tmpdir}
        mv ${tmpdir}$filename ${tmpdir}inputfile.vcf.gz
        gunzip ${tmpdir}inputfile.vcf.gz
else
        cp $input ${tmpdir}
        mv ${tmpdir}$filename ${tmpdir}inputfile.vcf
fi

fixedinput=${tmpdir}inputfile.vcf


#Here the header is stored in a file to be added back to the file in later stages.
echo "Analyzing input.."
grep '#' $fixedinput > ${tmpdir}header.txt

#Here the all the data without the header is taken from the file and "chr" is added infront of the chromosome number as it is required by certain databases.

grep -v '#'  $fixedinput | grep -v 'protein_coding' | tr ' ' '_' > ${tmpdir}startfile.txt
grep -v '#' $fixedinput | grep '5_prime_UTR'  | tr ' ' '_' >> ${tmpdir}startfile.txt

grep 'chr' $fixedinput > ${tmpdir}test.tmp
if [ -s ${tmpdir}test.tmp ];
then
        :
else

        sed -i 's/^/chr/' ${tmpdir}startfile.txt
fi


#sed -i 's/^/chr/' ${tmpdir}startfile.txt

#This takes the total number of columns in the file to ensure no problems when files have different number of columns.
awk '{print NF}' ${tmpdir}startfile.txt | sort -nu | tail -n 1 > ${tmpdir}numberofcolumns.tmp
numberofcolumns=$(cat ${tmpdir}numberofcolumns.tmp)

#This gives a default output folder if non is chosen.
#if [[ ${output} == 0 ]];
#then
#        echo "No output folder given, choosing default: /groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/output/"
#        output=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/output/
#else
#        echo "Your output folder is ${output}."
#fi


#This loads bedtools in the enviroment for multiple steps and starts the first analysis.
module load BEDTools
echo "Starting analysis.."

Analysis_script ${tmpdir}startfile.txt  ${tmpdir}header.txt $numberofcolumns ${tmpdir}

#This loads the validationscript
echo "Validating analysis.."
Validation_script ${tmpdir}output.txt ${tmpdir}header.txt $numberofcolumns ${tmpdir}

#This checks if there are any tissues indicated and will run the tissues script for additional annotation.
if [ -z "$tissue" ]
then
        echo "No tissues assigned.."
else
        echo "Adding tissue data.."
        Tissue_types $tissue ${tmpdir}output.txt $numberofcolumns ${tmpdir}
fi


#This starts the duplication fix script to remove additional lines.
echo "Fixing duplicated lines.."
Duplication_fix ${tmpdir}output.txt $numberofcolumns ${tmpdir}


#This recombines the found variants with the input file and stores it in a new output file.
echo "Adding info to input file.."
Output_share ${tmpdir}output.txt $fixedinput ${tmpdir}header.txt $output $filename ${tmpdir}

#This clears the tmp directory and removes any leftover .txt files it also shows the modules currently loaded.
rm ${tmpdir}/*
module list

echo "Analysis completed, have a nice day!"

###############################################################################################################

#######################################FINISHING UP############################################################

## deleting files created in the program.

rm -rf $tmpdir
###############################################################################################################
