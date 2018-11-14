#!bin/sh

input=$2
tissuetypes=$1
numberofcolumns=$3
tmpdir=$4

#This code gets the expression values for genes in the tissues that are available based on the dataset.
#Required datasets
genexprss=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct
tissuelist=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/tissuelist.txt

columnwithgenename=$((numberofcolumns + 12))
columnwithTFgenename=$((numberofcolumns + 10))

#This goes over the list of numbers that are given as input and pulls the result from the expression datafile. 
for number in $( echo $tissuetypes | sed 's/,/ /g' ); do
	#This corrects for the first 2 information columns
	numb=$(($number+2))
	tissue=$(awk -v var=$numb 'NR==var' $tissuelist)
	#This loop takes the gene and grabs the expression value from the file.
	for gene in $(cut -d$'\t' -f${columnwithgenename} $input); do
        	grep  $'\t'${gene}$'\t' $genexprss | cut -d$'\t' -f${numb} > ${tmpdir}geneexpressionvalue.tmp
		sed "s/^/${tissue}=/" ${tmpdir}geneexpressionvalue.tmp >> ${tmpdir}geneexprssval.tmp
	done
	#This statement ensures that every gene gets its own column to keep the information clear. 
	if [ -f ${tmpdir}Geneexprssval.tmp ];
	then 	
		paste ${tmpdir}Geneexprssval.tmp ${tmpdir}geneexprssval.tmp > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}Geneexprssval.tmp
	else
		paste ${tmpdir}geneexprssval.tmp > ${tmpdir}Geneexprssval.tmp
	fi
	#This is a repeat of the code for the gene to do the same with the transscription factor.
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

#This combines the output and then adds it to the file.
paste ${tmpdir}Genebindexprssval.tmp ${tmpdir}Geneexprssval.tmp > ${tmpdir}combinedexprssval.tmp
paste $input ${tmpdir}combinedexprssval.tmp > ${tmpdir}tmp && mv ${tmpdir}tmp ${tmpdir}output.txt

