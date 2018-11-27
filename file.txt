# IwaN1.1
Integrated Workload Analysis Non-coding
This program annotates variants found in transcription factor binding sites.
It also gives a score based on the predicted effect of the mutation.

Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.
Download the "IWAN1.1.sh" and the "config" files from https://github.com/IJHidding/SCA-pipeline 
Add both files in the same place, then follow instructions below to install the relevant databases.
Link the locations of the databases on your server in the config file.

Encode:
From this link : http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/ 
Get the wgEncodeRegTfbsClusteredV3.bed file and put it on the preferred location. 

GTEx:
Follow this link: https://gtexportal.org/home/datasets
and then download "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz"
and unpack in the preffered location.

Tissuelist from GTEx:
The tissuelist was generated based on the columns in the GTEx database. The tissuelist for GTEX_analysis version 2016-01-15 van be found on the github page name "tissuelist.txt"

Jaspardatafile:
from this link: http://jaspar.genereg.net/downloads/ go to other data and download the BED files for genomic coodrinates of the sequences. 
Unpack this file in its own folder on the preferred location. 
Only the .bed files are required. 

Jasparmatrixfile:
The Matrixfile can be downloaded from http://jaspar.genereg.net/downloads/ , then under the JASPAR collections (PFMs) select the individual PFMs (zip) and download the Non-redundant CORE collection.
Unpack this file in its own folder on the preferred location 

Then run this script on all the downloaded files.
"
#! bin/bash

input=$1

for file in $input; do
        sed 's/\[*\]*//g' $file | tail -4 > tmp1
        head -1 $file | cut -d$'\t' -f2 > tmp2
        cat tmp1 tmp2 > $input
	rm tmp*
done
"

Ensemble:
This dataset can be found in the archives of the Ensemble website: http://feb2014.archive.ensembl.org/Homo_sapiens/Info/Index
Then on top click the BioMart option, choose Ensembl Genes 75, then Homo sapiens genes (GRCh37.p13), for filters-REGION select Chromosome 1-Y.
Under attributes select structures untag the tagged options, select in order: GENE:Chromosome Name, EXON: Exon Chr Start (bp), EXON: Chr End (bp), GENE:Strand, GENE: Associated Gene Name.
Press results and download the file as cvs. Then put it in the preferred location. 
Should this dataset become unavailable, any version should work as long as it uses GRCh37. 
Keep in mind this could lead to slightly different results based on the version. 



Prerequisites
This version of IWAN1.1 only requires BEDTools v2.25.0


Installing
After downloading the databases and the IWAN1.1 and the config files, make sure that the databases are linked in the config file and that the config and IWAN1.1 files are in the same location. 

To test if the installation is correct download the test.vcf from the github page. Then run the program with the command $: sh IWAN1.1.sh -i test.vcf -t 1
This should generate a new folder called IWAN_output in your current directory. This folder will contain the test_output.vcf file, which when opened should contain a header and 3 variant lines. 
Only one of the lines should be annotated with: hg19_chr18:12377521-12377531(-)|MA0079.3|A99:C6703:G0:T1932|-768.46|SP1|12376947|AFG3L2|Validated|TFAdipose-Subcutaneous=40.84|Adipose-Subcutaneous=35.055

Versioning
This is version 1.1

Authors
Iwan J. Hidding - Initial work

##License
License information is available on the github page (https://github.com/IJHidding/SCA-pipeline) under LICENSE, COPYING and COPYING_LESSER.


Acknowledgments
Helpfunction code was used from the CramConversion.sh made by RoanKanninga: https://github.com/molgenis/ngs-utils/blob/master/CramConversion.sh
Lennart Johansson for help with specific code and general guidance.
Cleo van Diemen for general guidance. 

