#!bin/sh

#This is a currently disabled optional path of the pipeline that finds the gnomad data corresponding with each found variant. It is currently very slow in implementation and therefor not active. 
#This data is useful to the pipeline, because it shows how common the found variants are, which helps in determining harmful mutations.

#database for the gnomAD files.
gnomAD=/apps/data/gnomAD/release-170228/vcf/genomes/r2.0.2/
#Input is just the output file of the previous steps.
input=$1

module load parallel
export LC_ALL=C

#parallel --pipe --block 2M zgrep foo < bigfile



#This loop takes the position from the found variants and the chromosome and finds the corresponding file and then searches for the line containing the variant. Potential slowdown caused by the line not being in the file,
#which, while prefered, would mean the program goes through the whole 1-8GB file. 
for pos in $(cut -d$'\t' -f2 $input | sort -m | uniq -u); do
	grep $pos $input > line.tmp
	chr=$(head -1 line.tmp | cut -d$'\t' -f1)
	echo "precheck"
	parallel --pipe --block 2M zgrep $pos < ${gnomAD}gnomad.genomes.r2.0.2.sites.${chr}.normalized.vcf.gz > greppedline.tmp
	echo "aftercheck"
	gnomADval=$(cut -d$'\t' -f8 greppedline.tmp | cut -d ';' -f1,2)
	#gnomADval=$(zgrep $pos ${gnomAD}gnomad.genomes.r2.0.2.sites.${chr}.normalized.vcf.gz | cut -d$'\t' -f8 | cut -d ';' -f1,2)
	echo $gnomADval
	if [ -z "$gnomADval" ];
	then
		gnomADval="AC=0;AF=0"
	fi
	echo $gnomADval > gnomad.tmp
	paste line.tmp gnomad.tmp >> output.tmp
done



for pos in $(cut -d$'\t' -f2 $input | sort -m | uniq -d); do
        grep $pos $input > line.tmp
        chr=$(head -1 line.tmp | cut -d$'\t' -f1)
	parallel --pipe --block 2M zgrep $pos < ${gnomAD}gnomad.genomes.r2.0.2.sites.${chr}.normalized.vcf.gz > greppedline.tmp
	gnomADval=$(cut -d$'\t' -f8 greppedline.tmp | cut -d ';' -f1,2)
	#gnomADval=$(zgrep $pos ${gnomAD}gnomad.genomes.r2.0.2.sites.${chr}.normalized.vcf.gz | cut -d$'\t' -f8 | cut -d ';' -f1,2)
	if [ -z "$gnomADval" ];
	then
                gnomADval="AC=0;AF=0"
        fi
        echo $gnomADval > gnomad.tmp
	paste line.tmp gnomad.tmp >> output.tmp
done

#This ensures the current output is saved on the file. 
mv output.tmp gnomadtestoutput.txt

rm *.tmp
