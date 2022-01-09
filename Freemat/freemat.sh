##### copyright 2022 yyzou. All rights reserved.
##### date : 2022/1/7

input_bam=$1
input_bed=$2

out_name=$3
init_chroms_size=$4  ### for examplt : hg19.len, chrom_ID length
juicer_path=$5   ### for example : xxxxx/juicer_tools.jar

########################################################################## step 1: subdata ###########################################################################################
bedtools intersect -a ${input_bam} -b ${input_bed} >${out_name}.bam ;

samtools view -q 30 -h  ${out_name}.bam  |  pairtools parse -c ${init_chroms_size}   -o ${out_name}.pairsam.gz --assembly hg19;

if [ ! -d "tmp/" ];
then
	 mkdir tmp
fi

pairtools sort -o ${out_name}.sorted.pairsam.gz --memory 20G --compress-program gzip --tmpdir tmp ${out_name}.pairsam.gz;

pairtools select --chrom-subset /home/yyzou/hint/ref/hint-hg19/hg19.len -o ${out_name}_valid.sorted.pairsam.gz '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' ${out_name}.sorted.pairsam.gz; 

pairtools dedup --output  ${out_name}_valid.sorted.deduped.pairsam.gz ${out_name}_valid.sorted.pairsam.gz;

pairtools split --output-pairs  ${out_name}_merged_valid.pairs   ${out_name}_valid.sorted.deduped.pairsam.gz;

bgzip -f ${out_name}_merged_valid.pairs
pairix -f ${out_name}_merged_valid.pairs.gz

###
python pairix_prepare.py  ${input_bed} ${out_name} ;### generate : ${out_name}_pairix_region.txt ${out_name}.size

### python generate the pairix_region file
pairix -L  ${out_name}_merged_valid.pairs.gz  ${out_name}_pairix_region.txt  > ${out_name}_pairix_region_pairs.txt

############################################################### step 2. rename and position update  #####################################################################################################
python updata_position.py ${out_name}_pairix_region_pairs.txt   ${input_bed} ${out_name}_pairix_region_updated_pairs.txt    ### generate : ${out_name}_pairix_region_updated_pairs.txt
 
#################################################################### step 3. .hic file  generate   #################################################################################################
### add header
cat <(sed 's/^/#chromsize:\t/' ${out_name}.size ) <(echo "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"|sed 's/ /\t/g') <(cat ${out_name}_pairix_region_updated_pairs.txt) > ${out_name}_pairs_final.txt

pairtools sort -o ${out_name}_juicer.txt.gz --memory 20G --compress-program gzip --tmpdir tmp  ${out_name}_pairs_final.txt ;

java -jar ${juicer_path}  pre -n  ${out_name}_juicer.txt.gz   ${out_name}.hic ${out_name}.size

### if need to balance
#java -jar ${juicer_path}  addNorm -d -w 5000 -F  ${out_name}.hic 

rm ${out_name}.sorted.pairsam.gz ${out_name}_valid.sorted.pairsam.gz ${out_name}_valid.sorted.deduped.pairsam.gz ${out_name}_merged_valid.pairs.gz ${out_name}_merged_valid.pairs.gz.px2
rm ${out_name}_pairix_region_updated_pairs.txt ${out_name}_pairs_final.txt   ${out_name}_pairix_region.txt 

