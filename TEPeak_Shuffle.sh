module load bedtools/2.29.2
sample=$1
scriptFolder=$2
annotationFolder=$3
TEBedFile=$4
chromSizeFile=$5
familyCountFile=$6

##### shuffle peaks
i=1
mkdir $sample".shuffle/"

#mkdir shuffled_beds
echo "step1"
for i in {1..1000}               ##### loop
do

# non-desert
echo "shuffle:"$i
bedtools shuffle -i $sample"_summit.5UTROnly.bed" -g $annotationFolder${chromSizeFile} -incl $annotationFolder"5UTR.bed" >$sample"_summit.5UTROnly.shuffle.tmp.bed"
bedtools intersect -wa -u -a $annotationFolder${TEBedFile} -b $sample"_summit.5UTROnly.shuffle.tmp.bed" > $sample"_summit.5UTROnly.shuffle.TE.tmp.bed"
bedtools shuffle -i $sample"_summit.distalOnly.bed" -g $annotationFolder${chromSizeFile} -incl $annotationFolder"distal.bed" >$sample"_summit.distalOnly.shuffle.tmp.bed"
bedtools intersect -wa -u -a $annotationFolder${TEBedFile} -b $sample"_summit.distalOnly.shuffle.tmp.bed" > $sample"_summit.distalOnly.shuffle.TE.tmp.bed"
bedtools shuffle -i $sample"_summit.exonsOnly.bed" -g $annotationFolder${chromSizeFile} -incl $annotationFolder"exons.bed" >$sample"_summit.exonsOnly.shuffle.tmp.bed"
bedtools intersect -wa -u -a $annotationFolder${TEBedFile} -b $sample"_summit.exonsOnly.shuffle.tmp.bed" > $sample"_summit.exonsOnly.shuffle.TE.tmp.bed"
bedtools shuffle -i $sample"_summit.intronsOnly.bed" -g $annotationFolder${chromSizeFile} -incl $annotationFolder"introns.bed" >$sample"_summit.intronsOnly.shuffle.tmp.bed"
bedtools intersect -wa -u -a $annotationFolder${TEBedFile} -b $sample"_summit.intronsOnly.shuffle.tmp.bed" > $sample"_summit.intronsOnly.shuffle.TE.tmp.bed"
bedtools shuffle -i $sample"_summit.promotersOnly.bed" -g $annotationFolder${chromSizeFile} -incl $annotationFolder"promoters.bed" >$sample"_summit.promotersOnly.shuffle.tmp.bed"
bedtools intersect -wa -u -a $annotationFolder${TEBedFile} -b $sample"_summit.promotersOnly.shuffle.tmp.bed" > $sample"_summit.promotersOnly.shuffle.TE.tmp.bed"
bedtools shuffle -i $sample"_summit.proximalOnly.bed" -g $annotationFolder${chromSizeFile} -incl $annotationFolder"proximal.bed" >$sample"_summit.proximalOnly.shuffle.tmp.bed"
bedtools intersect -wa -u -a $annotationFolder${TEBedFile} -b $sample"_summit.proximalOnly.shuffle.tmp.bed" > $sample"_summit.proximalOnly.shuffle.TE.tmp.bed"
bedtools shuffle -i $sample"_summit.proximalOnly.bed" -g $annotationFolder${chromSizeFile} -incl $annotationFolder"tss.bed" >$sample"_summit.tssOnly.shuffle.tmp.bed"
bedtools intersect -wa -u -a $annotationFolder${TEBedFile} -b $sample"_summit.proximalOnly.shuffle.tmp.bed" > $sample"_summit.proximalOnly.shuffle.TE.tmp.bed"
# desert
bedtools shuffle -i $sample"_summit.desertOnly.bed" -g $annotationFolder${chromSizeFile} -excl $annotationFolder"non.desert.bed" >$sample"_summit.desertOnly.shuffle.tmp.bed"
bedtools intersect -wa -u -a $annotationFolder${TEBedFile} -b $sample"_summit.desertOnly.shuffle.tmp.bed" > $sample"_summit.desertOnly.shuffle.TE.tmp.bed"
# count
grep '[a-z0-9]' $sample"_summit."*".shuffle.TE.tmp.bed" | sed 's/:/ /g' |awk '{print$1,$5,$7,$8}'| sort |uniq -c | awk '{print$2,$3,$4,$5,$1}' >$sample".shuffle/"$sample"_summit.shuffle.TE."$i".counts.anno"
grep '[a-z0-9]' $sample"_summit."*".shuffle.tmp.bed" | sed 's/:/ /g' |awk '{print$1}' | sort | uniq -c >$sample".shuffle/"$sample"_summit.shuffle.TE."$i".counts.all"

truncate -s 0 $sample"_summit."*".shuffle."*".tmp.bed" 
done
rm $sample"_summit."*".shuffle"*".tmp.bed"

###### combined into the count table
echo "step2"
cp ${sample}_summit.TE.counts.anno ${sample}.shuffle/
cd ${sample}.shuffle/
ls | grep counts.anno >../${sample}.shuffle.fileList
python ${scriptFolder}Combined_TEenrichmentByTEfamily.py -i ../${sample}.shuffle.fileList -l ${annotationFolder}${familyCountFile} >../${sample}.counts

###### double check the shuffling process
echo "step3"
touch ../$sample".counts.all_check"
touch ../$sample".counts.anno_check"
for i in {1..1000}               ##### loop
do
	tmp1=${sample}.shuffle.${i}
	awk -v var=$tmp1 '{sum+=$1;} END{print var,"\t",sum;}' $sample"_summit.shuffle.TE."$i".counts.all" >>../$sample".counts.all_check"
	awk -v var=$tmp1 '{sum+=$5;} END{print var,"\t",sum;}' $sample"_summit.shuffle.TE."$i".counts.anno" >>../$sample".counts.anno_check"
done
cd ../

tar czvf $sample".shuffle.tar.gz" $sample".shuffle/"
rm -rf $sample".shuffle/"
rm $sample"_summit."*"Only.bed"
