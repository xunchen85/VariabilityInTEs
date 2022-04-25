#!/bin/bash
module load bedtools/2.29.2

###     Author:     Xun Chen
###     Date:       2020/7/10
###     Contact:    xunchen85@gmail.com

sample=$1
inputFolder=$2
ScriptFolder=$3
annotationFolder=$4
Type=$5
TEbedFileName=$6

#### 1. extract median as the summit
python ${ScriptFolder}Convert_to_PeakSummit.py -i $inputFolder$sample".bed" >$sample"_summit.bed"

#### 2. Get peaks within TEs
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.bed" >$sample"_summit.rmsk.bed"
bedtools intersect -wa -v -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.bed" >$sample"_summit.rmskNo.bed"
cat $sample"_summit.rmsk.bed" | sed 's/:/ /g' |awk '{print$4,$6,$7}' |uniq -c | awk '{print$2,$3,$4,$1}' >$sample"_summit.rmsk.counts.all"

echo "step1"
#### 3-1. Obtain percentage within each annotation file
# non-desert
bedtools intersect -u -a $sample"_summit.bed" -b $annotationFolder"5UTR.bed" >$sample"_summit.5UTR.bed"
bedtools intersect -u -a $sample"_summit.bed" -b $annotationFolder"3UTR.bed" >$sample"_summit.3UTR.bed"
bedtools intersect -u -a $sample"_summit.bed" -b $annotationFolder"distal.bed" >$sample"_summit.distal.bed"
bedtools intersect -u -a $sample"_summit.bed" -b $annotationFolder"exons.bed" >$sample"_summit.exons.bed"
bedtools intersect -u -a $sample"_summit.bed" -b $annotationFolder"introns.bed" >$sample"_summit.introns.bed"
bedtools intersect -u -a $sample"_summit.bed" -b $annotationFolder"non.desert.bed" >$sample"_summit.non.desert.bed"
bedtools intersect -u -a $sample"_summit.bed" -b $annotationFolder"promoters.bed" >$sample"_summit.promoters.bed"
bedtools intersect -u -a $sample"_summit.bed" -b $annotationFolder"proximal.bed" >$sample"_summit.proximal.bed"
bedtools intersect -u -a $sample"_summit.bed" -b $annotationFolder"tss.bed" >$sample"_summit.tss.bed"

# desert
bedtools intersect -v -a $sample"_summit.bed" -b $annotationFolder"non.desert.bed" >$sample"_summit.desert.bed"

#### 3-2. Keep each annotation only peaks
echo "step2"
# desert
cat $sample"_summit.5UTR.bed" $sample"_summit.3UTR.bed" $sample"_summit.distal.bed" $sample"_summit.exons.bed" $sample"_summit.introns.bed" $sample"_summit.promoters.bed" $sample"_summit.proximal.bed" $sample"_summit.tss.bed" | sort | uniq >$sample"_tmp.bed"
bedtools subtract -a $sample"_summit.desert.bed" -b $sample"_tmp.bed" >$sample"_summit.desertOnly.bed"
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.desertOnly.bed" >$sample"_summit.desertOnly.TE.bed"
bedtools intersect -wa -wb -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.desertOnly.bed" >$sample"_summit.desertOnly.TE2.bed"

# distal
cat $sample"_summit.5UTR.bed" $sample"_summit.3UTR.bed" $sample"_summit.exons.bed" $sample"_summit.introns.bed" $sample"_summit.promoters.bed" $sample"_summit.proximal.bed" $sample"_summit.tss.bed" | sort | uniq >$sample"_tmp.bed"
bedtools subtract -a $sample"_summit.distal.bed" -b $sample"_tmp.bed" >$sample"_summit.distalOnly.bed"
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.distalOnly.bed" >$sample"_summit.distalOnly.TE.bed"
bedtools intersect -wa -wb -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.distalOnly.bed" >$sample"_summit.distalOnly.TE2.bed"

# proximal
cat $sample"_summit.5UTR.bed" $sample"_summit.3UTR.bed" $sample"_summit.exons.bed" $sample"_summit.introns.bed" $sample"_summit.promoters.bed" $sample"_summit.tss.bed" | sort | uniq >$sample"_tmp.bed"
bedtools subtract -a $sample"_summit.proximal.bed" -b $sample"_tmp.bed" >$sample"_summit.proximalOnly.bed"
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.proximalOnly.bed" >$sample"_summit.proximalOnly.TE.bed"
bedtools intersect -wa -wb -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.proximalOnly.bed" >$sample"_summit.proximalOnly.TE2.bed"

# introns
cat $sample"_summit.5UTR.bed" $sample"_summit.3UTR.bed" $sample"_summit.exons.bed" $sample"_summit.promoters.bed" $sample"_summit.tss.bed" | sort | uniq >$sample"_tmp.bed"
bedtools subtract -a $sample"_summit.introns.bed" -b $sample"_tmp.bed" >$sample"_summit.intronsOnly.bed"
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.intronsOnly.bed" >$sample"_summit.intronsOnly.TE.bed"
bedtools intersect -wa -wb -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.intronsOnly.bed" >$sample"_summit.intronsOnly.TE2.bed"

# promoters
cat $sample"_summit.5UTR.bed" $sample"_summit.3UTR.bed" $sample"_summit.exons.bed" $sample"_summit.tss.bed" | sort | uniq >$sample"_tmp.bed"
bedtools subtract -a $sample"_summit.promoters.bed" -b $sample"_tmp.bed" >$sample"_summit.promotersOnly.bed"
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.promotersOnly.bed" >$sample"_summit.promotersOnly.TE.bed"
bedtools intersect -wa -wb -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.promotersOnly.bed" >$sample"_summit.promotersOnly.TE2.bed"

### original order: tss -> exons -> 5UTR ->3UTR before 2021/1/12
# tss
rm $sample"_tmp.bed"
touch $sample"_tmp.bed"
bedtools subtract -a $sample"_summit.tss.bed" -b $sample"_tmp.bed" >$sample"_summit.tssOnly.bed"
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.tssOnly.bed" >$sample"_summit.tssOnly.TE.bed"
bedtools intersect -wa -wb -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.tssOnly.bed" >$sample"_summit.tssOnly.TE2.bed"

# 3UTR
cat $sample"_summit.tss.bed" $sample"_summit.5UTR.bed" | sort | uniq >$sample"_tmp.bed"
bedtools subtract -a $sample"_summit.3UTR.bed" -b $sample"_tmp.bed" >$sample"_summit.3UTROnly.bed"
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.3UTROnly.bed" >$sample"_summit.3UTROnly.TE.bed"
bedtools intersect -wa -wb -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.3UTROnly.bed" >$sample"_summit.3UTROnly.TE2.bed"

# 5UTR
cat $sample"_summit.tss.bed" | sort | uniq >$sample"_tmp.bed"
bedtools subtract -a $sample"_summit.5UTR.bed" -b $sample"_tmp.bed" >$sample"_summit.5UTROnly.bed"
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.5UTROnly.bed" >$sample"_summit.5UTROnly.TE.bed"
bedtools intersect -wa -wb -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.5UTROnly.bed" >$sample"_summit.5UTROnly.TE2.bed"

# exons
cat $sample"_summit.5UTR.bed" $sample"_summit.3UTR.bed" $sample"_summit.tss.bed" | sort | uniq >$sample"_tmp.bed"
bedtools subtract -a $sample"_summit.exons.bed" -b $sample"_tmp.bed" >$sample"_summit.exonsOnly.bed"
bedtools intersect -wa -u -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.exonsOnly.bed" >$sample"_summit.exonsOnly.TE.bed"
bedtools intersect -wa -wb -a ${annotationFolder}${TEbedFileName} -b $sample"_summit.exonsOnly.bed" >$sample"_summit.exonsOnly.TE2.bed"

rm $sample"_tmp.bed"
rm $sample"_summit.5UTR.bed" $sample"_summit.3UTR.bed" $sample"_summit.distal.bed" $sample"_summit.exons.bed" $sample"_summit.introns.bed" $sample"_summit.non.desert.bed" $sample"_summit.promoters.bed" $sample"_summit.proximal.bed" $sample"_summit.tss.bed" $sample"_summit.desert.bed"
# move files
mkdir TEpeaks
grep '[a-z0-9]' $sample"_summit."*"Only.TE.bed" | sed 's/:/ /g' |awk '{print$1,$5,$7,$8}'| sort |uniq -c | awk '{print$2,$3,$4,$5,$1}' >$sample"_summit.TE.counts.anno"
grep '[a-zA-Z]' $sample"_summit."*"Only.TE2.bed" |sed 's/:/ /' >$sample"_summit.TE2.bed"

mv $sample"_summit."*".TE2.bed" TEpeaks/
mv $sample"_summit."*".TE.bed" TEpeaks/
mv $sample"_summit.rmsk.bed" TEpeaks/
mv $sample"_summit.rmskNo.bed" TEpeaks/

