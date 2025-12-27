#!/bin/bash
# Homework 2 – Sequence Alignment & Variant Calling
# Ornit Bhonkar

####################
# Prerequisites
####################
qsub -I -l select=1:ncpus=2:mem=4gb -l walltime=01:00:00 /bin/bash
cp /home/stavnaky/HW1_Input_Files_2025_2026/* $HOME

####################
# Q1 – Alignment
####################
REF=/home/stavnaky/ref_grch38/GRCh38.d1.vd1.fa

bwa mem $REF PDX_P0.chr17.R1.fastq.gz PDX_P0.chr17.R2.fastq.gz > PDX_P0.sam
bwa mem $REF PDX_P2.chr17.R1.fastq.gz PDX_P2.chr17.R2.fastq.gz > PDX_P2.sam

####################
# Q2 – SAM to BAM
####################
samtools view -b -o PDX_P0.bam PDX_P0.sam
samtools view -b -o PDX_P2.bam PDX_P2.sam

ls -lh PDX_P0.sam PDX_P0.bam
ls -lh PDX_P2.sam PDX_P2.bam

####################
# Q3 – Sort BAM files
####################
samtools sort -o PDX_P0.sorted.bam PDX_P0.bam
samtools sort -o PDX_P2.sorted.bam PDX_P2.bam

####################
# Q4 – Index BAM files
####################
samtools index PDX_P0.sorted.bam
samtools index PDX_P2.sorted.bam

####################
# Q5 – mpileup
####################
bcftools mpileup -f "$REF" PDX_P0.sorted.bam PDX_P2.sorted.bam > PDX.P0_P2.pileup.vcf

####################
# Q6 – Variant calling (SNPs only)
####################
bcftools call -mv -Ov -o PDX.variants.vcf PDX.P0_P2.pileup.vcf
bcftools view -v snps -Ov -o PDX.snps.vcf PDX.variants.vcf

####################
# Q7 – Germline variants
####################
bcftools view -i 'GT[0]!="0/0" && GT[0]!="./."' -Ov -o PDX.germline.vcf PDX.snps.vcf

grep -vc '^#' PDX.snps.vcf
grep -vc '^#' PDX.germline.vcf

####################
# Q8 – Quality filtering + statistics
####################
bcftools view -i 'INFO/BQBZ==0 && INFO/DP>8 && INFO/MQ>=60 && INFO/AN==4' \
  -Ov -o PDX.Q8.filtered.vcf PDX.germline.vcf

bcftools stats PDX.germline.vcf > stats.Q8.before.txt
bcftools stats PDX.Q8.filtered.vcf > stats.Q8.after.txt

grep '^SN' stats.Q8.before.txt | head -n 30
grep '^SN' stats.Q8.after.txt  | head -n 30

bcftools view -h PDX.germline.vcf | grep 'BQBZ'

bcftools query -f '%BQBZ\n' PDX.Q8.filtered.vcf | sort | uniq -c

grep -A10 '^DP' stats.Q8.before.txt | head
grep -A10 '^DP' stats.Q8.after.txt  | head

bcftools query -f '%MQ\n' PDX.Q8.filtered.vcf | sort | uniq -c
bcftools query -f '%MQ\n' PDX.germline.vcf | sort | uniq -c

bcftools query -f '%AN\n' PDX.germline.vcf | sort | uniq -c
bcftools query -f '%AN\n' PDX.Q8.filtered.vcf | sort | uniq -c

####################
# Q9 – Loss of Heterozygosity (LOH)
####################
bcftools view -i '(GT[0]="0/1" || GT[0]="1/0") && (GT[1]="0/0" || GT[1]="1/1")' \
  -Ov -o PDX.Q9.LOH.vcf PDX.Q8.filtered.vcf

bcftools view -i '(GT[0]="0/1" || GT[0]="1/0") && GT[1]="0/0"' \
  -Ov -o PDX.Q9.LOH_ref.vcf PDX.Q8.filtered.vcf

bcftools view -i '(GT[0]="0/1" || GT[0]="1/0") && GT[1]="1/1"' \
  -Ov -o PDX.Q9.LOH_alt.vcf PDX.Q8.filtered.vcf

echo "LOH total:" $(grep -vc '^#' PDX.Q9.LOH.vcf)
echo "LOH to ref:" $(grep -vc '^#' PDX.Q9.LOH_ref.vcf)
echo "LOH to alt:" $(grep -vc '^#' PDX.Q9.LOH_alt.vcf)

####################
# Q10 – Specific mutation
####################
bcftools view -H -i 'CHROM="chr17" && POS=78392462' PDX.variants.vcf
bcftools view -H -i 'CHROM="chr17" && POS=78392462' PDX.Q9.LOH.vcf

####################
# Q11 – Substitution statistics
####################
bcftools stats PDX.snps.vcf > stats.Q11.snps.txt
grep '^ST' stats.Q11.snps.txt

bcftools stats PDX.Q9.LOH.vcf > stats.Q11.loh.txt
grep '^ST' stats.Q11.loh.txt
