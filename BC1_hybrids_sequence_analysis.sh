
## BC1_hybrids_sequence_analysis.sh
## The Script originates from q21012104.sh.
## Needs path to bwa, samtools, bcftools
##
## Requires
## 1. A pair of fastq file of pair-end illumina short reads of WGS
## (e.g. SH8RSG001_R1.fastq and SH8RSG001_R2.fastq)
## 2. P. pacificus reference genome (El_Paco_genome.fa)
## 3. A chromosome length file of the reference (El_Paco_genome_chr_len.txt)
## 4. The region file for SNP calling that points the position of
##    P. exspectatus unique SNPs
##    (Pexspure0_snps_strsel_regionfile.txt)

## "...RSG" is an identifier of the QTL analysis
## 8RSG - BC1 hybrid male
## 10RSG - BC1 hybrid female
## 11RSG - BC1 hybrid hermaphrodite
## Here I show the example of 8RSG.

# Define ID number i
i=1
fnum=$i; printf -v fnum "%03d" $fnum

# P. pacificus Reference sequence data
REF="Dataset/El_Paco_genome.fa"

# Mapping
bwa mem $REF \
Dataset/seqreads/SH8RSG${fnum}_R1.fastq Dataset/seqreads/SH8RSG${fnum}_R2.fastq \
> Dataset/SH8RSG${fnum}.sam

samtools view -bS Dataset/SH8RSG${fnum}.sam > Dataset/SH8RSG${fnum}.bam

samtools sort -T Dataset/SH8RSG${fnum}_sorted -o Dataset/SH8RSG${fnum}_sorted.bam \
Dataset/SH8RSG${fnum}.bam

rm Dataset/SH8RSG${fnum}.sam
rm Dataset/SH8RSG${fnum}.bam

samtools index Dataset/SH8RSG${fnum}_sorted.bam

samtools idxstats Dataset/SH8RSG${fnum}_sorted.bam > Dataset/SH8RSG${fnum}_idxstat.txt

# SNP calling
samtools mpileup -f $REF -g Dataset/SH8RSG${fnum}_sorted.bam > Dataset/SH8RSG${fnum}.bcf

bcftools index Dataset/SH8RSG${fnum}.bcf

# Call SNPs at the sites where P. exspectatus pure species has specific SNPs

bcftools call -R Dataset/Pexspure0_snps_strsel_regionfile.txt \
-m Dataset/SH8RSG${fnum}.bcf > Dataset/SH8RSG${fnum}_RS.vcf

samtools depth -a \
Dataset/SH8RSG${fnum}_sorted.bam \
> Dataset/SH8RSG${fnum}_depth.txt

# The perl script makes 1kb non-overlapping sliding window data
# of the count of read depth. Each base pair was counted.

perl Depth_count_sliding.pl \
Dataset/SH8RSG${fnum}_depth.txt \
Dataset/El_Paco_genome_chr_len.txt \
> Dataset/SW_Depth_PpacEP_SH8RSG${fnum}.txt

# The perl script sums $ARGV[0] kb length of sliding window that
# shifting by $ARGV[1] kb. $ARGV[2] decides
# if the end of the contig is truncated (0) or summed as another window (1).

perl 1kb_sliding_sum3_shell.pl \
Dataset/SW_Depth_PpacEP_SH8RSG${fnum}.txt 100 100 0

rm Dataset/SH8RSG${fnum}_sorted.bam
rm Dataset/SH8RSG${fnum}_sorted.bam.bai
rm Dataset/SH8RSG${fnum}.bcf
rm Dataset/SH8RSG${fnum}_depth.txt
