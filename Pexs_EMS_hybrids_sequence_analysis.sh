
## Pexs_EMS_hybrids_sequence_analysis.sh
## The Script originates from q21051002.sh.
## Needs path to bwa, samtools, bcftools,
##
## Requires
## 1. A pair of fastq file of pair-end illumina short reads of WGS
## (e.g. PexsLM4G001_R1.fastq and PexsLM4G001_R2.fastq)
## 2. P. exspectatus reference genome (RS5522B_assembly_ver5_210505.fa)
## 3. A chromosome length file of the reference (RS5522B_assembly_ver5_210505_Chr_len.txt)
## 4. The region file for SNP calling that points the position of
##    unique SNPs of the grand-dam (D) or grand-sire (S)
##    (e.g. PexsLM_A3_region_file.txt)
## 5. Hush database file binding SNP position => "D" or "S" (e.g. snp_db_PexsLM2G_A3)

# Define ID number i
i=1
fnum=$i; printf -v fnum "%03d" $fnum

# P. exspectatus Reference sequence data
REF="Dataset/RS5522B_assembly_ver5_210505.fa"

# 4th Generation (analyzed hybrids)
gen="4G"

# Corresponding group to 103 individuals of 4th G
marray=("A" "A" "B" "B" "B" "C" "F" "F" "A" "B" "B" "B" "B" "B" "F" "F" "G" "A" "B" "B" "B" "B" "B" "B" "B" "B" "B" "C" "C" "C" "C" "C" "F" "F" "F" "F" "G" "A" "A" "A" "B" "B" "B" "B" "B" "B" "B" "B" "C" "C" "C" "C" "C" "C" "C" "C" "D" "F" "F" "F" "G" "A" "B" "B" "B" "B" "B" "B" "C" "C" "C" "F" "F" "F" "F" "G" "J" "J" "J" "B" "B" "F" "F" "F" "F" "F" "G" "J" "J" "J" "J" "A" "A" "B" "C" "C" "C" "C" "F" "F" "F" "F" "F")

((m=i-1))
mark="${marray[$m]}"

# Sample ID is "PexsLM${gen}${fnum}"
# Start from a pair of fastq files of the sample generated in illumina short read sequencer

# Mapping

bwa mem $REF \
Dataset/PexsLM${gen}${fnum}_R1.fastq \
Dataset/PexsLM${gen}${fnum}_R2.fastq > \
Dataset/PexsLM${gen}${fnum}.sam

samtools view -bS \
Dataset/PexsLM${gen}${fnum}.sam > \
Dataset/PexsLM${gen}${fnum}.bam

samtools sort -T \
Dataset/PexsLM${gen}${fnum}_sorted \
-o Dataset/PexsLM${gen}${fnum}_sorted.bam \
Dataset/PexsLM${gen}${fnum}.bam

samtools index \
Dataset/PexsLM${gen}${fnum}_sorted.bam

# SNP call
samtools mpileup -f $REF -g \
Dataset/PexsLM${gen}${fnum}_sorted.bam \
> Dataset/PexsLM${gen}${fnum}.bcf

bcftools index \
Dataset/PexsLM${gen}${fnum}.bcf

# The region where their grand-dam(D) or grand-sire(S) have unique SNPs
# were investigated.

bcftools call \
-R Dataset/PexsLM_${mark}3_region_file.txt \
-m Dataset/PexsLM${gen}${fnum}.bcf \
> Dataset/PexsLM${gen}${fnum}_region${mark}3.vcf

# The perl script generates count data of 1kb non-overlapping sliding window
# The coordinate is 0-base kb number.(e.g. 1-1000 => 0)

perl SW_DS_SNP_counting_PexsLM.pl \
Dataset/PexsLM${gen}${fnum}_region${mark}3.vcf \
Dataset/RS5522B_assembly_ver5_210505_Chr_len.txt \
Dataset/snp_db_PexsLM2G_${mark}3 \
> Dataset/SW_DS_SNP_count_PexsLM${gen}${fnum}_PexsCon_ver5.txt

# The perl script sums $ARGV[0] length of sliding window that
# shifting by $ARGV[1]. $ARGV[2] decides if the end of the contig is truncated (0) or summed as another window (1)

perl 1kb_sliding_sum3_shell.pl \
Dataset/SW_DS_SNP_count_PexsLM${gen}${fnum}_PexsCon_ver5.txt 500 500 0

rm PexsLM${gen}${fnum}.sam
rm PexsLM${gen}${fnum}.bam
rm PexsLM${gen}${fnum}_sorted.bam
rm PexsLM${gen}${fnum}_sorted.bam.bai
rm PexsLM${gen}${fnum}.bcf
rm PexsLM${gen}${fnum}.bcf.csi
