
## synteny_analysis_data_preparation.sh
## The Script originates from q21050503.sh and q21051101.sh
## Needs path to latz
##
## Requires
## 1. A pair of fasta file of the assembly
## (El_Paco_genome.fa & RS5522B_assembly_ver5_210505.fa)
## 2. A chromosome length file of one reference (RS5522B_assembly_ver5_210505_Chr_len.txt)

## Synteny analysis with lastz
## The general format is used for the susequent analyses.

lastz \
Dataset/El_Paco_genome.fa[multiple] \
Dataset/RS5522B_assembly_ver5_210505.fa \
--notransition --step=20 --nogapped \
--format=general:name1,start1,end1,name2,strand2,start2+,end2+,score \
> Dataset/lastz_general_PpacEP_vs_PexsCon_ver5.txt

## For Dot plot
## perl script generates single nucleotide plot between X and Y
## but the number is downsized.
## $ARGV: 0, lastz output file; 1, denominator of downsize
## Here it reduces by 1/100.
## Mutiple-hit regions are removed.

perl dotplot_GA_dupdel3.pl \
Dataset/lastz_general_PpacEP_vs_PexsCon_ver5.txt 100 \
> Dataset/dotplot_dupdel_PpacEP_vs_PexsCon_ver5.txt

## For Circos plot
## perl script generates count data of homologous nucleotide
## between 1kb window of X and Y
## Until1 means this allows only single hit regions. (single in X)

perl dotplot_GA_seq_A_until1.pl \
Dataset/lastz_general_PpacEP_vs_PexsCon_ver5.txt \
Dataset/RS5522B_assembly_ver5_210505_Chr_len.txt \
> Dataset/dot_plot_1kbtable_PpacEP_vs_PexsCon_ver5.txt

## perl script sums the counts of $ARGV[1] kb window combination between X and Y
## $ARGV[2] base pair is the minimum number. Otherwise just deleted.

perl dotplot_GA_seq_B.pl \
Dataset/dot_plot_1kbtable_PpacEP_vs_PexsCon_ver5.txt 100 10000 \
> Dataset/dot_plot_100kbtable_PpacEP_vs_PexsCon_ver5.txt
