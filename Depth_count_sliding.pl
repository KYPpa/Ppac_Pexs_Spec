#Update 190306
#SNP_count_sliding.pl-> Depth_count_sliding.pl

#! /usr/bin/perl

use strict;
use warnings;

my $snptfilename=$ARGV[0];#Depth file which was produced in the samtools depth -a
my $chrlenfilename=$ARGV[1];#Chromosome length file which was produced by fasta_count.pl

open(my $lenh,$chrlenfilename) or die;

my@length_list=<$lenh> or die;
chomp @length_list;
my$reshash={};
reshash_initialing(\@length_list,$reshash,1);

open(my $snpth,"$snptfilename") or die;
	
my $line;
while($line=<$snpth>){
    chomp $line;
    if($line=~/^#/){
        next;
    }
    my@info=split("\t",$line);
    my$pos=int(($info[1]-1)/1000);
    $reshash->{$info[0]}[$pos][0]+=$info[2];
}

my@chr=keys %$reshash;
@chr=sort(@chr);
print "Chromosome\tPosition\tCount\n";
	
for (my $ded=0; $ded<@chr; $ded++){
    for (my $dad=0; $dad<@{$reshash->{$chr[$ded]}}; $dad++){
        local $"="\t";
        my$position=$dad;#To reduce the file complexity I used this 0-based number of kb coordinate
        print "$chr[$ded]\t";
        print "$position\t";
        print "@{$reshash->{$chr[$ded]}[$dad]}\n";
    }
}


sub reshash_initialing{
	#Newlines should be chomped in $len_list before input.
	my($len_list,$res,$colnum)=@_;
	for(my$ded=0; $ded<@{$len_list}; $ded++){
		$len_list->[$ded]=~/([^\t]+)\t([^\t]+)/;
		my$chr=$1;
		my$len=$2;
		my$end=int(($len-1)/1000);
		for(my$dad=0; $dad<=$end; $dad++){
			for(my$dud=0; $dud<$colnum; $dud++){
				$res->{$chr}[$dad][$dud]=0;
			}
		}
	}
}
