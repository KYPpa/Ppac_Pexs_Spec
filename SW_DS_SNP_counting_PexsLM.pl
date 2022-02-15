#Wrote 210105
#SNP_count_sliding_shell_general.pl -> SW_DS_SNP_counting_PexsLM.pl


use strict;
use warnings;

my $file=$ARGV[0]; #whole path, vcf file name
my $chrlenfilename=$ARGV[1];
my $dbname=$ARGV[2];
open my $snpth, "$file" or die;
open my $lenh, $chrlenfilename or die;
my%snpdb=();
dbmopen(%snpdb,"$dbname",0644) or die "$!";

my@length_list=<$lenh> or die;
chomp @length_list;

my$reshash={};
reshash_initialing(\@length_list,$reshash,3);

my $line;
my $dcount=0;
my $scount=0;
my $rcount=0;
my $ecount=0;
while($line=<$snpth>){
	chomp $line;
	if($line=~/^#/){
		next;
	}
	$line=~/INDEL/ and next;
	my@info=split("\t", $line);
	if(defined $snpdb{"$info[0]:$info[1]"}){
		my $altnuc="";
		my $ds="";
		my $pos=int(($info[1]-1)/1000);
		if($snpdb{"$info[0]:$info[1]"}=~/^([DS])([ATGC])$/){
			$ds=$1;
			$altnuc=$2;
			if($altnuc eq $info[4]){
				if($ds eq "D"){
					$reshash->{$info[0]}[$pos][1]++;
					$dcount++;
				}else{
					$reshash->{$info[0]}[$pos][2]++;
					$scount++;
				}
			}else{
				$reshash->{$info[0]}[$pos][0]++;
				$rcount++;
			}
		}elsif($snpdb{"$info[0]:$info[1]"}=~/^D([ATGC])S([ATGC])$/){
			if($info[4] eq $1){
				$reshash->{$info[0]}[$pos][1]++;
			}elsif($info[4] eq $2){
				$reshash->{$info[0]}[$pos][2]++;
			}else{
				$reshash->{$info[0]}[$pos][0]++;
			}
		}
	}else{
		$ecount++;
	}
}

dbmclose(%snpdb);

my@chr=keys %$reshash;
@chr=sort(@chr);

print "Chromosome\tPosition\tRef\tD\tS\n";

for (my $ded=0; $ded<@chr; $ded++){
	for (my $dad=0; $dad<@{$reshash->{$chr[$ded]}}; $dad++){
		local $"="\t";
		my$position=$dad;#To reduce the file complexity I used this 0-start number of kb coordinate
		print "$chr[$ded]\t";
		print "$position\t";
		print "@{$reshash->{$chr[$ded]}[$dad]}\n";
	}
}

print STDERR "Finished properly\n$file\n$chrlenfilename\n$dbname\nR, $rcount; D, $dcount; S, $scount\nError, $ecount\n";

sub reshash_initialing{
	#Newlines should be chomped in $len_list before input.
	my($len_list,$res,$colnum)=@_;
	for(my$ded=1; $ded<@{$len_list}; $ded++){
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
