use strict;
use warnings;

my($latzfile,$yclenfile)=($ARGV[0],$ARGV[1]);
open(my$ih,$latzfile) or die;
open(my$sh,$yclenfile) or die;

my $line;
$line=<$ih>;
print "xchr\txpos\tychr\typos\tcount\n";
$line=<$ih>;
chomp $line;
my@info=split("\t",$line);
my$previousfileposition=tell($ih);
my$precon=$info[3];

# Produces string of 0 with the same length of Y Chr

my$yrefseq=yrefseq_preparation($precon,$sh);
my$headline=$line;

# Take the part of string corresponding to coordiates of latz hit 

my $ystring=substr($yrefseq,$info[5]-1,$info[6]-$info[5]+1);

# Count the number of the hit on the string part but higher number is not counted because only the single hit is important.

$ystring=~tr/0123/1233/;

# Put back the part of the string

substr($yrefseq,$info[5]-1,$info[6]-$info[5]+1,$ystring);

while($line=<$ih>){
	chomp $line;
	@info=split("\t",$line);
	if($precon ne $info[3]){
	
		# function to print count of the single hits of 1kb window combinations between genome X and Y using the string indicating the hit number at each site
		$headline=dotting_process_sw_seq_A($ih,$yrefseq,$previousfileposition,$headline,$precon);

		$previousfileposition=tell($ih);
		$precon=$info[3];
		$yrefseq=yrefseq_preparation($precon,$sh);
	}
	$ystring=substr($yrefseq,$info[5]-1,$info[6]-$info[5]+1);
	$ystring=~tr/0123/1233/;
	substr($yrefseq,$info[5]-1,$info[6]-$info[5]+1,$ystring);
}
dotting_process_sw_seq_A($ih,$yrefseq,$previousfileposition,$headline,$precon);

sub yrefseq_preparation{
	my($precon,$sh)=@_;
	seek($sh,0,0);
	my($line);
	while($line=<$sh>){
		chomp $line;
		$line=~/([^\t]+)\t([^\t]+)/;
		my $name=$1;
		my $length=$2;
		if($name eq $precon){
			my $string="";
			for my $c (1..$length){
				$string .= "0";
			}
			return($string);
			last;
		}
	}
}

sub dotting_process_sw_seq_A{
	my($ih,$seq,$fileposition,$line,$ychr)=@_;
	my $res={};
	process_to_resfile($res,$line,$seq);
	
	# Go to beginning of the previous contig. i.e. every line is read twice.
	seek($ih,$fileposition,0);
	my$headline="";
	while($line=<$ih>){
		chomp $line;
		my@info=split("\t",$line);
		if($ychr ne $info[3]){
			$headline=$line;
			last;
		}
		process_to_resfile($res,$line,$seq);
	}
	my @xchrs=keys(%{$res});
	@xchrs = sort { lc($a) cmp lc($b) } @xchrs;
	for my $xchr (@xchrs){
		my @xcos=keys(%{$res->{$xchr}});
		@xcos = sort { $a <=> $b } @xcos;
		for my $xco (@xcos){
			my @ycos=keys(%{$res->{$xchr}{$xco}});
			@ycos = sort { $a <=> $b } @ycos;
			for my $yco (@ycos){
				print "$xchr\t$xco\t$ychr\t$yco\t$res->{$xchr}{$xco}{$yco}\n";
			}
		}
	}
	return($headline);
}

sub process_to_resfile{
	# For each Y Chr, a hush, $res, is made.
	# This has the information of X Chr 1kb window more than one basepair hit with one 1kb window of Y chr
	# The hush is used because there are so many combination of coordinates between X and Y genome.
	my($res,$line,$seq)=@_;
	chomp $line;
	my @info=split("\t",$line);
	my $xchr=$info[0];
	my $xs=$info[1];
	my $xe=$info[2];
	my $ys=$info[5];
	my $ye=$info[6];
	for my $yco ($ys..$ye){
		my$dupno=substr($seq,$yco-1,1);
		# Here, the number of duplication is reflected.
		# <2 means only positions with the single hit are allowed.
		if($dupno<2){
			my$xco=0;
			if($info[4] eq "+"){
				$xco=($xe-$xs)*($yco-$ys)/($ye-$ys)+$xs;
			}else{
				$xco=($xs-$xe)*($yco-$ys)/($ye-$ys)+$xe;
			}
			$xco=int(($xco-1)/1000);
			$yco=int(($yco-1)/1000);
			if(defined $res->{$xchr}{$xco}{$yco}){
				$res->{$xchr}{$xco}{$yco}++;
			}else{
				$res->{$xchr}{$xco}{$yco}=1;
			}
		}
	}
}



