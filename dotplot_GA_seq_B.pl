use strict;
use warnings;

my($dotplotfile,$length)=($ARGV[0],$ARGV[1]);
my($minimum)=0;
if(defined $ARGV[2]){
	$minimum=$ARGV[2];
}
open(my$ih,$dotplotfile) or die;

my $line;
$line=<$ih>;
print "xchr\txpos\tychr\typos\tcount\n";
$line=<$ih>;
chomp $line;
my @info=split("\t",$line);
my $ychr=$info[2];
my $res={};
summing_count($res,\@info,$length);

while($line=<$ih>){
	chomp $line;
	@info=split("\t",$line);
	if($ychr ne $info[2]){
		print_sum_count($res,$ychr,$minimum);
		$ychr=$info[2];
		$res={};
	}
	summing_count($res,\@info,$length);
}
print_sum_count($res,$ychr,$minimum);

sub summing_count{
	my($res,$info,$length)=@_;
	my$xpos=int($info->[1]/$length);
	my$ypos=int($info->[3]/$length);
	if(defined $res->{$info->[0]}{$xpos}{$ypos}){
		$res->{$info->[0]}{$xpos}{$ypos}+=$info[4];
	}else{
		$res->{$info->[0]}{$xpos}{$ypos}=$info[4];
	}
}

sub print_sum_count{
	my($res,$ychr,$minimum)=@_;
	my @xchrs=keys(%{$res});
	@xchrs = sort { lc($a) cmp lc($b) } @xchrs;
	for my $xchr (@xchrs){
		my @xcos=keys(%{$res->{$xchr}});
		@xcos = sort { $a <=> $b } @xcos;
		for my $xco (@xcos){
			my @ycos=keys(%{$res->{$xchr}{$xco}});
			@ycos = sort { $a <=> $b } @ycos;
			for my $yco (@ycos){
				if($res->{$xchr}{$xco}{$yco}>=$minimum){
					print "$xchr\t$xco\t$ychr\t$yco\t$res->{$xchr}{$xco}{$yco}\n";
				}
			}
		}
	}
}
