use strict;
use warnings;

my($file,$interval)=($ARGV[0],$ARGV[1]);
open(my$ih,$file) or die;

my$newfile=$file."_count_dupdel.txt";
open(my$oh,"> $newfile") or die;

my $line;
$line=<$ih>;
print "xchr\txstart\tychr\tystart\n";
print $oh "Contigs\tChrI\tChrII\tChrIII\tChrIV\tChrV\tChrX\n";

my$preline;
my$precon="XXX";
my$data={};
my$matrix=[];
while($line=<$ih>){
	chomp $line;
	my@info=split("\t",$line);
	if($precon ne "XXX" and $precon ne $info[3]){
		dotting_process_dupdel($data,$matrix,$interval,$oh);
        $data={};
		$matrix=[];
    }
    $precon=$info[3];
	my$ind=chr_to_ind($info[0]);
	if($ind>5){next;}
	for my $nuc ($info[5]..$info[6]){
		if(defined $data->{$nuc}){
			$data->{$nuc}++;
		}else{
			$data->{$nuc}=1;
		}
	}
	push(@$matrix,(\@info));
}
dotting_process_dupdel($data,$matrix,$interval,$oh);

sub dotting_process_dupdel{
	my($data,$matrix,$interval,$oh)=@_;
	my$fdig=$matrix->[0][3];
	my@count=(0,0,0,0,0,0);
	for my $ded (0..(scalar(@$matrix)-1)){
		my $xs=$matrix->[$ded][1];
		my $xe=$matrix->[$ded][2];
		my $ys=$matrix->[$ded][5];
		my $ye=$matrix->[$ded][6];
		my$ind=chr_to_ind($matrix->[$ded][0]);
		if($ind>5){next;}
		
		for my $yco ($ys..$ye){
			my$dupno=$data->{$yco};
			if($dupno==1){
				$count[$ind]++;
                if(bernoulli_choose(1/$interval)){
                    my$xco=0;
                    if($matrix->[$ded][4] eq "+"){
                        $xco=($xe-$xs)*($yco-$ys)/($ye-$ys)+$xs;
                    }else{
                        $xco=($xs-$xe)*($yco-$ys)/($ye-$ys)+$xe;
                    }
                    print "$ind\t$xco\t$fdig\t$yco\n";
                }
			}
		}
	}
	print $oh "$fdig\t$count[0]\t$count[1]\t$count[2]\t$count[3]\t$count[4]\t$count[5]\n";
}

sub bernoulli_choose{
	my($prob)=@_;
	if(rand(1/$prob)<=1){
		return 1;
	}else{
		return 0;
	}
}

sub chr_to_ind{
	my($chr)=@_;
	my$ind=6;
	if($chr eq "ChrI"){
		$ind=0;
	}elsif($chr eq "ChrII"){
		$ind=1;
	}elsif($chr eq "ChrIII"){
		$ind=2;
	}elsif($chr eq "ChrIV"){
		$ind=3;
	}elsif($chr eq "ChrV"){
		$ind=4;
	}elsif($chr eq "ChrX"){
		$ind=5;
    }
	return($ind);
}
