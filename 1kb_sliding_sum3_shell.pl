#! /usr/bin/perl
#2018.6.25 1kb_sliding_sum3_dir.pl
#Change for header information
#and new coordinate system 

use strict;
use warnings;

my $slidename=$ARGV[0];
my $win=$ARGV[1]; #eg 100
my $shift=$ARGV[2]; #eg 100
my $remain=$ARGV[3]; #eg 0

open(my $snpth,"$slidename") or die;
my$newfile=add_to_filename($slidename,"${win}w-${shift}s");
open(my($outputh),"> $newfile") or die;

my $line;
my $slide_data=[];

#The first chromosome
$line=<$snpth>; #The first line is recognized as the header.
chomp $line;
my @info=split("\t",$line);
splice(@info,1,1,("Start","End"));
$line=join("\t",@info);
print $outputh "$line\n";

$line=<$snpth>;
chomp $line;
@info=split("\t",$line);
my $chromosome=$info[0];
my $log_num=1;
my $catnumber=@info-2;

my $dad=0;
while($dad<$catnumber){
	$slide_data->[0][$dad]=$info[2+$dad];
	$dad++;
}
my $ded=1;
my $lastposition=$info[1];

#Loop
while($line=<$snpth>){
	chomp $line;
	@info=split("\t",$line);
	if($info[0] ne $chromosome){
		$lastposition=sliding_window_summing2($slide_data,$lastposition,$win,$shift,$catnumber,$remain);
		sliding_window_writing_for_hetero($outputh,$chromosome,$slide_data,$lastposition,$win,$shift,$catnumber);
		$chromosome=$info[0];
		$log_num=1;
		$ded=0;
		$slide_data=[];
	}
	$dad=0;
	while($dad<$catnumber){
		$slide_data->[$ded][$dad]=$info[2+$dad];
		$dad++;
	}
	$ded++;
	$lastposition=$info[1];
}
$lastposition=sliding_window_summing2($slide_data,$lastposition,$win,$shift,$catnumber,$remain);
sliding_window_writing_for_hetero($outputh,$chromosome,$slide_data,$lastposition,$win,$shift,$catnumber);

sub sliding_window_summing2{
	my($slide,$last,$wind,$shif,$cat,$rem)=@_;
	my $newslide=[];
	my $framenumber=$last+1; #framenumber = last position of new coordinate system
	if($rem==0 and $framenumber<$win){
		#frameが一つ目のwindowも満たさない場合
		$last=0;
	}else{
		#frameから最初のwinを除いたとき、どれくらい動けるかによって数を決める。
		#$remain=1のときは切り上げ、0のときは切り下げ
		#$remain=1のとき、int((動けるフレーム-1)/動く幅)+1だけ動く。
		#$remain=0のとき、int((動けるフレーム)/動く幅))
		$last=int(($framenumber-$win-$rem)/$shif)+$rem+1; #last: 1 basenumber
		if($rem==1 and $framenumber<=$win){$last=1;}
		my $ded=0;
		while($ded<$last){
			my $pos=$shif*$ded;
			#remain=1のときにslideのデータが切れたら即終了
			unless(defined $slide->[$pos][0] and $slide->[$pos][0] ne ""){
				last;
			}
			my $dad=0;
			while($dad<$wind){
				my$dud=0;
				unless(defined $slide->[$pos+$dad][0] and $slide->[$pos+$dad][0] ne ""){
					last;
				}
				while($dud<$cat){
					if($dad==0){
						$newslide->[$ded][$dud]=$slide->[$pos][$dud];
					}else{
						$newslide->[$ded][$dud]+=$slide->[$pos+$dad][$dud];
					}
					$dud++;
				}
				$dad++;
			}
			$ded++;
		}
	}
	@$slide=@$newslide;
	return $last;
}


#done
sub sliding_window_writing_for_hetero{
	my($outh,$chromo,$slide,$last,$wind,$shif,$cat)=@_;
	my($ded)=0;
	while($ded<$last){
		my$start=$shif*$ded*1000+1;
		my$end=$shif*$ded*1000+$win*1000;
		print $outh "$chromo\t$start\t$end";
		my $dad=0;
		while($dad<$cat){
			print $outh "\t$slide->[$ded][$dad]";
			$dad++;
		}
		print $outh "\n";
		$ded++;
	}
}

sub add_to_filename{
	my($file,$added)=@_;
	$file=~/^(.+)(\.[^\.]+)$/ or die;
	$file=$1."_".$added.$2;
	return $file;
}

