#! /usr/bin/perl

use strict;
use warnings;

my($dirname)=$ARGV[0];
my $slidefilenames=[];
directorysearch_name($slidefilenames,$dirname,[".txt"]);

my $slidename=$slidefilenames->[0];
open(my $snpth,"$dirname/$slidename") or die;
$slidename=~s/\.([a-z]+)$//;
print STDERR "$slidename\n";


my$line=<$snpth>; #The first line is recognized as the header.
chomp $line;
my @info=split("\t",$line);
splice(@info,0,3,("File","Chr","Start","End"));
$line=join("\t",@info);
print "$line\n";

while($line=<$snpth>){
	chomp $line;
	print "$slidename\t$line\n";
}

my $dud=1;
while($dud<@$slidefilenames){
	my $slidename=$slidefilenames->[$dud];
	open(my $snpth,"$dirname/$slidename") or die;
	$slidename=~s/\.([a-z]+)$//;
	print STDERR "$slidename\n";

	my $line;

	#The first chromosome
	$line=<$snpth>; #The first line is recognized as the header.

	while($line=<$snpth>){
		chomp $line;
		print "$slidename\t$line\n";
	}
	$dud++;
}

sub directorysearch_name{
	my($file_name_list,$directory_name,$fileextend)=@_;
	my($directory_handle);
	opendir($directory_handle,$directory_name);
	my(@filelist)=readdir($directory_handle);
	my($number_a)=0;
	my($number_c)=0;
	while($number_a<@filelist){
		my($number_b)=0;
		while($number_b<@$fileextend){
			if($filelist[$number_a]=~/$fileextend->[$number_b]/){
				print STDERR "Found $filelist[$number_a].\n";
				$file_name_list->[$number_c]=$filelist[$number_a];
				$number_c++;
			}
			$number_b++;
		}
	$number_a++;
	}
}
