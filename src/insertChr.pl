#!/usr/bin/perl -w
use strict;

my ($file,$out,$chr) = @ARGV;

open IN,"<",$file or die "can not read the file $file\n";
open SE,">",$out or die "can not write the file $out\n";
while(<IN>){
  chomp($_);
  my @tmp = split(/\t/,$_);
  print SE $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$chr."\t".$tmp[5]."\t".$tmp[6]."\t".$tmp[7]."\n";
}
close IN;
close SE;