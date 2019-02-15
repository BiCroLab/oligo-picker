#!/usr/bin/perl -w
use strict;
use Log::Log4perl qw(:easy);
use IPC::System::Simple qw(capturex);
use File::Basename;

my ($vmfa, $vm_map, $vmFout, $arguments, $script, $log, $snictmp, $dir_out) = @ARGV;

require("$script/FileHandling.pm");
Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

my @arg  = split(/::/,$arguments);
my $kmer = $arg[2];

INFO "\tGet the oligonucleotides without multiple matches in the first run\n";
my %vmCoord       = FileHandling::process_vmatch_output($vm_map,$kmer,$log);
my ($decision,$x) = getSeq($vmfa,$vmFout,\%vmCoord,$log);
my %files         = %{$x};

if($decision eq "yes"){# there are files, which contains unique oligos
  INFO "\tCopy the input files for the second vmatch run to $dir_out\n";
  foreach my $el (keys %files){
    capturex("rsync",($el,$dir_out));
  }
}


sub getSeq{
  my $fa_file    = $_[0];
  my $fa_out     = $_[1];
  my %coord      = %{$_[2]};
  my $obj        = Bio::SeqIO->new(-file => $fa_file, -format => 'fasta');
  my @list       = ();
  my $flag       = "yes";
  my $path       = dirname($fa_out);
  my $base       = basename($fa_out);
  my ($no,$pid)  = ($1,$2) if($base =~ /.*_(\d+)_(\d+)\..*/);
  my ($i,$snr)   = (1,1);
  my ($size,$bp) = (250000,0);
  my %output     = ();

 INFO "working on input $fa_file and output $fa_out\n";     

  # extract only those sequences which are really unique
  while(my $seq = $obj->next_seq){
    my $nuc   = $seq->seq;
    my $id    = $seq->id;
    my $desc  = $seq->desc;
    
    if(exists($coord{$id})){
      my $str = pack("Z*Z*Z*",$id,$desc,$nuc);
      push(@list,$str);
    }
  }
  
  if(@list){
    foreach my $str (@list){
      my ($id,$des,$nuc) = unpack ("Z*Z*Z*",$str);
      my $batch          = $path."/batch_".$i."_".$no."_".$pid.".fa";
      $output{$batch}    = 1;
      $bp               += length($nuc);
      $snr++;
    
      if($bp <= $size){
        open WH,">>",$batch or LOGDIE "Cannot write the 2ndVmatchRun batch file $batch\n";
          print WH ">$id $des\n$nuc\n";
        close WH;
      }elsif($bp > $size and $snr == 1){
        open WH,">>",$batch or LOGDIE "Cannot write the 2ndVmatchRun batch file $batch\n";
          print WH ">$id $des\n$nuc\n";
        close WH;
      }else{
        open WH,">>",$batch or LOGDIE "Cannot write the 2ndVmatchRun batch file $batch\n";
          print WH ">$id $des\n$nuc\n";
        close WH;

        $bp  = length($nuc);
        $snr = 1;
        $i++;   
      }
    }
  }else{# array is empty
    $flag = "no";
  }
  
  return ($flag,\%output);
}
