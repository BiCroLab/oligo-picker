#!/usr/bin/perl -w
use strict;
use File::Basename;
use Log::Log4perl qw(:easy);
use Bio::SeqIO;
use IPC::System::Simple qw(capturex);

my ($fa, $coord, $output, $tm_min, $tm_max, $arguments, $pwd, $log, $snictmp, $result) = @ARGV;

require("$pwd/FileHandling.pm");
Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

my @arg      = split(/::/,$arguments);
my $degreeC  = $arg[10];
my $salt     = $arg[8];
my $kmer     = $arg[2];

INFO "\tGet the coordinates and the sequences\n";
my %vmCoord  = FileHandling::process_vmatch_output($coord,$kmer,$log);

INFO "\tCheck the melting temperature for $fa\n";
my ($ssfa,$pref) = filter_tm($fa,$output);

INFO "\tCalculate minimal free energy for $ssfa\n";
capturex("hybrid-ss-min",("-o",$pref,"-n","DNA","-t",$degreeC,"-T",$degreeC,"-N",$salt,$ssfa));

INFO "\tFilter oligos with stable sec structure out\n";
my $dg_out       = $pref.".ct";
my %dgInfo       = get_dg($dg_out);
filter_dg(\%dgInfo,\%vmCoord,$ssfa,$output);

INFO "\tCopy the output $output to $result\n";
capturex("rsync",($output,$result));



sub filter_dg{
  my %dg_val = %{$_[0]};
  my %vmcoor = %{$_[1]};
  my $nucseq = $_[2];
  my $out    = $_[3];

  my $obj = Bio::SeqIO->new(-file => $nucseq, -format => 'fasta');

  open DF,">",$out or LOGDIE "\tCannot write the file $out\n";
  print DF "#id\tseq\tGC\tTm\tDg\tstart\tstop\tchr\n";

  while(my $seq = $obj->next_seq){
    my $nuc      = $seq->seq;
    my $id       = $seq->id;
    my $desc     = $seq->desc;
    my ($gc,$tm) = ($1,$2) if($desc =~ /gc=([0-9-_.]+);tm=([0-9-_.]+)/);

    if(exists($dg_val{$id}) and $dg_val{$id} > 0){# a free energy was calculated and it is > 0
      if(exists($vmcoor{$id})){# and a genomic coordinate exists
        my ($chr,$start,$stop) = unpack("Z*ii",$vmcoor{$id});
        print DF "$id\t$nuc\t$gc\t$tm\t".$dg_val{$id}."\t$start\t$stop\t$chr\n";
      }
    }
  }
  close DF;
}


sub filter_tm{
  my $file = $_[0];
  my $out  = $_[1];
  my $obj  = Bio::SeqIO->new(-file => $file, -format => 'fasta');

  my $base   = basename($out); 
  my $ss_fa  = $snictmp."/".$base."_ss.fa";
  my $prefix = $snictmp."/".$base."_ss_dg";

  open FL,">",$ss_fa or LOGDIE "\tCannot write the file $ss_fa\n";
  while(my $seq = $obj->next_seq){
    my $nuc      = $seq->seq;
    my $id       = $seq->id;
    my $desc     = $seq->desc;
    my ($gc,$tm) = ($1,$2) if($desc =~ /gc=([0-9-_.]+);tm=([0-9-_.]+)/);

    next if($tm < $tm_min or $tm_max < $tm);                                    # filter tm, which lie outside the interval
    print FL ">$id gc=$gc;tm=$tm\n";
    print FL "$nuc\n";
  }
  close(FL);

  return ($ss_fa,$prefix);
}


sub get_dg{
  my $file = $_[0];
  my %list = ();

  open NH,"<",$file or LOGDIE "\tCannot read the file $file\n";
  while(<NH>){
    chomp($_);
    next if(/^$/);
    
    if(/dG\s+=\s+([0-9-.]+)\t(candidate\_[0-9]+\_[0-9]+)/){
      $list{$2} = $1;
    }else{
      next;
    }
  }
  close NH;

  return %list;
}
