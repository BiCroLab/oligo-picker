#!/usr/bin/perl -w
use strict;
use Log::Log4perl qw(:easy);
use IPC::System::Simple qw(capturex);
use File::Basename;

my ($res_dir, $arguments, $pid, $script, $tmp_dir, $out_dir, $log, $snictmp, $vmin, $vmout, $filterdir, $ref_count) = @ARGV;

require("$script/FileHandling.pm");
require("$script/DbHandling.pm");
Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

my $vmsh = $tmp_dir."/vmRun/vmatch_sh_".$pid;
INFO "\tDelete the vmatch sh, input and output directories for the pid $pid: $vmin, $vmout, $vmsh!\n";
#capturex("rm",("-r",$vmin));
#capturex("rm",("-r",$vmout));
#capturex("rm",("-r",$vmsh));

my @arg      = split(/::/,$arguments);
my $dist     = $arg[9];
my $kmer     = $arg[2];
my $gcmin    = $arg[5];
my $gcmax    = $arg[6];
my $tm_delta = $arg[7];
my $driver   = $arg[0];
my $dbname   = $arg[3];
my $homopol  = $arg[11];
my $hflag    = "yes";                                                  # set $hflag on "yes" if homopolymer filtering activated and on "no" otherwise
$hflag       = "no" if($homopol eq "na");
my @list;


my @results  = FileHandling::read_dir($res_dir,$log);

foreach my $file (@results){
  next if(-z $file);                                                    # next if file exists and has zero size
  
  open HS,"<",$file or LOGDIE "\tCannot read the file $file\n";
  while(<HS>){
    chomp($_);
    next if(/^#id.*/);
    my ($id,$seq,$gc,$tm,$dg,$start,$stop,$chr) = split(/\t/,$_);
    my $new_id = $1."_".$chr if($id =~ /(candidate_\d+)_\d+/);          # change the pid to chr id

    my $str = pack("Z*Z*iiZ*fff",$new_id,$seq,$start,$stop,$chr,$gc,$tm,$dg);
    push(@list,$str);
  }
  close HS;  
}

my @filtered = FileHandling::filter_overlapping_pos(\@list,$dist,$pid,$log);                  # find non-overlapping kmers
LOGDIE "There are no oligonucleotides left over after all filtering steps! Please check your parameters. Maybe they are to stringent for process id $pid\n" if(!@filtered);

# check if the DB loading and enquiry can be started
my $fin  = $tmp_dir."/afterFiltering";                                                                                                     # create the output directory for the filtered files, if it does not exist
capturex("mkdir", ($fin))   if(!(-d $fin));
my $lock = $tmp_dir."/lock";                                                                                                                # lock file

while(-e $lock){}                                                                                                                          # do nothing if the lock file exists. Another process use it, wait until the process finishes
LOGDIE "\tThat's strange! The lock file should not exist any more.\n" if(-e $lock);
capturex("touch", ($lock));

# write the filtered sequences out on the SNIC_TMP and move them to the output dir afterwards
my $fin_seq  = $snictmp."/oligos_pid_".$pid.".tab";
my $flag     = FileHandling::write_tmp_file(\@filtered,$fin_seq,$log,$pid);
LOGDIE "\tThe output file for PID $pid does not exist\n" if(!(-e $fin_seq));
capturex("rsync",($fin_seq,$fin));

INFO "\tCheck if all chromosomes were processed\n";
# check if all files have been written
my $all_done = 1;
for(my $i = 0; $i < $ref_count; $i++){                                                             # 22 chromosomes, X and Y chromosomes, and the MT genome
  $all_done = (FileHandling::check($i,$fin) and $all_done);
}

capturex("rm", ($lock));
capturex("rm",("-r",$filterdir));                                                          # remove the filtered results (sh and tab files)

INFO "\t$all_done\n";

# start the last script which filter out the data and load it to the DB if all chromosomes were processed
if($all_done == 1){
  my @unique = FileHandling::mergeChromosomes($fin,$log);
  
  DbHandling::createDB($driver,$dbname,$kmer,$tm_delta,$gcmin,$gcmax,$hflag,$log);
  my $tmpdb = $snictmp . "/" . basename($dbname);
  capturex("rsync",($dbname,$tmpdb));

  DbHandling::loadData2Db($driver,$tmpdb,$kmer,$tm_delta,$gcmin,$gcmax,\@unique,$hflag,$log);
  capturex("rsync",($tmpdb,$dbname));
  
  if(scalar(@arg) == 12){
    INFO "\t".scalar(@unique)." unique oligonucleotides with a length of $kmer bases were successfully inserted into the database $dbname.\n";
    INFO "\tThere is no DB enquiry requested!\n";
  }elsif(scalar(@arg) == 15){                                                                                    # create output directory
    capturex("mkdir", ($out_dir)) if(!(-d $out_dir));
    
    my ($chr,$start,$stop) = ($arg[12],$arg[13],$arg[14]);
    my $result = DbHandling::enquire_db($chr,$start,$stop,$kmer,$tm_delta,$gcmin,$gcmax,$driver,$dbname,$hflag,$log);                            # extract the data from the database
    my $out    = $out_dir."/requestedRegion_".$chr."_".$start."_".$stop."_kmer".$kmer."_gcMin".$gcmin."_gcMax".$gcmax."_dTm".$tm_delta."_hpol".$hflag.".bed";
 
    FileHandling::write_DBresult($result,$out,$log);                                                                                      # write the output of the database to a bed file
  }else{
    LOGDIE "unknown number of elements: ".scalar(@arg)."\n";
  }
  
  capturex("rm",("-r",$tmp_dir));                                                                                                         # delete the tmp directory
  INFO "\t....All jobs finished....Done\n";
}else{
  INFO "\tJob $pid done! File $fin_seq written to $fin!\n";
}
