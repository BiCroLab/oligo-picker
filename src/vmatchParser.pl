#!/usr/bin/perl -w
use strict;
use Log::Log4perl qw(:easy);
use IPC::System::Simple qw(capturex);
use File::Basename;

my ($vm_out,$vm_in,$arguments,$pid,$tm,$workdir,$script,$tmpdir,$log,$logdir,$outdir,$ref_count,$mail,$snictmp,$idx,$jpt,$preScript) = @ARGV;

require("$script/FileHandling.pm");
Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

my @arg                       = split(/::/,$arguments);
my ($project,$kmer,$tm_delta) = ($arg[1], $arg[2], $arg[7]);
my $user                      = $arg[13];
my ($todo,$rerun_dir)         = check_jobs($jpt,$tmpdir,$pid);


if($todo eq "yes"){
  my $cronjob = FileHandling::createCronjob($rerun_dir,$log,$tmpdir,$pid,$project,$user,$mail,$preScript,"vm1st");
  my $id      = FileHandling::submit_job($cronjob,0,0,$log,"cron");
}else{
  capturex("rm","-f","-r",($rerun_dir));

  my @seq_files                 = FileHandling::read_dir($vm_in,$log);                                             # collect the vm_in (fasta) and vm_out (tab) files
  my @coord_files               = FileHandling::read_dir($vm_out,$log);                                            # read the fasta files from vm_in                                                                                                                 
  LOGDIE "\tThe number of Vmatch input fasta files does not correspond with the number of Vmatch output files\n" if(scalar(@seq_files) != scalar(@coord_files));
  my @list                      = FileHandling::combine_seq_coord(\@seq_files,\@coord_files,$log);                 # combine the fasta files with the coord files
  my $tm_min                    = $tm - $tm_delta;
  my $tm_max                    = $tm + $tm_delta;

  my ($g,$dir,$res)             = FileHandling::createSubmissionScript_filtering(\@list,$log,$pid,$tmpdir,$mail,$project,$arguments,$tm_min,$tm_max,$script);
  my $sh_ctrl                   = FileHandling::createSubmissionScript_ctrl2ndVmRun($pid,$arguments,$mail,$tmpdir,$log,$script,$vm_in,$vm_out,$ref_count,$outdir,$res,$dir,$project,$idx,$logdir,$workdir);

  # create cronjob to run the filtering step after the first vmatch run
  my $cronjob     = FileHandling::createCronjob($g,$log,$tmpdir,$pid,$project,$user,$mail,$sh_ctrl,"filter");  
  my $id          = FileHandling::submit_job($cronjob,0,0,$log,"cron");
}


sub check_jobs{
  my $jfile = $_[0];
  my $tmp   = $_[1];
  my $pjd   = $_[2];
  my $rerun = $tmp."/rerun_".$pjd;
  my $flag  = "no";
  
  capturex("mkdir","-p", ($rerun)) if(!(-d $rerun));
  if(-f $jfile){
    open my $fs,"<",$jfile or die "can't read the file $jfile\n";
    while(<$fs>){
      chomp($_);
      next if(/^$/);
      my ($id,$dat) = split(/\t/,$_);
      
      my $cmd  = 'sacct -X --format=State --jobs='.$id." \| grep -v -E \"State\|-\"";    # construct command
      my $out  = `$cmd`;
      ($out =~ /([A-Z]+)/) or LOGDIE "\tError executing $cmd\n";
      my $stat = $1;
      
      if($stat =~ /TIMEOUT/){# restart the job after increasing the time to 200 hours
        my $file    = $dat;
        my $cmd_sed = 'sed -i \'\/#SBATCH -t \/c\\#SBATCH -t 100:00:00\' '.$file;
        my $done    = `$cmd_sed`;
        $flag       = "yes";
        capturex("rsync",($file,$rerun));      
        INFO "\tThe job $file will be restarted due to the state \"timeout\". The new time has been set to 100 hours\n";
      }elsif($stat =~ /CANCELLED/){
        my $file     = $dat;
        my $cmd_core = 'perl -pi -e \'s/#SBATCH -n (\d+)/"#SBATCH -n " . ($1 < 16 ? $1+1 : LOGDIE "\tJobid:".$jids[$i]." failed with the state \"$stat\". Please restart it manually\n")/ge\' '.$file;
        my $core     = $cmd_core;
        $flag        = "yes";
        capturex("rsync",($file,$rerun));
      }elsif($stat =~ /NODE_FAIL|FAILED/){
        my $file = $dat;
        $flag    = "yes";
        capturex("rsync",($file,$rerun));
        INFO "\tThe job $file will be restarted due to the state \"$stat\"\n";
      }elsif($stat =~ /COMPLETED/){# nothing to do - this is the ideal case
      }else{
        LOGDIE "\tJobid:".$id." failed with the state $stat. Please restart it manually\n";
      }
    }
  }
  
  return ($flag,$rerun);
}
