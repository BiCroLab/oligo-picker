#!/usr/bin/perl -w
use strict;
use Log::Log4perl qw(:easy);
use IPC::System::Simple qw(capturex);
use File::Basename;

my ($vm_out, $vm_in, $pid, $arguments, $log, $tmpdir, $workdir, $script, $mail, $ref_count, $idx, $outdir, $logdir, $snictmp, $jpt, $preScript) = @ARGV;
require("$script/FileHandling.pm");
require("$script/DbHandling.pm");

Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

my @arg                                 = split(/::/,$arguments);
my $project                             = $arg[1];
my ($dist,$kmer,$gcmin,$gcmax)          = ($arg[9],$arg[2],$arg[5],$arg[6]);
my ($tm_delta,$driver,$dbname,$homopol) = ($arg[7],$arg[0],$arg[3],$arg[12]);
my ($user,$refname, $overlap)           = ($arg[13],$arg[14],$arg[15]);
my $hflag                               = "yes";                                                                             # set $hflag on "yes" if homopolymer filtering activated and on "no" otherwise
$hflag                                  = "no" if($homopol eq "na");
my ($todo,$rerun_dir)                   = check_jobs($jpt,$tmpdir,$pid);


if($todo eq "yes"){# there are jobs which needs to be rerun -> prepare cronjob
  my $cronjob     = FileHandling::createCronjob($rerun_dir,$log,$tmpdir,$pid,$project,$user,$mail,$preScript,"vm2nd");            # 0 = flag to show that there are no jobids yet
  my $id          = FileHandling::submit_job($cronjob,0,0,$log,"cron");
}else{# no jobs to rerun due to timeout
  capturex("rm","-f","-r",($rerun_dir));                                                                  # delete the directory since it should be empty

  # collect the vm_in (fasta) and vm_out (tab) files
  my @seq_files   = FileHandling::read_dir($vm_in,$log);                                                                       # read the fasta files from vm_in
  my @coord_files = FileHandling::read_dir($vm_out,$log);                                                                      # read the vmatch coordinates from vm_out
  LOGDIE "\tThe number of Vmatch input fasta files does not correspond with the number of Vmatch output files\n" if(scalar(@seq_files) != scalar(@coord_files));
  my @list        = FileHandling::combine_seq_coord(\@seq_files,\@coord_files,$log);                                           # combine the fasta files with the coord files

  # filter out those kmers with multiple matches
  my @filtered    = FileHandling::filter_2ndVmRun(\@list,$log,$pid,$tmpdir);
  LOGDIE "There are no oligonucleotides left over after all filtering steps! Please check your parameters. Maybe they are to stringent for process id $pid\n" if(!@filtered);
  INFO "Number of oligos for $pid: ".scalar(@filtered)."\n";

  # check if the DB loading and enquiry can be started
  my $final_dir = $tmpdir."/afterFiltering";
  capturex("mkdir", ($final_dir)) if(!(-d $final_dir));
  my $lock = $tmpdir."/lock";                                                                                                  # lock the file

  while(-e $lock){}                                                                                                            # do nothing if the lock file exists. Another process use it, wait until the process finishes
  LOGDIE "\tThat's strange! The lock file should not exist any more.\n" if(-e $lock);
  capturex("touch",($lock));

  # write the filtered sequences out ont the SNIC_TMP and move them to the output dir afterwards
  my $final_file = $snictmp."/oligos_pid_".$pid.".tab";
  my $flag       = FileHandling::write_tmp_file(\@filtered,$final_file, $log);
  LOGDIE "\tThe output file for PID $pid does not exist\n" if(!(-e $final_file));
  capturex("rsync",($final_file,$final_dir));

  INFO "\tCheck if all chromosomes were processed\n";
  # check if all files have been written
  my $all_done = 1;
  for(my $i = 0; $i < $ref_count; $i++){                                                                                       # 21 chromosomes, X and Y chromosomes, and the MT genome
    $all_done = (FileHandling::check($i,$final_dir) and $all_done);
  }

  capturex("rm","-f",($lock));
  INFO "\tDelete the vmatch sh, input and output directories for the pid $pid!\n";
  my $vmsh = $tmpdir."/vm_2ndRun/vm_sh_".$pid;
  capturex("rm","-f","-r",($vm_in));
  capturex("rm","-f","-r",($vm_out));
  capturex("rm","-f","-r",($vmsh));
  INFO "\t$all_done\n";

  # start the last script which filter out the data and load it to the DB if all chromosomes were processed
  if($all_done == 1){
    my @unique = FileHandling::mergeChromosomes($final_dir,$log);

    DbHandling::createDB($driver,$dbname,$kmer,$tm_delta,$gcmin,$gcmax,$hflag,$refname,$overlap,$log);
    my $tmpdb = $snictmp . "/" . basename($dbname);
    capturex("rsync",($dbname,$tmpdb));

    DbHandling::loadData2Db($driver,$tmpdb,$kmer,$tm_delta,$gcmin,$gcmax,\@unique,$hflag,$refname,$overlap,$log);
    capturex("rsync",($tmpdb,$dbname));
  
    if(scalar(@arg) == 16){
      INFO "\t".scalar(@unique)." unique oligonucleotides with a length of $kmer bases were successfully inserted into the database $dbname.\n";
      INFO "\tThere is no DB enquiry requested!\n";
    }elsif(scalar(@arg) == 19){                                                                                                                    # create output directory
      capturex("mkdir", ($outdir)) if(!(-d $outdir));
    
      my ($chr,$start,$stop) = ($arg[16],$arg[17],$arg[18]);
      my $result = DbHandling::enquire_db($chr,$start,$stop,$kmer,$tm_delta,$gcmin,$gcmax,$driver,$dbname,$hflag,$refname,$overlap,$log);                            # extract the data from the database
      my $out    = $outdir."/requestedRegion_".$chr."_".$start."_".$stop."_kmer".$kmer."_gcMin".$gcmin."_gcMax".$gcmax."_dTm".$tm_delta."_hpol".$hflag."_".$refname."_overlap".$overlap.".bed";
 
      FileHandling::write_DBresult($result,$out,$log);                                                                                             # write the output of the database to a bed file
    }else{
      LOGDIE "unknown number of elements: ".scalar(@arg)."\n";
    }
  
    capturex("rm","-f","-r",($tmpdir));                                                                                                               # delete the tmp directory
    INFO "\t....All jobs finished....Done\n";
  }else{
    INFO "\tJob $pid done! File $final_file written to $final_dir!\n";
  }
}


sub check_jobs{
  my $jfile = $_[0];
  my $tmp   = $_[1];
  my $pjd   = $_[2];
  my $rerun = $tmp."/rerun_".$pjd;
  my $flag  = "no";
  
  capturex("mkdir","-p", ($rerun)) if(!(-d $rerun));
  if(-f $jfile){
    open my $fs,"<",$jfile or LOGDIE "\tcan't read the file $jfile\n";
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
      }elsif($stat =~ /NODE_FAIL/){
        my $file = $dat;
        $flag    = "yes";
        capturex("rsync",($file,$rerun));
        INFO "\tThe job $file will be restarted due to the state \"node_fail\"\n";
      }elsif($stat =~ /COMPLETED/){# nothing to do - this is the ideal case
      }else{
        LOGDIE "\tJobid:".$id." failed with the state $stat. Please restart it manually\n";
      }
    }
    close $fs;
  }
  
  return ($flag,$rerun);
}
