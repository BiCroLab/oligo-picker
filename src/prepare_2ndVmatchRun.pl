#!/usr/bin/perl -w
use strict;
use Log::Log4perl qw(:easy);
use IPC::System::Simple qw(capturex);
use File::Basename;

my ($fa_dir, $pid, $arguments, $mail, $script, $tmpdir, $idx, $workdir, $log, $ref_count, $outdir, $logdir, $jpt, $preScript) = @ARGV;
  
require("$script/FileHandling.pm");
Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

my @arg                             = split(/::/,$arguments);
my ($project,$kmer,$dist,$mismatch) = ($arg[1],$arg[2],$arg[9],$arg[11]);
my $user                            = $arg[13];
my ($todo,$rerun_dir)               = check_jobs($jpt,$tmpdir,$pid);

if($todo eq "yes"){
  my $cronjob     = FileHandling::createCronjob($rerun_dir,$log,$tmpdir,$pid,$project,$user,$mail,$preScript,"filter");                  # 0 = flag to show that there are no jobids yet
  my $id          = FileHandling::submit_job($cronjob,0,0,$log,"cron");
}else{# no jobs to rerun due to timeout
  capturex("rm","-f","-r",($rerun_dir));                                                                  # delete the directory since it should be empty

  # create the sh, vm_in and vm_out folder for the second run
  my $run_dir = $tmpdir."/vm_2ndRun";
  my $vm_sh   = $run_dir."/vm_sh_".$pid;
  my $vm_in   = $run_dir."/vm_in_".$pid;
  my $vm_out  = $run_dir."/vm_out_".$pid;
  capturex("mkdir","-p",($run_dir)) if(!(-d $run_dir));
  capturex("mkdir","-p",($vm_sh))   if(!(-d $vm_sh));
  capturex("mkdir","-p",($vm_in))   if(!(-d $vm_in));
  capturex("mkdir","-p",($vm_out))  if(!(-d $vm_out));

  # collect the tm, sec. struc., and gc filtered sequences and save them in an hash
  my @fa_filtered = FileHandling::read_dir($fa_dir,$log);
  my @fa_list;

  foreach my $file (@fa_filtered){
    next if(-z $file);                        # next if file exists and has 0 size
    open HS,"<",$file or LOGDIE "\tCannot read the file $file\n";
    while(<HS>){
      chomp($_);
      next if(/^$/);
      next if(/^#id.*/);
      my ($id,$seq,$gc,$tm,$dg,$start,$stop,$chr) = split(/\t/,$_);
      my $new_id = $1."_".$chr if($id =~ /(candidate_\d+)_\d+/);             # change the pid to chr id
      my $str    = pack("Z*Z*iiZ*fff",$new_id,$seq,$start,$stop,$chr,$gc,$tm,$dg);
      push(@fa_list,$str);
    }
    close(HS);
  }

  # filter out overlapping kmers
  my @kmer_noOverlap  = FileHandling::filter_overlapping_pos(\@fa_list, $dist, $pid, $log);
  LOGDIE "There are no oligonucleotides left over after all filtering steps! Please check your parameters. Maybe they are to stringent for process id $pid\n" if(!@kmer_noOverlap);

  # split the files in smaller batches and prepare the second vmatch run
  my %batches         = create_batches(\@kmer_noOverlap,$vm_in);
  my ($sh_vm,$out_vm) = FileHandling::prepare_vmatch_2ndRun(\%batches,$tmpdir,$vm_out,$vm_sh,$log,$pid,$kmer,$project,$mail,$idx,$mismatch);
  my $sh_vmParser     = FileHandling::createSubmissionScript_vmParser_2ndRun($out_vm,$vm_in,$pid,$arguments,$log,$tmpdir,$workdir,$script,$mail,$ref_count,$idx,$outdir,$logdir); 

  # create cronjob to run the filtering step after the first vmatch run
  my $cronjob         = FileHandling::createCronjob($sh_vm,$log,$tmpdir,$pid,$project,$user,$mail,$sh_vmParser,"vm2nd");  
  my $id              = FileHandling::submit_job($cronjob,0,0,$log,"cron");

  # delete the results of the first vmatch run, since they are not needed any more
  my $run1st = $tmpdir."/vm_1stRun/*_".$pid."/";
  system("rm -r $run1st");
}


# METHODS
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


sub create_batches{
  my @list               = @{$_[0]};
  my $path               = $_[1];
  my ($size,$bp,$i,$snr) = (5000,0,1,1);
  my %out_files          = ();
  
  foreach my $row (@list){
    my ($name,$seq,$start,$stop,$chr,$gc,$tm,$dg) = unpack("Z*Z*iiZ*fff", $row);
    my $desc           = $start."::".$stop."::".$chr."::".$gc."::".$tm."::".$dg;
    my $batch          = $path."/batch_".$i."_".$pid.".fa";
    $out_files{$batch} = 1;
    $bp               += length($seq);
    $snr++;
    
    if($bp <= $size){
      open WH,">>",$batch or LOGDIE "Cannot write the 2ndVmatchRun batch file $batch\n";
        print WH ">$name $desc\n$seq\n";
      close(WH);
    }elsif($bp > $size and $snr == 1){
      open WH,">>",$batch or LOGDIE "Cannot write the 2ndVmatchRun batch file $batch\n";
        print WH ">$name $desc\n$seq\n";
      close(WH);
    }else{
      open WH,">>",$batch or LOGDIE "Cannot write the 2ndVmatchRun batch file $batch\n";
        print WH ">$name $desc\n$seq\n";
      close(WH);
      $bp  = length($seq);
      $snr = 1;
      $i++;
    }
  }

  return %out_files;
}

