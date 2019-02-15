package FileHandling;

use strict;
use warnings;
use IPC::System::Simple qw(capturex);
use Log::Log4perl qw(:easy);
use Bio::SeqIO;
use File::Basename;
use POSIX;
require Exporter;

our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(&createCronjob &createSubmissionScript_vmParser_2ndRun &filter_2ndVmRun &createSubmissionScript_filtering &combine_seq_coord &submit_job &filter_overlapping_pos &get_log10 &write_DBresult &process_jellyfish_output &prepare_vmatch_1stRun &createSubmissionScript_jellyfish &createSubmissionScript_jellyfishParser &read_dir &process_vmatch_output &check &write_tmp_file &mergeChromosomes &createSubmissionScript_vmParser &createSubmissionScript_vmatchIdx &prepare_vmatch_2ndRun &createSubmissionScript_ctrl2ndVmRun &createSubmissionScript_ctrlLastStep &write_tmp_aftervm1);


#############################
# Preloaded methods go here #
#############################
sub createCronjob{
  my ($vm_path,$log,$tmpdir,$pid,$project,$user,$mail,$parser,$flag) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  # create cronjob for the submission of vmatch jobs
  my $cj_dir = $tmpdir."/cronjobs/";
  capturex("mkdir","-p", ($cj_dir)) if(!(-d $cj_dir));
  my $cj       = $cj_dir."cronjob_".$pid."_".$flag.".sh";
  my $myscript = $cj_dir."check_".$pid."_".$flag.".sh";
  my $fjid     = $cj_dir."jobids_paths_".$pid."_".$flag.".tab";
  my $jobids   = $cj_dir."jobids_".$pid."_".$flag.".tab";
 
  open HJ,">",$myscript or LOGDIE "\tCannot write the help cronjob script\n";
    print HJ "#!/bin/bash -l\n\n";
    print HJ "out=\$1\n";
    print HJ "jobids=\$2\n";
    print HJ "file=\$3\n";
    print HJ "\n";
    print HJ "jobnr=`squeue -u $user | grep -v \"JOBID\" | wc -l`\n";
    print HJ "submit=\$file.toSubmit\n\n";
    print HJ "# if file extension \"toSubmit\" exists, then submit the job\n";
    print HJ "if [[ -e \"\$submit\" ]]\n";
    print HJ "  then\n";
    print HJ "    if (( \"\$jobnr\" < 300 ))\n";
    print HJ "      then\n";
    print HJ "        # submit the job\n";
    print HJ "        line=`sbatch \$file 2>&1`\n";
    print HJ "        err=`echo \$line | grep -ic \"error\"`\n";
    print HJ "        if (( \"\$err\" == 0 ))   # no error \n";
    print HJ "          then\n";
    print HJ "            id=`echo \$line | grep -o \'[0-9]*\'`\n\n";
    print HJ "            # add jobid to the other ids\n";
    print HJ "            echo -e \"\$id\\t\$file\" >> \$out\n";
    print HJ "            echo \$id >> \$jobids\n\n";
    print HJ "            usleep 100000\n";
    print HJ "            # delete processed file\n";
    print HJ "            rm \$submit\n";
    print HJ "        fi\n";
    print HJ "    fi\n";
    print HJ "fi\n";
  close HJ;
  capturex("chmod","a+wx", ($myscript));

  open CJ,">",$cj or LOGDIE "\tCannot write the cronjob script to organise the submission of vmatch/filtering jobs\n";
    print CJ "#!/bin/bash -l\n";
    print CJ "#SBATCH -A $project\n";
    print CJ "#SBATCH --qos=interact\n";
    print CJ "#SBATCH -p core -n 1\n";
    print CJ "#SBATCH -t 10:00:00\n";
    print CJ "#SBATCH -J cronJob\n";
    print CJ "#SBATCH --open-mode=append\n";
    print CJ "#SBATCH -o $log\n";
    print CJ "#SBATCH -e $log\n";
    print CJ "#SBATCH --mail-type=FAIL\n";
    print CJ "#SBATCH --mail-user=$mail\n\n";
    print CJ "path=\"$vm_path\"\n";
    print CJ "parser=\"$parser\"\n";
    print CJ "fid=\"$fjid\"\n";
    print CJ "jid=\"$jobids\"\n";
    print CJ "\n\n";
    print CJ "# go through all sh jobs and submit them to the queue if the number of submitted jobs does not succeed 300\n";
    print CJ "/bin/find \$path -type f -name \'*.sh\' -print0 | /usr/bin/xargs -n 1 -0 $myscript \$fid \$jid\n";
    print CJ "\n";
    print CJ "# check if there are still jobs, which need to be submitted\n";
    print CJ "check=`/bin/find \$path/ -type f -name '*.toSubmit' 2>/dev/null | /usr/bin/wc -l`\n";
    print CJ "if [ \$? -ne 0 ]\n";
    print CJ "  then\n";
    print CJ "    echo \"unable to check how many files still need to be submitted\"\n";
    print CJ "    exit 1\n";
    print CJ "fi\n\n";
    print CJ "if [[ \$check != 0 ]]\n";
    print CJ "  then\n";
    print CJ "    # resubmit this script with --begin set to run the job every 15 minutes\n";
    print CJ "    sbatch --quiet --begin=now+20hour $cj\n";
    print CJ "  else\n";
    print CJ "    sbatch --dependency=afterok:\$(/bin/cat \$jid \| /usr/bin/tr '\\n' '\:' \| /bin/sed 's/\:\$//') \$parser \$fid \$parser\n";
    print CJ "fi\n";
  close(CJ);
  capturex("chmod","a+wx", ($cj)); 
  return $cj;
}


sub write_tmp_aftervm1{
  my @list = @{$_[0]};
  my $file = $_[1];
  my $log  = $_[2];
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  open my $fh,">",$file or LOGDIE "Cannot write the file $file\n";
  foreach my $elm (@list){
    my($sid,$seq,$start,$stop,$chr,$gc,$tm,$dg) = unpack("Z*Z*iiZ*fff");
    print $fh "$sid\t$seq\t$start\t$stop\t$chr\t$gc\t$tm\t$dg\t2\n";                       # '2' specifies that only the first vmatch was executed, the second was not requested
  }
  close($fh);
  
  my $flag = "done";
  return $flag;
}


sub write_tmp_file{
  my @list = @{$_[0]};
  my $file = $_[1];
  my $log  = $_[2];
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  open my $fh,">",$file or LOGDIE "Cannot write the file $file\n";
  foreach my $elm (@list){
    my ($sid,$desc,$nuc,$flag)         = unpack("Z*Z*Z*i",$elm);
    my ($start,$stop,$chr,$gc,$tm,$dg) = split(/::/,$desc);
    print $fh "$sid\t$nuc\t$start\t$stop\t$chr\t$gc\t$tm\t$dg\t$flag\n";
  }
  close($fh);

  my $flag = "done";
  
  return $flag;
}


sub filter_2ndVmRun{
  my @list   = @{$_[0]};
  my $log    = $_[1];
  my $pid    = $_[2];
  my $tmpdir = $_[3];
  my @uniq = ();
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  foreach my $elm (@list){
    my @files = unpack("Z*Z*Z*",$elm);
    my $fa    = $files[0];
    my $vm    = $files[1];
    my %tmp_vm;

    open my $vh,"<",$vm or LOGDIE "\tCannot read the 2nd Vmatch output: $vm\n";
    while(<$vh>){
      chomp($_);
      next if(/^$/);                            # drop empty lines
      next if(/^#/);                            # drop description lines
      my ($chr,$pos,$name) = split(/\t/,$_);
      $tmp_vm{$name} = $chr."::".$pos;
    }
    close($vh);
 
    # select process_vmatch_output
    my $obj = Bio::SeqIO->new(-file => $fa, -format => 'fasta');
    while(my $o = $obj->next_seq){
      my $head = $o->id;
      if(exists($tmp_vm{$head})){# 2nd vm run successful -> discarded = false = 0
        my $str = pack("Z*Z*Z*i",$head,$o->desc,$o->seq,0);
        push(@uniq,$str);
      }else{# filtered out in the 2nd vm run -> discarded = true = 1
        my $str = pack("Z*Z*Z*i",$head,$o->desc,$o->seq,1);
        push(@uniq,$str);
      }
    }
  }

  return @uniq;
}


sub prepare_vmatch_2ndRun{
  my %vm_in                        = %{$_[0]};
  my ($tmpdir,$vm_out,$vm_sh,$log) = ($_[1],$_[2],$_[3],$_[4]);
  my ($pid,$kmer,$proj,$mail,$idx) = ($_[5],$_[6],$_[7],$_[8],$_[9]);
  my $mismatch                     = $_[10];
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  INFO "\tCreate 2nd Vmatch run submission scripts for process id $pid\n";
  my $mism       = ceil(($mismatch/$kmer)*100);
  my @list       = keys(%vm_in);
  my $id         = 0;
  my $vm2process = 48;
  my $idxn       = basename($idx);
  
  # create the vmatch job itself
  my $run = $vm_sh."/run_vm2nd_".$pid.".sh";
  open my $hw,">",$run or LOGDIE "\tcan't write the file $run\n";
    print $hw "#!/bin/bash -l\n\n";
    print $hw "file=\$(basename \$1)\n";
    print $hw "pid=\$(basename \$1 .fa \| cut -d \"_\" -f 2,3)\n";    
    print $hw "dir=$vm_out\n";
    print $hw "out=\"coordinates_uniqueOligos_\"\$pid\".tab\"\n";
    print $hw "MISMATCH=$mism\n\n";
    print $hw "# check if the result of this sample is already created\n";
    print $hw "if [ -e \$dir/\$out ]\n";
    print $hw "  then exit 0\n";
    print $hw "fi\n\n";
    print $hw "# copy sample to node\n";
    print $hw "rsync \$1 \$SNIC_TMP\n\n";
    print $hw "# run vmatch\n";
    print $hw "vmatch -q \$SNIC_TMP/\$file -complete -p -d -h \${MISMATCH}b -selfun mydb.so -showdesc 100 \$SNIC_TMP/$idxn > \$SNIC_TMP/\$out\n\n";
    print $hw "# copy the results back to network storage\n";
    print $hw "rsync \$SNIC_TMP/\$out \$dir/\$out\.incomplete\n\n";
    print $hw "# rename file once the copying is done\n";
    print $hw "mv \$dir/\$out\.incomplete \$dir/\$out\n";
    print $hw "rm \$SNIC_TMP/\$out\n";
  close $hw;
  capturex("chmod","a+wx", ($run));

  # create the wrapper which manages 48 batches in one job, running 16 at one time
  for(my $i = 0; $i < scalar(@list); $i+=$vm2process){
    my $tmp_ls = $vm_sh."/seqList_".$pid."_".$id.".txt";                   # write a tmp text file with the links to the files
    
    open my $wh,">",$tmp_ls or LOGDIE "\tcan't write the file $tmp_ls\n";
    for(my $j = $i; $j < ($i + $vm2process) and $j < scalar(@list); $j++){
      print $wh $list[$j]."\n";
    }
    close $wh;
    
    # create the shell script and a dummy file to control the submission
    my $sh    = $vm_sh."/vmatch_".$pid."_".$id.".sh";
    my $dummy = $sh.".toSubmit";
    my $jname = "vm2nd_".$pid."_".$id;
    my $dir   = dirname($idx);
    capturex("touch",($dummy));
    
    open my $vm, ">", $sh or LOGDIE "\tCannot write the shell script for the 2nd Vmatch run for batch file $sh (pid $pid)!\n";
      print $vm "#!/bin/bash -l\n";
      print $vm "#SBATCH -A $proj\n";
      print $vm "#SBATCH -p node\n";
      print $vm "#SBATCH -n 1\n";
      print $vm "#SBATCH -t 170:00:00\n";
      print $vm "#SBATCH -J $jname\n";
      print $vm "#SBATCH --open-mode=append\n";
      print $vm "#SBATCH -o $log\n";
      print $vm "#SBATCH -e $log\n";
      print $vm "#SBATCH --mail-type=FAIL\n";
      print $vm "#SBATCH --mail-user=$mail\n\n";
      print $vm "REF=$idx\n\n";
      print $vm "rsync -r $dir/ \$SNIC_TMP\n";
      print $vm "wait\n\n";
      print $vm "threads=\$SLURM_CPUS_ON_NODE\n";
      print $vm "cat $tmp_ls \| xargs -i --max-procs=\$threads bash -c \"bash $run {}\"\n";
    close $vm;
    capturex("chmod","a+wx", ($sh));   
    $id++;
  }
 
  INFO $id." jobs submitted for the second vmatch run\n";
  return ($vm_sh,$vm_out);
}


sub createSubmissionScript_ctrl2ndVmRun{
  my ($pid,$arguments,$mail,$tmpdir,$log,$script,$vm_in,$vm_out,$ref_count,$out_dir,$res,$dir,$project,$index,$logdir,$workdir) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
     
  my $sh_name = $tmpdir."/my_prepare2ndVmatchRun_".$pid.".sh";
  open OS,">",$sh_name or LOGDIE "\tCannot write the shell script which organises the second vmatch run\n";
    print OS "#!/bin/bash -l\n";
    print OS "#SBATCH -A $project\n";
    print OS "#SBATCH --qos=interact\n";
    print OS "#SBATCH -p core\n";
    print OS "#SBATCH -n 5\n";
    print OS "#SBATCH -t 11:59:00\n";
    print OS "#SBATCH -J prepare_2ndVmRun\n";
    print OS "#SBATCH --open-mode=append\n";
    print OS "#SBATCH -o $log\n";
    print OS "#SBATCH -e $log\n";
    print OS "#SBATCH --mail-type=FAIL\n";
    print OS "#SBATCH --mail-user=$mail\n\n";
    print OS "ARG=\'$arguments\'\n";
    print OS "JPT=\$1\n";
    print OS "NAME=\$2\n";
    print OS "perl $script/prepare_2ndVmatchRun.pl $res $pid \$ARG $mail $script $tmpdir $index $workdir $log $ref_count $out_dir $logdir \$JPT \$NAME\n";
  close(OS);
  capturex("chmod","a+wx", ($sh_name));

  return $sh_name;
}


sub createSubmissionScript_ctrlLastStep{
  my ($pid,$arguments,$mail,$tmpdir,$log,$script,$vm_in,$vm_out,$ref_count,$out_dir,$res,$dir,$project,$index,$logdir,$workdir) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
     
  my $sh_name = $tmpdir."/my_prepareLastStep_".$pid.".sh";
  open OS,">",$sh_name or LOGDIE "\tCannot write the shell script which organises the last step after running only one vmatch run\n";
    print OS "#!/bin/bash -l\n";
    print OS "#SBATCH -A $project\n";
    print OS "#SBATCH --qos=interact\n";
    print OS "#SBATCH -p core\n";
    print OS "#SBATCH -n 5\n";
    print OS "#SBATCH -t 11:59:00\n";
    print OS "#SBATCH -J lastStep\n";
    print OS "#SBATCH --open-mode=append\n";
    print OS "#SBATCH -o $log\n";
    print OS "#SBATCH -e $log\n";
    print OS "#SBATCH --mail-type=FAIL\n";
    print OS "#SBATCH --mail-user=$mail\n\n";
    print OS "ARG=\'$arguments\'\n";
    print OS "JPT=\$1\n";
    print OS "NAME=\$2\n";
    print OS "perl $script/finish_afterOnlyOneVmatchRun.pl $res $pid \$ARG $mail $script $tmpdir $index $workdir $log $ref_count $out_dir $logdir \$JPT \$NAME \$SNIC_TMP\n";
  close(OS);
  capturex("chmod","a+wx", ($sh_name));

  return $sh_name;
}


sub get_log10{
  my $n = shift;
  return log($n)/log(10);
}

sub combine_seq_coord{
  my @seq  = @{$_[0]};
  my @coor = @{$_[1]};
  my $log  = $_[2];
  my @list = ();
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  for(my $i = 0; $i < scalar(@seq); $i++){
    my $id = $1 if($seq[$i] =~ /batch_(\d+.*\d+)\.fa/);
  
    for(my $j = 0; $j < scalar(@coor);$j++){
      my $cid = $1 if($coor[$j] =~ /coordinates_uniqueOligos_(\d+_\d+)\.tab/);

      # not only the ids should overlap, but also the files should not be empty
      if(($cid eq $id) and !(-z $coor[$j])){
        my $out = "filtered_uniqueOligos_".$id.".tab";
        my $str = pack("Z*Z*Z*",$seq[$i],$coor[$j],$out);
        push(@list,$str);
      }
    }
  }
  
  return @list;
}


sub process_jellyfish_output{
  my ($in,$pid,$tmpdir,$gcmin,$gcmax,$kmer,$salt,$formamid,$homopoly,$log,$snictmp) = @_;
  my ($loc,$nuc_count,$bad_gc,$bp,$id,$sum,$poly) = (0,0,0,0,0,0,0,0);
  my ($i,$seq_nr)                                 = (1,1);
  my $batch_size                                  = 5000000;                                       # size of a chunck: 5 Mb
  my %collect                                     = ();
  
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  # create vmatch input directory
  my $dir = $tmpdir."/vm_1stRun/vmatch_in_".$pid;
  capturex("mkdir","-p", ($dir)) if(!(-d $dir));
  
  # copy the jellyfish output to SNIC_TMP
  my $tmpin = $snictmp."/jellyfishOut_".$pid.".fa";
  capturex("rsync", ($in,$tmpin));
  
  # check and prepare for homopolymer filtering
  if($homopoly ne "na"){
    $homopoly = $1 if($homopoly =~ /(.*);$/);                                    # remove ";" if it exists at the end of the string
    $homopoly =~ tr/;/|/ if($homopoly =~ /;/);                                   # replace ";" with "|"
    
    INFO "\tProcess the jellyfish output for pid $pid: filter homopolymers, calculate tm and GC, filter the oligos without a correct GC, and split the file in 5 Mb chunks.\n";
  }else{
    INFO "\tProcess the jellyfish output for pid $pid: calculate tm and GC, filter the oligos without a correct GC value, and split the file in 5 Mb chunks.\n";
  }
  
  my $obj = Bio::SeqIO->new(-file => $tmpin, -format => 'fasta');                # read the fasta file using Bio::SeqIO
  
  while(my $seq = $obj->next_seq){
    my $nuc = $seq->seq;
    my $gc  = sprintf("%.1f", ( (($nuc =~ tr/C//) + ($nuc =~ tr/G//)) / (($nuc =~ tr/A//) + ($nuc =~ tr/T//) + ($nuc =~ tr/C//) + ($nuc =~ tr/G//)) * 100));
    $nuc_count++;

    if($gc < $gcmin or $gcmax < $gc){# remove oligos with an GC-content that lies outside of the interval
      $bad_gc++;
      next;
    }
    
    if($homopoly ne "na" and $nuc =~ /$homopoly/){# remove homopolymers if they exists and the flag not "na"
      $poly++;
      next;
    }
    
    my $tm = 81.5 + (0.41 * $gc) + (16.6 * get_log10($salt)) - (500 / $kmer) - (0.62 * $formamid);                    # calculate tm temperature
    $sum  += $tm;
    
    my $batch_file        = $snictmp."/batch_".$i."_".$pid.".fa";
    $collect{$batch_file} = 1;
    $bp                  += length($nuc);
    $seq_nr++;
    $id++;
    
    if($bp <= $batch_size){
      open SP,'>>',$batch_file or LOGDIE "Cannot write the Vmatch batch file $batch_file\n";
        print SP ">candidate_".$id."_".$pid." gc=$gc;tm=$tm\n$nuc\n";
      close SP;    
    }elsif($bp > $batch_size and $seq_nr == 1){
      open SP,'>>',$batch_file or LOGDIE "Cannot write the Vmatch batch file $batch_file\n";
        print SP ">candidate_".$id."_".$pid." gc=$gc;tm=$tm\n$nuc\n";
      close SP;    
    }else{
      open SP,'>>',$batch_file or LOGDIE "Cannot write the Vmatch batch file $batch_file\n";
        print SP ">candidate_".$id."_".$pid." gc=$gc;tm=$tm\n$nuc\n";
      close SP;
      
      # set the counters back
      $bp     = length($nuc);
      $seq_nr = 1;
      $i++;
    }
  }
  
  my $tm_avg = sprintf("%.3f",($sum/$id));

  INFO "\tunique jellyfish oligos\t$nuc_count\n\tGC condition fulfilled\t$id\nGC rejected\t$bad_gc\n\taverage melting temperature\t$tm_avg\n";
  INFO "\tremoved oligos due to homopolymers\t$poly\n" if($homopoly ne "na");
  INFO "\t$i batch files generated\n";
    
  foreach my $elm (keys %collect){                      # copy the batch files from snictmp to tmp
    capturex("rsync", ($elm,$dir));
  }
  capturex("rm", ($in));                                # remove the jellyfish input file, since it is no longer necessary
  
  return ($dir,$tm_avg);
}


sub mergeChromosomes{
  my ($dir,$log) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  my @list  = read_dir($dir,$log);
  my $flag  = "yes";                                             # first file in the list
  my %collect;
  my $count = 0;
  
  INFO "\tProcess the ".scalar(@list)." chromosomes by removing redundant oligonucleotides\n";
  # identify redundant oligonucleotides
  foreach my $file (@list){
    open(GU,$file) or LOGDIE "Cannot read the file $file\n";
      while(<GU>){
        chomp($_);
        my @tmp = split(/\t/,$_);
        $count++;
        
        if(exists($collect{$tmp[1]})){
          my @elm = unpack("Z*iiZ*fffii",$collect{$tmp[1]});
          $elm[8] += 1;
          my $str = pack("Z*iiZ*fffii",$elm[0],$elm[1],$elm[2],$elm[3],$elm[4],$elm[5],$elm[6],$elm[7],$elm[8]);
          $collect{$tmp[1]} = $str;
        }else{
          my $x   = 1;
          my $str = pack("Z*iiZ*fffii",$tmp[0],$tmp[2],$tmp[3],$tmp[4],$tmp[5],$tmp[6],$tmp[7],$tmp[8],$x);
          $collect{$tmp[1]} = $str;
        }
      }
    close(GU);
  }
  
  # remove redundant oligonucleotides
  my @uniq;
  
  foreach my $oligo (keys %collect){
    my @tmp = unpack("Z*iiZ*fffii",$collect{$oligo});
    if($tmp[8] == 1){
      my $str = pack("Z*Z*iiZ*fffi",$tmp[0],$oligo,$tmp[1],$tmp[2],$tmp[3],$tmp[4],$tmp[5],$tmp[6],$tmp[7]);
      push(@uniq,$str);
    }
  }
  
  INFO "\t".scalar(@uniq)." oligos (out of $count) are unique in the whole human genome\n";
  
  return @uniq;
}


sub submit_job{
  my ($script_file,$is_chain,$job_id,$log,$flag) = @_;                                             # SLURM script to be run; 0: not chained, 1: chained; job ID variable: 0 if no chaining, > 0 if chained; log file
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
              
  # construct call to sbatch
  my $cmd = 'sbatch';
  $cmd   .= " --dependency=afterOK:$job_id" if($job_id !~ /^0$/ && $is_chain == 1);          # if chained add dependency
  $cmd   .= " $script_file";
  $cmd   .= " none none none" if($flag eq "cron");
  my $out = `$cmd`;                                                                          # sumbit job
  ($out   =~ /^Submitted batch job (\d+)/) or LOGDIE "\tError executing $cmd\n";             # abort if no job ID is returned
  $job_id = $1;
  
  return $job_id;
}


sub filter_overlapping_pos{
  my @list             = @{$_[0]};
  my ($dist,$pid,$log) = ($_[1],$_[2],$_[3]);
  my $count = 0;
  my @uniq  = ();
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  @list = sort { my ($a_id,$a_seq,$a_start,$a_stop,$a_chr,$a_gc,$a_tm,$a_dg) = unpack("Z*Z*iiZ*fff", $a);
                 my ($b_id,$b_seq,$b_start,$b_stop,$b_chr,$b_gc,$b_tm,$b_dg) = unpack("Z*Z*iiZ*fff", $b);
                 $a_chr cmp $b_chr || $a_start <=> $b_start; } @list;
  
  INFO "\tCheck for pid $pid if the coordinates of the ".scalar(@list)." oligonucleotides overlap\n";
  my ($cmp_id,$cmp_seq,$cmp_start,$cmp_stop,$cmp_chr,$cmp_gc,$cmp_tm,$cmp_dg) = unpack("Z*Z*iiZ*fff",$list[0]);                                 # initialise values
  
  for(my $i = 1; $i < scalar(@list); $i++){
    my ($c_id,$c_seq,$c_start,$c_stop,$c_chr,$c_gc,$c_tm,$c_dg) = unpack("Z*Z*iiZ*fff",$list[$i]);    

    if($i == (scalar(@list) - 1)){                                                                                                     # current element is the LAST ELEMENT in the loop: no further comparisons possible
      if($cmp_chr eq $c_chr){                                                                                                          # SAME CHROMOSOME
        if((($c_start - $cmp_stop) - 1) >= $dist){                                                                                     # distance greater than $dist: keep both
          my $str_cmp = pack("Z*Z*iiZ*fff",$cmp_id,$cmp_seq,$cmp_start,$cmp_stop,$cmp_chr,$cmp_gc,$cmp_tm,$cmp_dg);
          my $str_c   = pack("Z*Z*iiZ*fff",$c_id,$c_seq,$c_start,$c_stop,$c_chr,$c_gc,$c_tm,$c_dg);
          push(@uniq,$str_cmp);
          push(@uniq,$str_c);
        }else{                                                                                                                         # overapping coordinates
          if($cmp_gc >= $c_gc){                                                                                                        # gc of the compared value is better than that of the current value: keep compared value
            my $str_cmp = pack("Z*Z*iiZ*fff",$cmp_id,$cmp_seq,$cmp_start,$cmp_stop,$cmp_chr,$cmp_gc,$cmp_tm,$cmp_dg);
            push(@uniq,$str_cmp);
          }else{                                                                                                                       # gc of the current value is better: keep current value
            my $str_c   = pack("Z*Z*iiZ*fff",$c_id,$c_seq,$c_start,$c_stop,$c_chr,$c_gc,$c_tm,$c_dg);
            push(@uniq,$str_c);
          }
        }
      }else{                                                                                                                           # DIFFERENT CHROMOSOMES
        my $str_cmp = pack("Z*Z*iiZ*fff",$cmp_id,$cmp_seq,$cmp_start,$cmp_stop,$cmp_chr,$cmp_gc,$cmp_tm,$cmp_dg);
        my $str_c   = pack("Z*Z*iiZ*fff",$c_id,$c_seq,$c_start,$c_stop,$c_chr,$c_gc,$c_tm,$c_dg);
        push(@uniq,$str_cmp);                                                                                                          # save the values of the compared object
        push(@uniq,$str_c);                                                                                                            # save the values of the current element too, since they lie on different chromosomes
      }
    }else{                                                                                                                             # NOT LAST
      if($cmp_chr eq $c_chr){                                                                                                          # SAME CHROMOSOME
        if((($c_start - $cmp_stop) - 1) >= $dist){                                                                                     # distance between the two oligos is greater than $dist: keep the compared value as unique and set $cmp on the current element
          my $str_cmp = pack("Z*Z*iiZ*fff",$cmp_id,$cmp_seq,$cmp_start,$cmp_stop,$cmp_chr,$cmp_gc,$cmp_tm,$cmp_dg);
          push(@uniq,$str_cmp);
          ($cmp_id,$cmp_seq,$cmp_start,$cmp_stop,$cmp_chr,$cmp_gc,$cmp_tm,$cmp_dg) = ($c_id,$c_seq,$c_start,$c_stop,$c_chr,$c_gc,$c_tm,$c_dg);
        }else{                                                                                                                         # overlapping coordinates
          if($cmp_gc >= $c_gc){                                                                                                        # gc of the compared value is better than that of the current value: go to next element, ignore current
            next;            
          }else{                                                                                                                       # current gc better: set $cmp to current element
            ($cmp_id,$cmp_seq,$cmp_start,$cmp_stop,$cmp_chr,$cmp_gc,$cmp_tm,$cmp_dg) = ($c_id,$c_seq,$c_start,$c_stop,$c_chr,$c_gc,$c_tm,$c_dg);
          }
        }
      }else{                                                                                                                           # DIFFERENT CHROMOSOMES
        my $str_cmp = pack("Z*Z*iiZ*fff",$cmp_id,$cmp_seq,$cmp_start,$cmp_stop,$cmp_chr,$cmp_gc,$cmp_tm,$cmp_dg);
        push(@uniq,$str_cmp);                                                                                                          # save the value of the compared object: it is the last on that chromosome
        ($cmp_id,$cmp_seq,$cmp_start,$cmp_stop,$cmp_chr,$cmp_gc,$cmp_tm,$cmp_dg) = ($c_id,$c_seq,$c_start,$c_stop,$c_chr,$c_gc,$c_tm,$c_dg);         # set values of the current object as values of the compared object
      }
    }
  }  
  INFO "\tPID:$pid -> ".scalar(@uniq)." oligonucleotides without overlapping coordinates.\n";

  return @uniq;
}


sub process_vmatch_output{
  my ($vmout,$kmer,$log) = @_;
  my $count              = 0;
  my %uniq;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  open IN,"<",$vmout or LOGDIE "\tCan't read the Vmatch output: $vmout!\n";
  while(<IN>){
    chomp($_);
    next if(/^$/);                         # drop empty lines
    next if(/^#/);                         # drop description line
    my ($chr,$start,$name) = split(/\t/,$_);
    $start += 1;                           # vmatch coordinates are 0-based
    my $stop = ($start + $kmer) - 1;
    my $str  = pack("Z*ii",$chr,$start,$stop);
    $uniq{$name} = $str;
  }
  close IN;
   
  INFO "\t\t".keys(%uniq)." oligonucleotides mapped uniq on the reference chromosome for: $vmout\n";
  
  return %uniq;
}


sub read_dir{
  my ($path,$log) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  my @list = ();
  my $file = undef;
  
  opendir(DIR,$path) or LOGDIE "\tCan't open directory $path: $!\n";
  while(defined($file = readdir(DIR))){
    next if $file =~ /^\.\.?$/;
    my $in = $path."/".$file;
    push(@list,$in);
  }
  closedir(DIR);
  
  return @list;
}


sub get_reference{
  my ($path,$log) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  my @list = ();
  my $file = undef;
  
  opendir(DIR,$path) or LOGDIE "\tCan't open directory $path: $!\n";
  while(defined($file = readdir(DIR))){
    next if $file =~ /^\.\.?$/;
    next if $file !~ /\.fa$|\.fasta$|\.fas$/;
    push(@list,$file);
  }
  closedir(DIR);
  
  return @list;
}


sub write_DBresult{
  my ($result,$output,$log) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  INFO "\tStart writing the result of the database enquiry\n";
  
  open my $bh,">",$output or LOGDIE "\tCannot write the file containing the oligonucleotides in the requested region: $output $!\n";
  print $bh "chromosom\tstart\tstop\tk-mer\ttm\tGC\tfree energy (dG)\n";  
  foreach my $row (@{$result}){
    my ($valChr,$valStart,$valStop,$valSeq,$valTm,$valGc,$valDg) = @{$row};
    
    print $bh "$valChr\t$valStart\t$valStop\t$valSeq\t$valTm\t$valGc\t$valDg\n";
  }
  close $bh;
  
  INFO "\tNumber of unique oligonucleotides found in this region: ", 0 + @{$result},"\n";
  INFO "\tThe selected oligonucleotides have been written to the file: $output\n";
}


sub check{
  my ($i,$dir) = @_;
  
  my $flag = 0;
  my $file = $dir."/oligos_pid_".$i.".tab";
  $flag    = 1 if(-e $file);                       # file exists: set flag to 1
  
  return $flag;
}


sub createSubmissionScript_vmatchIdx{
  my ($wgs,$tmp_dir,$log,$project,$mail) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  my $dir_idx = $tmp_dir."/vmatch_idx";
  capturex("mkdir", ($dir_idx)) if(!(-d $dir_idx));

  INFO "\tGenerate the Vmatch-DB index for: $wgs\n";
  my $ref_idx = "hg19_wgs";
  my $symb    = $dir_idx."/".$ref_idx;
  my $idx_sh  = $dir_idx."/vmatch_indexTree.sh";
  my $ref     = basename($wgs);
  
  open my $fh,'>', $idx_sh or LOGDIE "\tCannot write the shell script to indices the ref. genome $wgs: file $idx_sh\n";
    print $fh "#!/bin/bash -l\n";
    print $fh "#SBATCH -A $project\n";
    print $fh "#SBATCH -p core\n";
    print $fh "#SBATCH -n 5\n";
    print $fh "#SBATCH -t 04:00:00\n";
    print $fh "#SBATCH --open-mode=append\n";
    print $fh "#SBATCH -o $log\n";
    print $fh "#SBATCH -e $log\n";
    print $fh "#SBATCH --mail-type=FAIL\n";
    print $fh "#SBATCH --mail-user=$mail\n\n";  
    print $fh "IN=$wgs\n";
    print $fh "rsync \$IN \$SNIC_TMP\n\n";
    print $fh "DB=\$SNIC_TMP/$ref\n";
    print $fh "OUT=\$SNIC_TMP/$ref_idx\n\n";    
    print $fh "mkvtree -db \$DB -indexname \$OUT -dna -allout -pl\n\n";
    print $fh "rsync \$SNIC_TMP/$ref_idx* $dir_idx\n";
  close $fh;
  capturex("chmod","a+wx", ($idx_sh));
 
  return ($idx_sh,$symb);
}


sub createSubmissionScript_vmParser{
  my ($vm_out,$vm_in,$pid,$arguments,$tm,$work,$script,$tmpdir,$log,$logdir,$out_dir,$ref_count,$mail,$idx) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  my $name    = $tmpdir."/my_VmParserFiltering_".$pid.".sh";
  my @tmp     = split(/::/,$arguments);
  my $account = $tmp[1];
  
  INFO "\tCreate the Vmatch parser and filtering submission script for pid $pid!\n";
  open my $fh,">",$name or LOGDIE "\tCannot write the shell script for the Vmatch parser $name (pid $pid)\n";
    print $fh "#!/bin/bash -l\n";
    print $fh "#SBATCH -A $account\n";
    print $fh "#SBATCH -p core\n";
    print $fh "#SBATCH -n 1\n";
    print $fh "#SBATCH -t 100:00:00\n";
    print $fh "#SBATCH -J vm1stParserFiltering\n";
    print $fh "#SBATCH --mail-type=FAIL\n";
    print $fh "#SBATCH --mail-user=$mail\n";
    print $fh "#SBATCH --open-mode=append\n";
    print $fh "#SBATCH -o $log\n";
    print $fh "#SBATCH -e $log\n";
    print $fh "JPT=\$1\n";
    print $fh "NAME=\$2\n";
    print $fh "VM_OUT=$vm_out\n";
    print $fh "VM_IN=$vm_in\n";
    print $fh "ARG=\'$arguments\'\n\n";
    print $fh "perl $script/vmatchParser.pl \$VM_OUT \$VM_IN \$ARG $pid $tm $work $script $tmpdir $log $logdir $out_dir $ref_count $mail \$SNIC_TMP $idx \$JPT \$NAME\n";    
  close $fh;
  capturex("chmod","a+wx", ($name));
  
  return $name;
}


sub createSubmissionScript_vmParser_2ndRun{
  my ($vm_out,$vm_in,$pid,$arguments,$log,$tmpdir,$workdir,$script,$mail,$ref_count,$idx,$outdir,$logdir) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  my $name    = $tmpdir."/my_2ndVmParserFiltering_".$pid.".sh";
  my @tmp     = split(/::/,$arguments);
  my $project = $tmp[1];

  INFO "\tCreate the 2nd Vmatch parser and filtering submission script for pid $pid!\n";
  open my $fh,">",$name or LOGDIE "\tCannot write the shell script for the 2nd Vmatch parser $name (pid $pid)\n";
    print $fh "#!/bin/bash -l\n";
    print $fh "#SBATCH -A $project\n";
    print $fh "#SBATCH -p core\n";
    print $fh "#SBATCH -n 3\n";
    print $fh "#SBATCH -t 25:00:00\n";
    print $fh "#SBATCH -J vm2ndParserFiltering\n";
    print $fh "#SBATCH --mail-type=FAIL\n";
    print $fh "#SBATCH --mail-user=$mail\n";
    print $fh "#SBATCH --open-mode=append\n";
    print $fh "#SBATCH -o $log\n";
    print $fh "#SBATCH -e $log\n";
    print $fh "JPT=\$1\n";
    print $fh "NAME=\$2\n";
    print $fh "VM_OUT=$vm_out\n";
    print $fh "VM_IN=$vm_in\n";
    print $fh "ARG=\'$arguments\'\n\n";
    print $fh "perl $script/vmatchParser_2ndRun.pl \$VM_OUT \$VM_IN $pid \$ARG $log $tmpdir $workdir $script $mail $ref_count $idx  $outdir $logdir \$SNIC_TMP \$JPT \$NAME\n";
  close $fh;
  capturex("chmod","a+wx", ($name));

  return $name;
}


sub createSubmissionScript_jellyfishParser{
  my ($infile,$ref,$pid,$project,$arguments,$work_dir,$script,$tmp_dir,$log_dir,$pid_log,$dbRes_dir,$vmIndex,$nr_ref,$mail,$max_overlap) = @_;
  
  Log::Log4perl->easy_init({level => $INFO, file => ">> $pid_log"});
  my $name = $tmp_dir."/my_jellyfishParser_".$pid.".sh";
  
  INFO "\tParse Jellyfish output for pid $pid -> $infile";
  open my $fh,">",$name or LOGDIE "\tCannot write the shell script for the jellyfish parser $name\n";
    print $fh "#!/bin/bash -l\n";
    print $fh "#SBATCH -A $project\n";
    print $fh "#SBATCH -p core\n";
    print $fh "#SBATCH -n 1\n";
    print $fh "#SBATCH -t 10:00:00\n";
    print $fh "#SBATCH -J jellyfishParser_$pid\n";
    print $fh "#SBATCH --open-mode=append\n";
    print $fh "#SBATCH -o $pid_log\n";
    print $fh "#SBATCH -e $pid_log\n";
    print $fh "#SBATCH --mail-type=FAIL\n";
    print $fh "#SBATCH --mail-user=$mail\n\n";
    print $fh "ARG=\'$arguments\'\n\n";
    print $fh "perl $script/jellyfishParser.pl $infile $ref $pid \$ARG $work_dir $script $tmp_dir $log_dir $pid_log $dbRes_dir $vmIndex $nr_ref $mail $max_overlap \$SNIC_TMP\n";
  close $fh;
  capturex("chmod","a+wx", ($name));
 
  return $name;
}


sub createSubmissionScript_jellyfish{
  my ($indir,$project,$kmer,$tmp_dir,$log,$mail) = @_;  
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  my %tmp_out;
  
  # create the output directories and output files
  my $name    = $tmp_dir."/my_jellyfish.sh";
  my $dir_out = $tmp_dir."/jellyfish";
  my @list    = get_reference($indir,$log);
  my $nr_ref  = scalar(@list);
  capturex("mkdir",($dir_out)) if(!(-d $dir_out));
  
  open my $fh,">",$name or LOGDIE "\tCannot write the shell script for the jellyfish invoice $name\n";
    print $fh "#!/bin/bash -l\n";
    print $fh "#SBATCH -A $project\n";
    print $fh "#SBATCH -p core\n";
    print $fh "#SBATCH -n 10\n";
    print $fh "#SBATCH -t 03:30:00\n";
    print $fh "#SBATCH -J jellyfish\n";
    print $fh "#SBATCH --open-mode=append\n";
    print $fh "#SBATCH -o $log\n";
    print $fh "#SBATCH -e $log\n";
    print $fh "#SBATCH --mail-type=FAIL\n";
    print $fh "#SBATCH --mail-user=$mail\n\n";
    print $fh "module load bioinfo-tools jellyfish &>/dev/null\n\n";
    
    foreach my $dat (@list){
      my $tout       = "merCount_".$1."jf" if $dat =~ /(.*)\.fa/;                       # binary file
      my $output     = "uniqueOligos_kmer".$kmer."_".$dat;                              # fasta file     
      my $ref        = $indir."/".$dat;                                                 # output files
      my $fin        = $dir_out."/".$output;
      $tmp_out{$fin} = $ref;
      
      print $fh "\n";
      print $fh "jellyfish count -o \$SNIC_TMP/$tout -m $kmer -s 2G -t 10 -U 1 $ref\n";
      print $fh "jellyfish dump \$SNIC_TMP/$tout -o \$SNIC_TMP/$output\n\n";
    }    
    print $fh "scp \$SNIC_TMP/uniqueOligos_kmer* $dir_out\n";
    print $fh "module unload bioinfo-tools jellyfish &>/dev/null\n";
  close $fh;
  capturex("chmod","a+wx", ($name));

  return ($name,$nr_ref,\%tmp_out);
}


sub prepare_vmatch_1stRun{
  my ($tmpdir,$pid,$kmer,$idx,$vm_inp,$project,$log,$logdir,$mail,$max_over) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  # vmatch sh-scripts and output directories
  my $sh_dir = $tmpdir."/vm_1stRun/vmatch_sh_".$pid;
  my $vm_out = $tmpdir."/vm_1stRun/vmatch_out_".$pid;
  capturex("mkdir", "-p", ($vm_out)) if(!(-d $vm_out));
  capturex("mkdir", "-p", ($sh_dir)) if(!(-d $sh_dir));

  # create the shell scripts to run vmatch in parallel. Each submission script should run 1000 vmatch jobs
  INFO "\tCreate Vmatch submission scripts for process id $pid\n";
  my @list = read_dir($vm_inp,$log);
  LOGDIE "\tThe Vmatch directory does not contain any batch files for process id $pid. Check if any files were written out.\n" if(!@list);

  my $vm2process = 1000;                                                   # process 1000 batch files at one time
  my $id         = 0;
  my $idxn       = basename($idx);
  
  # create the vmatch job itself
  my $run = $sh_dir."/run_vmatch_".$pid.".sh";
  open my $hw,">",$run or LOGDIE "\tcan't write the file $run\n";
    print $hw "#!/bin/bash -l\n\n";
    print $hw "file=\$(basename \$1)\n";
    print $hw "pid=\$(basename \$1 .fa \| cut -d \"_\" -f 2,3)\n";    
    print $hw "dir=$vm_out\n";
    print $hw "out=\"coordinates_uniqueOligos_\"\$pid\".tab\"\n";
    print $hw "KMER=$max_over\n\n";
    print $hw "# check if the result of this sample is already created\n";
    print $hw "if [ -e \$dir/\$out ]\n";
    print $hw "  then exit 0\n";
    print $hw "fi\n\n";
    print $hw "# copy sample to node\n";
    print $hw "rsync \$1 \$SNIC_TMP\n\n";
    print $hw "# run vmatch\n";
    print $hw "vmatch -q \$SNIC_TMP/\$file -l \$KMER -p -d -selfun mydb.so -showdesc 100 \$SNIC_TMP/$idxn > \$SNIC_TMP/\$out\n\n";
    print $hw "# copy the results back to network storage\n";
    print $hw "rsync \$SNIC_TMP/\$out \$dir/\$out\.incomplete\n\n";
    print $hw "# rename file once the copying is done\n";
    print $hw "mv \$dir/\$out\.incomplete \$dir/\$out\n";
  close $hw;
  capturex("chmod","a+wx", ($run));
 
  # create the wrapper which manages 1000 batches in one job, running 16 at one time
  for(my $i = 0; $i < scalar(@list); $i+=$vm2process){
    my $tmp_ls = $sh_dir."/seqList_".$pid."_".$id.".txt";                   # write a tmp text file with the links to the files
    
    open my $wh,">",$tmp_ls or LOGDIE "\tcan't write the file $tmp_ls\n";
    for(my $j = $i; $j < ($i + $vm2process) and $j < scalar(@list); $j++){
      print $wh $list[$j]."\n";
    }
    close $wh;

    # create the shell script and a dummy file to control the submission
    my $sh    = $sh_dir."/vmatch_".$pid."_".$id.".sh";
    my $dummy = $sh.".toSubmit";
    my $jname = "vm1st_".$pid."_".$id;
    my $dir   = dirname($idx);
    capturex("touch",($dummy));
    
    open my $vm, ">", $sh or LOGDIE "\tCannot write the shell script for the 1st Vmatch run for batch file $sh (pid $pid)!\n";
      print $vm "#!/bin/bash -l\n";
      print $vm "#SBATCH -A $project\n";
      print $vm "#SBATCH -p node\n";
      print $vm "#SBATCH -n 1\n";
      print $vm "#SBATCH -t 240:00:00\n";
      print $vm "#SBATCH -J $jname\n";
      print $vm "#SBATCH --open-mode=append\n";
      print $vm "#SBATCH -o $log\n";
      print $vm "#SBATCH -e $log\n";
      print $vm "#SBATCH --mail-type=FAIL\n";
      print $vm "#SBATCH --mail-user=$mail\n\n";
      print $vm "REF=$idx\n\n";
      print $vm "rsync -r $dir/ \$SNIC_TMP\n";
      print $vm "wait\n\n";
      print $vm "threads=\$SLURM_CPUS_ON_NODE\n";
      print $vm "cat $tmp_ls \| xargs -i --max-procs=\$threads bash -c \"bash $run {}\"\n";
    close $vm;
    capturex("chmod","a+wx", ($sh));   
    $id++;
  }

  INFO $id." jobs submitted for the first vmatch run\n";
  return ($sh_dir,$vm_out);          # path to all VMATCH shell scripts and output directories
}


sub createSubmissionScript_filtering{
  my @list  = @{$_[0]};
  my ($log,$pid,$tmp,$mail,$proj,$arg,$tm_min,$tm_max,$script) = ($_[1],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7],$_[8],$_[9]);
  my $count = 0;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  INFO "\tCheck the melting temperature, calculate the minimal free energy, and filter out those with stable sec. structure\n";
  my $dir   = $tmp."/filter_".$pid;
  my $shdir = $dir."/sh";
  my $res   = $dir."/results";
  capturex("mkdir", ($dir))   if(!(-d $dir));
  capturex("mkdir", ($shdir)) if(!(-d $shdir));
  capturex("mkdir", ($res))   if(!(-d $res));

  # set the number of cores for each shell script
  my $core  = 16;
  $core     = scalar(@list) if(scalar(@list) <= $core);
  
  for(my $x = 0; $x < scalar(@list); $x+=$core){
    my $sh_filter = $shdir."/filter_".$pid."_".$count.".sh";
    my $dummy     = $sh_filter.".toSubmit";
    my @commands  = ();
    my @copies    = ();
    $core         = (scalar(@list) - $x) if(($x + $core) > scalar(@list));                      # adjust the core number if less than 16 jobs remaining
  
    # create dummy file for the shell script which should be submitted
    capturex("touch",($dummy));

    open my $fh,">",$sh_filter or LOGDIE "\tCannot write the shell script for the filtering step $sh_filter\n";
      print $fh "#!/bin/bash -l\n";
      print $fh "#SBATCH -A $proj\n";
      print $fh "#SBATCH -p core\n";
      print $fh "#SBATCH -n $core\n";
      print $fh "#SBATCH -t 50:00:00\n";
      print $fh "#SBATCH -J filtering\n";
      print $fh "#SBATCH --open-mode=append\n";
      print $fh "#SBATCH -o $log\n";
      print $fh "#SBATCH -e $log\n";
      print $fh "#SBATCH --mail-type=FAIL\n";
      print $fh "#SBATCH --mail-user=$mail\n\n";
      print $fh "ARG=\'$arg\'\n\n";
    
      # get the individual filter calls
      for(my $y = $x; $y < ($x+$core) and $y < scalar(@list); $y++){
        my @tmp   = unpack("Z*Z*Z*",$list[$y]);
        my $fa    = basename($tmp[0]);
        my $cord  = basename($tmp[1]);
        my $out   = $tmp[2];
      
        my $bef = "rsync ".$tmp[0]." ".$tmp[1]." \$SNIC_TMP &";
        push(@copies,$bef);

        my $str = "perl $script/filteringParser.pl \$SNIC_TMP/$fa \$SNIC_TMP/$cord \$SNIC_TMP/$out $tm_min $tm_max \$ARG $script $log \$SNIC_TMP $res &";
        push(@commands,$str);
      }
    
      # write out the commands to copy the files to snic_tmp
      foreach my $lem (@copies){
        print $fh $lem."\n";
      }
      print $fh "wait\n\n";

      # write out the commands to run the filtering step
      foreach my $row (@commands){
        print $fh $row."\n";
      }
      print $fh "wait\n\n";
    close $fh;
    capturex("chmod","a+wx", ($sh_filter)); 

    $count++;
  }
  
  return ($shdir,$dir,$res);
}
1;
__END__

