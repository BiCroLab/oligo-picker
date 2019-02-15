package FileHandling;

use strict;
use warnings;
use IPC::System::Simple qw(capturex);
use Log::Log4perl qw(:easy);
use Bio::SeqIO;
use File::Basename;
require Exporter;

our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw(&createSubmissionScript_checkAndCombine &createSubmissionScript_filtering &combine_seq_coord &submit_job &filter_overlapping_pos &calculate_filter_tm &get_gc &get_dg &log10 &write_DBresult &process_jellyfish_output &prepare_vmatch_run &createSubmissionScript_jellyfish &createSubmissionScript_jellyfishParser &read_dir &createSubmissionScript_vmatchParser &process_vmatch_output &check &write_tmp_file &mergeChromosomes &createSubmissionScript_vmatchIdx);


#############################
# Preloaded methods go here #
#############################
sub combine_seq_coord{
  my @seq  = @{$_[0]};
  my @coor = @{$_[1]};
  my @list = ();
  
  for(my $i = 0; $i < scalar(@seq); $i++){
    my $id = $1 if($seq[$i] =~ /batch_(\d+_\d+)\.fa/);
  
    for(my $j = 0; $j < scalar(@coor);$j++){
      my $cid = $1 if($coor[$j] =~ /coordinates_uniqueOligos_(\d+_\d+)\.tab/);

      if($cid eq $id){
        my $out = "filtered_uniqueOligos_".$id.".tab";
        my $str = pack("Z*Z*Z*",$seq[$i],$coor[$j],$out);
        push(@list,$str);
      }
    }
  }
  
  return @list;
}


sub process_jellyfish_output{
  my ($in,$ref,$pid,$tmpdir,$gcmin,$gcmax,$kmer,$salt,$formamid,$homopoly,$log,$snictmp) = @_;
  my ($loc,$nuc_count,$bad_gc,$bp,$id,$sum,$poly) = (0,0,0,0,0,0,0,0);
  my ($i,$seq_nr)                                 = (1,1);
  my $batch_size                                  = 40000000;                                       # size of a chunck: 40 Mb
  my %collect                                     = ();
  
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  # create vmatch input directory
  my $dir = $tmpdir."/vmatch_in_".$pid;
  capturex("mkdir", ($dir)) if(!(-d $dir));
  
  # copy the jellyfish output to SNIC_TMP
  my $tmpin = $snictmp."/jellyfishOut_".$pid.".fa";
  capturex("rsync", ($in,$tmpin));
  
  # check and prepare for homopolymer filtering
  if($homopoly ne "na"){
    $homopoly = $1 if($homopoly =~ /(.*);$/);                                    # remove ";" if it exists at the end of the string
    $homopoly =~ tr/;/|/ if($homopoly =~ /;/);                                   # replace ";" with "|"
    
    INFO "\tProcess the jellyfish output for pid $pid: filter homopolymers, calculate tm and GC, filter the oligos without a correct GC, and split the file in 40 Mb chunks.\n";
  }else{
    INFO "\tProcess the jellyfish output for pid $pid: calculate tm and GC, filter the oligos without a correct GC value, and split the file in 40 Mb chunks.\n";
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
    
    my $tm = 81.5 + (0.41 * $gc) + (16.6 * log10($salt)) - (500 / $kmer) - (0.62 * $formamid);                    # calculate tm temperature
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

  INFO "\t$nuc_count unique oligos found. $id sequences fulfiled the GC condition ($bad_gc removed). $i batch files generated.\n";
  INFO "\t$poly sequences removed due to homopolymers filtering\n" if($homopoly ne "na");
  INFO "\tAverage melting temperature: $tm_avg\n";
  
  INFO "\tCopy the batch files to $dir";
  foreach my $elm (keys %collect){                      # copy the batch files from snictmp to tmp
    capturex("rsync", ($elm,$dir));
  }
  
  # remove the jellyfish input file, since it is no longer necessary
  INFO "\tRemove the jellyfish original output, since it is no longer requested: $in\n";
  capturex("rm", ($in));
  
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
          my @elm = unpack("Z*iiZ*fffi",$collect{$tmp[1]});
          $elm[7] += 1;
          my $str = pack("Z*iiZ*fffi",$elm[0],$elm[1],$elm[2],$elm[3],$elm[4],$elm[5],$elm[6],$elm[7]);
          $collect{$tmp[1]} = $str;
        }else{
          my $x   = 1;
          my $str = pack("Z*iiZ*fffi",$tmp[0],$tmp[2],$tmp[3],$tmp[4],$tmp[5],$tmp[6],$tmp[7],$x);
          $collect{$tmp[1]} = $str;
        }
      }
    close(GU);
  }
  
  # remove redundant oligonucleotides
  my @uniq;
  
  foreach my $oligo (keys %collect){
    my @tmp = unpack("Z*iiZ*fffi",$collect{$oligo});
    if($tmp[7] == 1){
      my $str = pack("Z*Z*iiZ*fff",$tmp[0],$oligo,$tmp[1],$tmp[2],$tmp[3],$tmp[4],$tmp[5],$tmp[6]);
      push(@uniq,$str);
    }
  }
  
  INFO "\t".scalar(@uniq)." oligos (out of $count) are unique in the whole human genome\n";
  
  return @uniq;
}


sub submit_job{
  my ($script_file,$is_chain,$job_id,$log) = @_;                                             # SLURM script to be run; 0: not chained, 1: chained; job ID variable: 0 if no chaining, > 0 if chained; log file
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
              
  # construct call to sbatch
  my $cmd = 'sbatch';
  $cmd   .= " --dependency=afterOK:$job_id" if($job_id !~ /^0$/ && $is_chain == 1);          # if chained add dependency
  $cmd   .= " $script_file";  
  my $out = `$cmd`;                                                                          # sumbit job
  ($out   =~ /^Submitted batch job (\d+)/) or LOGDIE "\tError executing $cmd\n";             # abort if no job ID is returned
  $job_id = $1;
  
  INFO "\tJob submitted to the queue: id=$job_id, name=$script_file\n\n";
  
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


sub log10{
  my $n = shift;
  return log($n)/log(10);
}


sub process_vmatch_output{
  my ($vmout,$log) = @_;
  my $count        = 0;
  my (%list,%uniq);
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  open IN,"<",$vmout or LOGDIE "\tCan't read the Vmatch output: $vmout!\n";
  while(<IN>){
    chomp($_);
    next if(/^$/);                         # drop empty lines
    next if(/^#/);                         # drop description line
    my $line = $_;
    $line    =~ s/^\s+|\s+$//g;            # remove white space at the beginning and at the end of the line
    my @tmp  = split(/\s+/,$line);
    my $start = $tmp[2] + 1;  
    my $stop = ($start + $tmp[4]) - 1;
    my $chr  = $1 if($tmp[1] =~ /^(MT|Y|X|\d+)_.*/);
    my $name = $1 if($tmp[5] =~ /(candidate_\d+_\d+)_gc/);
      
    if(exists($list{$name})){# multiple occurrance in the reference
      my ($t_beg,$t_end,$t_chr,$flag) = unpack("iiZ*Z*",$list{$name});
      $flag   = "no";                                                       # oligonucleotide is not unique: flag = no
      my $str = pack("iiZ*Z*",$t_beg,$t_end,$t_chr,$flag);
      $list{$name} = $str;
    }else{
      my $flag = "yes";
      my $str  = pack("iiZ*Z*",$start,$stop,$chr,$flag);                   # start, stop, chr, flag=yes -> unique oligo
      $list{$name} = $str;
    }
  }
  close IN;
 
  # check if all oligonucleotides have a unique match to the reference sequence
  foreach my $qid (keys %list){
    my @tmp  = unpack("iiZ*Z*",$list{$qid});
    my $flag = $tmp[3];
    
    if($flag eq "no"){# oligo has multiple matches to the reference genome
      $count++;
    }else{# oligo is unique, flag = yes
      my $str = pack("Z*ii",($tmp[2]+1),$tmp[0],$tmp[1]);
      $uniq{$qid} = $str;
    }
  }
  
  my $percent = sprintf("%.1f",($count/keys(%list)*100));
  INFO "\t\t".keys(%list)." oligonucleotides mapped on the reference chromosome: $vmout\n";
  INFO "\t\t$count ($percent%) oligonucleotides mapped to multiple locations on the reference genome (i.e. reverse complemented matches, palindroms)\n";
  
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


sub write_tmp_file{
  my $list = $_[0];
  my $file = $_[1];
  my $log  = $_[2];
  my $pid  = $_[3];
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  my $flag = "yes";
  
  INFO "\tWrite out the temporary file (pid:$pid): $file\n";
  open(TW,">$file") or LOGDIE "\tCannot write the file $file for pid $pid\n";
  foreach my $i (@{$list}){
    my ($id,$seq,$start,$stop,$chr,$gc,$tm,$dg) = unpack("Z*Z*iiZ*fff",$i);
    print TW "$id\t$seq\t$start\t$stop\t$chr\t$gc\t$tm\t$dg\n";
  }
  close(TW);
  
  return $flag;
}


sub createSubmissionScript_vmatchIdx{
  my ($ref_chr,$tmpdir,$log,$pid,$project,$mail) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  
  INFO "\tCreate vmatch index and sh directories";
  my $dir_idx = $tmpdir."/vmatch_idx";
  my $dir_sh  = $tmpdir."/vmatch_sh_".$pid;
  capturex("mkdir", ($dir_idx)) if(!(-d $dir_idx));
  capturex("mkdir", ($dir_sh))  if(!(-d $dir_sh));

  INFO "\tCreate shell script to generate the Vmatch-DB index for $pid: $ref_chr\n";
  my $ref_idx = "hg19_pid".$pid;
  my $symb    = $dir_idx."/hg19_pid".$pid;
  my $idx     = $dir_sh."/vmatch_indexTree.sh";
  my $ref     = basename($ref_chr);
  
  open my $fh,'>', $idx or LOGDIE "\tCannot write the shell script to indices the ref. chromosome $ref_chr, pid $pid: file $idx\n";
    print $fh "#!/bin/bash -l\n";
    print $fh "#SBATCH -A $project\n";
    print $fh "#SBATCH -p core\n";
    print $fh "#SBATCH -n 1\n";
    print $fh "#SBATCH -t 02:00:00\n";
    print $fh "#SBATCH --open-mode=append\n";
    print $fh "#SBATCH -o $log\n";
    print $fh "#SBATCH -e $log\n";
    print $fh "#SBATCH --mail-type=FAIL\n";
    print $fh "#SBATCH --mail-user=$mail\n\n";  
    print $fh "IN=$ref_chr\n";
    print $fh "rsync \$IN \$SNIC_TMP\n\n";
    print $fh "DB=\$SNIC_TMP/$ref\n";
    print $fh "OUT=\$SNIC_TMP/$ref_idx\n\n";    
    print $fh "mkvtree -db \$DB -indexname \$OUT -dna -allout -pl\n\n";
    print $fh "rsync \$SNIC_TMP/$ref_idx* $dir_idx\n";
  close $fh;
  
  return ($idx,$symb);
}


sub createSubmissionScript_vmatchParser{
  my ($vm_out,$vm_in,$pid,$arguments,$tm,$work,$script,$tmpdir,$log,$logdir,$out_dir,$ref_count,$mail) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  my $name    = $tmpdir."/my_vmatchParserFiltering_".$pid.".sh";
  my @tmp     = split(/::/,$arguments);
  my $account = $tmp[1];
  
  INFO "\tCreate the Vmatch parser and filtering submission script for pid $pid!\n";
  open my $fh,">",$name or LOGDIE "\tCannot write the shell script for the Vmatch parser $name (pid $pid)\n";
    print $fh "#!/bin/bash -l\n";
    print $fh "#SBATCH -A $account\n";
    print $fh "#SBATCH -p core\n";
    print $fh "#SBATCH -n 1\n";
    print $fh "#SBATCH -t 05:00:00\n";
    print $fh "#SBATCH -J vmatchParserFiltering\n";
    print $fh "#SBATCH --mail-type=FAIL\n";
    print $fh "#SBATCH --mail-user=$mail\n";
    print $fh "#SBATCH --open-mode=append\n";
    print $fh "#SBATCH -o $log\n";
    print $fh "#SBATCH -e $log\n";
    print $fh "VM_OUT=$vm_out\n";
    print $fh "VM_IN=$vm_in\n";
    print $fh "ARG=\'$arguments\'\n\n";
    print $fh "perl $script/vmatchParser.pl \$VM_OUT \$VM_IN \$ARG $pid $tm $work $script $tmpdir $log $logdir $out_dir $ref_count $mail \$SNIC_TMP\n";
  close $fh;
  
  return $name;
}


sub createSubmissionScript_jellyfishParser{
  my ($infile,$ref,$pid,$account,$arguments,$work,$script,$tmp,$logdir,$log,$out_dir,$idx,$ref_count,$mail) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  my $name = $tmp."/my_jellyfishParser_".$pid.".sh";
  
  INFO "\tCreate the Jellyfish parser submission script for pid $pid -> $infile";
  open my $fh,">",$name or LOGDIE "\tCannot write the shell script for the jellyfish parser $name\n";
    print $fh "#!/bin/bash -l\n";
    print $fh "#SBATCH -A $account\n";
    print $fh "#SBATCH -p core\n";
    print $fh "#SBATCH -n 1\n";
    print $fh "#SBATCH -t 10:00:00\n";
    print $fh "#SBATCH -J jellyfishParser_$pid\n";
    print $fh "#SBATCH --open-mode=append\n";
    print $fh "#SBATCH -o $log\n";
    print $fh "#SBATCH -e $log\n";
    print $fh "#SBATCH --mail-type=FAIL\n";
    print $fh "#SBATCH --mail-user=$mail\n\n";
    print $fh "ARG=\'$arguments\'\n\n";
    print $fh "perl $script/jellyfishParser.pl $infile $ref $pid \$ARG $work $script $tmp $logdir $log $out_dir $idx $ref_count $mail \$SNIC_TMP\n";
  close $fh;
  
  return $name;
}


sub createSubmissionScript_jellyfish{
  my ($indir,$account,$kmer,$tmp,$log,$mail) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});
  my %tmp_out;
  INFO "\tCreate the Jellyfish submission script. Call Jellyfish for each chromosome individually!";  

  # create the output directories and output files
  my $name      = $tmp."/my_jellyfish.sh";
  my $dir_out   = $tmp."/jellyfish";
  my @list      = get_reference($indir,$log);
  my $ref_count = scalar(@list);
  capturex("mkdir",($dir_out)) if(!(-d $dir_out));
  
  open my $fh,">",$name or LOGDIE "\tCannot write the shell script for the jellyfish invoice $name\n";
    print $fh "#!/bin/bash -l\n";
    print $fh "#SBATCH -A $account\n";
    print $fh "#SBATCH -p core\n";
    print $fh "#SBATCH -n 10\n";
    print $fh "#SBATCH -t 02:30:00\n";
    print $fh "#SBATCH -J jellyfish\n";
    print $fh "#SBATCH --open-mode=append\n";
    print $fh "#SBATCH -o $log\n";
    print $fh "#SBATCH -e $log\n";
    print $fh "#SBATCH --mail-type=FAIL\n";
    print $fh "#SBATCH --mail-user=$mail\n\n";
    print $fh "module load bioinfo-tools jellyfish\n\n";
    
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
    print $fh "module unload bioinfo-tools jellyfish\n";
  close $fh;

  return ($name,$ref_count,\%tmp_out);
}


sub prepare_vmatch_run{
  my ($tmpdir,$pid,$kmer,$idx,$vm_inp,$project,$log,$logdir,$mail) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  # vmatch sh-scripts and output directories
  my $sh_dir = $tmpdir."/vmatch_sh_".$pid;
  my $vm_out = $tmpdir."/vmatch_out_".$pid;
  capturex("mkdir", ($vm_out)) if(!(-d $vm_out));
  
  # create the shell scripts to run vmatch in parallel. Each submission script contain 16 vmatch calls
  INFO "\tCreate Vmatch submission scripts for process id $pid\n";
  my $id   = 0;
  my @list = read_dir($vm_inp,$log);
  LOGDIE "\tThe Vmatch directory does not contain any batch files for process id $pid. Check if any files were written out.\n" if(!@list);
  
  # set the number of cores for each shell script
  my $core  = 16;
  $core     = scalar(@list) if(scalar(@list) <= 16);
  
  for(my $i = 0; $i < scalar(@list); $i+=$core){
    my $sh   = $sh_dir."/vmatch_".$pid."_".$id.".sh";
    my @copy = ();
    my @comm = ();
    my @copi = ();
    
    # adjust the core number if less than 16 jobs remaining
    $core = (scalar(@list) - $i) if(($i + $core) > scalar(@list));
  
    open my $vm,">",$sh or LOGDIE "\tCannot write the shell script to run Vmatch for batch file $sh (pid $pid)!\n";
      print $vm "#!/bin/bash -l\n";
      print $vm "#SBATCH -A $project\n";
      print $vm "#SBATCH -p core\n";
      print $vm "#SBATCH -n $core\n";
      print $vm "#SBATCH -t 02:00:00\n";
      print $vm "#SBATCH -J vmatch\n";
      print $vm "#SBATCH --open-mode=append\n";
      print $vm "#SBATCH -o $log\n";
      print $vm "#SBATCH -e $log\n";
      print $vm "#SBATCH --mail-type=FAIL\n";
      print $vm "#SBATCH --mail-user=$mail\n\n";
      print $vm "KMER=$kmer\n";
      print $vm "REF=$idx\n\n";
      
      # copy the reference to the SNIC_TMP
      print $vm "rsync $idx* \$SNIC_TMP\n";
      print $vm "wait\n\n";
      
      # get the individual vmatch calls for the current submission script
      for(my $j   = $i; $j < ($i+$core) and $j < scalar(@list); $j++){
        my $tid   = $1 if($list[$j] =~ /batch_(\d+_\d+)\.fa/);                               # contain pid and unique number
        my $out   ="coordinates_uniqueOligos_".$tid.".tab";
        my $in    = $list[$j];
        my $name  = basename($list[$j]);
        my $idxn  = basename($idx);

        my $moc = "srun -c 1 -n 1 --exclusive rsync ".$list[$j]." \$SNIC_TMP &";
        push(@copy,$moc);

        my $com = "srun -c 1 -n 1 --exclusive vmatch -q \$SNIC_TMP/$name -l \$KMER -p -d -showdesc 100 \$SNIC_TMP/$idxn > \$SNIC_TMP/$out &";
        push(@comm,$com);

        my $str = "srun -c 1 -n 1 --exclusive rsync \$SNIC_TMP/$out $vm_out &";
        push(@copi,$str);
      }
      
      # copy the files to the SNIC_TMP
      foreach my $elm (@copy){
        print $vm $elm."\n";
      }
      print $vm "wait\n\n";
      
      # run vmatch
      foreach my $lem (@comm){
        print $vm $lem."\n";
      }
      print $vm "wait\n\n";
      
      # go through the generated output files and copy them back to the login node
      foreach my $elm (@copi){
        print $vm "$elm\n";
      }
      print $vm "wait\n";
    close $vm;
    $id++;
  }
   
  return ($sh_dir,$vm_out);          # path to all VMATCH shell scripts and output directories
}


sub createSubmissionScript_checkAndCombine{
  my ($pid,$vm_in,$vm_out,$dir,$ref_count,$mail,$log,$outdir,$tmpdir,$script,$arguments,$res,$project) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  my $overlap_sh = $tmpdir."/my_overlapCheck_DBload_".$pid.".sh";
  open OS,">",$overlap_sh or LOGDIE "\tCannot write the shell script for the filtering step $overlap_sh\n";
    print OS "#!/bin/bash -l\n";
    print OS "#SBATCH -A $project\n";
    print OS "#SBATCH -p core\n";
    print OS "#SBATCH -n 5\n";
    print OS "#SBATCH -t 10:00:00\n";
    print OS "#SBATCH -J overlapCheck_DBload\n";
    print OS "#SBATCH --open-mode=append\n";
    print OS "#SBATCH -o $log\n";
    print OS "#SBATCH -e $log\n";
    print OS "#SBATCH --mail-type=ALL\n";
    print OS "#SBATCH --mail-user=$mail\n\n";
    print OS "ARG=\'$arguments\'\n";
    print OS "perl $script/checkAndCombine.pl $res \$ARG $pid $script $tmpdir $outdir $log \$SNIC_TMP $vm_in $vm_out $dir $ref_count\n";
  close(OS);

  return $overlap_sh;
}


sub createSubmissionScript_filtering{
  my @list  = @{$_[0]};
  my $log   = $_[1];
  my $pid   = $_[2];
  my $tmp   = $_[3];
  my $mail  = $_[4];
  my $proj  = $_[5];
  my $arg   = $_[6];
  my $tm_min= $_[7];
  my $tm_max= $_[8];
  my $script= $_[9];
  my $count = 0;
  my @collect;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  my $dir   = $tmp."/filter_".$pid;
  my $shdir = $dir."/sh";
  my $res   = $dir."/results";
  capturex("mkdir", ($dir))   if(!(-d $dir));
  capturex("mkdir", ($shdir)) if(!(-d $shdir));
  capturex("mkdir", ($res))   if(!(-d $res));

  # set the number of cores for each shell script
  my $core  = 16;
  $core     = scalar(@list) if(scalar(@list) <= 16);
  
  for(my $x = 0; $x < scalar(@list); $x+=$core){
    my $sh_filter = $shdir."/filter_".$pid."_".$count.".sh";
    my @commands  = ();
    my @copies    = ();
    $core         = (scalar(@list) - $x) if(($x + $core) > scalar(@list));                      # adjust the core number if less than 16 jobs remaining
  
    open my $fh,">",$sh_filter or LOGDIE "\tCannot write the shell script for the filtering step $sh_filter\n";
      print $fh "#!/bin/bash -l\n";
      print $fh "#SBATCH -A $proj\n";
      print $fh "#SBATCH -p core\n";
      print $fh "#SBATCH -n $core\n";
      print $fh "#SBATCH -t 15:00:00\n";
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
      
        my $bef = "srun -c 1 -n 1 --exclusive rsync ".$tmp[0]." ".$tmp[1]." \$SNIC_TMP &";
        push(@copies,$bef);

        my $str = "srun -c 1 -n 1 --exclusive perl $script/filteringParser.pl \$SNIC_TMP/$fa \$SNIC_TMP/$cord \$SNIC_TMP/$out $tm_min $tm_max \$ARG $script $log \$SNIC_TMP $res &";
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
  
    $count++;
    push(@collect,$sh_filter);
  }
  
  return (\@collect,$dir,$res);
}
1;
__END__
