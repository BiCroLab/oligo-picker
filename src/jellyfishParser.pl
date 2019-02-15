#!/usr/bin/perl -w
use strict;
use Log::Log4perl qw(:easy);

my ($infile,$ref,$pid,$arguments,$work_dir,$script,$tmp_dir,$log_dir,$log,$out_dir,$idx,$ref_count,$mail,$max_over,$snictmp) = @ARGV;           #J output, refChr, processID, arguments (12 or 15), workDir, scriptDir, tmpDir, logDir, log, outDir, vmatchIDX, maximal overlap in bp,snictmp

require("$script/FileHandling.pm");
Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

my @tmp_arg   = split(/::/,$arguments);
my $proj      = $tmp_arg[1];
my $kmer      = $tmp_arg[2];
my $gcmin     = $tmp_arg[5];
my $gcmax     = $tmp_arg[6];
my $salt      = $tmp_arg[8];
my $formamid  = $tmp_arg[4];
my $homopoly  = $tmp_arg[11];

my ($vm_in_dir,$tm) = FileHandling::process_jellyfish_output($infile,$pid,$tmp_dir,$gcmin,$gcmax,$kmer,$salt,$formamid,$homopoly,$log,$snictmp);                                              # process the jellyfish output: set unique IDs, filter GC

# prepare the first vmatch run using the following options: vmatch -p -d -l $max_over -showdesc 100
my ($vm_sh,$vm_out) = FileHandling::prepare_vmatch_1stRun($tmp_dir,$pid,$kmer,$idx,$vm_in_dir,$proj,$log,$log_dir,$mail,$max_over);                                                           # path to vmatch shell scripts, vmatch index, vmatch output and input, maximal overlap
my $sh_vmParser     = FileHandling::createSubmissionScript_vmParser($vm_out,$vm_in_dir,$pid,$arguments,$tm,$work_dir,$script,$tmp_dir,$log,$log_dir,$out_dir,$ref_count,$mail,$idx);   # submission script for vmatch output processing submit the vmatch index job (slurm script, 0:not chained/1:chained, job id: 0 = no chaining, >0 if chained)
my @vmsh            = FileHandling::read_dir($vm_sh,$log);                                          # submit the regular vmatch jobs
my $job_ids         = "";

foreach my $file (@vmsh){
  next if($file =~ /vmatch_indexTree.sh/);                                                          # index file already called
  my $id = FileHandling::submit_job($file,0,0,$log);                                                # run only after the index job finished
  
  if($job_ids eq ""){                                                                               # collect all job ids
    $job_ids = $id;   
  }else{
    $job_ids = $job_ids.":".$id;
  }
}

my $vpF_id = FileHandling::submit_job($sh_vmParser,1,$job_ids,$log);                                # submit the vmatch parser -> run only after the vmatch jobs finished
