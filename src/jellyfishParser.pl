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
my $homopoly  = $tmp_arg[12];
my $user      = $tmp_arg[13];

my ($vm_in_dir,$tm) = FileHandling::process_jellyfish_output($infile,$pid,$tmp_dir,$gcmin,$gcmax,$kmer,$salt,$formamid,$homopoly,$log,$snictmp);                                              # process the jellyfish output: set unique IDs, filter GC

# prepare the first vmatch run using the following options: vmatch -p -d -l $max_over -showdesc 100
my ($vm_sh,$vm_out) = FileHandling::prepare_vmatch_1stRun($tmp_dir,$pid,$kmer,$idx,$vm_in_dir,$proj,$log,$log_dir,$mail,$max_over);                                                           # path to vmatch shell scripts, vmatch index, vmatch output and input, maximal overlap
my $sh_vmParser     = FileHandling::createSubmissionScript_vmParser($vm_out,$vm_in_dir,$pid,$arguments,$tm,$work_dir,$script,$tmp_dir,$log,$log_dir,$out_dir,$ref_count,$mail,$idx);   # submission script for vmatch output processing submit the vmatch index job (slurm script, 0:not chained/1:chained, job id: 0 = no chaining, >0 if chained)

# prepare and submit the cronjob to submit batches of vmatch calls and the vmatch parser after all other jobs finished
my $cronjob         = FileHandling::createCronjob($vm_sh,$log,$tmp_dir,$pid,$proj,$user,$mail,$sh_vmParser,"vm1st");                                                                                # 0 = flag to show that there are no jobids yet
my $id              = FileHandling::submit_job($cronjob,0,0,$log,"cron");

