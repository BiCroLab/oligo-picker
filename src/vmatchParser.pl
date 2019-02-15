#!/usr/bin/perl -w
use strict;
use Log::Log4perl qw(:easy);
use IPC::System::Simple qw(capturex);
use File::Basename;

my ($vm_out,$vm_in,$arguments,$pid,$tm,$workdir,$script,$tmpdir,$log,$logdir,$outdir,$ref_count,$mail,$snictmp) = @ARGV;

require("$script/FileHandling.pm");
Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

my @arg      = split(/::/,$arguments);
my $project  = $arg[1];
my $tm_delta = $arg[7];

my $tm_min   = $tm - $tm_delta;                                                                                    # calculate tm interval
my $tm_max   = $tm + $tm_delta;

# collect the vm_in (fasta) and vm_out (tab) files
my @seq_files   = FileHandling::read_dir($vm_in,$log);
my @coord_files = FileHandling::read_dir($vm_out,$log);
LOGDIE "\tThe number of Vmatch input fasta files does not correspond with the number of Vmatch output files\n" if(scalar(@seq_files) != scalar(@coord_files));

my @list          = FileHandling::combine_seq_coord(\@seq_files,\@coord_files);                                                                                      # combine the fasta files with the coord files
my ($g,$dir,$res) = FileHandling::createSubmissionScript_filtering(\@list,$log,$pid,$tmpdir,$mail,$project,$arguments,$tm_min,$tm_max,$script);                      # generate the submission of the filtering scripts
my @collect       = @{$g};

# submit the filtering jobs
my $jobid = "";
foreach my $shfile (@collect){
  my $jid = FileHandling::submit_job($shfile,0,0,$log);
  
  if($jobid eq ""){
    $jobid = $jid;
  }else{
    $jobid = $jobid.":".$jid;
  }
}

my $overlap_sh = FileHandling::createSubmissionScript_checkAndCombine($pid,$vm_in,$vm_out,$dir,$ref_count,$mail,$log,$outdir,$tmpdir,$script,$arguments,$res,$project);        # create submission script to check overlappings and load the data to the DB and submit the job after the previous ones are finished
my $oid        = FileHandling::submit_job($overlap_sh,1,$jobid,$log);
