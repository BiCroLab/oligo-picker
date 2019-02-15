#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Log::Log4perl qw(:easy);
use IPC::System::Simple qw(capturex);
use Math::Round qw(round);

# SET DEFAULT VARIABLES
my $mail        = "gabriele.girelli\@scilifelab.se";
my $driver      = "SQLite";
my $proj        = "b2015233";
my $work_dir    = "/proj/b2015233/nobackup/uniqueOligoPipeline_v0.2";
my $script      = "/proj/b2015233/bils/kmer_pipeline_v0.2_apr2016";
my $wgs         = "/sw/data/uppnex/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa";
my $ref         = "/sw/data/uppnex/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Chromosomes";
my $db_name     = $work_dir."/unique_oligonucleotides_2016.db";
my $formamid    = 50;
my $salt        = 0.42;
my $gc_min      = 35;
my $gc_max      = 80;
my $tmDelta     = 10;
my $dist        = 10;
my $sstmp       = 65;
my $hpol        = "na";
my $max_over_p  = 60;                                       # maximal overlap between two oligomers in percentage

require("$script/DbHandling.pm");
require("$script/FileHandling.pm");

# TIME STEMP
my @months            = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my ($s,$m,$h,$d,$mon) = localtime();
my $time              = $h."h_".$m."m_".$s."s_".$d."".$months[$mon];

# CREATE WORKING DIRECTORIES: RESULTS, LOG, TMP
capturex("mkdir",($work_dir)) if(!(-d $work_dir));
my $out_dir = $work_dir."/db_enquiries";
my $log_dir = $work_dir."/log_".$time;
my $tmp_dir = $work_dir."/tmp_".$time;
capturex("mkdir",($log_dir)) if(!(-d $log_dir));
capturex("mkdir",($tmp_dir)) if(!(-d $tmp_dir));

# INITIATE LOG
my $log = $log_dir."/log_".$time.".log";
Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

# CHECK THE PASSED PARAMETERS
LOGDIE "\tIncorrect number of arguments! Use option -h or --help for usage!\n" if(@ARGV == 0);     # no parameters found
GetOptions( 'p|project=s'     => \$proj,
            'r|reference=s'   => \$ref,
            'W|wgs=s'         => \$wgs,
            'd|db_name=s'     => \$db_name,
            'f|formamid=f'    => \$formamid, 
            'g|gc_min=f'      => \$gc_min,
            'G|gc_max=f'      => \$gc_max,
            't|tm_delta=f'    => \$tmDelta,
            's|salt=f'        => \$salt,
            'S|ss_temp=f'     => \$sstmp,
            'D|distance=i'    => \$dist,
            'w|workdir=s'     => \$work_dir,
            'v|maxoverlap=i'  => \$max_over_p,
            'a|all=s'         => \ my $db_flag,
            'k|kmer=i'        => \ my $kmer,
            'c|chromosom=s'   => \ my $chr,
            'b|begin_pos=i'   => \ my $begin,
            'e|end_pos=i'     => \ my $end,
            'o|homopolymer=s' => \ $hpol,
            'F|hpol_flag=s'   => \ my $hflag,
            'm|mail=s'        => \$mail,
            'h|help'          => \&display_help);

my $flow = -1;
if($db_flag and $kmer){
  LOGDIE "\tWrong argument for the option '-a'! This option accepts only 'yes' or 'no' as arguments. Use option -h or --help for usage!\n" if(lc($db_flag) ne "yes" and lc($db_flag) ne "no");
  $flow  = check_arguments($db_flag);
}else{
  LOGDIE "\n\tThe options '-a' and '-k' are mandatory. Use -h for usage.\n"  if(!$db_flag and !$kmer);
  LOGDIE "\n\tThe option '-a' is mandatory. Use -h for usage.\n"             if(!$db_flag);
  LOGDIE "\n\tThe option '-k' is mandatory. Use -h for usage.\n"             if(!$kmer);
}

# START THE PIPELINE WORKFLOW
main($flow);


###########
# METHODS #
###########
sub main{
  my $option = $_[0];
  my $str    = "";
  
  # calculate the maximal overlap between 2 oligos in base pairs
  my $max_over = round(($kmer * $max_over_p) / 100);
  
  if($option == 2){# option = 2: DB enquiry for a specific region.
    LOGDIE "\tThe file $db_name does not exists\n" if(!(-e $db_name));                                                                    # check if database exists
    capturex("mkdir",($out_dir)) if(!(-d $out_dir));                                                                                      # create output directory

    my $result = DbHandling::enquire_db($chr, $begin, $end, $kmer, $tmDelta, $gc_min, $gc_max, $driver, $db_name, $hflag, $log);          # enquire the database
    my $out    = $out_dir."/requestedRegion_".$chr."_".$begin."_".$end."_kmer".$kmer."_gcMin".$gc_min."_gcMax".$gc_max."_dTm".$tmDelta."_hpol".$hflag.".bed";
    FileHandling::write_DBresult($result,$out,$log);                                                                                      # write the result to a file
    
    capturex("rm",("-r",$tmp_dir)) if(-d $tmp_dir);                                                                                       # remove the tmp directory
    INFO "\tThe temporary directory successfully removed\n";
  }else{
    $str = $driver."::".$proj."::".$kmer."::".$db_name."::".$formamid."::".$gc_min."::".$gc_max."::".$tmDelta."::".$salt."::".$dist."::".$sstmp."::".$hpol if($option == 0);                                         # option = 0: identify all unique oligos, no DB enquiry.
    $str = $driver."::".$proj."::".$kmer."::".$db_name."::".$formamid."::".$gc_min."::".$gc_max."::".$tmDelta."::".$salt."::".$dist."::".$sstmp."::".$hpol."::".$chr."::".$begin."::".$end if($option == 1);         # option = 1: identify all unique oligos + DB enquiry for a specific region.    
    LOGDIE "\tI'm confused. Which part of the pipeline do you want to run? Please check your parameters! Use -h for usage\n" if($option != 0 and $option != 1);

    # JELLYFISH BATCH SCRIPT: IDENTIFY UNIQUE KMERS
    my ($sh_j,$ref_count,$x) = FileHandling::createSubmissionScript_jellyfish($ref,$proj,$kmer,$tmp_dir,$log,$mail);
    my %j_out                = %{$x};
    my $id_j                 = FileHandling::submit_job($sh_j,0,0,$log);                                                                  # slurm script, 0:not chained/1:chained, job id: 0 = no chaining, >0 if chained

    # VMATCH INDEX TREE
    my ($sh_vmidx,$idx)      = FileHandling::createSubmissionScript_vmatchIdx($wgs,$tmp_dir,$log,$proj,$mail);                            # run Vmatch index for the entire genome
    my $id_vIdx              = FileHandling::submit_job($sh_vmidx,0,0,$log);
    my $chain_id             = $id_j.":".$id_vIdx;                                                                                        # start the parser only after Vmatch index and jellyfish have finished

    # JELLYFISH PARSER
    my $pid = 0;
    foreach my $in (keys %j_out){
      my $plog               = $log_dir."/log_".$time."_".$pid.".log";
      my $sh_jParser         = FileHandling::createSubmissionScript_jellyfishParser($in,$j_out{$in},$pid,$proj,$str,$work_dir,$script,$tmp_dir,$log_dir,$plog,$out_dir,$idx,$ref_count,$mail,$max_over);
      my $id_jP              = FileHandling::submit_job($sh_jParser,1,$chain_id,$plog);                                                   # slurm script, 0:not chained/1:chained, job id: 0 = no chaining, >0 if chained    
      $pid++;
    }
  }
}
  

sub check_arguments{
  my $flag   = lc($_[0]);
  my $option = -1;
  
  LOGDIE "\tOption '-k' or '--kmer' requires a positive number as argument\n" if(!$kmer);
  
  if($flag eq "no"){# only database enquiry: kmer, gc_min, gc_max, tm_delta, db name, chromosom, start and end position options are required
    LOGDIE "\tOption '-c' or '--chromosom' requires an argument (e.x. '1', 'MT', or 'X').\n" if(!$chr);
    LOGDIE "\tOption '-b' or '--begin_pos' requires a positive number as argument.\n"        if(!$begin);
    LOGDIE "\tOption '-e' or '--end_pos' requires a positive number as argument.\n"          if(!$end);
    LOGDIE "\tOption '-F' or '--hpol_flag' requires an argument. Please use 'yes' if the table you want to enquire contained homopolymers filtered data, otherwise 'no'.\n"  if(!$hflag);
    LOGDIE "\tThe file $db_name does not exists. Please check your parameters.\n"            if(!(-e $db_name));
    
    $option = 2;                                                                             # database enquiry
    INFO "\tTO DO: enquire \"$db_name\" for unique oligonucleotides using the following parameters: -k $kmer -g $gc_max -G $gc_min -t $tmDelta -c $chr -b $begin -e $end\n";
  }else{
    if(!$chr and !$begin and !$end){
      $option = 0;                                                                           # identify unique kmer and load them to the db without searching for a specific region
      INFO "\tTO DO: identify unique oligonucleotides and insert them into the database. PARAMETERS: -p $proj -r $ref -d $db_name -k $kmer -g $gc_min -G $gc_max -t $tmDelta -D $dist -f $formamid -s $salt -S $sstmp -w $work_dir\n";
      
      if($hpol eq "na"){# no homopolymers filtering
        INFO "\tNo homopolymer filtering activated\n";
      }else{
        $hpol = $1 if($hpol =~ /\"|\'(.*)\"|\'/);
        INFO "\tHomopolymer filtering activated: -o $hpol\n";
      }
    }else{
      LOGDIE "\tThe options -c, -b, and -e are mandatory if you want to enquire the database $db_name. Use -h for usage\n"  if(!$chr and !$begin and !$end);
      LOGDIE "\tThe option '-c'/'--chromosom' requires an argument (e.x. 'chr1' or 'chrX').\n"                              if(!$chr);
      LOGDIE "\tThe option '-b'/'--begin_pos' requires a positive number as argument.\n"                                    if(!$begin);
      LOGDIE "\tThe option '-e'/'--end_pos' requires a positive number as argument.\n"                                      if(!$end);
      
      $option = 1;                                          # identify unique kmers, load them to the db, and search for a specific region
      INFO "\tTO DO: identify unique oligonucleotides, insert them into the database and extract oligos from a specific region.\n\tPARAMETERS: -p $proj -r $ref -d $db_name -k $kmer -g $gc_min -G $gc_max -t $tmDelta -D $dist -f $formamid -s $salt -S $sstmp -c $chr -b $begin -e $end -w $work_dir\n";
      
      if($hpol eq "na"){
        INFO "\tNo homopolymer filtering activated\n";
      }else{
        INFO "\tHomopolymer filtering activated: -o $hpol\n";
      }
    }
  }

  return $option;
}


sub display_help{
  print "\n\tUsage:\n\t\tperl $0 [options]\n";
  print "\tOptions:\n";
  print "\t\t-h, --help\t\tShow this help message and exit.\n";
  print "\t\t-w, --workdir\t\tThe absolute path to the working directory. [default: $work_dir]\n";
  print "\t\t-a, --all\t\tSTR\tThis option accepts only 'yes' or 'no' as arguments (mandatory).\n\t\t\t\t\t\t'yes': generate a new database for a given kmer size,\n\t\t\t\t\t\t'no': enquire existing database for a given kmer size.\n";
  print "\t\t-k, --kmer\t\tINT\tThe length of the unique oligonucleotide (e.x. 40, mandatory).\n";
  print "\t\t-p, --project\t\tSTR\tUppmax project name. [default: b2015233]\n";
  print "\t\t-r, --reference\t\tSTR\tThe path to the individual chromosomes of the genome of interest. [default: /sw/data/uppnex/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Chromosomes]\n";
  print "\t\t-W, --wgs\t\tSTR\tThe path to the whole genome sequence of interest. [default: /sw/data/uppnex/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa]\n";
  print "\t\t-v, --maxoverlap\t\tINT\tThe maximal overlap between 2 oligomers in percentage\n";
  print "\t\t-d, --db_name\t\tSTR\tThe path and name of the SQLite database. [default: unique_oligonucleotides.db]\n";
  print "\t\t-f, --formamid\t\tINT\tThe percentage of the formamid concentration. [default: 50]\n";
  print "\t\t-s, --salt\t\tFLOAT\tThe salt concentration [Na+] in M. [default: 0.42]\n";
  print "\t\t-g, --gc_min\t\tINT\tThe minimum oligonucleotide GC content in percentage. A positive integer between 0 and 100 is expected. [default: 35]\n";
  print "\t\t-G, --gc_max\t\tINT\tThe maximum oligonucleotide GC content in percentage. A positive integer between 0 and 100 is expected. [default: 80]\n";
  print "\t\t-t, --tm_delta\t\tINT\tThe degree of which the average melting temperature should be increased and decreased to form an interval. A positive integer is expected. [default: 10]\n";
  print "\t\t-D, --distance\t\tINT\tThe minimum distance between two unique oligonucleotides. [default: 10]\n";
  print "\t\t-S, --ss_temp\t\tFLOAT\tThe threshold temperature for the secondary structure prediction [default: 65]\n";
  print "\t\t-o, --homopolymer\tSTR\tA list of homopolymers separate by semi-colon, which should not be contained in the oligonucleotide sequences (optional, e.x. \'GGGG;CCCC;AAAAA;TTTTT\'). Please use quotation marks (\'\'). [default: na]\n";
  print "\t\t-c, --chromosom\t\tSTR\tThe chromosome of the region of interest (e.x. '1', 'MT' or 'X').\n";
  print "\t\t-b, --begin_pos\t\tINT\tThe start position of the region of interest.\n";
  print "\t\t-e, --end_pos\t\tINT\tThe stop position of the region of interest.\n";
  print "\t\t-F, --hpol_flag\t\tSTR\tThis option accepts only 'yes' or 'no' as arguments. 'yes' means that homopolymers were filtered in the table you request, 'no' means there was no homopolymer filtering done for the table of interest. It is mandatory if you want to run only the database enquiry\n";
  print "\t\t-m, --mail\t\tSTR\tThe email address to which the SLURM system should report if the jobs fail or complete. Please use quotation marks (\'\') around the email address. [default: 'mihaela.martis\@bils.se']\n";
  exit;
}
