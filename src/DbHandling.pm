package DbHandling;

use strict;
use warnings;
use DBI;
use Log::Log4perl qw(:easy);
require Exporter;

our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw() ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw( &enquire_db &connect2Db &loadData2Db &createDB &table_exists);


############################
# Preloaded methods go here#
############################
sub connect2Db{
  my ($driver,$db,$log) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  my $dsn          = "DBI:$driver:dbname=$db";
  my ($user,$pass) = ("","");
  
  my $dbh = DBI->connect($dsn,$user,$pass, {RaiseError=>1, PrintError=>0,ShowErrorStatement=>1});
  
  return $dbh;
}


sub createDB{
  my ($driver,$db,$kmer,$dtm,$gcmin,$gcmax,$hflag,$refname,$overlap,$log) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  my $table = "kmer".$kmer."_dtm".$dtm."_gcmin".$gcmin."_gcmax".$gcmax."_hpol".$hflag."_".$refname."_overlap".$overlap;  
  my $dbh   = connect2Db($driver,$db,$log);
  my $flag  = table_exists($dbh,$table);

  # delete the table if it already exists
  if($flag == 1){
    INFO "\tThe table $table already exists. It's content will be dropped and a new table will be created!\n";
    
    my $sql  = qq(DROP TABLE IF EXISTS $table);
    my $sth  = $dbh->do($sql) or LOGDIE $DBI::errstr; 
  }
  
  # create the table
  my $sqlt = qq(CREATE TABLE $table
               (ID INTEGER PRIMARY KEY  AUTOINCREMENT,
                SEQ            TEXT NOT NULL,
                GC             REAL NOT NULL,
                TM             REAL NOT NULL,
                DG             REAL NOT NULL,
                START          INT  NOT NULL,
                STOP           INT  NOT NULL,
                NAME           CHAR(70),
                CHR            CHAR(50));
             );
  my $stht = $dbh->do($sqlt);
  
  if($stht < 0){
    LOGDIE $DBI::errstr."\n";
  }else{
    INFO "\tThe table $table has been successfully created.\n";
  }
  
  $dbh->disconnect();
}


sub table_exists{
  my ($dbh,$table_name) = @_;
  my @tables            = $dbh->tables('','','','TABLE');
  my $flag              = 0;
  
  if(@tables){# non empty array: test if table is present
    for(@tables){
      next unless $_;
      $flag = 1 if $_ eq $table_name;
    }
  }else{# array is empty (DB is empty or DBMS fails to report). To be on the safe side try a SELECT statement with a false WHERE clause
    # if it returns an error: table does not exist
    eval{
      local $dbh->{PrintError} = 0;
      local $dbh->{RaiseError} = 1;
      $dbh->do(qq{SELECT * FROM $table_name WHERE 1 = 0});
    };
    return 1 unless $@;
  }

  return $flag;
}


sub enquire_db{
  my ($chr,$beg,$end,$kmer,$dtm,$gcmin,$gcmax,$driver,$name,$hflag,$refname,$overlap,$log) = @_;
  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"});

  my $tab = "kmer".$kmer."_dtm".$dtm."_gcmin".$gcmin."_gcmax".$gcmax."_hpol".$hflag."_".$refname."_overlap".$overlap;
  my $dbh = connect2Db($driver,$name,$log);
  my $all;
  
  INFO "\tStart DB enquiry for table $tab, chromosome: $chr, start: $beg, stop: $end, kmer: $kmer, gcmin: $gcmin, gcmax: $gcmax, dTm: $dtm\n";
  # check if the table exists in the database
  if(table_exists($dbh,$tab)){
    my $sql = qq(SELECT CHR,START,STOP,SEQ,TM,GC,DG FROM $tab WHERE CHR = \'$chr\' AND $beg <= START AND START <= $end ORDER BY START);        # get all oligos with the start position in between the requested region (even if the stop is beyond the requested border)
    my $sth = $dbh->prepare($sql) or LOGDIE $DBI::errstr."\n";
    
    $sth->{RaiseError} = 1;                                                                                                                    # any error will cause the DBI module to 'die' with an appropriate message
    $sth->execute();
    $all = $sth->fetchall_arrayref();                                                                                                          # fetch all data from the select statement at one step
 
    $sth->finish();
    $dbh->disconnect();
  }else{
    LOGDIE "\tThe table $tab (kmer size $kmer) does not exist in the database $name. Please check your parameters.\n";
  }
  INFO "\t....Done\n";

  return $all;
}


sub loadData2Db{
  my $driver = $_[0];
  my $db     = $_[1];
  my $kmer   = $_[2];
  my $tdelta = $_[3];
  my $gcmin  = $_[4];
  my $gcmax  = $_[5];
  my @list   = @{$_[6]};
  my $hflag  = $_[7];
  my $refname= $_[8];
  my $overlap= $_[9];
  my $log    = $_[10];

  Log::Log4perl->easy_init({level => $INFO, file => ">> $log"}); 
  my $tab    = "kmer".$kmer."_dtm".$tdelta."_gcmin".$gcmin."_gcmax".$gcmax."_hpol".$hflag."_".$refname."_overlap".$overlap;
  my $dbh    = connect2Db($driver,$db,$log);
  
  INFO "\tStart loading records to the database\n";
  @list = sort {
            my ($a_id,$a_seq,$a_start,$a_stop,$a_chr,$a_gc,$a_tm,$a_dg) = unpack("Z*Z*iiZ*fff", $a);
            my ($b_id,$b_seq,$b_start,$b_stop,$b_chr,$b_gc,$b_tm,$b_dg) = unpack("Z*Z*iiZ*fff", $b);
            $a_chr cmp $b_chr || $a_start <=> $b_start;
          } @list;
  
  #$dbh->begin_work;
  my $insert_count   = 0;
  $dbh->{AutoCommit} = 0;
  $dbh->do("PRAGMA synchronous=OFF");

  my $sql = qq(INSERT INTO $tab (CHR,START,STOP,SEQ,NAME,TM,GC,DG) VALUES (?,?,?,?,?,?,?,?));                          # ? - placeholders
  my $sth = $dbh->prepare($sql);
  
  for(my $i = 0; $i < scalar(@list); $i++){
    $insert_count++;
    my ($id,$seq,$start,$stop,$chr,$gc,$tm,$dg) = unpack("Z*Z*iiZ*fff",$list[$i]);
    $sth->bind_param(1,$chr);
    $sth->bind_param(2,$start);
    $sth->bind_param(3,$stop);
    $sth->bind_param(4,$seq);
    $sth->bind_param(5,$id);
    $sth->bind_param(6,$tm);
    $sth->bind_param(7,$gc);
    $sth->bind_param(8,$dg);
    $sth->execute();
    LOGDIE "\t$DBI::errstr\n" if($dbh->err());
    
    if($insert_count == 10000){
      $dbh->commit;
    }
  }
 
  $dbh->commit;
  $dbh->{AutoCommit} = 1;
  $dbh->disconnect();
}

1;
__END__
