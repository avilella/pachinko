#!/usr/local/ensembl/bin/perl -w

use strict;
use DBI;
use Getopt::Long;
use Time::HiRes qw(time gettimeofday tv_interval);
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::Production::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::GenomeDB;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::URLFactory;
use Bio::EnsEMBL::Registry;

my $conf_file;
my @speciesList = ();
my @setList = ();
my %hive_params ;
my %dnds_params;
my %sitewise_dnds_params;
my %pachinko_params;

my %compara_conf = ();
$compara_conf{'-port'} = 3306;

my ($help, $host, $user, $pass, $dbname, $port, $compara_conf, $adaptor, $ensembl_genomes);
my ($subset_id, $genome_db_id, $prefix, $fastadir, $verbose);
$verbose = 0;
GetOptions('help'            => \$help,
           'conf=s'          => \$conf_file,
           'ensembl_genomes' => \$ensembl_genomes,
           'v' => \$verbose,
          );

$verbose = 1 if ($verbose);
if ($help) { usage(); }

# ok this is a hack, but I'm going to pretend I've got an object here
# by creating a blessed hash ref and passing it around like an object
# this is to avoid using global variables in functions, and to consolidate
# the globals into a nice '$self' package
my $self = bless {};
my $starttime = time();

Bio::EnsEMBL::Registry->no_version_check(1);
parse_conf($conf_file);

$host   =   $compara_conf{'-host'} unless (defined $host);     $compara_conf{'-host'} = $host  ;
$port   =   $compara_conf{'-port'} unless (defined $port);     $compara_conf{'-port'} = $port  ;
$dbname = $compara_conf{'-dbname'} unless (defined $dbname); $compara_conf{'-dbname'} = $dbname;
$user   =   $compara_conf{'-user'} unless (defined $user);     $compara_conf{'-user'} = $user  ;
$pass   =   $compara_conf{'-pass'} unless (defined $pass);     $compara_conf{'-pass'} = $pass  ;

unless(defined($host) and defined($user) and defined($dbname)) {
  print "\nERROR : must specify host, user, and database to connect to compara\n\n";
  usage();
}

# Create qdirs
throw("Couldnt find qdirs parameter. Check config file\n") if (!defined($pachinko_params{qdirs}));
foreach my $qdir (split(":",$pachinko_params{'qdirs'})) {
  print STDERR "Creating qdirs directory: $qdir\n";
  system ("mkdir -p $qdir");
  my $lustre_cmd = "lfs setstripe $qdir -c -1";
  eval { system ("$lustre_cmd 1>/dev/null 2>/dev/null"); } ;
  throw "directory not available: $qdir [check your config file]\n" if ((not(-d $qdir)));
  push @{$self->{qdirs_array}}, $qdir;
}
throw("Couldnt create qdirs. Check config file\n") if (!defined($self->{qdirs_array}));

# Check the runnables
print STDERR "Checking runnables exist\n";
foreach my $runnable_key (keys %pachinko_params) {
  next unless ($runnable_key =~ /_exe$/);
  my $runnable = $pachinko_params{$runnable_key};
  throw("Couldnt find executable [$runnable]. Check your config file") unless (!-z $runnable);
}

# We create the new database here
my $dsn = "DBI:mysql:host=$host;port=$port";
print STDERR "Connecting to server $dsn $user:$pass\n";
my $dbh = DBI->connect( $dsn, $user, $pass, {RaiseError => 1})
  or die "Could not connect to database host : " . DBI->errstr;
print STDERR "Attempting to drop old database $dbname\n";
eval { $dbh->do("DROP DATABASE IF EXISTS $dbname"); };
print STDERR "Creating database $dbname\n";
$dbh->do( "CREATE DATABASE " . $dbname )
  or die "Could not create database $dbname: " . $dbh->errstr;
$dbh->disconnect();
# Load all tables
print STDERR "Loading all sql tables from sqlfile parameters\n";
foreach my $sqlfile_key (keys %pachinko_params) {
  next unless ($sqlfile_key =~ /^sqlfile_/);
  my $sqlfile = $pachinko_params{$sqlfile_key};
  throw("Couldnt find sql file [$sqlfile]. Check your config file") unless (-e $sqlfile);
  my $cmd = "mysql -h$host -P$port -u$user -p$pass $dbname < $sqlfile";
  system($cmd) == 0 or die "Failed to $cmd : $? . Check your config file\n";
}
# Make sure all tables are InnoDB
my $cmd = "mysql -h$host -P$port -u$user -p$pass $dbname -N -e \"show tables\"";
my @tables = split("\n",`$cmd`);
print STDERR "Converting tables to InnoDB\n";
foreach my $table (@tables) {
  my $innodb_cmd = "mysql -h$host -P$port -u$user -p$pass $dbname -N -e \"ALTER TABLE $table engine=InnoDB\"";
  system($innodb_cmd) == 0 or die "Failed to $innodb_cmd: $!";
}

# Connecting to the newly created database
$self->{'comparaDBA'} = new Bio::EnsEMBL::Compara::Production::DBSQL::DBAdaptor(%compara_conf);
$self->{gdba}         = $self->{'comparaDBA'}->get_GenomeDBAdaptor();
$self->{'hiveDBA'}    = new Bio::EnsEMBL::Hive::DBSQL::DBAdaptor(-DBCONN => $self->{'comparaDBA'}->dbc);

if (%hive_params) {
  if (defined($hive_params{'hive_output_dir'})) {
    die("\nERROR!! hive_output_dir doesn't exist, can't configure\n  ", $hive_params{'hive_output_dir'} , "\n")
      if(($hive_params{'hive_output_dir'} ne "") and !(-d $hive_params{'hive_output_dir'}));
    $self->{'comparaDBA'}->get_MetaContainer->delete_key('hive_output_dir');
    $self->{'comparaDBA'}->get_MetaContainer->store_key_value('hive_output_dir', $hive_params{'hive_output_dir'});
  }
  if (defined($hive_params{'name'})) {
    $self->{'comparaDBA'}->get_MetaContainer->delete_key('name');
    $self->{'comparaDBA'}->get_MetaContainer->store_key_value('name', $hive_params{'name'});
  }
}

# Connecting to the reuse_db
throw "Couldnt find reuse_db parameters. Check your conf file\n" if (!defined($pachinko_params{reuse_db}));
my $reuse_url = $pachinko_params{reuse_db};
$self->{comparaDBA_reuse} = Bio::EnsEMBL::Hive::URLFactory->fetch($reuse_url, 'compara');

my $db_version; $self->{comparaDBA_reuse}->dbc->dbname =~ /\_(\d+)$/; $db_version = $1;
Bio::EnsEMBL::Registry->load_registry_from_db
  (-host => $self->{comparaDBA_reuse}->dbc->host,
   -port => $self->{comparaDBA_reuse}->dbc->port,
   -user => $self->{comparaDBA_reuse}->dbc->username,
   -pass => $self->{comparaDBA_reuse}->dbc->password,
   -species_suffix => 'reuse',
   -db_version=> $db_version,
   -verbose => $verbose);

# We load the ncbi_taxonomy tables from the reuse_db to the production db
my $reuse_username = $self->{comparaDBA_reuse}->dbc->username;
my $reuse_password = $self->{comparaDBA_reuse}->dbc->password;
my $reuse_pass = ''; $reuse_pass = "-p$reuse_password " if ($reuse_password);
my $reuse_host = $self->{comparaDBA_reuse}->dbc->host;
my $reuse_port = $self->{comparaDBA_reuse}->dbc->port;
my $reuse_dbname = $self->{comparaDBA_reuse}->dbc->dbname;
my $mysqldumpcmd = "mysqldump --compress --lock_tables=FALSE --add-drop-table=TRUE ";
my @ncbi_tables = ('ncbi_taxa_name','ncbi_taxa_node');
for my $tbl_name (@ncbi_tables) {
  printf(STDERR "  %1.3f secs\n", (time()-$starttime));
  print STDERR "Importing taxonomy information: $tbl_name (this may take a while)\n";
  $cmd = "$mysqldumpcmd -u $reuse_username $reuse_pass -h $reuse_host -P$reuse_port $reuse_dbname $tbl_name";
  $cmd .= " | perl -pe 's/MyISAM/InnoDB/' | mysql -u $user -p$pass -h $host -P$port $dbname";
  my $ret = system($cmd);
  if (0 != $ret) {
    throw("Error importing $tbl_name: $ret\n");
  }
}
my $analysisDBA = $self->{'hiveDBA'}->get_AnalysisAdaptor;
my $submitGenome = $analysisDBA->fetch_by_logic_name('SubmitGenome');
unless($submitGenome) {
  $submitGenome = Bio::EnsEMBL::Analysis->new
    (-db_version      => '1',
     -logic_name      => 'SubmitGenome',
     -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy');
  $analysisDBA->store($submitGenome);
  my $stats = $submitGenome->stats;
  $stats->batch_size(100);
  $stats->hive_capacity(-1);
  $stats->update();
  $self->{submitGenome_analysis} = $submitGenome;
}

foreach my $genome (@{$self->{comparaDBA_reuse}->get_GenomeDBAdaptor->fetch_all}) {
  my $this_genome_name = $genome->name;
  my $core_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($this_genome_name.'reuse', 'core');
  unless (defined($core_dba)) {
    warn("Couldnt connect to [$this_genome_name]. Skipping...\n");
    next;
  }
  my $new_locator = "Bio::EnsEMBL::DBSQL::DBAdaptor/user=" . 
    $self->{comparaDBA_reuse}->dbc->username . ";species=" . $genome->name . ";port=" .
      $self->{comparaDBA_reuse}->dbc->port . ";disconnect_when_inactive=1;dbname=" 
        . $core_dba->dbc->dbname . ";host=" . $self->{comparaDBA_reuse}->dbc->host;
  my $new_genome = Bio::EnsEMBL::Compara::GenomeDB->new();
  $new_genome->taxon_id($genome->taxon_id);
  $new_genome->name($genome->name);
  $new_genome->assembly($genome->assembly);
  my $new_genebuild = $genome->genebuild; if (0 >= length($new_genebuild)) {$new_genebuild = 'none'};
  $new_genome->genebuild($new_genebuild);
  $new_genome->locator($new_locator);
  $new_genome->dbID($genome->dbID);

  next if ($new_genome->name =~ /ancestral seq/i);

  print STDERR "Loading " . $new_genome->dbID . " " . $new_genome->name . "\n";
  # Manually inserting values
  my $sql = qq(
           INSERT into genome_db (genome_db_id,name,assembly,taxon_id,assembly_default,genebuild,locator)
           VALUES (?,?,?,?,?,?,?)
        );
  my $gdb_sth = $self->{'comparaDBA'}->dbc->prepare($sql);
  eval { $gdb_sth->execute
           ($new_genome->dbID,
            $new_genome->name,
            $new_genome->assembly,
            $new_genome->taxon_id,
            $new_genome->assembly_default,
            $new_genome->genebuild,
            $new_genome->locator
           ); };
  if ($@) {
    # Trying with the store method
    $self->{'comparaDBA'}->get_GenomeDBAdaptor->store($new_genome);
  }

  my $genomeHash = {};
  $genomeHash->{'gdb'} = $genome->dbID;
  my $input_id = encode_hash($genomeHash);
  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob 
      (-input_id       => $input_id,
       -analysis       => $submitGenome,
       -input_job_id   => 0);
}

# FIXME: new non-reused genomes would go here
# foreach my $speciesPtr (@speciesList) {
#   $self->submitGenome($speciesPtr);
# }

foreach my $setPtr (@setList) {
  $self->submitSet($setPtr);
}


# build the pachinko analyses
$self->build_PachinkoSystem();

exit(0);


#######################
#
# subroutines
#
#######################

sub usage {
  print "pachinko_load.pl [options]\n";
  print "  -help                  : print this help\n";
  print "  -conf <path>           : config file describing compara, templates, and external genome databases\n";
  print "Execute like this:\n";
  print "pachinko_load.pl -conf /home/username/pachinko_dbs/db_genome1/p.conf\n";

  exit(1);
}


sub parse_conf {
  my($conf_file) = shift;

  if($conf_file and (-e $conf_file)) {
    #read configuration file from disk
    my @conf_list = @{do $conf_file};

    foreach my $confPtr (@conf_list) {
      my $type = $confPtr->{TYPE};
      delete $confPtr->{TYPE};
      print("HANDLE type $type\n") if($verbose);
      if($type eq 'COMPARA') {
        %compara_conf = %{$confPtr};
      }
      elsif($type eq 'dNdS') {
        %dnds_params = %{$confPtr};
      }
      elsif($type eq 'sitewise_dNdS') {
        %sitewise_dnds_params = %{$confPtr};
      }
      elsif($type eq 'PACHINKO') {
        %pachinko_params = %{$confPtr};
      }
      elsif($type eq 'SPECIES') {
        push @speciesList, $confPtr;
      }
      elsif($type eq 'SET') {
        push @setList, $confPtr;
      }
      elsif($type eq 'HIVE') {
        %hive_params = %{$confPtr};
      }
    }
  }
}


sub submitGenome
{
  my $self     = shift;
  my $species  = shift;  #hash reference

  print("SubmitGenome for ".$species->{abrev}."\n") if($verbose);

  #
  # connect to external genome database
  #
  my $genomeDBA = undef;
  my $locator = $species->{dblocator};

  unless($locator) {
    print("  dblocator not specified, building one\n")  if($verbose);
    $locator = $species->{module}."/host=".$species->{host};
    $species->{port}   && ($locator .= ";port=".$species->{port});
    $species->{user}   && ($locator .= ";user=".$species->{user});
    $species->{pass}   && ($locator .= ";pass=".$species->{pass});
    $species->{dbname} && ($locator .= ";dbname=".$species->{dbname});
    $species->{species} && ($locator .= ";species=".$species->{species});
    $species->{species_id} && ($locator .=";species_id=".$species->{species_id});
  }
  $locator .= ";disconnect_when_inactive=1";
  print("    locator = $locator\n")  if($verbose);

  eval {
    $genomeDBA = Bio::EnsEMBL::DBLoader->new($locator);
  };

  unless($genomeDBA) {
    print("ERROR: unable to connect to genome database $locator\n\n");
    return;
  }

  my $meta = $genomeDBA->get_MetaContainer;
  my $taxon_id = $meta->get_taxonomy_id;

#   # If we are in E_G then we need to look for a taxon in meta by 'NAME.species.taxonomy_id'
#   if($ensembl_genomes) {
#     if(!defined $taxon_id or $taxon_id == 1) {
#       # We make the same call as in the MetaContainer code, but with the NAME appendage
#       my $key = $species->{eg_name}.'.'.'species.taxonomy_id';
#       my $arrRef = $meta->list_value_by_key($key);
#       if( @$arrRef ) {
#         $taxon_id = $arrRef->[0];
#         print "Found taxonid ${taxon_id}\n" if $verbose;
#       }
#       else {
#         warning("Please insert meta_key '${key}' in meta table at core db.\n");
#       }
#     }
#   }

  my $ncbi_taxon = $self->{'comparaDBA'}->get_NCBITaxonAdaptor->fetch_node_by_taxon_id($taxon_id);
  throw("Couldnt find taxonomy entry for $taxon_id. Check meta container\n") if (!defined($ncbi_taxon));
  my $genome_name = $self->_get_name($species, $meta, $ncbi_taxon);

  my ($cs) = @{$genomeDBA->get_CoordSystemAdaptor->fetch_all()};
  my $assembly = $cs->version;
  my $genebuild = ($meta->get_genebuild or "");

  if($species->{taxon_id} && ($taxon_id ne $species->{taxon_id})) {
    throw("$genome_name taxon_id=$taxon_id not as expected ". $species->{taxon_id});
  }

  my $genome = Bio::EnsEMBL::Compara::GenomeDB->new();
  $genome->taxon_id($taxon_id);
  $genome->name($genome_name);
  $genome->assembly($assembly);
  $genome->genebuild($genebuild);
  $genome->locator($locator);
  $genome->dbID($species->{'genome_db_id'}) if(defined($species->{'genome_db_id'}));

 if($verbose) {
    print("  about to store genomeDB\n");
    print("    taxon_id = '".$genome->taxon_id."'\n");
    print("    name = '".$genome->name."'\n");
    print("    assembly = '".$genome->assembly."'\n");
		print("    genebuild = '".$genome->genebuild."'\n");
    print("    genome_db id=".$genome->dbID."\n");
  }

  $self->{'comparaDBA'}->get_GenomeDBAdaptor->store($genome);
  $species->{'genome_db'} = $genome;
  print "  ", $genome->name, " STORED as genome_db id = ", $genome->dbID, "\n";

  #
  # now fill table genome_db_extra
  #
  eval {
    my ($sth, $sql);
    $sth = $self->{'comparaDBA'}->dbc->prepare("SELECT genome_db_id FROM genome_db_extn
        WHERE genome_db_id = ".$genome->dbID);
    $sth->execute;
    my $dbID = $sth->fetchrow_array();
    $sth->finish();

    if($dbID) {
      $sql = "UPDATE genome_db_extn SET " .
                "phylum='" . $species->{phylum}."'".
                ",locator='".$locator."'".
                " WHERE genome_db_id=". $genome->dbID;
    }
    else {
      $sql = "INSERT INTO genome_db_extn SET " .
                " genome_db_id=". $genome->dbID.
                ",phylum='" . $species->{phylum}."'".
                ",locator='".$locator."'";
    }
    print("$sql\n") if($verbose);
    $sth = $self->{'comparaDBA'}->dbc->prepare( $sql );
    $sth->execute();
    $sth->finish();
    print("done SQL\n") if($verbose);
  };

  #
  # now configure the input_id_analysis table with the genome_db_id
  #
  my $analysisDBA = $self->{'hiveDBA'}->get_AnalysisAdaptor;
  my $submitGenome = $analysisDBA->fetch_by_logic_name('SubmitGenome');

  unless($submitGenome) {
    $submitGenome = Bio::EnsEMBL::Analysis->new(
        -db_version      => '1',
        -logic_name      => 'SubmitGenome',
        -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy'
      );
    $analysisDBA->store($submitGenome);
    my $stats = $submitGenome->stats;
    $stats->batch_size(100);
    $stats->hive_capacity(-1);
    $stats->update();
  }

  my $genomeHash = {};
  $genomeHash->{'gdb'} = $genome->dbID;
  if(defined($species->{'pseudo_stableID_prefix'})) {
    $genomeHash->{'pseudo_stableID_prefix'} = $species->{'pseudo_stableID_prefix'};
  }
  my $input_id = encode_hash($genomeHash);

  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob (
        -input_id       => $input_id,
        -analysis       => $submitGenome,
        -input_job_id   => 0
        );
}

sub submitSet
{
  my $self     = shift;
  my $set  = shift;  #hash reference

  print STDERR "Loading Sets:\n";
  print("SubmitSet for ".$set->{set_id}."\n") if($verbose);

  #
  # connect to external genome database
  #
  my $genomeDBA = undef;
  my $locator = undef;
  my $url = $set->{url};
  my $set_id = $set->{set_id};

#   # If we are in E_G then we need to look for a taxon in meta by 'NAME.species.taxonomy_id'
#   if($ensembl_genomes) {
#     if(!defined $taxon_id or $taxon_id == 1) {
#       # We make the same call as in the MetaContainer code, but with the NAME appendage
#       my $key = $set->{eg_name}.'.'.'species.taxonomy_id';
#       my $arrRef = $meta->list_value_by_key($key);
#       if( @$arrRef ) {
#         $taxon_id = $arrRef->[0];
#         print "Found taxonid ${taxon_id}\n" if $verbose;
#       }
#       else {
#         warning("Please insert meta_key '${key}' in meta table at core db.\n");
#       }
#     }
#   }


  my $taxon_id = $set->{taxon_id};
  my $ncbi_taxon = $self->{'comparaDBA'}->get_NCBITaxonAdaptor->fetch_node_by_taxon_id($taxon_id);
  my $genome_name;
  eval { $genome_name = $ncbi_taxon->name;};
  if ($@) {
    throw ("Couldnt find $taxon_id in ncbi_taxonomy tables. Check config file\n");
  }

  if($set->{taxon_id} && ($taxon_id ne $set->{taxon_id})) {
    throw("$genome_name taxon_id=$taxon_id not as expected ". $set->{taxon_id});
  }

  my $genome = Bio::EnsEMBL::Compara::GenomeDB->new();
  $genome->taxon_id($taxon_id);
  $genome->name($genome_name);
  $genome->assembly($set_id);
  # $genome->genebuild("$taxon_id.$set_id"); # Leave undef
  my $this_locator;

  # Bio::EnsEMBL::DBSQL::DBAdaptor/host=ens-staging2;port=3306;user=ensro;dbname=tarsius_syrichta_core_57_1c;species=Tarsius syrichta;disconnect_when_inactive=1

  my $format = 'fastq'; # FIXME implement other formats
  $this_locator = "Bio::EnsEMBL::DBSQL::DBAdaptor/";
  $this_locator .= "format=$format;";
  my $readlen = $set->{readlen} || 'short';
  $this_locator .= "readlen=$readlen;";
  $this_locator .= "type=ftp;" if ($url =~ /^ftp/);
  $this_locator .= "type=file;" unless ($url =~ /^ftp/);
  my $this_qdir = shift @{$self->{qdirs_array}}; push @{$self->{qdirs_array}}, $this_qdir;
  $this_locator .= "qdir=$this_qdir;";
  $this_locator .= "url=$url;";
  print STDERR "$genome_name - $this_locator\n";
  $genome->locator($this_locator);
  $genome->dbID($set->{'genome_db_id'}) if(defined($set->{'genome_db_id'}));

  if($verbose) {
    print("  about to store genomeDB\n");
    print("    taxon_id = '".$genome->taxon_id."'\n");
    print("    name = '".$genome->name."'\n");
    print("    assembly = '".$genome->assembly."'\n");
    # print("    genebuild = '".$genome->genebuild."'\n");
    print("    genome_db id=".$genome->dbID."\n") if(defined($set->{'genome_db_id'}));
  }

  eval { $self->{'comparaDBA'}->get_GenomeDBAdaptor->store($genome);};
  if ($@) {
    $self->{'comparaDBA'}->get_GenomeDBAdaptor->store($genome);
  }
  $set->{'genome_db'} = $genome;
  print "  ", $genome->name, " STORED as genome_db id = ", $genome->dbID, "\n";

  #
  # now fill table genome_db_extra
  #
  eval {
    my ($sth, $sql);
    $sth = $self->{'comparaDBA'}->dbc->prepare("SELECT genome_db_id FROM genome_db_extn
        WHERE genome_db_id = ".$genome->dbID);
    $sth->execute;
    my $dbID = $sth->fetchrow_array();
    $sth->finish();

    $locator = '' if (!defined $locator);
    if($dbID) {
      $sql = "UPDATE genome_db_extn SET " .
                "phylum='" . $set->{phylum}."'".
                ",locator='".$locator."'".
                " WHERE genome_db_id=". $genome->dbID;
    }
    else {
      $sql = "INSERT INTO genome_db_extn SET " .
        " genome_db_id=". $genome->dbID.
          ",phylum='" . $set->{taxon_id}."'".
            ",locator='".$locator."'";
    }
    print("$sql\n") if($verbose);
    $sth = $self->{'comparaDBA'}->dbc->prepare( $sql );
    $sth->execute();
    $sth->finish();
    print("done SQL\n") if($verbose);
  };

  #
  # now configure the input_id_analysis table with the genome_db_id
  #
  my $analysisDBA = $self->{'hiveDBA'}->get_AnalysisAdaptor;
  my $submitset = $analysisDBA->fetch_by_logic_name('SubmitSet');

  unless($submitset) {
    $submitset = Bio::EnsEMBL::Analysis->new(
        -db_version      => '1',
        -logic_name      => 'SubmitSet',
        -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy'
      );
    $analysisDBA->store($submitset);
    my $stats = $submitset->stats;
    $stats->batch_size(100);
    $stats->hive_capacity(-1);
    $stats->update();
    $self->{submitSet_analysis} = $submitset;
  }

  my $genomeHash = {};
  $genomeHash->{'gdb'} = $genome->dbID;
  if(defined($set->{'pseudo_stableID_prefix'})) {
    $genomeHash->{'pseudo_stableID_prefix'} = $set->{'pseudo_stableID_prefix'};
  }
  my $input_id = encode_hash($genomeHash);

  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob (
        -input_id       => $input_id,
        -analysis       => $submitset,
        -input_job_id   => 0
        );
}



sub _get_name {
	my ($self, $species, $meta, $ncbi_taxon) = @_;

	my $genome_name;

	#If we have both taxon & species then check rules & assign
	if(defined $ncbi_taxon && defined $species->{species}) {
		print "Taxon & Species given; checking for equality\n" if $verbose;
		if($ncbi_taxon->binomial() ne $species->{species}) {
			if($ensembl_genomes) {
				print "Taxon bionmial & species name differ; using species name as we are in ensembl_genomes mode\n" if $verbose;
				$genome_name = $species->{species};
			}
			else {
				print "Taxon bionmial & species name differ; using taxon binomial as we are in compara mode\n" if $verbose;
				$genome_name = $ncbi_taxon->binomial();
			}
		}
		else {
			print "Taxon binomial & species agree with name\n";
			$genome_name = $ncbi_taxon->binomial();
		}
  }
  #If just had taxon then taxon wins
  elsif(defined $ncbi_taxon) {
  	print "Taxon found but no species name given; using taxon binomal\n" if $verbose;
		$genome_name = $ncbi_taxon->binomial();
  }
  #Otherwise best guess from meta
  else {
  	$verbose && print"  Cannot get binomial from NCBITaxon, try Meta...\n";
  	if (defined($species->{species})) {
      $genome_name = $species->{species};
    }
    else {
    	$genome_name = (defined $meta->get_Species()) ? $meta->get_Species()->binomial() : $species->{species};
    }
  }

  return $genome_name;
}

sub build_PachinkoSystem
{
  #yes this should be done with a config file and a loop, but...
  my $self = shift;
  print STDERR "Loading Pachinko System:\n";
  my $dataflowRuleDBA  = $self->{'hiveDBA'}->get_DataflowRuleAdaptor;
  my $ctrlRuleDBA      = $self->{'hiveDBA'}->get_AnalysisCtrlRuleAdaptor;
  my $analysisStatsDBA = $self->{'hiveDBA'}->get_AnalysisStatsAdaptor;
  my $analysisDBA      = $self->{'hiveDBA'}->get_AnalysisAdaptor;
  my $stats;

#   #
#   # SubmitGenome
#   print STDERR "SubmitGenome...\n";
#   #
#   my $submitgenome = Bio::EnsEMBL::Analysis->new(
#       -db_version      => '1',
#       -logic_name      => 'SubmitGenome',
#       -input_id_type   => 'genome_db_id',
#       -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy'
#     );
#   $analysisDBA->store($submitgenome);
#   $stats = $analysisStatsDBA->fetch_by_analysis_id($submitgenome->dbID);
#   $stats->batch_size(100);
#   $stats->hive_capacity(-1);
#   $stats->update();

#   return $submitgenome unless(defined $pachinko_params{fasta_dir} && defined $pachinko_params{fastq_dir});

#   #
#   # SubmitSet
#   print STDERR "SubmitSet...\n";
#   #
#   my $submitset = Bio::EnsEMBL::Analysis->new(
#       -db_version      => '1',
#       -logic_name      => 'SubmitSet',
#       -input_id_type   => 'set_id',
#       -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy'
#     );
#   $analysisDBA->store($submitset);
#   $stats = $analysisStatsDBA->fetch_by_analysis_id($submitset->dbID);
#   $stats->batch_size(100);
#   $stats->hive_capacity(-1);
#   $stats->update();

#   return $submitset unless(defined $pachinko_params{fasta_dir} && defined $pachinko_params{fastq_dir});

  my $submitgenome = $self->{submitGenome_analysis};
  my $submitset = $self->{submitSet_analysis};

  #
  # modelsearch_template
  #
  # create an unlinked analysis called modelsearch_template
  # it will not have rules so it will never execute
  # used to store module,parameters... to be used as template for
  # the dynamic creation of the analyses

  my $pachinko_parameters_string = "{" . join(", ", map(    "\'$_\'=>\'$pachinko_params{$_}\'", keys %pachinko_params) ) . "}";
  $self->{cpunum} = $pachinko_params{cpunum} || 200;
  my $modelsearch_template = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'Modelsearch_template',
      -parameters      => $pachinko_parameters_string
    );
  my $modelsearch_template_analysis_data_id =
    $self->{'hiveDBA'}->get_AnalysisDataAdaptor->store_if_needed($modelsearch_template->parameters);
  my $parameters = undef;
  if (defined $modelsearch_template_analysis_data_id) {
    $parameters = "{'analysis_data_id'=>'$modelsearch_template_analysis_data_id'}";
    $modelsearch_template->parameters($parameters);
  }
  eval { $analysisDBA->store($modelsearch_template); };

  #
  # PachinkoGeneTreeReuse
  print STDERR "PachinkoGeneTreeReuse...\n";
  #
  my $pachinkogenetreereuse = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'PachinkoGeneTreeReuse',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::PachinkoGeneTreeReuse',
      -parameters      => $modelsearch_template->parameters
    );
  $analysisDBA->store($pachinkogenetreereuse);
  $stats = $analysisStatsDBA->fetch_by_analysis_id($pachinkogenetreereuse->dbID);
  $stats->batch_size(1);
  $stats->hive_capacity(20);
  $stats->update();

  my @reuse_tbl = qw
 (subset_member
  member
  sequence
  sequence_cds
  subset
  protein_tree_member
  protein_tree_node
  protein_tree_tag
  protein_tree_hmmprofile
  super_protein_tree_node
  super_protein_tree_member
  super_protein_tree_tag
  protein_tree_stable_id
  mapping_session
  sitewise_aln);
;
#  sequence_exon_bounded
#  stable_id_history

  foreach my $tbl_name (@reuse_tbl) {
    Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob
        (-input_id       => $tbl_name,
         -analysis       => $pachinkogenetreereuse,);
  }

  #
  # GenomeLoadMembers
  print STDERR "GenomeLoadMembers...\n";
  #
  # Uses GenomeLoadReuseMembers but will not execute run section if it's a reused genome
  my $genomeloadmembers = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'GenomeLoadMembers',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::GenomeLoadReuseMembers',
      -parameters      => $modelsearch_template->parameters
    );
  $analysisDBA->store($genomeloadmembers);
  $stats = $analysisStatsDBA->fetch_by_analysis_id($genomeloadmembers->dbID);
  $stats->batch_size(1);
  $stats->hive_capacity(-1); #unlimited
  $stats->failed_job_tolerance(5); # Some of the biggest clusters can fail and go through the other options
  $stats->update();

  $dataflowRuleDBA->create_rule($submitgenome, $genomeloadmembers);

  #
  # BuildSet
  print STDERR "BuildSet...\n";
  #
  my $buildset_exe = $pachinko_params{'buildset_exe'} || '/nfs/users/nfs_a/avilella/src/maq/maq-0.7.1/scripts/fq_all2std.pl';
  my $buildset = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'BuildSet',
      -input_id_type   => 'genome_db_id',
      -program_file    => $buildset_exe,
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::BuildSet',
      -parameters      => $modelsearch_template->parameters,
    );
  $analysisDBA->store($buildset);
  $stats = $analysisStatsDBA->fetch_by_analysis_id($buildset->dbID);
  $stats->batch_size(1);
  $stats->hive_capacity(25);
  $stats->update();

  $dataflowRuleDBA->create_rule($submitset, $buildset);

  #
  # IndexSet
  print STDERR "IndexSet...\n";
  #
  my $indexset_exe = $pachinko_params{'indexset_exe'} || '/nfs/users/nfs_a/avilella/src/cdbfasta/latest/cdbfasta/cdbfasta';

  my $indexset = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'IndexSet',
      -input_id_type   => 'genome_db_id',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::IndexSet',
      -program_file    => $indexset_exe,
      -parameters      => $modelsearch_template->parameters,
    );
  $analysisDBA->store($indexset);
  $stats = $analysisStatsDBA->fetch_by_analysis_id($indexset->dbID);
  $stats->batch_size(1);
  $stats->hive_capacity(25);
  $stats->update();

  $dataflowRuleDBA->create_rule($buildset, $indexset);

  #
  # CreateStoreSeqCDSJobs
  print STDERR "CreateStoreSeqCDSJobs...\n";
  #
  $parameters = $pachinko_params{'cluster_params'};
  my $createStoreSeqCDSJobs = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'CreateStoreSeqCDSJobs',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::CreateStoreSeqCDSJobs',
      -parameters      => $modelsearch_template->parameters,
#      -parameters      => $parameters
  );
  $analysisDBA->store($createStoreSeqCDSJobs);

  if(defined($self->{'hiveDBA'})) {
    $stats = $analysisStatsDBA->fetch_by_analysis_id($createStoreSeqCDSJobs->dbID);
    $stats->batch_size(1);
    $stats->hive_capacity(-1);
    $stats->status('READY');
    $stats->update();
  }

  #
  # StoreSeqCDS
  print STDERR "StoreSeqCDS...\n";
  #
  my $storeseqcds = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'StoreSeqCDS',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::StoreSeqCDS',
      -parameters      => $modelsearch_template->parameters
    );
  $analysisDBA->store($storeseqcds);
  $stats = $analysisStatsDBA->fetch_by_analysis_id($storeseqcds->dbID);
  $stats->batch_size(1);
  $stats->hive_capacity($self->{cpunum});
  $stats->update();

  #
  # UsearchPrepare
  print STDERR "UsearchPrepare...\n";
  #
  $parameters = '{' . "qdirs=>'" . $pachinko_params{qdirs} . "'}";
  my $usearchprepare = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'UsearchPrepare',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::UsearchPrepare',
      -parameters      => $modelsearch_template->parameters,
#      -parameters      => $parameters
    );
  $analysisDBA->store($usearchprepare);
  $stats = $analysisStatsDBA->fetch_by_analysis_id($usearchprepare->dbID);
  $stats->batch_size(1);
  $stats->hive_capacity($self->{cpunum});
  $stats->update();

  #
  # Usearch
  print STDERR "Usearch...\n";
  #
  my $uclust_exe = $pachinko_params{'uclust_exe'} || '/nfs/users/nfs_a/avilella/src/uclust/latest/uclust';
  my $exonerate_exe = $pachinko_params{'exonerate_exe'} || '/nfs/users/nfs_a/avilella/src/exonerate/latest/exonerate-2.2.0-x86_64/bin/exonerate';
  $parameters = $modelsearch_template->parameters;
  $parameters =~ s/\A{//;
  $parameters =~ s/}\Z//;
  $parameters = '{' . $parameters;
  if($pachinko_params{exonerate_exe}) {
    $parameters .= qq|, exonerate_exe=>'$pachinko_params{exonerate_exe}'|;
  }
  $parameters .=  "}";
  my $usearch = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'Usearch',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::Usearch',
      -program_file    => $uclust_exe,
      -parameters      => $modelsearch_template->parameters,
#      -parameters      => $parameters
    );
  $analysisDBA->store($usearch);
  $stats = $analysisStatsDBA->fetch_by_analysis_id($usearch->dbID);
  $stats->batch_size(1);
  $stats->failed_job_tolerance(5);
  $stats->hive_capacity($self->{cpunum});
  $stats->update();

  #
  # Createmodelsearchjobs
  print STDERR "Createmodelsearchjobs...\n";
  #
  # FIXME
  my $createmodelsearchjobs = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'Createmodelsearchjobs',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::Createmodelsearchjobs',
      -parameters      => $modelsearch_template->parameters,
  );
  $analysisDBA->store($createmodelsearchjobs);

  if(defined($self->{'hiveDBA'})) {
    $stats = $analysisStatsDBA->fetch_by_analysis_id($createmodelsearchjobs->dbID);
    $stats->batch_size(1);
    $stats->hive_capacity(-1);
    $stats->status('READY');
    $stats->update();
  }

  $ctrlRuleDBA->create_rule($submitgenome, $pachinkogenetreereuse);
  $ctrlRuleDBA->create_rule($pachinkogenetreereuse, $genomeloadmembers);
  $ctrlRuleDBA->create_rule($pachinkogenetreereuse, $createmodelsearchjobs);

  #
  # BuildHMMaa
  print STDERR "BuildHMMaa...\n";
  #
  my $buildhmm_program = $pachinko_params{'buildhmm_exe'} || '/software/ensembl/compara/hmmer3/hmmer-3.0/src/hmmbuild';
  my $sreformat_program = $pachinko_params{'sreformat_exe'} || '/usr/local/ensembl/bin/sreformat';
  my $buildhmm_hive_capacity = $hive_params{buildhmm_hive_capacity} || 200;
  my $buildhmm_batch_size = $hive_params{buildhmm_batch_size} || 1;

  my $buildhmmaa = Bio::EnsEMBL::Analysis->new
    (
     -db_version      => '1',
     -logic_name      => 'BuildHMMaa',
     -module          => 'Bio::EnsEMBL::Compara::RunnableDB::BuildHMM',
     -program_file    => $buildhmm_program,
      -parameters      => $modelsearch_template->parameters,
#     -parameters      => "{sreformat => '${sreformat_program}'}"
    );
  $analysisDBA->store($buildhmmaa);
  $stats = $buildhmmaa->stats;
  $stats->batch_size($buildhmm_batch_size);
  $stats->hive_capacity($buildhmm_hive_capacity);
  $stats->failed_job_tolerance(5);
  $stats->status('READY');
  $stats->update();

  #
  # BuildHMMcds
  print STDERR "BuildHMMcds...\n";
  #
  $parameters = '';
  $parameters = "{cdna=>1, sreformat => '${sreformat_program}'}";
  my $buildhmmcds = Bio::EnsEMBL::Analysis->new
    (
     -db_version      => '1',
     -logic_name      => 'BuildHMMcds',
     -module          => 'Bio::EnsEMBL::Compara::RunnableDB::BuildHMM',
     -program_file    => $buildhmm_program,
     -parameters      => $modelsearch_template->parameters,
#     -parameters      => $parameters,
    );
  $analysisDBA->store($buildhmmcds);
  $stats = $buildhmmcds->stats;
  $stats->batch_size($buildhmm_batch_size);
  $stats->hive_capacity($self->{cpunum});
  $stats->failed_job_tolerance(5);
  $stats->status('READY');
  $stats->update();

  #
  # Createspeciestree
  print STDERR "Createspeciestree...\n";
  #

  my $createspeciestree = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'Createspeciestree',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::Createspeciestree',
  );
  $analysisDBA->store($createspeciestree);

  if(defined($self->{'hiveDBA'})) {
    $stats = $analysisStatsDBA->fetch_by_analysis_id($createspeciestree->dbID);
    $stats->batch_size(1);
    $stats->hive_capacity(-1);
    $stats->status('READY');
    $stats->update();
  }

  $ctrlRuleDBA->create_rule($createmodelsearchjobs, $createspeciestree);

  #
  # Modelsearch
  print STDERR "Modelsearch...\n";
  #
  my $hmmsearch_exe = $pachinko_params{'hmmsearch_exe'} || '/software/ensembl/compara/hmmer3/hmmer-3.0/src/hmmsearch';
  my $modelsearch = Bio::EnsEMBL::Analysis->new(
      -logic_name      => 'Modelsearch',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::Modelsearch',
      -program_file    => $hmmsearch_exe,
      -parameters      => $modelsearch_template->parameters,
    );
  $analysisDBA->store($modelsearch);
  $stats = $modelsearch->stats;
  $stats->batch_size(1);
  $stats->hive_capacity($self->{cpunum});
  $stats->failed_job_tolerance(5);
  $stats->status('READY');
  $stats->update();

  $dataflowRuleDBA->create_rule($createmodelsearchjobs, $modelsearch, 1);
  $ctrlRuleDBA->create_rule($createmodelsearchjobs, $modelsearch, 1);
  $dataflowRuleDBA->create_rule($createmodelsearchjobs, $buildhmmcds, 2);
  $dataflowRuleDBA->create_rule($createmodelsearchjobs, $buildhmmaa, 3);
  $ctrlRuleDBA->create_rule($buildhmmcds,$modelsearch);
  $ctrlRuleDBA->create_rule($buildhmmaa,$modelsearch);
  $dataflowRuleDBA->create_rule($buildhmmcds,$modelsearch);
  $dataflowRuleDBA->create_rule($buildhmmaa,$modelsearch);

  #
  # Pachinkoalign
  print STDERR "Pachinkoalign...\n";
  #

  my $pachinkoalign_exe = $pachinko_params{'pachinkoalign_exe'} || '/nfs/users/nfs_a/avilella/src/papaya/latest/papaya/pagan';
  my $pachinkoalign = Bio::EnsEMBL::Analysis->new(
      -logic_name      => 'Pachinkoalign',
      -program_file    => $pachinkoalign_exe,
      -parameters      => $modelsearch_template->parameters,
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::Pachinkoalign'
    );

  $analysisDBA->store($pachinkoalign);
  $stats = $pachinkoalign->stats;
  $stats->batch_size(1);
  $stats->hive_capacity($self->{cpunum});
  $stats->failed_job_tolerance(5);
  $stats->update();

  #
  # TreeBeST
  print STDERR "TreeBeSTMagicMirror...\n";
  #
  my $tree_best_program = $pachinko_params{'treebest'} || '/nfs/users/nfs_a/avilella/src/treesoft/trunk/treebest/treebest';
  my $treebestmagicmirror = Bio::EnsEMBL::Analysis->new(
      -logic_name      => 'TreeBeSTMagicMirror',
      -program_file    => $tree_best_program,
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::TreeBeSTMagicMirror',
      -parameters      => $modelsearch_template->parameters,
    );
  $analysisDBA->store($treebestmagicmirror);
  $stats = $treebestmagicmirror->stats;
  $stats->batch_size(1);
  $stats->failed_job_tolerance(5);
  my $njtree_hive_capacity = $hive_params{'njtree_hive_capacity'};
  $njtree_hive_capacity = 400 unless defined $njtree_hive_capacity;
  $stats->hive_capacity($self->{cpunum});
  $stats->update();

  #
  # Change $dnds_params{'method_link_type'} to {'method_link_types'}
  #
  if( defined( $dnds_params{'method_link_type'} ) ){
    warn('[WARN] dNdS => method_link_type is deprecated. '
         .'Use method_link_types instead');
    $dnds_params{'method_link_types'}
      ||= '['. $dnds_params{'method_link_type'} . ']';
  }

  #
  # OrthoTree
  print STDERR "OrthoTree...\n";
  #
  my $with_options_orthotree = 0;
  my $ortho_params = '';
  if (defined $pachinko_params{'honeycomb_dir'}) {
    $ortho_params = "'honeycomb_dir'=>'".$pachinko_params{'honeycomb_dir'}."'";
    $with_options_orthotree = 1;
  }
  if (defined $dnds_params{'species_sets'}) {
    $ortho_params .= ',species_sets=>' . $dnds_params{'species_sets'};
    if( defined $dnds_params{'method_link_types'} ){
      $ortho_params .= ',method_link_types=>'
          . $dnds_params{'method_link_types'};
    }
    $with_options_orthotree = 1;
  }
  if(defined $pachinko_params{'species_tree_file'}) {
    my $tree_file = $pachinko_params{'species_tree_file'};
    $ortho_params .= ",'species_tree_file'=>'${tree_file}'";
    $with_options_orthotree = 1;
  }

  $ortho_params .= ", 'use_genomedb_id'=>1" if defined $pachinko_params{use_genomedb_id};

  #EDIT Originally created a anon hash which caused problems with OrthoTree when using eval
  if($with_options_orthotree) {
    $parameters =~ s/\A{//;
    $parameters =~ s/}\Z//;
    $parameters = '{' . $parameters . ',' .  $ortho_params . '}'
  }

  my $analysis_data_id = $self->{'hiveDBA'}->get_AnalysisDataAdaptor->store_if_needed($parameters);
  if (defined $analysis_data_id) {
    $parameters = "{'analysis_data_id'=>'$analysis_data_id'}";
  }

  my $orthotree = Bio::EnsEMBL::Analysis->new(
      -logic_name      => 'OrthoTree',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::OrthoTree',
      -parameters      => $modelsearch_template->parameters,
      );
  $analysisDBA->store($orthotree);
  $stats = $orthotree->stats;
  $stats->batch_size(1);
  $stats->failed_job_tolerance(5); 
  my $ortho_tree_hive_capacity = $hive_params{'ortho_tree_hive_capacity'};
  $ortho_tree_hive_capacity = 200 unless defined $ortho_tree_hive_capacity;
  $stats->hive_capacity($self->{cpunum});

  $stats->update();

  #
  # QuickTreeBreak
  #
  print STDERR "QuickTreeBreak...\n";
  #
   $parameters = "{max_gene_count=>".$pachinko_params{'max_gene_count'};
  if($pachinko_params{sreformat_exe}) {
  	$parameters .= qq|, sreformat_exe=>'$pachinko_params{sreformat_exe}'|;
  }
  $parameters .= '}';
  my $quicktreebreak = Bio::EnsEMBL::Analysis->new(
      -logic_name      => 'QuickTreeBreak',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::QuickTreeBreak',
      -parameters      => $modelsearch_template->parameters,
#      -parameters      => $parameters
      );
  $quicktreebreak->program_file($pachinko_params{quicktree}) if $pachinko_params{quicktree};
  $analysisDBA->store($quicktreebreak);
  $stats = $quicktreebreak->stats;
  $stats->batch_size(1);
  my $quicktreebreak_hive_capacity = 1; # Some deletes in OrthoTree can be hard on the mysql server
  $stats->hive_capacity($self->{cpunum});
  $stats->update();


  # turn these two on if you need dnds from the old homology system
  #
  # CreateHomology_dNdSJob
  print STDERR "CreateHomology_dNdSJob...\n";
  #
  my $createhomology_dndsjob = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'CreateHomology_dNdSJob',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::CreateHomology_dNdSJobs',
  );
  $analysisDBA->store($createhomology_dndsjob);

  if(defined($self->{'hiveDBA'})) {
    $stats = $analysisStatsDBA->fetch_by_analysis_id($createhomology_dndsjob->dbID);
    $stats->batch_size(1);
    $stats->hive_capacity(-1);
    $stats->status('BLOCKED');
    $stats->update();
    $ctrlRuleDBA->create_rule($orthotree,$createhomology_dndsjob);
    $ctrlRuleDBA->create_rule($treebestmagicmirror,$createhomology_dndsjob);
#    $ctrlRuleDBA->create_rule($BreakPAFCluster,$createhomology_dndsjob);
  }
  if (defined $dnds_params{'species_sets'}) {
    $self->{'hiveDBA'}->get_AnalysisJobAdaptor->CreateNewJob
        (
         -input_id       => ( '{species_sets=>'
                              . $dnds_params{'species_sets'}
#                               . ',method_link_types=>'
#                               . $dnds_params{'method_link_types'}
                              . '}' ),
         -analysis       => $createhomology_dndsjob,
        );
  }

  #
  # Homology_dNdS
  print STDERR "Homology_dNdS...\n";
  #
  my $homology_dNdS = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'Homology_dNdS',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::Homology_dNdS'
  );
  $self->store_codeml_parameters(\%dnds_params);
  if (defined $dnds_params{'dNdS_analysis_data_id'}) {
    $homology_dNdS->parameters('{dNdS_analysis_data_id=>' . $dnds_params{'dNdS_analysis_data_id'} . '}');
  }
  $analysisDBA->store($homology_dNdS);
  if(defined($self->{'hiveDBA'})) {
    $stats = $analysisStatsDBA->fetch_by_analysis_id($homology_dNdS->dbID);
    $stats->batch_size(1);
    my $homology_dnds_hive_capacity = $hive_params{homology_dnds_hive_capacity};
  	$homology_dnds_hive_capacity = 200 unless defined $homology_dnds_hive_capacity;
  	$stats->hive_capacity($self->{cpunum});
    $stats->failed_job_tolerance(2);
    $stats->status('BLOCKED');
    $stats->update();
    $ctrlRuleDBA->create_rule($createhomology_dndsjob,$homology_dNdS);
  }

  #
  # OtherParalogs
  print STDERR "OtherParalogs...\n";
  #
  my $otherparalogs = Bio::EnsEMBL::Analysis->new(
      -logic_name      => 'OtherParalogs',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::OtherParalogs',
      # -parameters      => $parameters
      );
  $analysisDBA->store($otherparalogs);
  $stats = $otherparalogs->stats;
  $stats->batch_size(1);
  my $otherparalogs_hive_capacity = 50;
  $stats->hive_capacity($self->{cpunum});
  $stats->update();

  #
  # Threshold_on_dS
  print STDERR "Threshold_on_dS...\n";
  #
  my $threshold_on_dS = Bio::EnsEMBL::Analysis->new(
      -db_version      => '1',
      -logic_name      => 'Threshold_on_dS',
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::Threshold_on_dS'
  );
  $analysisDBA->store($threshold_on_dS);

  if(defined($self->{'hiveDBA'})) {
    $stats = $analysisStatsDBA->fetch_by_analysis_id($threshold_on_dS->dbID);
    $stats->batch_size(1);
    $stats->hive_capacity(-1);
    $stats->status('BLOCKED');
    $stats->update();
    $ctrlRuleDBA->create_rule($homology_dNdS,$threshold_on_dS);
  }
  if (defined $dnds_params{'species_sets'}) {
    $self->{'hiveDBA'}->get_AnalysisJobAdaptor->CreateNewJob
        (
         -input_id       => ( '{species_sets=>'
                              . $dnds_params{'species_sets'}
#                               . ',method_link_types=>'
#                               . $dnds_params{'method_link_types'}
                              . '}' ),
         -analysis       => $threshold_on_dS,
        );
  }

  #
  # Sitewise_dNdS
  print STDERR "Sitewise_dNdS...\n";
  #

  if (defined $sitewise_dnds_params{'saturated'}) {
    $parameters = '';
    my $with_options_sitewise_dnds = 0;
    if (defined $pachinko_params{'honeycomb_dir'}) {
      $parameters = "'honeycomb_dir'=>'".$pachinko_params{'honeycomb_dir'}."',";
      $with_options_sitewise_dnds = 1;
    }
    if (defined $sitewise_dnds_params{'saturated'}) {
      $parameters .= "'saturated'=>" . $sitewise_dnds_params{'saturated'};
      $with_options_sitewise_dnds = 1;
    }
    if (defined $sitewise_dnds_params{gblocks}) {
    	$parameters .= q{'gblocks_exe'=>} . $sitewise_dnds_params{'gblocks_exe'};
    	$with_options_sitewise_dnds = 1;
    }
    $parameters = '{' . $parameters .'}' if (1==$with_options_sitewise_dnds);

    my $Sitewise_dNdS = Bio::EnsEMBL::Analysis->new
      (
       -db_version      => '1',
       -logic_name      => 'Sitewise_dNdS',
       -module          => 'Bio::EnsEMBL::Compara::RunnableDB::Sitewise_dNdS',
       -program_file    => $sitewise_dnds_params{'program_file'} || ''
      );

    #If params exceed 254 then use the analysis_data table.
    if(length($parameters) > 254) {
      my $ad_dba =  $self->{'hiveDBA'}->get_AnalysisDataAdaptor();
      my $adi = $ad_dba->store_if_needed($parameters);
      $parameters = "{'analysis_data_id'=>${adi}}";
    }
    $Sitewise_dNdS->parameters($parameters);

    $analysisDBA->store($Sitewise_dNdS);

    if(defined($self->{'hiveDBA'})) {
      $stats = $analysisStatsDBA->fetch_by_analysis_id($Sitewise_dNdS->dbID);
      $stats->batch_size(1);
      $stats->hive_capacity($self->{cpunum});
      $stats->failed_job_tolerance(5); 
      $stats->status('BLOCKED');
      $stats->update();
      $ctrlRuleDBA->create_rule($orthotree,$Sitewise_dNdS);
    }

    # Only start sitewise dnds after the pairwise dnds and threshold on dS has finished
    $ctrlRuleDBA->create_rule($threshold_on_dS,$Sitewise_dNdS);

    # When a Sitewise_dNdS job is saturated, we reincorporate the
    # subtrees in the analysis to rerun them again
    $dataflowRuleDBA->create_rule($orthotree, $Sitewise_dNdS, 1);
    # Saturated jobs create new jobs of the same kind
    $dataflowRuleDBA->create_rule($Sitewise_dNdS, $Sitewise_dNdS, 2);
  }
  #
  # build graph of control and dataflow rules
  #

  $ctrlRuleDBA->create_rule($pachinkogenetreereuse,$createStoreSeqCDSJobs);
  $ctrlRuleDBA->create_rule($createStoreSeqCDSJobs, $storeseqcds);
  $ctrlRuleDBA->create_rule($genomeloadmembers,$usearchprepare);
  $ctrlRuleDBA->create_rule($storeseqcds,$usearchprepare);
  $ctrlRuleDBA->create_rule($submitgenome, $createmodelsearchjobs);
  $ctrlRuleDBA->create_rule($submitset, $createmodelsearchjobs);

  $ctrlRuleDBA->create_rule($buildset, $modelsearch);
  $ctrlRuleDBA->create_rule($indexset, $pachinkoalign);
  $dataflowRuleDBA->create_rule($modelsearch, $pachinkoalign, 2);

  $ctrlRuleDBA->create_rule($usearchprepare, $usearch);
  $dataflowRuleDBA->create_rule($usearchprepare, $usearch, 1);
  $ctrlRuleDBA->create_rule($usearch,     $pachinkoalign);
  $ctrlRuleDBA->create_rule($usearch,     $modelsearch);
  $ctrlRuleDBA->create_rule($modelsearch, $pachinkoalign);
  $dataflowRuleDBA->create_rule($pachinkoalign, $treebestmagicmirror);

  $dataflowRuleDBA->create_rule($treebestmagicmirror, $treebestmagicmirror, 2);
  $dataflowRuleDBA->create_rule($treebestmagicmirror, $quicktreebreak, 3);
  $dataflowRuleDBA->create_rule($treebestmagicmirror, $otherparalogs, 3);

  $dataflowRuleDBA->create_rule($orthotree, $buildhmmaa, 1);
  $dataflowRuleDBA->create_rule($orthotree, $buildhmmcds, 1);
  $dataflowRuleDBA->create_rule($orthotree, $quicktreebreak, 2);
  $dataflowRuleDBA->create_rule($orthotree, $otherparalogs, 2);

  # OtherParalogs are calculated for every QuickTreeBreak, but only
  # after all clusters are analysed, to avoid descriptions clashing
  # between OrthoTree and OtherParalogs.
  $ctrlRuleDBA->create_rule($orthotree,   $otherparalogs);
  $ctrlRuleDBA->create_rule($treebestmagicmirror,$otherparalogs);

  $dataflowRuleDBA->create_rule($treebestmagicmirror, $orthotree, 1);

  #
  # create initial job
  print STDERR "create initial job...\n";
  #

  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob
    (
     -input_id       => 1,
     -analysis       => $createmodelsearchjobs,
    );

  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob
    (
     -input_id       => 1,
     -analysis       => $modelsearch,
    );


  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob
    (
     -input_id       => 1,
     -analysis       => $createspeciestree,
    );

  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob
    (
     -input_id       => 1,
     -analysis       => $createStoreSeqCDSJobs,
    );

  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob
    (
     -input_id       => 1,
     -analysis       => $storeseqcds,
    );

  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob
    (
     -input_id       => 1,
     -analysis       => $usearchprepare,
    );

  # at least one so that it doesnt block further analysis
  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob
      (
       -input_id       => 1,
       -analysis       => $buildhmmaa,
      );


  # at least one so that it doesnt block further analysis
  Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor->CreateNewJob
    (
     -input_id       => 1,
     -analysis       => $buildhmmcds,
    );

  print STDERR "Done.\n";
  return 1;
}

sub store_codeml_parameters
{
  my $self = shift;
  my $dNdS_Conf = shift;

  my $options_hash_ref = $dNdS_Conf->{'codeml_parameters'};
  return unless($options_hash_ref);

  my @keys = keys %{$options_hash_ref};
  my $options_string = "{\n";
  foreach my $key (@keys) {
    $options_string .= "'$key'=>'" . $options_hash_ref->{$key} . "',\n";
  }
  $options_string .= "}";

  $dNdS_Conf->{'dNdS_analysis_data_id'} =
         $self->{'hiveDBA'}->get_AnalysisDataAdaptor->store_if_needed($options_string);

  $dNdS_Conf->{'codeml_parameters'} = undef;
}


