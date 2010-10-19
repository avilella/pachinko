#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::IndexSet

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $indexset = Bio::EnsEMBL::Compara::RunnableDB::IndexSet->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$indexset->fetch_input(); #reads from DB
$indexset->run();
$indexset->output();
$indexset->write_output(); #writes to DB

=cut


=head1 DESCRIPTION

This Analysis will take the sequences from a cluster, the cm from
nc_profile and run a profiled alignment, storing the results as
cigar_lines for each sequence.

=cut


=head1 CONTACT

  Contact Albert Vilella on module implementation/design detail: avilella@ebi.ac.uk
  Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=cut


=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Compara::RunnableDB::IndexSet;

use strict;
use Getopt::Long;
use Time::HiRes qw(time gettimeofday tv_interval);
use Bio::SeqIO;
use Bio::Index::Fasta;
use Bio::Index::Fastq;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use LWP::Simple;

use Bio::EnsEMBL::Hive;
our @ISA = qw(Bio::EnsEMBL::Hive::Process);


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data from the database
    Returns :   none
    Args    :   none

=cut


sub fetch_input {
  my( $self) = @_;

  $self->{'clusterset_id'} = 1;

  #create a Compara::DBAdaptor which shares the same DBI handle
  #with the Pipeline::DBAdaptor that is based into this runnable
  $self->{'comparaDBA'} = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
    (
     -DBCONN=>$self->db->dbc
    );

  # Get the needed adaptors here
  $self->{gdba} = $self->{'comparaDBA'}->get_GenomeDBAdaptor;

  $self->get_params($self->parameters);
  $self->get_params($self->input_id);

# For long parameters, look at analysis_data
  if($self->{analysis_data_id}) {
    my $analysis_data_id = $self->{analysis_data_id};
    my $analysis_data_params = $self->db->get_AnalysisDataAdaptor->fetch_by_dbID($analysis_data_id);
    $self->get_params($analysis_data_params);
  }

  my $genome_db_id = $self->{gdb};

  #get the Compara::GenomeDB object for the genome_db_id
  $self->{genome_db} = $self->{gdba}->fetch_by_dbID($genome_db_id);
  throw("Couldnt find a genome_db for gdb_id $genome_db_id") unless (defined($self->{genome_db}));

  return 1;
}


sub get_params {
  my $self         = shift;
  my $param_string = shift;

  return unless($param_string);
  print("parsing parameter string : ",$param_string,"\n") if($self->debug);

  my $params = eval($param_string);
  return unless($params);

  if($self->debug) {
    foreach my $key (keys %$params) {
      print("  $key : ", $params->{$key}, "\n");
    }
  }

#   foreach my $key (qw[gdb fastq_dir fasta_dir analysis_data_id]) {
  foreach my $key (qw[gdb analysis_data_id]) {
    my $value = $params->{$key};
    $self->{$key} = $value if defined $value;
  }

  return;
}


=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   runs something
    Returns :   none
    Args    :   none

=cut

sub run {
  my $self = shift;

# FIXME: disabled
#  $self->run_indexset;
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   stores something
    Returns :   none
    Args    :   none

=cut


sub write_output {
  my $self = shift;

}


##########################################
#
# internal methods
#
##########################################

sub run_indexset {
  my $self = shift;

  $self->load_set;
  return 1;
}

sub load_set {
  my $self = shift;

  my $locator = $self->{genome_db}->locator;
  throw("Cannot find locator $locator") unless (defined($locator));

  my $locator_hash;
  my $string = $locator;
  $string =~ /(\S+?)\/([\S+\s*]+)/ || die "Could not parse [$string] as a locator. Needs database_module/params";
  my $module = $1;
  my $param = $2;
  my @param = split(/;/,$param);
  foreach my $keyvalue ( @param ) {
    $keyvalue =~ /(\S+?)=([\S*\s*]*)/ 
      || do { warn("In loading $keyvalue, could not split into keyvalue. Ignoring"); next; };
    my $key = $1;
    my $value = $2;
    $locator_hash->{$key} = $value;
  }

  throw("could not find format for $string") unless (defined($locator_hash->{format}));

  my $type = $locator_hash->{type};
  my $url  = $locator_hash->{url};
  throw ("error no $url") unless (defined($url));

  next unless ($type =~ /ftp/); #FIXME implement other types

  my $starttime = time();

  my $fastq_directory = $self->{fastq_dir};
  my $fasta_directory = $self->{fasta_dir};
  my $compression_ext = "gz"; #FIXME other methods

  my $assembly = $self->{genome_db}->assembly;
  throw ("error assembly id for $assembly and $url") unless ($url =~ /$assembly/);

  # Loading locator info
  my $locator = $self->{genome_db}->locator;
  throw("Cannot find locator $locator") unless (defined($locator));

  my $locator_hash;
  my $string = $locator;
  $string =~ /(\S+?)\/([\S+\s*]+)/ || die "Could not parse [$string] as a locator. Needs database_module/params";
  my $module = $1;
  my $param = $2;
  my @param = split(/;/,$param);
  foreach my $keyvalue ( @param ) {
    $keyvalue =~ /(\S+?)=([\S*\s*]*)/ 
      || do { warn("In loading $keyvalue, could not split into keyvalue. Ignoring"); next; };
    my $key = $1;
    my $value = $2;
    $locator_hash->{$key} = $value;
  }

  throw("could not find format for $string") unless (defined($locator_hash->{format}));

  my $fasta_out_filename = $assembly . "." . 'fastq';
  my $qdir  = $locator_hash->{qdir};
  my $fasta_out_file = $qdir . "/" . $fasta_out_filename;
  my $index_out_file = $qdir . "/" . $fasta_out_filename . ".cidx";

  print STDERR "Creating $fasta_out_filename index file...\n" if ($self->debug);

  my $cmd = $self->analysis->program_file || "/nfs/users/nfs_a/avilella/src/cdbfasta/latest/cdbfasta/cdbfasta";
  $cmd   .= " -Q";
  $cmd   .= " $fasta_out_file -o $index_out_file";

  unless(system("cd $fastq_directory; $cmd") == 0) {
    print("$cmd\n") if ($self->debug);
    throw("[$cmd] error indexing $fasta_out_file $!\n");
  }


  return 1;
}

1;
