#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::BuildSet

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $buildset = Bio::EnsEMBL::Compara::RunnableDB::BuildSet->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$buildset->fetch_input(); #reads from DB
$buildset->run();
$buildset->output();
$buildset->write_output(); #writes to DB

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


package Bio::EnsEMBL::Compara::RunnableDB::BuildSet;

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

  foreach my $key (qw[gdb tag indexset_exe cdbyank_exe analysis_data_id]) {
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

  $self->run_buildset;
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

sub run_buildset {
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

#   my ($key,$url) = split("=",$locator);
#   my ($loader,$type) = split("/",$key);
#   my $format = $loader; $format =~ s/.+\:(\S+)/$1/;

  my $type = $locator_hash->{type};
  my $url  = $locator_hash->{url};
  my $qdir  = $locator_hash->{qdir};
  throw ("error no $url") unless (defined($url));

  my $starttime = time();

  my $assembly = $self->{genome_db}->assembly;
  throw ("error assembly id for $assembly and $url") unless ($url =~ /$assembly/);
  my $tmp_filename = $assembly . "." . $locator_hash->{format};
  my $tmp_file = $self->worker_temp_directory . $tmp_filename;
  my $final_file = $qdir . "/" . $tmp_filename;
  my $source_file = $url;
  # $DB::single=1;1;
  my $compression_ext = "gz"; #FIXME other methods
  $compression_ext = "bz2" if ($source_file =~ /\.bz2$/);
  my $compressed_tmp_file = $tmp_file . ".$compression_ext";
  my $cmd;

  $cmd = "rm -f $compressed_tmp_file";
  if (-e $compressed_tmp_file) {
    unless(system("$cmd") == 0) {
      print("$cmd\n") if ($self->debug);
      throw("error deleting previously downloaded file $!\n");
    }
  }
  printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);

  ## This is not needed since cdbfasta can overwrite previous files
  #   $cmd = "rm -f $final_file";
  #   if (-e $final_file) {
  #     unless(system("$cmd") == 0) {
  #       print("$cmd\n") if ($self->debug);
  #       throw("error deleting previously downloaded file $!\n");
  #     }
  #   }

  if ($type =~ /ftp/) {
    print STDERR "fetching $source_file to $compressed_tmp_file\n" if ($self->debug);
    my $status = getstore($source_file, $compressed_tmp_file);
    die "load_set: error $status on $source_file" unless is_success($status);
  } elsif ($type =~ /file/) {
    print STDERR "copying $source_file to $compressed_tmp_file\n" if ($self->debug);
    my $status = system("cp $source_file $compressed_tmp_file");
    die "load_set: error $status on $source_file" unless (0 == $status);
  }
  my $decompress_method = 'gunzip';
  $decompress_method = 'bunzip2' if ($compression_ext eq 'bz2');
#   $cmd = "gunzip -c $compressed_tmp_file > $final_file";
  $cmd = "$decompress_method $compressed_tmp_file";
  unless(system("$cmd") == 0) {
    print("# $cmd\n") if ($self->debug);
    throw("error expanding file: $!\n");
  }
  printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);

  print STDERR "Creating fasta file...\n" if ($self->debug);
  my $fasta_out_filename = $assembly . "." . 'fasta';
  my $fasta_out_file = $qdir . "/" . $fasta_out_filename;

  my $buildset_executable = $self->analysis->program_file;
  unless (-e $buildset_executable) {
    $buildset_executable = "/nfs/users/nfs_a/avilella/src/maq/maq-0.7.1/scripts/fq_all2std.pl";
  }

  $cmd = "perl $buildset_executable fq2fa $tmp_file > $fasta_out_file";
  unless(system("$cmd") == 0) {
    print("$cmd\n") if ($self->debug);
    throw("error creating fastq to fasta [$tmp_file] to [$fasta_out_file] $!\n");
  }

  open TMPFILE, "$tmp_file" or die $!;
  my $file_size = 0; my $chunk_id = 0;
  my $this_fh = $self->create_new_chunk($tmp_file,$chunk_id);
  my $count = 0;
  while (<TMPFILE>) {
    # Read in 4 lines
    my $entry = $_;
    $entry .= <TMPFILE>;
    $entry .= <TMPFILE>;
    $entry .= <TMPFILE>;

    # next if all nucleotides are Ns
    next if ($entry =~ /^n+$/i);

    $file_size += length($entry);

    if ($file_size > 2000000000) {
      close $this_fh;
      $chunk_id++;
      $this_fh = $self->create_new_chunk($tmp_file,$chunk_id);
      $file_size = 0;

    }

    print $this_fh $entry;
    $count++;
  }
  close $this_fh;
  close TMPFILE;
  printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);

  print STDERR "$count entries in file\n" if ($self->debug);
  eval { $self->_store_tagvalue($assembly,1,'set_numreads',$count);};

  my $format = 'fastq';
  foreach my $id (0 .. $chunk_id) {
    my $new_name = $tmp_file . "_" . sprintf("%07d",$id);
    # $DB::single=1;1;
    print STDERR "Creating compressed indexed file $new_name...\n" if ($self->debug);
    my $fastq_out_filename = 
    my $compressed_indexed_out_file = $qdir . "/" . $assembly . ".$format.z_" . sprintf("%07d",$id);
    my $tag = $self->{tag} || 'HWUSI';
    my $fastq_tag = '@' . $tag;
    my $cmd = $self->{indexset_exe} || "/nfs/users/nfs_a/avilella/src/cdbfasta/latest/cdbfasta/cdbfasta";
    $cmd   .= " -r $fastq_tag -z $compressed_indexed_out_file $new_name";

    unless(system("$cmd") == 0) {
      print("$cmd\n") if ($self->debug);
      throw("[$cmd] error indexing $compressed_indexed_out_file $!\n");
    }
  }
  printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);

  return 1;
}

sub create_new_chunk {
  my $self = shift;
  my $orig_filename = shift;
  my $chunk_id = shift;
  my $fh;

  my $new_name = $orig_filename . "_" . sprintf("%07d",$chunk_id);
  open($fh, '>', "$new_name") or die $!;

  return $fh;
}

sub _store_tagvalue {
  my $self = shift;
  my $stable_id = shift;
  my $node_id = shift;
  my $tag = shift;
  my $value = shift;

  $value="" unless(defined($value));

  my $sql = "INSERT ignore into dna_hmm_tag (stable_id, node_id,tag) values (\"$stable_id\",$node_id,\"$tag\")";
  #print("$sql\n");
  $self->dbc->do($sql);

  $sql = "UPDATE dna_hmm_tag set value=\"$value\" where stable_id=\"$stable_id\" and node_id=$node_id and tag=\"$tag\"";
  #print("$sql\n");
  $self->dbc->do($sql);
}


1;
