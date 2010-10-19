#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::UsearchPrepare

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $usearchprepare = Bio::EnsEMBL::Compara::RunnableDB::UsearchPrepare->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$usearchprepare->fetch_input(); #reads from DB
$usearchprepare->run();
$usearchprepare->output();
$usearchprepare->write_output(); #writes to DB

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


package Bio::EnsEMBL::Compara::RunnableDB::UsearchPrepare;

use strict;
use Getopt::Long;
use Time::HiRes qw(time gettimeofday tv_interval);

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive;
use Bio::Seq;
use Bio::SeqIO;
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

  # This is the length of the sequences given to uclust,
  # if longer, they are substring'ed with an overlap

  # Length of the uclust reference sequences to be split (uclust is fast for ~1000)
  $self->{uclust_length} = 1000;
  # Overlap for the subseq splitting
  $self->{uclust_overlap} = 100;
  # Amount of bps for each reference file (~100Mbp is 1.5GB RAM in uclust)
  $self->{uclust_reference_size} = 100000000;
  # Number of blocks in which each fastq file is split per job
  $self->{uclust_block_num} = 100;

  #create a Compara::DBAdaptor which shares the same DBI handle
  #with the Pipeline::DBAdaptor that is based into this runnable
  $self->{'comparaDBA'} = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
    (
     -DBCONN=>$self->db->dbc
    );

  # Get the needed adaptors here
  # $self->{silly_adaptor} = $self->{'comparaDBA'}->get_SillyAdaptor;
  $self->{gdba} = $self->{'comparaDBA'}->get_GenomeDBAdaptor;

  # list of genome_db sets
  foreach my $gdb (@{$self->{gdba}->fetch_all}) {
    next if (defined($gdb->genebuild) && $gdb->genebuild ne '');
    $self->{gdb_sets}{$gdb->dbID} = $gdb;
  }


  $self->get_params($self->parameters);
  # $self->get_params($self->input_id);

# For long parameters, look at analysis_data
  if($self->{analysis_data_id}) {
    my $analysis_data_id = $self->{analysis_data_id};
    my $analysis_data_params = $self->db->get_AnalysisDataAdaptor->fetch_by_dbID($analysis_data_id);
    $self->get_params($analysis_data_params);
  }

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

  foreach my $key (qw[qdirs uclust_block_num uclust_reference_size uclust_overlap uclust_length analysis_data_id]) {
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

  $self->run_usearchprepare;
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

  $self->dataflow_usearch_jobs;
}


##########################################
#
# internal methods
#
##########################################

sub run_usearchprepare {
  my $self = shift;
  my $length  = $self->{uclust_length};
  my $overlap = $self->{uclust_overlap};

  my $starttime = time();
  throw("undefined qdirs. Check your config file\n") if (!defined $self->{qdirs});
  my @qdirs = split(":",$self->{qdirs});
  my $qdir = shift @qdirs;
  my $orig_filename = $qdir . "/" . "reference.fasta";

  my $sql1 = "SELECT ptn.root_id, ptn.parent_id from protein_tree_node ptn where ptn.node_id=ptn.root_id";
  my $sql2 = "SELECT ptm.member_id, s.sequence_cds, ptm.root_id from protein_tree_member ptm, sequence_cds s WHERE ptm.member_id=s.member_id and ptm.root_id=?";
  my $sth1 = $self->{comparaDBA}->dbc->prepare($sql1);
  my $sth2 = $self->{comparaDBA}->dbc->prepare($sql2);
  $sth1->execute;
  my $file_size = 0; my $chunk_id = 0;
  my $this_seqio = $self->create_new_seqio_chunk($orig_filename,$chunk_id);
  while (my $result1 = $sth1->fetchrow_hashref) {
    next if ($result1->{parent_id} < 1);
    my $root_id = $result1->{root_id};
    print STDERR "Fetching reference sequences: $root_id\n";
    $sth2->execute($root_id);
    printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
    my $sequence_hash;
    while (my $result2 = $sth2->fetchrow_hashref) {
      my $seq_string = $result2->{sequence_cds};
      my $id         = $result2->{root_id} . "_" . $result2->{member_id};
      next if (length($seq_string) < 10);
      my $seq = Bio::Seq->new
        (-display_id => $id,
         -seq => $seq_string);
      my $seq_length = $seq->length();
      # Loops through the sequence
      for (my $i = 1; $i < $seq_length; $i += $length) {
        my $end = (($i + $length) > $seq_length) ? ($seq_length + 1) : ($i + $length);
        my $seq_range = (($i + $length) > $seq_length) ? "$i-".($end - 1) : "$i-$end";
        my $id = $seq->id();
        $id .= "_$seq_range";
        my $trunc_seq = $seq->trunc($i, $end - 1);
        $trunc_seq->id($id);
        # Decreases the $i value with the offset value
        $i -= $overlap;

        # We only output identical sequences once
        if (defined $sequence_hash->{$trunc_seq->seq}) {
          next;
        }
        $sequence_hash->{$trunc_seq->seq} = 1;

        $this_seqio->write_seq($trunc_seq);
        $file_size += $trunc_seq->length;
      }

      # about 100Mbp per file
      if ($file_size > $self->{uclust_reference_size}) {
        $this_seqio->close;
        $chunk_id++;
        $this_seqio = $self->create_new_seqio_chunk($orig_filename,$chunk_id);
        $file_size = 0;
      }
    }
  }
  $this_seqio->close;

  $sth1->finish;
  $sth2->finish;
  return 1;
}

sub create_new_seqio_chunk {
  my $self = shift;
  my $orig_filename = shift;
  my $chunk_id = shift;
  my $fh;

  my $new_name = $orig_filename . "_" . sprintf("%07d",$chunk_id);
  my $this_seqio = Bio::SeqIO->new
    (-format => 'fasta',
     -file => ">$new_name");

  print STDERR "Reference file $new_name\n" if ($self->debug);

  $self->{reference_files}{$new_name} = 1;
  return $this_seqio;
}

sub dataflow_usearch_jobs {
  my $self = shift;
  my $block_num = $self->{uclust_block_num} || 100;

  my $count = 0;
  foreach my $reference_file (keys %{$self->{reference_files}}) {
    foreach my $genome_db_id (keys %{$self->{gdb_sets}}) {
      my $genome_db = $self->{gdb_sets}{$genome_db_id};
      my $assembly = $genome_db->assembly;
      my $locator = $genome_db->locator;
      my $format = 'fasta';
      my $qdir = $self->get_qdir($locator);
      my $read_file = $qdir . "/" . "$assembly.$format";
      foreach my $block (0 .. $block_num) {
        my $input_id = "{\'reference_file\'=>\'$reference_file\',\'read_file\'=>\'$read_file\',\'block\'=>\'$block/$block_num\'}";
        $self->dataflow_output_id($input_id, 1);
        $count++;
      }
    }
  }
  print STDERR "Number of jobs created: $count\n" if ($self->debug);

  return;
}

sub get_qdir {
  my $self = shift;
  my $locator = shift;

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

  throw("could not find qdir for $string") unless (defined($locator_hash->{qdir}));

  return $locator_hash->{qdir};
}


1;
