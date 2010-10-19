#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::Createmodelsearchjobs

=cut

=head1 SYNOPSIS

my $aa = $sdba->get_AnalysisAdaptor;
my $analysis = $aa->fetch_by_logic_name('Createmodelsearchjobs');
my $rdb = new Bio::EnsEMBL::Compara::RunnableDB::Createmodelsearchjobs();
$rdb->fetch_input;
$rdb->run;

=cut

=head1 DESCRIPTION

This is a compara specific runnableDB, that based on an input
of arrayrefs of genome_db_ids, creates StoreSeqCDS jobs in the hive 
analysis_job table.

=cut

=head1 CONTACT

avilella@ebi.ac.uk

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Compara::RunnableDB::Createmodelsearchjobs;

use strict;
use Time::HiRes qw(time gettimeofday tv_interval);
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive;

use Bio::EnsEMBL::Hive::Process;
our @ISA = qw(Bio::EnsEMBL::Hive::Process);

sub fetch_input {
  my( $self) = @_;

  $self->{'reuse_clusterset_id'} = undef;
  $self->throw("No input_id") unless defined($self->input_id);

  #create a Compara::DBAdaptor which shares the same DBI handle
  #with the pipeline DBAdaptor that is based into this runnable
  $self->{'comparaDBA'} = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $self->{memberDBA} = $self->{'comparaDBA'}->get_MemberAdaptor;
  $self->{gdbDBA} = $self->{'comparaDBA'}->get_GenomeDBAdaptor;
  $self->{treeDBA} = $self->{'comparaDBA'}->get_ProteinTreeAdaptor;

  $self->get_params($self->parameters);

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
  print("parsing parameter string : ",$param_string,"\n");

  my $params = eval($param_string);
  return unless($params);

  foreach my $key (keys %$params) {
    print("  $key : ", $params->{$key}, "\n");
  }

  foreach my $key (qw[reuse_clusterset_id analysis_data_id]) {
    my $value = $params->{$key};
    $self->{$key} = $value if defined $value;
  }

  if (defined $params->{'reuse_clusterset_id'}) {
    $self->{'reuse_clusterset_id'} = $params->{'reuse_clusterset_id'};
  }

#   if (defined $params->{'species_set'}) {
#     $self->{'species_set'} = $params->{'species_set'};
#   }

  print("parameters...\n");
  printf("  reuse_clusterset_id    : %s\n", $self->{'reuse_clusterset_id'});

  return;
}

sub run
{
  my $self = shift;

  return 1 unless($self->{'reuse_clusterset_id'});
  $self->create_analysis_jobs($self->{'reuse_clusterset_id'});

  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

##########################################
#
# internal methods
#
##########################################

sub create_analysis_jobs {
  my $self = shift;
  my $starttime = time();

  my @dna_gdbs;
  foreach my $genome_db (@{$self->{gdbDBA}->fetch_all}) {
    my $genebuild = $genome_db->genebuild;
    if (!defined($genebuild) || $genebuild eq '') {
      push @dna_gdbs, $genome_db->dbID;
    }
  }

  my @aa_gdbs;
  @aa_gdbs = @{$self->{species_set}} if (defined($self->{species_set}));


  my $reuse_clusterset_id = $self->{reuse_clusterset_id};
  print STDERR "Fetching clusterset->children...\n" if ($self->debug);
  my $clusterset = $self->{treeDBA}->fetch_node_by_node_id($reuse_clusterset_id);
  printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
  my $counter = 0;
  my $ranked_ids;
  foreach my $root (@{$clusterset->children}) {
    my $root_id = $root->node_id;
    my $aln_percent_identity = $root->get_tagvalue('aln_percent_identity');
    my @leaves = @{$root->get_all_leaves};
    my $num_leaves = scalar @leaves;
    my $species_hash;
    foreach my $leaf (@leaves) {$species_hash->{$leaf->taxon_id} = 1;};
    my $num_species = scalar keys %$species_hash;
    my $total_score = 1;
    $total_score += ($aln_percent_identity-70)*10 if ($aln_percent_identity ne '' && $aln_percent_identity > 70 && $aln_percent_identity < 100 && $num_species > 2);
    $total_score += ($num_species*3);
    $total_score += $num_leaves;
    $ranked_ids->{$total_score}{$root_id} = 1;
    $root->release_tree;
    if($counter % 100 == 0) { 
      printf("%10d clusters\n", $counter) if ($self->debug);
      printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
    }
    $counter++;
  }

  $counter = 0;
  foreach my $score (sort {$b<=>$a} keys %{$ranked_ids}) {
    foreach my $root_id (keys %{$ranked_ids->{$score}}) {
      # $DB::single=1;1;
      if (0 < scalar(@dna_gdbs)) {
        my $type = 'dna';

        # Check if hmmprofile exists
        my $query = "SELECT hmmprofile FROM protein_tree_hmmprofile WHERE type=\"type\" AND node_id=$root_id";
        # print STDERR "$query\n" if ($self->debug);
        my $sth = $self->{comparaDBA}->dbc->prepare($query);
        $sth->execute();
        my $result = $sth->fetchrow_hashref;
        $sth->finish;
        my $output_id = sprintf("{'protein_tree_id'=>%d, 'clusterset_id'=>%d, 'type'=>\'%s\', 'cdna'=>'1'}", $root_id, $reuse_clusterset_id, $type);
        # print STDERR "Dataflow # $output_id\n" if ($self->debug);
        if (!defined($result->{hmmprofile})) {
          $self->dataflow_output_id($output_id, 2);
        } else {
          $self->dataflow_output_id($output_id, 1);
        }
      }
      if (0 < scalar(@aa_gdbs)) {
        my $type = 'aa';
        my $query = "SELECT hmmprofile FROM protein_tree_hmmprofile WHERE type=\"type\" AND node_id=$root_id";
        #      print STDERR "$query\n" if ($self->debug);
        my $sth = $self->{comparaDBA}->dbc->prepare($query);
        $sth->execute();
        my $result = $sth->fetchrow_hashref;
        $sth->finish;
        my $output_id = sprintf("{'protein_tree_id'=>%d, 'clusterset_id'=>%d, 'type'=>\'%s\'}", $root_id, $reuse_clusterset_id, $type);
        # print STDERR "Dataflow # $output_id\n" if ($self->debug);
        if (!defined($result->{hmmprofile})) {
          $self->dataflow_output_id($output_id, 3);
        } else {
          $self->dataflow_output_id($output_id, 1);
        }
      }
      if($counter % 100 == 0) { 
        printf("%10d clusters\n", $counter) if ($self->debug);
      }
      $counter++;
    }
  }

  return 1;
}

1;
