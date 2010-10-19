#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::Createspeciestree

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $createspeciestree = Bio::EnsEMBL::Compara::RunnableDB::Createspeciestree->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$createspeciestree->fetch_input(); #reads from DB
$createspeciestree->run();
$createspeciestree->output();
$createspeciestree->write_output(); #writes to DB

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


package Bio::EnsEMBL::Compara::RunnableDB::Createspeciestree;

use strict;
use Getopt::Long;
use Time::HiRes qw(time gettimeofday tv_interval);

use Bio::EnsEMBL::Compara::Graph::NewickParser;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
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
  $self->{treeDBA} = $self->{'comparaDBA'}->get_ProteinTreeAdaptor;
  $self->{taxonDBA} = $self->{'comparaDBA'}->get_NCBITaxonAdaptor;

#   $self->get_params($self->parameters);
#   $self->get_params($self->input_id);

# # For long parameters, look at analysis_data
#   if($self->{analysis_data_id}) {
#     my $analysis_data_id = $self->{analysis_data_id};
#     my $analysis_data_params = $self->db->get_AnalysisDataAdaptor->fetch_by_dbID($analysis_data_id);
#     $self->get_params($analysis_data_params);
#   }

  foreach my $gdb (@{$self->{gdba}->fetch_all}) {
    $self->{species_set_hash}{$gdb->taxon_id} = 1;
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

  foreach my $key (qw[param1 param2 param3 analysis_data_id]) {
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

  $self->run_createspeciestree;
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

sub run_createspeciestree {
  my $self = shift;
  my $root;

  my $errors = 0;
  my $clusterset = $self->{treeDBA}->fetch_node_by_node_id(1);
  foreach my $taxon_id (keys %{$self->{species_set_hash}}) {
    my $taxon = $self->{taxonDBA}->fetch_node_by_taxon_id($taxon_id);
    if (!defined ($taxon)) {
      $clusterset->store_tag('undefined_taxon',$taxon_id);
      $errors = 1;
    }
    print STDERR "  [$taxon_id]\n";
    $taxon->release_children;
    $root = $taxon->root unless($root);
    $root->merge_node_via_shared_ancestor($taxon);
  }

  throw("Undefined taxa, check 'undefined_taxon' tag in protein_tree_tag table") if ($errors);

  $root = $root->minimize_tree;
  my $species_tree_string = $root->newick_format('njtree');

  my $old_species_tree_string = $clusterset->get_tagvalue('species_tree_string');
  if (0 < length($old_species_tree_string)) {
    $clusterset->store_tag('old_species_tree_string',$old_species_tree_string);
  }

  my $eval_species_tree;
  eval {
    $eval_species_tree = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($species_tree_string);
    my @leaves = @{$eval_species_tree->get_all_leaves};
  };

  unless ($@) {
    $clusterset->store_tag('species_tree_string',$species_tree_string);
  } else {
    throw("Error creating new species_tree_string: \n $species_tree_string \n");
  }

  return 1;
}

1;
