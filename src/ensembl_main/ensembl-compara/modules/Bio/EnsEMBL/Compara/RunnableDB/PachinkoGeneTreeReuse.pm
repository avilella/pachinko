#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::PachinkoGeneTreeReuse

=cut

=head1 SYNOPSIS

my $aa = $sdba->get_AnalysisAdaptor;
my $analysis = $aa->fetch_by_logic_name('PachinkoGeneTreeReuse');
my $rdb = new Bio::EnsEMBL::Compara::RunnableDB::PachinkoGeneTreeReuse(
                         -input_id   => "{'species_set'=>[1,2,3,14]}",
                         -analysis   => $analysis);

$rdb->fetch_input
$rdb->run;

=cut

=head1 DESCRIPTION

This analysis will import peptide_align_feature tables from previous
compara release using a mysql-dependent mysqldump pipe mysql command.

=cut

=head1 CONTACT

  Contact Albert Vilella on module implemetation/design detail: avilella@ebi.ac.uk
  Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Compara::RunnableDB::PachinkoGeneTreeReuse;

use strict;
use Switch;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::URLFactory;               # Blast_reuse
use Bio::EnsEMBL::Compara::MethodLinkSpeciesSet;
use Time::HiRes qw(time gettimeofday tv_interval);

our @ISA = qw(Bio::EnsEMBL::Hive::Process);

sub fetch_input {
  my( $self) = @_;

  $self->throw("No input_id") unless defined($self->input_id);
  $self->{tbl_name} = $self->input_id;

  #create a Compara::DBAdaptor which shares the same DBI handle
  #with the Pipeline::DBAdaptor that is based into this runnable
  $self->{'comparaDBA'} = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $self->{gdba}         = $self->{'comparaDBA'}->get_GenomeDBAdaptor;

  my $p = eval($self->analysis->parameters);
  if (defined $p->{'analysis_data_id'}) {
    my $analysis_data_id = $p->{'analysis_data_id'};
    my $ada = $self->db->get_AnalysisDataAdaptor;
    my $new_params = eval($ada->fetch_by_dbID($analysis_data_id));
    if (defined $new_params) {
      $p = $new_params;
    }
  }
  $self->{p} = $p;
  $self->{null_cigar} = $p->{null_cigar} if (defined($p->{null_cigar}));

  if ($self->debug) {
    print("parameters...\n");
    foreach my $key (keys %{$self->{p}}) { print "  $key -- ", $self->{p}{$key}, "\n"; }
  }

#   # Check if this is one that we need to reuse
#   foreach my $reusable_gdb (@{$p->{reuse_gdb}}) {
#     $self->{reusable_gdb}{$reusable_gdb} = 1;
#   }

  return 1;
}

sub run
{
  my $self = shift;

  $self->import_genetree_data;
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

sub import_genetree_data {
  my $self = shift;
  my $starttime = time();

  $self->{comparaDBA_reuse} = Bio::EnsEMBL::Hive::URLFactory->fetch($self->{p}{reuse_db}, 'compara');
  my $reuse_username = $self->{comparaDBA_reuse}->dbc->username;
  my $reuse_password = $self->{comparaDBA_reuse}->dbc->password;
  my $pass = "-p$reuse_password " if ($reuse_password);
  my $reuse_host = $self->{comparaDBA_reuse}->dbc->host;
  my $reuse_port = $self->{comparaDBA_reuse}->dbc->port;
  my $reuse_dbname = $self->{comparaDBA_reuse}->dbc->dbname;

  my $dest_username = $self->dbc->username;
  my $dest_password = $self->dbc->password;
  my $dest_pass = "-p$dest_password" if ($dest_password);
  my $dest_host = $self->dbc->host;
  my $dest_port = $self->dbc->port;
  my $dest_dbname = $self->dbc->dbname;

#  my $genome_dbs = join (",",keys %{$self->{reusable_gdb}});
  my $tbl_name = $self->{tbl_name};

  my $mysqldumpcmd = "mysqldump --compress --lock_tables=FALSE --add-drop-table=TRUE ";

  print STDERR "mysqldump $tbl_name\n" if ($self->debug);

  my $cmd = "$mysqldumpcmd -u $reuse_username $pass -h $reuse_host -P$reuse_port $reuse_dbname $tbl_name";
  $cmd .= " | perl -pe 's/MyISAM/InnoDB/' | mysql -u $dest_username $dest_pass -h $dest_host -P$dest_port $dest_dbname";
#  $cmd .= " | mysql -u $dest_username $dest_pass -h $dest_host -P$dest_port $dest_dbname";

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  my $ret = system($cmd);
  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);
  printf("  %1.3f secs to mysqldump $tbl_name\n", (time()-$starttime));
  if (0 != $ret) { throw("Error importing $tbl_name: $ret\n"); }

  $self->test_table($self->{comparaDBA_reuse},$self->{comparaDBA},$tbl_name);

# member
# subset
# subset_member
# sequence
# sequence_cds
# sequence_exon_bounded
# protein_tree_member
# protein_tree_node
# protein_tree_tag
# protein_tree_hmmprofile
# super_protein_tree_node
# super_protein_tree_member
# super_protein_tree_tag
# protein_tree_stable_id
# mapping_session
# stable_id_history
# sitewise_aln
# # homology
# # homology_member

  # my $cmd = "mysqldump --compress --where=\"hgenome_db_id in ($hgenome_dbs)\" -u $reuse_username $pass -h $reuse_host -P$reuse_port $reuse_dbname $tbl_name";

  # my $mysqldumpcmd = "mysqldump --compress --lock_tables=FALSE --add-drop-table=FALSE --no-create-info";
  # my $mysqldumpcmd = "mysqldump --compress --lock_tables=FALSE --add-drop-table=TRUE ";


#   # Fill the member table
#   $tbl_name = 'member';
#   # $self->create_if_not_exists($tbl_name);
#   print STDERR "mysqldump $tbl_name\n" if ($self->debug);
#   my $cmd = "$mysqldumpcmd --where=\"genome_db_id in ($genome_dbs)\" -u $reuse_username $pass -h $reuse_host -P$reuse_port $reuse_dbname $tbl_name";
#   $cmd .= " | mysql -u $dest_username $dest_pass -h $dest_host -P$dest_port $dest_dbname";

#   my $ret = system($cmd);
#   printf("  %1.3f secs to mysqldump $tbl_name\n", (time()-$starttime));
#   if (0 != $ret) { throw("Error importing $tbl_name: $ret\n"); }

#   # Now we do all member_id FK tables
#   my @reuse_member_id_fk = ('subset_member');
#   foreach my $tbl_name (@reuse_member_id_fk) {
#     # $self->create_if_not_exists($tbl_name);
#     print STDERR "mysqldump $tbl_name\n" if ($self->debug);
#     my $cmd = "$mysqldumpcmd --where=\"member_id in (select member_id from member where genome_db_id in ($genome_dbs))\" -u $reuse_username $pass -h $reuse_host -P$reuse_port $reuse_dbname $tbl_name";
#     $cmd .= " | mysql -u $dest_username $dest_pass -h $dest_host -P$dest_port $dest_dbname";

#     my $ret = system($cmd);
#     printf("  %1.3f secs to mysqldump $tbl_name\n", (time()-$starttime));
#     if (0 != $ret) { throw("Error importing $tbl_name: $ret\n"); }
#   }

#   # sequence
#   # sequence_cds
#   # sequence_exon_bounded
#   my @sequence_id_fk = ('sequence','sequence_cds','sequence_exon_bounded');
#   foreach my $tbl_name (@sequence_id_fk) {
#     # $self->create_if_not_exists($tbl_name);
#     print STDERR "mysqldump $tbl_name\n" if ($self->debug);
#     my $field = $tbl_name . "_id";

#     my $cmd = "$mysqldumpcmd --where=\"$field in (select sequence_id from member where genome_db_id in ($genome_dbs))\" -u $reuse_username $pass -h $reuse_host -P$reuse_port $reuse_dbname $tbl_name";
#     $cmd .= " | mysql -u $dest_username $dest_pass -h $dest_host -P$dest_port $dest_dbname";

#     my $ret = system($cmd);
#     printf("  %1.3f secs to mysqldump $tbl_name\n", (time()-$starttime));
#     if (0 != $ret) { throw("Error importing $tbl_name: $ret\n"); }
#   }

  # The ones below should only be run once, for no specific genome_db_id
    # protein_tree_member
    # protein_tree_node
    # protein_tree_tag
    # protein_tree_hmmprofile
    # super_protein_tree_node
    # super_protein_tree_member
    # super_protein_tree_tag
    # protein_tree_stable_id
    # mapping_session
    # stable_id_history
    # sitewise_aln
#     my @all_data_in_table = ('subset_member','member','sequence','sequence_cds','sequence_exon_bounded','subset','protein_tree_member','protein_tree_node','protein_tree_tag','protein_tree_hmmprofile','super_protein_tree_node','super_protein_tree_member','super_protein_tree_tag','protein_tree_stable_id','mapping_session','stable_id_history','sitewise_aln');
#     foreach my $tbl_name (@all_data_in_table) {
#       print STDERR "mysqldump $tbl_name\n" if ($self->debug);
#       my $cmd = "$mysqldumpcmd -u $reuse_username $pass -h $reuse_host -P$reuse_port $reuse_dbname $tbl_name";
#       $cmd .= " | mysql -u $dest_username $dest_pass -h $dest_host -P$dest_port $dest_dbname";

#       my $ret = system($cmd);
#       printf("  %1.3f secs to mysqldump $tbl_name\n", (time()-$starttime));
#       if (0 != $ret) { throw("Error importing $tbl_name: $ret\n"); }
#     }

  return 1;
}


=head2 test_table

=cut

sub test_table {
  my ($self, $db1, $db2, $table_name) = @_;

  die "Cannot test table with no name\n" if (!$table_name);

  ## check the table is not empty
  my $count1 = $db1->dbc->db_handle->selectrow_array(
      "SELECT COUNT(*) FROM $table_name");

  ## check the table is not empty
  my $count2 = $db2->dbc->db_handle->selectrow_array(
      "SELECT COUNT(*) FROM $table_name");

  if ($count1 != $count2) {
    die("There is a different number of entries in the $table_name: $count1 != $count2\n");
  } else {
    print "Table $table_name OK.\n";
  }

}


1;
