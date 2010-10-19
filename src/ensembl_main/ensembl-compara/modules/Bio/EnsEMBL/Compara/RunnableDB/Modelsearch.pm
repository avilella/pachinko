#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::Modelsearch

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $search_hmm = Bio::EnsEMBL::Compara::RunnableDB::Modelsearch->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$search_hmm->fetch_input(); #reads from DB
$search_hmm->run();
$search_hmm->output();
$search_hmm->write_output(); #writes to DB

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


package Bio::EnsEMBL::Compara::RunnableDB::Modelsearch;

use strict;
use Getopt::Long;
use Time::HiRes qw(time gettimeofday tv_interval);

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
  $self->{max_evalue} = 0.05;

  if (1 == $self->input_id) {$self->{done} = 1; return;}

  #create a Compara::DBAdaptor which shares the same DBI handle
  #with the Pipeline::DBAdaptor that is based into this runnable
  $self->{'comparaDBA'} = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
    (
     -DBCONN=>$self->db->dbc
    );

  # Get the needed adaptors here
  $self->{treeDBA} = $self->{'comparaDBA'}->get_ProteinTreeAdaptor;
  $self->{gdba} = $self->{'comparaDBA'}->get_GenomeDBAdaptor;

  $self->get_params($self->parameters);
  $self->get_params($self->input_id);

# For long parameters, look at analysis_data
  if($self->{analysis_data_id}) {
    my $analysis_data_id = $self->{analysis_data_id};
    my $analysis_data_params = $self->db->get_AnalysisDataAdaptor->fetch_by_dbID($analysis_data_id);
    $self->get_params($analysis_data_params);
  }

  if(defined($self->{protein_tree_id})) {
    $self->{tree} = 
         $self->{treeDBA}->fetch_node_by_node_id($self->{protein_tree_id});
    printf("  protein_tree_id : %d\n", $self->{protein_tree_id}) if ($self->debug);
  }


  foreach my $genome_db (@{$self->{gdba}->fetch_all}) {
    my $genebuild = $genome_db->genebuild;
    if ((!defined($genebuild) || $genebuild eq '') && $self->{type} eq 'dna') {
      push @{$self->{dna_gdb}}, $genome_db;
    } elsif ((defined($genebuild) && $genebuild eq '') && $self->{type} eq 'aa') {
      # do something else
    }
  }

  if ($self->{search_mode} eq 'fast') {
    $self->{done} = 1;
    return undef;
  }
  unless ($self->{search_mode} =~ /hmm/) {
    $self->{done} = 1;
    return undef;
  }
  # Fetch hmm_profile
  $self->fetch_hmmprofile;

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

  foreach my $key (qw[search_mode qtaxon_id clusterset_id protein_tree_id type gdb cdna fasta_dir analysis_data_id]) {
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

  if ($self->{done}) { return; }
  $self->run_search_hmm;
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

  $self->search_hmm_store_hits unless ($self->{done});
  $self->dataflow_hit_clusters;
}


##########################################
#
# internal methods
#
##########################################


sub fetch_hmmprofile {
  my $self = shift;

  my $hmm_type = $self->{type} || 'aa';
  my $node_id = $self->{tree}->node_id;
  print STDERR "type = $hmm_type\n" if ($self->debug);

  my $query = "SELECT hmmprofile FROM protein_tree_hmmprofile WHERE type=\"$hmm_type\" AND node_id=$node_id";
  print STDERR "$query\n" if ($self->debug);
  my $sth = $self->{comparaDBA}->dbc->prepare($query);
  $sth->execute;
  my $result = $sth->fetchrow_hashref;
  $self->{hmmprofile} = $result->{hmmprofile} if (defined($result->{hmmprofile}));
  $sth->finish;

  return 1;
}

sub run_search_hmm {
  my $self = shift;

  my $node_id = $self->{tree}->node_id;
  my $type = $self->{type};
  my $hmmprofile = $self->{hmmprofile};
  my $fasta_dir = $self->{fasta_dir};

  my $species_list;
  foreach my $leaf (@{$self->{tree}->get_all_leaves}) {
    $species_list->{$leaf->taxon_id};
  }

  my $tempfilename = $self->worker_temp_directory . $node_id . "." . $type . ".hmm";
  open FILE, ">$tempfilename" or die "$!";
  print FILE $hmmprofile;
  close FILE;
  delete $self->{hmmprofile};

  my $search_hmm_executable = $self->analysis->program_file;
  unless (-e $search_hmm_executable) {
    $search_hmm_executable = "/software/ensembl/compara/hmmer3/hmmer-3.0/src/hmmsearch";
  }

  # Do the dna_gdb like this
  foreach my $genome_db (@{$self->{dna_gdb}}) {
    my $assembly = $genome_db->assembly;
    my $locator = $genome_db->locator;
    my $format = 'fasta';
    my $qtaxon_id = $genome_db->taxon_id;
    my $qdir = $self->get_qdir($locator);
    $DB::single=1;1;
    # Decide if we want to run it for this assembly
    my $readlen = 'short'; $locator =~ /readlen\=(\w+)\;/;
    $readlen = $1 if (defined $1);
    next if ($readlen != 'long');
    my $fastafile = $qdir . "/" . "$assembly.$format";
    if ($type eq 'aa') {
      my $name = $genome_db->name; $name =~ s/\ /\_/g;
      $fastafile = $qdir . "/" . $name . "_" . "$assembly.$format";

    }

    my $dtfile = $self->worker_temp_directory . $node_id . "." . $type . ".dt";

    my $noalioption = "--noali"; $noalioption = "" if ($self->debug);
    my $cmd = "$search_hmm_executable --cpu 1 --domtblout $dtfile $noalioption $tempfilename $fastafile";
    print STDERR "# $cmd\n" if ($self->debug);
    my $starttime = time();
    $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
    unless(system("$cmd") == 0) {
      print("$cmd\n");
      $self->throw("error running modelsearch $!\n");
    }
    $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

    open DTFILE,"$dtfile" or die "$!";
    while (<DTFILE>) {
      next if ($_ =~ /^#/);
      my ($target_name,$q_accession,$tlen,$query_name,$h_accession,$qlen,$full_E_value,$full_score,$full_bias,$number,$of,$c_Evalue,$i_Evalue,$domain_score,$domain_bias,$hmm_from,$hmm_to,$ali_from,$ali_to,$env_from,$env_to,$acc,$description) = split(" ",$_);
      # $DB::single=1;1;
      $self->{hits}{$target_name}{Score} = $domain_score;
      $self->{hits}{$target_name}{Evalue} = $i_Evalue;
      $self->{hits}{$target_name}{hmm_from} = $hmm_from;
      $self->{hits}{$target_name}{hmm_to} = $hmm_to;
    }
    close(DTFILE);

## #                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
## # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
## #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
## SRR034143.6280281    -             75 1929                 -            822   0.00014   28.9   0.0   1   1   4.1e-11   0.00014   28.9   0.0   464   536     2    74     1    75 0.94 -
## SRR034143.5856316    -             75 1929                 -            822       5.6   13.7   0.0   1   1   1.6e-06       5.6   13.7   0.0   547   588    33    74    12    75 0.93 -
## ~
##

    $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);
    print STDERR scalar (keys %{$self->{hits}}), " hits - ",(time()-$starttime)," secs...\n";

    $self->search_hmm_store_hits($qtaxon_id,$assembly);
    delete $self->{hits};
  }

  return 1;
}

sub search_hmm_store_hits {
  my $self = shift;
  my $qtaxon_id = shift;
  my $assembly = shift;
  my $type = $self->{type};
  my $node_id = $self->{tree}->node_id;

  my $sth = $self->{comparaDBA}->dbc->prepare
    ("INSERT INTO hmmsearch
       (stable_id,
        node_id,
        hmm_from,
        hmm_to,
        evalue,
        score,
        type,
        qtaxon_id,
        set_id) VALUES (?,?,?,?,?,?,?,?,?)");

  my $evalue_count = 0;
  foreach my $stable_id (keys %{$self->{hits}}) {
    my $evalue = $self->{hits}{$stable_id}{Evalue};
    next if ($evalue > 0.01); #FIXME do we want to store everything?
    my $score = $self->{hits}{$stable_id}{Score};
    my $hmm_from = $self->{hits}{$stable_id}{hmm_from};
    my $hmm_to   = $self->{hits}{$stable_id}{hmm_to};
    next unless (defined($stable_id) && $stable_id ne '');
    next unless (defined($score));
    next unless ($evalue < $self->{max_evalue});
    $evalue_count++;
    $sth->execute($stable_id,
                  $node_id,
                  $hmm_from,
                  $hmm_to,
                  $evalue,
                  $score,
                  $type,
                  $qtaxon_id,
                  $assembly);
    $self->{qtaxon_id_with_hits}{$qtaxon_id} = $type;
    printf("%10d hits stored\n", $evalue_count) if(($evalue_count> 0) && ($evalue_count % 10 == 0) && 0 < $self->debug);
  }
  printf("%10d hits stored\n", $evalue_count) if($self->debug);
  $sth->finish();
  return 1;
}

sub dataflow_hit_clusters {
  my $self = shift;
  my $node_id = $self->{tree}->node_id;
  my $clusterset_id = $self->{tree}->parent->node_id;
  my $type = $self->{type};

  foreach my $genome_db (@{$self->{dna_gdb}}) {
    my $qtaxon_id = $genome_db->taxon_id;
    next if (defined $self->{qtaxon_id_list}{$qtaxon_id});
    $self->{qtaxon_id_list}{$qtaxon_id} = 1;
    # Dataflow clusters
    # This will create a new MSA alignment job for each of the newly generated clusters
    my $output_id = sprintf("{'protein_tree_id'=>%d, 'clusterset_id'=>%d, 'qtaxon_id'=>%d, 'type'=>\'%s\'}", 
                            $node_id, $clusterset_id, $qtaxon_id, $type);
    print STDERR "Dataflow # $output_id\n" if ($self->debug);
    $self->dataflow_output_id($output_id, 2);
  }

#   foreach my $qtaxon_id (keys %{$self->{qtaxon_id_with_hits}}) {
#     my $type = $self->{qtaxon_id_with_hits}{$qtaxon_id};
#     # Dataflow clusters
#     # This will create a new MSA alignment job for each of the newly generated clusters
#     my $output_id = sprintf("{'protein_tree_id'=>%d, 'clusterset_id'=>%d, 'qtaxon_id'=>%d, 'type'=>\'%s\'}", 
#                             $node_id, $clusterset_id, $qtaxon_id, $type);
#     print STDERR "Dataflow # $output_id\n" if ($self->debug);
#     $self->dataflow_output_id($output_id, 2);
#   }

  return 1;
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
