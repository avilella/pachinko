#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::Usearch

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $search_u = Bio::EnsEMBL::Compara::RunnableDB::Usearch->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$search_u->fetch_input(); #reads from DB
$search_u->run();
$search_u->output();
$search_u->write_output(); #writes to DB

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


package Bio::EnsEMBL::Compara::RunnableDB::Usearch;

use strict;
use Getopt::Long;
use Time::HiRes qw(time gettimeofday tv_interval);

use File::Basename;

use Bio::SeqIO;
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
  $self->{uclust_id}  = 0.70;

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

  # Fetch u_profile
  $self->fetch_reads_block;

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

  foreach my $key (qw[exonerate_exe uclust_exe search_mode uclust_id reference_file read_file block qtaxon_id clusterset_id type analysis_data_id]) {
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

  $self->run_ublast_block;
  # $self->run_search_u;
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

  # $self->search_u_store_hits;
  # $self->dataflow_hit_clusters;
}


##########################################
#
# internal methods
#
##########################################

sub fetch_reads_block {
  my $self = shift;
  my $starttime = time();

  # /nfs/users/nfs_a/avilella/src/exonerate/latest/exonerate-2.2.0-x86_64/bin/fastasplit -f /lustre/scratch103/ensembl/bb2/hive/pa58m/bDB/61GFCAAXX_2_1_pf.fasta -o /tmp/fastasplit/ -c 100

  my ($block_numer,$block_denom) = split("/",$self->{block});
  my $read_file = $self->{read_file};
  my $worker_temp_directory = $self->worker_temp_directory;
  my ($infilebase,$path,$type) = fileparse($read_file);
  my $assembly = $infilebase; $assembly =~ s/\.fasta//;  $self->{this_assembly} = $assembly;

  my $keep_file = $worker_temp_directory . $infilebase . "_chunk_" . sprintf("%07d",$block_numer);

  my $exonerate_executable = $self->{exonerate_exe} || "/nfs/users/nfs_a/avilella/src/exonerate/latest/exonerate-2.2.0-x86_64/bin/exonerate";
  $exonerate_executable =~ s/exonerate$/fastasplit/;
  # $DB::single=1;1;
  my $cmd = "$exonerate_executable -f $read_file -o $worker_temp_directory -c $block_denom";
  print STDERR "Creating chunks $worker_temp_directory ...\n" if ($self->debug);
  unless(system("$cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running fastasplit $!\n");
  }

  print STDERR "Deleting chunks not for block $block_numer: $keep_file ...\n" if ($self->debug);
  foreach my $this_block (0 .. $block_denom) {
    my $this_file = $worker_temp_directory . $infilebase . "_chunk_" . sprintf("%07d",$this_block);
    next if ($this_file eq $keep_file);
    $cmd = "rm -f $this_file";
    unless(system("$cmd") == 0) {
      print("$cmd\n");
      $self->throw("error running rm -f $this_file $!\n");
    }
  }

  print STDERR "worker_temp_directory $worker_temp_directory\n" if ($self->debug);

  $self->{useqs} = $keep_file;
  print STDERR (time()-$starttime)," secs...\n" if ($self->debug);

  return 1;
}

sub fetch_useqs {
  my $self = shift;

  my $u_type = $self->{type} || 'aa';
  my $node_id = $self->{tree}->node_id;
  print STDERR "type = $u_type\n" if ($self->debug);

  my $seqs_file = $self->worker_temp_directory . "$node_id.fasta";
  open SEQS, ">$seqs_file" or die $!;
  my $clusterset = $self->{treeDBA}->fetch_node_by_node_id(1);
  my $count = 0;
  my $starttime = time();

  foreach my $root (@{$clusterset->children}) {
    foreach my $leaf (@{$root->get_all_leaves}) {
      my $sequence;
      eval {$sequence = $leaf->sequence_cds };
      next if (length($sequence) < 1);
      $count++;
      $sequence = $leaf->sequence if ($u_type eq 'aa');
      print SEQS ">" . $leaf->stable_id . "\n" . $sequence . "\n";
    }
    print STDERR "[$count]\n" if ($self->debug);
    last if ($count > 50000);
    $root->release_tree;
  }
  close SEQS;

#   foreach my $leaf (@{$self->{tree}->get_all_leaves}) {
#     my $sequence;
#     eval {$sequence = $leaf->sequence_cds };
#     next if (length($sequence) < 1);
#     $count++;
#     $sequence = $leaf->sequence if ($u_type eq 'aa');
#     print SEQS ">" . $leaf->stable_id . "\n" . $sequence . "\n";
#   }
#   close SEQS;

  print STDERR scalar "[$count] " ,(time()-$starttime)," secs...\n";

  $self->{useqs} = $seqs_file;
  return 1;
}
sub run_ublast_block {
  my $self = shift;
  my $uclust_id = $self->{uclust_id};

  my $reference_file = $self->{reference_file};
  my $useqsfile = $self->{useqs};
  my $assembly = $self->{this_assembly};
  throw("[$reference_file reference_file] [$useqsfile useqsfile] [$assembly assembly] need to be defined") unless (defined $reference_file && defined $useqsfile && defined $assembly);

  my $search_u_executable = $self->analysis->program_file;
  unless (-e $search_u_executable) {
    $search_u_executable = "/nfs/users/nfs_a/avilella/src/uclust/latest/uclust";
  }

  my $uclust_table = "usearch_" . $assembly;
  my $create = $self->{comparaDBA}->dbc->prepare
    ("CREATE TABLE IF NOT EXISTS $uclust_table like usearch");
  $create->execute;
  $create->finish;

  my $sth = $self->{comparaDBA}->dbc->prepare
    ("INSERT IGNORE INTO $uclust_table
       (stable_id,root_id) VALUES (?,?)");

  ## This is for references as db and reads as query
  my $cmd = "$search_u_executable --maxlen 1000 --db $reference_file --query $useqsfile --uc /dev/stdout --id $uclust_id";
  $cmd .= " 2>/dev/null" unless ($self->debug);
  print STDERR "# $cmd\n" if ($self->debug);
  my $starttime = time();
  my $fh;
  eval { open($fh, "$cmd |") || die $!; };
  if ($@) {
    warn("error running usearch $@ $!");
    return;
  }

  my $count = 0;
  while (<$fh>) {
    next if ($_ =~ /^#/);
    next unless ($_ =~ /^H/);
    chomp $_;
    my ($Type, $ClusterNr, $SeqLengthorClusterSize, $PctId, $Strand, $QueryStart, $SeedStart, $Alignment, $QueryLabel, $TargetLabel) = split("\t",$_);
    $DB::single=1;1;#??
    # We dont take hits that are smaller than 20 in size
    next if ($SeqLengthorClusterSize < 20);
    ## This is for references as db and reads as query
    my ($root_id,$member_id,$range) = split("_",$TargetLabel);
    print STDERR "$QueryLabel, [$root_id] [$member_id] [$range] ($count hits)\n" if ($count % 10 == 0 && $self->debug);
    $sth->execute($QueryLabel, $root_id);
    $count++;
    print STDERR "$count hits - ",(time()-$starttime)," secs...\n" if ($count % 100 ==0 && $self->debug);
  }

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);
  print STDERR "Finished.\n";
  print STDERR "$count hits - ",(time()-$starttime)," secs...\n";

  return 1;
}

sub search_u_store_hits {
  my $self = shift;
  my $qtaxon_id = shift;
  my $assembly = shift;
  my $type = $self->{type};
  my $node_id = $self->{tree}->node_id;

  my $sth = $self->{comparaDBA}->dbc->prepare
    ("INSERT INTO usearch
       (stable_id,
        node_id,
        u_from,
        u_to,
        evalue,
        score,
        type,
        qtaxon_id,
        set_id) VALUES (?,?,?,?,?,?,?,?,?)");

  my $evalue_count = 0;
  foreach my $stable_id (keys %{$self->{hits}}) {
    my $evalue = $self->{hits}{$stable_id}{Evalue};
    $DB::single=1;1;
    next if ($evalue > 0.01); #FIXME do we want to store everything?
    my $score = $self->{hits}{$stable_id}{Score};
    my $u_from = $self->{hits}{$stable_id}{u_from};
    my $u_to   = $self->{hits}{$stable_id}{u_to};
    next unless (defined($stable_id) && $stable_id ne '');
    next unless (defined($score));
    next unless ($evalue < $self->{max_evalue});
    $evalue_count++;
    $sth->execute($stable_id,
                  $node_id,
                  $u_from,
                  $u_to,
                  $evalue,
                  $score,
                  $type,
                  ($qtaxon_id*100),
                  $assembly);
    # FIXME: qtaxon_id*100 is to 
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

  foreach my $qtaxon_id (keys %{$self->{qtaxon_id_with_hits}}) {
    my $type = $self->{qtaxon_id_with_hits}{$qtaxon_id};
    # Dataflow clusters
    # This will create a new MSA alignment job for each of the newly generated clusters
    my $output_id = sprintf("{'protein_tree_id'=>%d, 'clusterset_id'=>%d, 'qtaxon_id'=>%d, 'type'=>\'%s\'}", 
                            $node_id, $clusterset_id, $qtaxon_id, $type);
    print STDERR "Dataflow # $output_id\n" if ($self->debug);
    $self->dataflow_output_id($output_id, 2);
  }

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

# sub run_search_u {
#   my $self = shift;

#   my $node_id = $self->{tree}->node_id;
#   my $type = $self->{type};
#   my $fasta_dir = $self->{fasta_dir};

# #   my $species_list;
# #   foreach my $leaf (@{$self->{tree}->get_all_leaves}) {
# #     $species_list->{$leaf->taxon_id};
# #   }
# #   if ($self->{search_mode} eq 'fast' && 2 > scalar keys %{$species_list}) {
# #     return undef;
# #   }

#   my $useqsfile = $self->{useqs};
#   my $resultsfile = $self->{useqs} . ".blast6out";

#   my $search_u_executable = $self->analysis->program_file;
#   unless (-e $search_u_executable) {
#     $search_u_executable = "/nfs/users/nfs_a/avilella/src/uclust/latest/uclust";
#   }

#   foreach my $genome_db (@{$self->{gdba}->fetch_all}) {
#     my $genebuild = $genome_db->genebuild;
#     if ((!defined($genebuild) || $genebuild eq '') && $type eq 'dna') {
#       push @{$self->{dna_gdb}}, $genome_db;
#     } elsif ((defined($genebuild) && $genebuild eq '') && $type eq 'aa') {
#       # FIXME something else
#     }
#   }

#   # Do the dna_gdb like this
#   foreach my $genome_db (@{$self->{dna_gdb}}) {
#     my $assembly = $genome_db->assembly;
#     my $locator = $genome_db->locator;
#     my $format = 'fasta';
#     my $qtaxon_id = $genome_db->taxon_id;
#     my $qdir = $self->get_qdir($locator);
#     my $fastafile = $qdir . "/" . "$assembly.$format";
#     if ($type eq 'aa') {
#       my $name = $genome_db->name; $name =~ s/\ /\_/g;
#       $fastafile = $qdir . "/" . $name . "_" . "$assembly.$format";
#     }

#     my $uclust_table = "usearch_" . $assembly;
#     my $create = $self->{comparaDBA}->dbc->prepare
#       ("CREATE TABLE IF NOT EXISTS $uclust_table like usearch");
#     $create->execute;
#     $create->finish;

#     my $sth = $self->{comparaDBA}->dbc->prepare
#       ("INSERT IGNORE INTO $uclust_table
#        (stable_id,member_id) VALUES (?,?)");

#     ## This is for references as query and reads as db
#     my $cmd = "$search_u_executable --maxlen 1000 --db $fastafile --query $useqsfile --uc /dev/stdout --id 0.80";
#     $cmd .= " 2>/dev/null" unless ($self->debug);
#     print STDERR "# $cmd\n" if ($self->debug);
#     my $starttime = time();
#     my $fh;
#     eval { open($fh, "$cmd |") || die $!; };
#     if ($@) {
#       warn("error running usearch $@ $!");
#       return;
#     }

#     while (<$fh>) {
#       next if ($_ =~ /^#/);
#       $DB::single=1;1;
#       next unless ($_ =~ /^H/);
#       my ($Type, $ClusterNr, $SeqLengthorClusterSize, $PctId, $Strand, $QueryStart, $SeedStart, $Alignment, $QueryLabel, $TargetLabel) = split("\t",$_);
#       # $sth->execute($QueryLabel, $TargetLabel);

#       ## This is for references as query and reads as db
#       ##     TargetLabel is the read_id (last column)
#       ##     QueryLabel is the reference sequence (ninth column)
#       ##     H       148260  957     81.3    +       0       0       520D75M362D     33231   HWUSI-EAS1599_0001:6:6:10267:4988#0/2
#       $sth->execute($TargetLabel,$QueryLabel);
#     }

# ## #                                                                            --- full sequence --- -------------- this domain -------------   u coord   ali coord   env coord
# ## # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
# ## #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
# ## SRR034143.6280281    -             75 1929                 -            822   0.00014   28.9   0.0   1   1   4.1e-11   0.00014   28.9   0.0   464   536     2    74     1    75 0.94 -
# ## SRR034143.5856316    -             75 1929                 -            822       5.6   13.7   0.0   1   1   1.6e-06       5.6   13.7   0.0   547   588    33    74    12    75 0.93 -
# ## ~
# ##

#     $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);
#     print STDERR scalar (keys %{$self->{hits}}), " hits - ",(time()-$starttime)," secs...\n";

#     $self->search_u_store_hits($qtaxon_id,$assembly);
#     delete $self->{hits};
#   }

#   return 1;
# }
