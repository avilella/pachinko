#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::Pachinkoalign

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $pachinkoalign = Bio::EnsEMBL::Compara::RunnableDB::Pachinkoalign->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$pachinkoalign->fetch_input(); #reads from DB
$pachinkoalign->run();
$pachinkoalign->output();
$pachinkoalign->write_output(); #writes to DB

=cut


=head1 DESCRIPTION

Does something, and it is something cool.

=cut


=head1 CONTACT

  Contact Albert Vilella on module implementation/design detail: avilella@ebi.ac.uk

=cut


=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Compara::RunnableDB::Pachinkoalign;

use strict;
use Getopt::Long;
use Time::HiRes qw(time gettimeofday tv_interval);
use Scalar::Util qw(looks_like_number);
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Compara::Graph::NewickParser;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive;

# #bwasw
# use Bio::DB::Sam;
# use Bio::Tools::Run::BWA;

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);

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
  $self->{ancestor_depth}  = 2;
  # Post august 2010 pagan deals better with Ns in 2x genomes
  $self->{no_N_edit} = 1;

  #create a Compara::DBAdaptor which shares the same DBI handle
  #with the Pipeline::DBAdaptor that is based into this runnable
  $self->{'comparaDBA'} = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
    (
     -DBCONN=>$self->db->dbc
    );

  # Get the needed adaptors here
  $self->{treeDBA} = $self->{'comparaDBA'}->get_ProteinTreeAdaptor;
  $self->{memberDBA} = $self->{'comparaDBA'}->get_MemberAdaptor;
  $self->{gdba} = $self->{'comparaDBA'}->get_GenomeDBAdaptor;
  $self->{mlssDBA} = $self->{'comparaDBA'}->get_MethodLinkSpeciesSetAdaptor;

  eval {$self->{mlssid} = @{$self->{mlssDBA}->fetch_all_by_method_link_type('PROTEIN_TREES')}->[0]->dbID;};
#  throw("no mlssid defined for PROTEIN_TREES") unless (0 < $self->{mlssid});
  $self->{mlssid} = 99999 unless (0 < $self->{mlssid});

  $self->get_params($self->parameters);
  $self->get_params($self->input_id);

# For long parameters, look at analysis_data
  if($self->{analysis_data_id}) {
    my $analysis_data_id = $self->{analysis_data_id};
    my $analysis_data_params = $self->db->get_AnalysisDataAdaptor->fetch_by_dbID($analysis_data_id);
    $self->get_params($analysis_data_params);
  }

  if(defined($self->{protein_tree_id})) {
    print STDERR "Fetching tree...\n" if ($self->debug);

    $self->{tree} = 
         $self->{treeDBA}->fetch_node_by_node_id($self->{protein_tree_id});
    printf("  protein_tree_id : %d\n", $self->{protein_tree_id}) if ($self->debug);
  }

  foreach my $gdb (@{$self->{gdba}->fetch_all}) {
    if (defined($gdb->genebuild) && $gdb->genebuild ne '') {
      $self->{species_set_hash}{$gdb->taxon_id} = 1;
    }
  }

  if ($self->{type} eq 'dna') {
    foreach my $genome_db (@{$self->{gdba}->fetch_all}) {
      my $genebuild = $genome_db->genebuild;
      if (!defined($genebuild) || $genebuild eq '') {
        my $assembly = $genome_db->assembly;
        $self->{qtaxon_list}->{$genome_db->taxon_id}{$assembly} = $genome_db->dbID;
      }
    }
  } elsif ($self->{type} eq 'aa') {
    foreach my $genome_db_id ( @{$self->{species_set}} ) {
      my $genome_db = $self->{gdba}->fetch_by_dbID($genome_db_id);
      my $genebuild = $genome_db->genebuild;
      if (0 < length($genebuild)) {
        my $assembly = $genome_db->assembly;
        $self->{qtaxon_list}->{$genome_db->taxon_id}{$assembly} = $genome_db->dbID;
      }
    }
  }
  if (defined ($self->{taxa}) && 0 < scalar @{$self->{taxa}}) {
    foreach my $taxon_id (@{$self->{taxa}}) {
      $self->{this_qtaxon_list}{$taxon_id} = 1;
    }
  } else {
    $self->{this_qtaxon_list} = $self->{qtaxon_list};
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

  foreach my $key (qw[search_mode ancestor_depth species_set gdb qtaxon_id fastq_dir tag velvet_exe bwa_exe abyss_exe samtools_exe prankancestor_exe cdbyank_exe exonerate_exe protein_tree_id clusterset_id type taxa analysis_data_id self_collapse]) {
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
  my $node_id = $self->{tree}->node_id;

  $self->run_pachinko;

  return;
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

  $self->dataflow_pachinko_tree;

  return;
}


##########################################
#
# internal methods
#
##########################################

sub run_pachinko {
  my $self = shift;

  print STDERR $self->worker_temp_directory, "\n" if ($self->debug);
  $self->{rectree_aa} = $self->run_prankancestor_pre('aa');
  $self->{rectree_cds} = $self->run_prankancestor_pre('cds');

  $self->run_pachinkoplacement;
  return undef if (1 == $self->{done});

  $self->define_read_placement_sets;
  foreach my $read_set_id (keys %{$self->{read_placement_sets}}) {
    my $pre_abyssfile = $self->run_preassemble_abyss($read_set_id);
    $self->run_bwasw('pre',$read_set_id,$self->{pre_abyss_dir},$pre_abyssfile) if (defined($pre_abyssfile));
    $self->run_pachinkoalign($read_set_id,$self->{consfile});
    return 1 unless($self->parse_alignment);
    $self->basic_stats;
    $self->indel_stats;
    return 1 unless defined($self->{pachinko_aln});
    $self->run_prankancestor;
    $self->{rectree} = $self->parse_prankancestor($self->{ancfile},'final');
    $self->run_this_columbus('post',$read_set_id);
  }

  return 0;
}

sub run_pachinkoplacement {
  my $self = shift;

  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';
  my $type = $self->{type};
  my $fastq_dir = $self->{fastq_dir};
  # FIXME fetch_reads_list_uclust is old system, use fetch_reads_list_usearch
  $self->fetch_reads_list_uclust      if ($self->{comparaDBA}->dbc->dbname eq 'avilella_p_58m');
  $self->fetch_reads_list_usearch unless ($self->{comparaDBA}->dbc->dbname eq 'avilella_p_58m');
  print STDERR "Fetching reads list...\n" if ($self->debug);
  # FIXME: disabled hmmsearch fetching
  $self->fetch_reads_list unless (defined($self->{search_mode}) && $self->{search_mode} =~ /hmm/);
  print STDERR "Fetching reads temp file...\n" if ($self->debug);
  my $readsfilename = $self->fetch_reads_temp_file;

  # No reads, nothing else to do
  if (-z $readsfilename) {
    print STDERR "No reads found for this tree. Finishing...\n" if ($self->debug);
    $self->{done} = 1;
    return undef;
  }
  $self->{readsfilename} = $readsfilename;

  # FETCHING THE REFERENCE CDS ALIGNMENT
  # It may have been re-run so we need to clean up existing attached nodes
  my @preclean;
  foreach my $leaf (@{$self->{tree}->get_all_leaves}) {
#     if ($leaf->stable_id =~ /$tag/ || $leaf->stable_id eq '' || $leaf->source_name eq 'EXTERNALCDS') {
    if (!defined $self->{species_set_hash}{$leaf->taxon_id}) {
      push @preclean, $leaf;
      # self_collapse option is only for testing purposes
    } elsif (1 == $self->{self_collapse} && defined($self->{qtaxon_list}{$leaf->taxon_id})) {
      push @preclean, $leaf;
    }
  }

  $self->{tree} = $self->{tree}->remove_nodes(\@preclean);
  foreach my $node (@preclean) {
    if ($node->taxon_id == $self->{qtaxon_id} || 0 == $node->taxon_id || $node->description == $self->{qtaxon_id}) {
      # Only delete the nodes for the same qtaxon_id as in this job
      $self->delete_old_pachinko_member($node);
    }
  }

  # We edit all NNNs to gaps so that they dont confuse the parsimony
  # process
  unless (defined($self->{no_N_edit})) {
    foreach my $leaf (@{$self->{tree}->get_all_leaves}) {
      my $cdna;
      eval { $cdna = $leaf->cdna_alignment_string; };
      if ($@) {
        my @clean; push @clean, $leaf;
        $self->{tree} = $self->{tree}->remove_nodes(\@clean);
        next;
        # throw("error fetching cdna sequence\n");
      }
      $cdna =~ s/N/\-/g;
      $leaf->{'cdna_alignment_string'} = $cdna;
      my $composed_id = $leaf->member_id . "_" . $leaf->taxon_id;
      $self->{_edited_refseqs}{$composed_id} = $cdna;
    }
  }

  print STDERR "fetching simplealign...\n" if ($self->debug);
  $self->{gtsa} = $self->{tree}->get_SimpleAlign
    (
     -id_type => 'MEMBER',
     -cdna => 1,
     -stop2x => 1,
     -append_taxon_id => 1,
    );
  $self->{gtsa}->set_displayname_flat(1);
  $self->{gtsa} = $self->{gtsa}->remove_gaps(undef,1);
  $self->{gtsa}->set_displayname_flat(1);

  my $cds_seqfile = $self->worker_temp_directory . "$node_id.cds";

  my $alignIO = Bio::AlignIO->new
    (-file => ">$cds_seqfile",
     -format => 'fasta'
    );
  $alignIO->write_aln($self->{gtsa});
  $alignIO->close;
  $self->{cds_seqfile} = $cds_seqfile;
  ########################################

  # FETCHING THE REFERENCE TREES
  my $cds_treefile = $self->worker_temp_directory . "$node_id.nh";
  my $newick_string = $self->{tree}->newick_format('member_id_taxon_id');
  open NEWICK, ">$cds_treefile" or die "$!";
  print NEWICK "$newick_string\n";
  close NEWICK;
  $self->{cds_treefile} = $cds_treefile;

  # Some alignments can have missing sequences if the core db has
  # changed place, so we trim the nhx accordingly
  my @alnclean;
  foreach my $leaf (@{$self->{tree}->get_all_leaves}) {
    my $id = $leaf->member_id . "_" . $leaf->taxon_id;
    my $alnleaf = $self->{gtsa}->get_seq_by_id($id);
    if (!defined($alnleaf)) {
      push @alnclean, $leaf;
    }
  }
  $self->{tree} = $self->{tree}->remove_nodes(\@alnclean);

  my $nhx_treefile = $self->worker_temp_directory . "$node_id.nhx";
  if (-e $nhx_treefile && !-z $nhx_treefile) {
    # we've already produced the nhx file once, no need to repeat
  } else {
    my $nhx_string = $self->{tree}->nhx_format('member_id_taxon_id');
    open NHX, ">$nhx_treefile" or die "$!";
    print NHX "$nhx_string\n";
    close NHX;
  }
  $self->{nhx_treefile} = $nhx_treefile;
  ########################################
  # PLACEMENT (1 STEP)
  my $pachinko_executable = $self->analysis->program_file;
  unless (-e $pachinko_executable) {
    $pachinko_executable = "/nfs/users/nfs_a/avilella/src/pagan/latest/pagan/pagan";
  }

  my $placement_file_prefix = $self->worker_temp_directory . $node_id . ".pagan.placement";
  my $placement_file = $placement_file_prefix . ".tsv";
  my $outfile_prefix = $self->worker_temp_directory . $node_id . ".pagan";
  my $cmd = "$pachinko_executable ";
  # This is PACHINKOPLACEMENT (1 STEP)
  # OPTIONS and INPUTFILES
  # $cmd .= " --trim-read-ends --rank-reads-for-nodes"; # we dont assert reads for placement
  $cmd .= " --pair-end ";  # PEs
  $cmd .= " --output-nhx-tree";
  $cmd .= " --454 " if (1 == $self->{454});  # longreads
  $cmd .= " --cds-seqfile $cds_seqfile";
  $cmd .= " --cds-treefile $nhx_treefile";
  $cmd .= " --readsfile $readsfilename";
  $cmd .= " --outfile $outfile_prefix";
  $cmd .= " --placement-only --placement-file $placement_file_prefix";
  $cmd .= " --silent" unless (2 == $self->debug);

  # $cmd .= " --cds-treefile $cds_treefile"; # This used to be nh, now we are using NHX to label TIDs
  # $cmd .= " --discard-overlapping-identical-reads"; # Be careful, overlapping will disappear from alignment and be reintroduced in sam file
  # $cmd .= " --min-reads-overlap 0.2 "; # default 0.1
  # $cmd .= " --reads-distance 0.5 ";  # default 0.5

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  my $starttime = time();
  print STDERR "pachinkoplacement...\n" if ($self->debug);
  my $worker_temp_directory = $self->worker_temp_directory;
  print STDERR "# $cmd\n" if ($self->debug);

  unless(system("cd $worker_temp_directory; $cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running pachinkoplacement pagan, $!\n");
  }

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);
  throw ("Couldnt find placement file: $placement_file") if (-z $placement_file);
  $self->{placement_file} = $placement_file;

  return 1;
}

sub define_read_placement_sets {
  my $self = shift;
  my $placement_file = $self->{placement_file};

  if (-z $placement_file) {
    $self->{read_placement_sets}{1} = 1;
    #my @tids = values %{$self->{tid_qtaxon_id}};
    #$self->{tree}->store_tag('TID',$tids[0]);
    $self->{to_root} = 1;
  } else {
    open PLACEMENTFILE, $placement_file or die $!;
    while (<PLACEMENTFILE>) {
      chomp $_;
      my ($stable_id,$read_set_id) = split(" ",$_);
      $self->{read_placement_sets}{$read_set_id}{$stable_id} = 1;
    }
    close PLACEMENTFILE;
  }
}

sub run_preassemble {
  my $self = shift;
  my $read_set_id = shift;

  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';
  my $type = $self->{type};
  my $fastq_dir = $self->{fastq_dir};

  my $cds_seqfile  =  $self->{cds_seqfile}  ;
  my $cds_treefile =  $self->{cds_treefile} ;
  my $nhx_treefile =  $self->{nhx_treefile} ;

  my $fullreadsfilename = $self->{readsfilename};

  # No reads, nothing else to do
  return undef if (-z $fullreadsfilename);

  # WRITE DOWN READ SET FROM MAIN READ FILE
  my $format = 'fastq';
  my $readsfilename = $self->worker_temp_directory . $node_id . ".$read_set_id" . ".reads";
  my $read_setout = Bio::SeqIO->new
    (-format => $format,
     -file => ">$readsfilename");
  my $fastqio = Bio::SeqIO->new(-file => "<$fullreadsfilename", -format => $format);
  while (my $seq  = $fastqio->next_seq) {
    my $stable_id = $seq->display_id;
    next if (!defined $self->{to_root} && !defined ($self->{read_placement_sets}{$read_set_id}{$stable_id}));
    $read_setout->write_seq($seq);
  }
  $read_setout->close;
  ########################################

  my $velvet_executable = $self->{velvet_exe};
  unless (-e $velvet_executable) {
    $velvet_executable = "/nfs/users/nfs_a/avilella/src/velvet/latest/velvet/velveth"
  }

  my $velvet_dir = $self->worker_temp_directory . "$read_set_id.velvet";
  my $cmd = "$velvet_executable $velvet_dir 63 -fastq $readsfilename";
  print("$cmd\n");
  unless(system("$cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running velveth, $!\n");
  }

  $velvet_executable =~ s/velveth/velvetg/;
  $cmd = "$velvet_executable $velvet_dir -exp_cov 1 -unused_reads yes";
  print("$cmd\n");

  unless(system("$cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running velvetg, $!\n");
  }

  my $velvetfile   = $velvet_dir . "/" . "contigs.fa";
  my $velvetunused = $velvet_dir . "/" . "UnusedReads.fa";

  my $final_velvet_file = $velvet_dir . "/" . "velvet_final.fa";
  my $final_velvet_io = Bio::SeqIO->new(-file => ">$final_velvet_file",
     -format => 'fasta');
  my $velvetfile_io = Bio::SeqIO->new(-file => "$velvetfile",
     -format => 'fasta');
  my $velvetunused_io = Bio::SeqIO->new(-file => "$velvetunused",
     -format => 'fasta');
  while (my $seq  = $velvetfile_io->next_seq) { $final_velvet_io->write_seq($seq); }
  while (my $seq  = $velvetunused_io->next_seq) { $final_velvet_io->write_seq($seq); }
  if (-z $final_velvet_file) {
    print STDERR "no contigs produced\n";
      return undef;
  } else {
    $self->{velvetfile}   = $final_velvet_file;
    $self->{pre_velvet_dir} = $velvet_dir;
    return $self->{velvetfile};
  }

  return 0;
}


sub run_preassemble_abyss {
  my $self = shift;
  my $read_set_id = shift;

  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';
  my $type = $self->{type};
  my $fastq_dir = $self->{fastq_dir};

  my $cds_seqfile  =  $self->{cds_seqfile}  ;
  my $cds_treefile =  $self->{cds_treefile} ;
  my $nhx_treefile =  $self->{nhx_treefile} ;

  my $fullreadsfilename = $self->{readsfilename};

  # No reads, nothing else to do
  return undef if (-z $fullreadsfilename);

  # WRITE DOWN READ SET FROM MAIN READ FILE
  my $format = 'fastq';
  my $readsfilename = $self->worker_temp_directory . $node_id . ".$read_set_id" . ".reads";
  my $read_setout = Bio::SeqIO->new
    (-format => $format,
     -file => ">$readsfilename");
  my $fastqio = Bio::SeqIO->new(-file => "<$fullreadsfilename", -format => $format);
  while (my $seq  = $fastqio->next_seq) {
    my $stable_id = $seq->display_id;
    next if (!defined $self->{to_root} && !defined ($self->{read_placement_sets}{$read_set_id}{$stable_id}));
    $read_setout->write_seq($seq);
  }
  $read_setout->close;
  ########################################

  my $abyss_executable = $self->{abyss_exe};
  unless (-e $abyss_executable) {
    $abyss_executable = "/nfs/users/nfs_a/avilella/src/abyss/latest/abyss-1.1.2/ABYSS/ABYSS";
  }

  #/nfs/users/nfs_a/avilella/src/abyss/latest/abyss-1.1.2/ABYSS/ABYSS -k$kmer -e0 -c0 $reads -o $contigs
  my $abyss_dir = $self->worker_temp_directory;
  my $abyssfile = $abyss_dir . "contigs.fa";
  my $cmd = "$abyss_executable -k23 -e0 -c0 $readsfilename -o $abyssfile";
  print("$cmd\n");
  unless(system("$cmd") == 0) {
    print("$cmd\n");
    warn("abyss returned an error, could be no contigs produced $!\n");
  }

  my $abyssunused = $readsfilename;

  my $final_abyss_file = $abyss_dir . "/" . "abyss_final.fa";
  my $final_abyss_io = Bio::SeqIO->new(-file => ">$final_abyss_file",
     -format => 'fasta');
  my $abyssfile_io = Bio::SeqIO->new(-file => "$abyssfile",
     -format => 'fasta');
  my $abyssunused_io = Bio::SeqIO->new(-file => "$abyssunused",
     -format => 'fastq');
  while (my $seq  = $abyssfile_io->next_seq) { 
    my $seqname = "NODE_ " . $seq->display_id . "_" . $seq->description;
    $seq->description('');
    $seqname =~ s/\ /\_/g; $seq->display_id($seqname);
    $final_abyss_io->write_seq($seq);
  }
  while (my $seq  = $abyssunused_io->next_seq) { $final_abyss_io->write_seq($seq); }
  if (-z $final_abyss_file) {
    print STDERR "no contigs produced\n";
      return undef;
  } else {
    $self->{abyssfile}   = $final_abyss_file;
    $self->{pre_abyss_dir} = $abyss_dir;
    return $self->{abyssfile};
  }

  return 0;
}

sub run_pachinkoalign {
  my $self = shift;
  my $read_set_id = shift;
  my $specific_readsfile = shift;

  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';
  my $type = $self->{type};
  my $fastq_dir = $self->{fastq_dir};

  my $cds_seqfile  =  $self->{cds_seqfile}  ;
  my $cds_treefile =  $self->{cds_treefile} ;
  my $nhx_treefile =  $self->{nhx_treefile} ;

  my $fullreadsfilename = $self->{readsfilename};

  # No reads, nothing else to do
  return undef if (-z $fullreadsfilename);

  # WRITE DOWN READ SET FROM MAIN READ FILE
  my $format = 'fastq';
  my $readsfilename;

  if (!defined $specific_readsfile) {
    $readsfilename = $self->worker_temp_directory . $node_id . ".$read_set_id" . ".reads";
    my $read_setout = Bio::SeqIO->new
      (-format => $format,
       -file => ">$readsfilename");
    my $fastqio = Bio::SeqIO->new(-file => "<$fullreadsfilename", -format => $format);
    while (my $seq  = $fastqio->next_seq) {
      my $stable_id = $seq->display_id;
      next if (!defined $self->{to_root} && !defined ($self->{read_placement_sets}{$read_set_id}{$stable_id}));
      my @values = values %{$self->{tid_qtaxon_id}};
      my $new_description =  $seq->description; $new_description .= " TID=" . $values[0] if (1 == scalar @values);
      $seq->description($new_description);
      $read_setout->write_seq($seq);
    }
    $read_setout->close;
  } else {
    $readsfilename = $self->worker_temp_directory . $node_id . ".$read_set_id" . ".reads";
    my $read_setout = Bio::SeqIO->new
      (-format => $format,
       -file => ">$readsfilename");
    my $fastqio = Bio::SeqIO->new(-file => "<$specific_readsfile", -format => $format);
    while (my $seq  = $fastqio->next_seq) {
      my $stable_id = $seq->display_id;

      my @values = values %{$self->{tid_qtaxon_id}};
      my $new_description =  $seq->description; $new_description .= " TID=" . $values[0] if (1 == scalar @values);
      $seq->description($new_description);
      $read_setout->write_seq($seq);
    }
    $read_setout->close;
  }
  ########################################

#   my $nhx_string = $self->{tree}->nhx_format('member_id_taxon_id');
#   open NHX, ">$nhx_treefile" or die "$!";
#   print NHX "$nhx_string\n";
#   close NHX;
#   $self->{nhx_treefile} = $nhx_treefile;
#   ########################################
  # PACHINKOALIGN (2 STEP)
  my $pachinko_executable = $self->analysis->program_file;
  unless (-e $pachinko_executable) {
    $pachinko_executable = "/nfs/users/nfs_a/avilella/src/pagan/latest/pagan/pagan";
  }

  # This is PACHINKOALIGN (2 STEP)
  my $outfile_prefix = $self->worker_temp_directory . $node_id . ".pagan";
  my $cmd = "$pachinko_executable ";
  # OPTIONS and INPUTFILES
  # specific to PACHINKO (2 STEP), we dont use trim-read-ends
  $cmd .= " --rank-reads-for-nodes";
#  $cmd .= " --output-ancestors";
  $cmd .= " --pair-end ";  # PEs
  $cmd .= " --output-nhx-tree";
  $cmd .= " --454 " if (1 == $self->{454});  # longreads
  $cmd .= " --cds-seqfile $cds_seqfile";
  $cmd .= " --cds-treefile $nhx_treefile";
  $cmd .= " --readsfile $readsfilename";
  $cmd .= " --outfile $outfile_prefix";
  $cmd .= " --silent" unless (2 == $self->debug);
  # $cmd .= " --cds-treefile $cds_treefile"; # This used to be nh, now we are using NHX to label TIDs
  # $cmd .= " --discard-overlapping-identical-reads"; # Be careful, overlapping will disappear from alignment and be reintroduced in sam file
  # $cmd .= " --min-reads-overlap 0.2 "; # default 0.1
  # $cmd .= " --reads-distance 0.5 ";  # default 0.5

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  my $starttime = time();
  print STDERR "pachinkoalign...\n" if ($self->debug);
  my $worker_temp_directory = $self->worker_temp_directory;
  print STDERR "# $cmd\n" if ($self->debug);

  unless(system("cd $worker_temp_directory; $cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running pachinkoalign pagan, $!\n");
  }

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

  $self->{outfile} = "$worker_temp_directory"."$node_id.pagan.fas";
  $self->{xmlfile} = "$worker_temp_directory"."$node_id.pagan.xml";
  $self->{nhx_tree_file} = "$worker_temp_directory"."$node_id.pagan.nhx_tree";

  return 1;
}

sub fetch_reads_list {
  my $self = shift;
  my $type = $self->{type};
  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';

  my $qtaxon_string = join(",",keys %{$self->{this_qtaxon_list}});
#  my $sql1 = "SELECT stable_id, hmm_from, hmm_to, set_id FROM hmmsearch WHERE type=\"$type\" AND node_id=$node_id AND qtaxon_id in ($qtaxon_string)";
  # FIXME
  # Get the whole lot instead, it will be faster if no index is in place
  # REPLACE is used to mapback
  my $sql1 = "SELECT stable_id, hmm_from, hmm_to, set_id, score, evalue FROM hmmsearch WHERE type=\"$type\" AND node_id=$node_id";
  # print STDERR "$sql1\n" if ($self->debug);
  my $sth1 = $self->{comparaDBA}->dbc->prepare($sql1);
  $sth1->execute;
  my $result1;
  my $sql2 = "SELECT node_id, qtaxon_id, score FROM hmmsearch WHERE type=\"$type\" AND stable_id=?";
  my $sth2 = $self->{comparaDBA}->dbc->prepare($sql2);

  while ($result1 = $sth1->fetchrow_hashref) {
    if ($result1->{set_id} =~ /\_([12])$/) {
      my $pe_stable_id = $result1->{stable_id} . "/$1";
      $pe_stable_id =~ s/\:/\@/g; #mapback

      $self->{reads_pe}{$pe_stable_id}{hmm_from} = $result1->{hmm_from};
      $self->{reads_pe}{$pe_stable_id}{hmm_to} = $result1->{hmm_to};
      $self->{reads_pe}{$pe_stable_id}{score} = $result1->{score};
      $self->{reads_pe}{$pe_stable_id}{evalue} = $result1->{evalue};
    }
    $sth2->execute($result1->{stable_id});
    my $result2;
    my $max_score = 0; my $max_node_id = undef;
    my $qtaxon_id = undef;
    while ($result2 = $sth2->fetchrow_hashref) {
      $qtaxon_id = $result2->{qtaxon_id};
      if ($result2->{score} > $max_score) {
        $max_score = $result2->{score}; $max_node_id = $result2->{node_id};
      }
    }
    $sth2->finish;
    my $mapback_stable_id = $result1->{stable_id}; $mapback_stable_id =~ s/\:/\@/g;
    $self->{reads_list}{$qtaxon_id}{$mapback_stable_id} = $max_score if ($max_node_id == $node_id);
    $self->{reads_taxon}{$mapback_stable_id} = $qtaxon_id;
    # $self->{reads_list}{$qtaxon_id}{$result1->{stable_id}} = 1 if ($max_node_id == $node_id);
  }
  $sth1->finish;

  # We do a last pass to fetch well-scored PEs
  foreach my $remaining_stable_id (keys %{$self->{reads_taxon}}) {
    my $qtaxon_id = $self->{reads_taxon}{$remaining_stable_id};
    next if (defined $self->{reads_list}{$qtaxon_id}{$remaining_stable_id});
    my $pe_stable_id1 = $remaining_stable_id . "/1";
    my $pe_stable_id2 = $remaining_stable_id . "/2";
#     if (
#         defined($self->{reads_pe}{$pe_stable_id1}) 
#         && defined($self->{reads_pe}{$pe_stable_id2}) 
#         && $self->{reads_pe}{$pe_stable_id1}{evalue} < 0.001
#         && $self->{reads_pe}{$pe_stable_id2}{evalue} < 0.001
#        ) {
    # no evalue filtering
    if (
        defined($self->{reads_pe}{$pe_stable_id1}) 
        && defined($self->{reads_pe}{$pe_stable_id2}) 
       ) {
      my $combined_avg_score = ($self->{reads_pe}{$pe_stable_id1}{score} + $self->{reads_pe}{$pe_stable_id2}{score})/2;
      $self->{reads_list}{$qtaxon_id}{$remaining_stable_id} = $combined_avg_score;
      print STDERR "$remaining_stable_id, $combined_avg_score - remaining_stable_id\n" if ($self->debug);
      $self->_store_tagvalue($remaining_stable_id,$node_id,'remaining_pe',$combined_avg_score);
    }
  }
  print STDERR "Reads: ", scalar keys %{$self->{reads_list}{$self->{qtaxon_id}}}, "\n";

  return 1;
}

sub fetch_reads_list_uclust {
  my $self = shift;
  my $type = $self->{type};
  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';
  my $starttime = time();

  foreach my $genome_db (@{$self->{gdba}->fetch_all}) {
    next unless (defined($self->{this_qtaxon_list}{$genome_db->taxon_id}));
    my $qtaxon_id = $genome_db->taxon_id;
    print STDERR "Assembly: ", $genome_db->assembly , " ", $genome_db->taxon_id, "\n";
    my $member_list = join(",",map {$_->member_id} @{$self->{tree}->get_all_leaves});

    my $assembly = $genome_db->assembly;
    my $table = "uclust_member_" . $assembly;
    my $sql1 = "SELECT stable_id, member_id FROM $table WHERE member_id in ($member_list)";
    print STDERR "$sql1\n" if ($self->debug);
    my $sth1 = $self->{comparaDBA}->dbc->prepare($sql1);
    eval { $sth1->execute; };
    printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
    next if ($@);
    my $result1;

    while ($result1 = $sth1->fetchrow_hashref) {
      my $pe_stable_id = $result1->{stable_id} . "/$1";
      if ($assembly =~ /\_([12])$/) {
        my $pe_stable_id = $result1->{stable_id} . "/$1";
      }
      $pe_stable_id =~ s/\:/\@/g; #mapback

      $self->{reads_pe}{$pe_stable_id}{hmm_from} = 0;
      $self->{reads_pe}{$pe_stable_id}{hmm_to} = 0;
      $self->{reads_pe}{$pe_stable_id}{score} = 0;
      $self->{reads_pe}{$pe_stable_id}{evalue} = 0;
      my $mapback_stable_id = $result1->{stable_id}; $mapback_stable_id =~ s/\:/\@/g;
      $self->{reads_list}{$self->{qtaxon_id}}{$mapback_stable_id} = 40;
      $self->{reads_taxon}{$mapback_stable_id} = $self->{qtaxon_id};
      # $self->{reads_list}{$qtaxon_id}{$result1->{stable_id}} = 1 if ($max_node_id == $node_id);
    }
    $sth1->finish;

    print STDERR "Reads after uclust: ", scalar keys %{$self->{reads_list}{$self->{qtaxon_id}}}, "\n";
    $self->{tree}->store_tag('num_uclust_reads',scalar keys %{$self->{reads_list}{$self->{qtaxon_id}}});
  }

  return 1;
}

sub fetch_reads_list_usearch {
  my $self = shift;
  my $type = $self->{type};
  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';
  my $starttime = time();

  foreach my $genome_db (@{$self->{gdba}->fetch_all}) {
    next unless (defined($self->{this_qtaxon_list}{$genome_db->taxon_id}));
    my $qtaxon_id = $genome_db->taxon_id;
    print STDERR "Assembly: ", $genome_db->assembly , " ", $genome_db->taxon_id, "\n";

    my $assembly = $genome_db->assembly;
    my $table = "usearch_" . $assembly;
    my $sql1 = "SELECT stable_id FROM $table WHERE root_id=$node_id";
    print STDERR "$sql1\n" if ($self->debug);
    my $sth1 = $self->{comparaDBA}->dbc->prepare($sql1);
    eval { $sth1->execute; };
    printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
    next if ($@);
    my $result1;

    while ($result1 = $sth1->fetchrow_hashref) {
      my $pe_stable_id = $result1->{stable_id} . "/$1";
      if ($assembly =~ /\_([12])$/) {
        my $pe_stable_id = $result1->{stable_id} . "/$1";
      }
      $pe_stable_id =~ s/\:/\@/g; #mapback

      $self->{reads_pe}{$pe_stable_id}{hmm_from} = 0;
      $self->{reads_pe}{$pe_stable_id}{hmm_to} = 0;
      $self->{reads_pe}{$pe_stable_id}{score} = 0;
      $self->{reads_pe}{$pe_stable_id}{evalue} = 0;
      my $mapback_stable_id = $result1->{stable_id}; $mapback_stable_id =~ s/\:/\@/g;
      $self->{reads_list}{$self->{qtaxon_id}}{$mapback_stable_id} = 40;
      $self->{reads_taxon}{$mapback_stable_id} = $self->{qtaxon_id};
      # $self->{reads_list}{$qtaxon_id}{$result1->{stable_id}} = 1 if ($max_node_id == $node_id);
    }
    $sth1->finish;

    print STDERR "Reads after uclust: ", scalar keys %{$self->{reads_list}{$self->{qtaxon_id}}}, "\n";
    $self->{tree}->store_tag('num_uclust_reads',scalar keys %{$self->{reads_list}{$self->{qtaxon_id}}});
  }

  return 1;
}


sub fetch_reads_temp_file {
  my $self = shift;

  my $starttime = time();
  print STDERR "fetch_reads_temp_file...\n" if ($self->debug);
  eval {
    require Bio::Index::Fasta;
    require Bio::Index::Fastq;
    require Bio::SeqIO;
      };
  $self->throw("No Bio::DB::Fasta") if ($@);

  my $sql1 = "select value from protein_tree_tag where tag='species_tree_string'";
  my $sth1 = $self->dbc->prepare($sql1);
  $sth1->execute;
  my $species_tree_string = $sth1->fetchrow_hashref;
  $sth1->finish;
  my $eval_species_tree;
  eval {
    $eval_species_tree = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($species_tree_string->{value});
    my @leaves = @{$eval_species_tree->get_all_leaves};
  };

  if($@) {
    unless(-e $self->{'species_tree_file'}) {
      throw("can't find species_tree\n");
    }
  } else {
    $self->{species_tree_string} = $species_tree_string->{value};
    my $spfilename = $self->worker_temp_directory . "spec_tax.nh";
    open SPECIESTREE, ">$spfilename" or die "$!";
    print SPECIESTREE $self->{species_tree_string};
    close SPECIESTREE;
    $self->{'species_tree_file'} = $spfilename;
  }

  my $taxon_tree = $self->load_species_tree_from_file;
  $self->{'taxon_tree'} = $taxon_tree;

  print STDERR $self->worker_temp_directory, "\n" if ($self->debug);

  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';
  my $gdbs = $self->{gdb};

  my $readsfilename = $self->worker_temp_directory . $node_id . ".reads";
  # my $format = 'fasta';
  my $format = 'fastq';
  my $out = Bio::SeqIO->new(-format => $format,
                            -file => ">$readsfilename");

  my $any_tid = undef;
  foreach my $genome_db (@{$self->{gdba}->fetch_all}) {
    next unless (defined($self->{this_qtaxon_list}{$genome_db->taxon_id}));
    my $qtaxon_id = $genome_db->taxon_id;
    print STDERR "Assembly: ", $genome_db->assembly , " ", $genome_db->taxon_id, "\n";
    # Find the closest internal node where this qtaxon could belong
    my $taxon_leaf = $self->{taxon_tree}->find_leaf_by_node_id($qtaxon_id);
    throw("species tree doesnt contain $qtaxon_id") if (!defined($taxon_leaf));
    # Traverse from the leaf to the root of the species tree, then try
    # to find the taxon level in the reference genetree
    my $taxon_id = $taxon_leaf->node_id;
    my $parent = $taxon_leaf->parent;
    my $found_tid = 0;
    while ($parent->has_parent) {
      last if (defined $self->{tid_qtaxon_id}{$qtaxon_id});
      last if ($qtaxon_id == $found_tid);
      my @ref_subnodes = @{$self->{tree}->get_all_nodes_by_tag_value('taxon_id',$parent->node_id)};
      $parent = $parent->parent;

      # We found a node in the reference genetree that is the smallest
      # taxon level in the species tree for the new taxon
      next unless (defined(@ref_subnodes) && 0 < scalar @ref_subnodes);
      foreach my $ref_subnode (@ref_subnodes) {
        my $tid_value = $ref_subnode->get_tagvalue('taxon_id');
        print STDERR "Adding tag: $tid_value / ", $ref_subnode->node_id, " ", $ref_subnode->get_tagvalue('taxon_name'), "\n" if ($self->debug);
        $found_tid = $qtaxon_id;
        $ref_subnode->store_tag('TID',$tid_value); $any_tid = 1;
        $self->{tid_qtaxon_id}{$qtaxon_id} = $tid_value;
      }
    }

    # Rank ordered by best scores first
    my @stable_ids = sort {$self->{reads_list}{$qtaxon_id}{$b} <=> $self->{reads_list}{$qtaxon_id}{$a}} keys %{$self->{reads_list}{$qtaxon_id}};
    next unless (0 < scalar @stable_ids);

    if ($self->{type} eq 'dna') {
      my $genebuild = $genome_db->genebuild;
      next if (0 < length($genebuild));
      my $assembly = $genome_db->assembly;
      my $locator = $genome_db->locator;

      my $format = 'fastq';
      my $qdir = $self->get_qdir($locator);

      my $fastafile = $qdir . "/" . "$assembly.$format";
#       my $indexfile = $qdir . "/" . "$assembly.$format.z.cidx";
#       throw ("unable to find $indexfile") unless (-e $indexfile);
#       print STDERR "  [$qtaxon_id]\t$indexfile\n" if ($self->debug);

      my $count = 0;
      my $worker_temp_directory = $self->worker_temp_directory;
      my $readlist = $worker_temp_directory . $assembly . ".readlist";
      open READLIST,">$readlist" or die ($!);
      foreach my $stable_id (@stable_ids) {
        $stable_id =~ s/\@/\:/g; # mapback
        my $discard = $self->illumina_same_set($assembly,$stable_id);

        # Delete the tag from the name so that it matches the "-r
        # @$tag" in the cidx index
        $stable_id =~ s/$tag//;
        print READLIST "$stable_id\n" unless (1 == $discard);
      }
      close READLIST;
      next if (-z $readlist);
      my $tmpfastq = $worker_temp_directory . $assembly . ".tmp";
      my $cmd;
      my $chunk_id = 0;
      my $indexfile = $qdir . "/" . $assembly . ".$format.z_" . sprintf("%07d",$chunk_id) . ".cidx";
      while (-e $indexfile) {
        throw ("unable to find $indexfile") unless (-e $indexfile);
        print STDERR "  [$qtaxon_id]\t$indexfile\n" if ($self->debug);
        $cmd = $self->{cdbyank_exe} || "/nfs/users/nfs_a/avilella/src/cdbfasta/latest/cdbfasta/cdbyank";
        $cmd   .= " $indexfile <$readlist >> $tmpfastq";

        print STDERR "$cmd\n" if ($self->debug);
        unless(system("$cmd") == 0) {
          print("$cmd\n");
          $self->throw("error running cdbyank $!\n");
        }
        printf(STDERR "  %1.3f secs\n", (time()-$starttime)) if ($self->debug);

        $chunk_id++;
        $indexfile = $qdir . "/" . $assembly . ".$format.z_" . sprintf("%07d",$chunk_id) . ".cidx";
      }
      my $fastqio = Bio::SeqIO->new(-file => "<$tmpfastq", -format => $format);
      while (my $seq  = $fastqio->next_seq) {
        my $stable_id = $seq->display_id;
        my $desc = $seq->desc;
        $stable_id =~ s/\:/\@/g; # mapback
        my $full_id = $stable_id . $desc;
        my $length = length($seq->seq);
        my $score = int($self->{reads_list}{$qtaxon_id}{$stable_id});

        $self->{coef_score}{$stable_id} = $score/$length;
        $desc .= " gtscore=$score";
        if (defined($self->{tid_qtaxon_id}{$qtaxon_id})) {
          $desc .= " TID=" . $self->{tid_qtaxon_id}{$qtaxon_id};
        }
        $seq->desc($desc);
        $self->{read_objs}{$score}{$full_id}{$assembly} = $seq;
        $count++;
        print STDERR "    [$count]\n" if (0 == $count % 10 && $self->debug);
      }
      $self->{tmpfastqs}{$tmpfastq} = 1;
      unlink($tmpfastq);
    } elsif ($self->{type} eq 'aa') {
      # Do something
    }
  }
  if (!defined $any_tid) { $self->{tree}->store_tag('TID', $self->{tree}->get_tagvalue('taxon_id')); }

  my $num_input_reads = 0;
  my $pairs;

  foreach my $score (sort {$b <=> $a} keys %{$self->{read_objs}}) {
    foreach my $full_id (keys %{$self->{read_objs}{$score}}) {
      foreach my $assembly (keys %{$self->{read_objs}{$score}{$full_id}}) {
        my $seq = $self->{read_objs}{$score}{$full_id}{$assembly};
        $num_input_reads++;
        my $flip = 'fwd_pair';
        if ($assembly =~ /\_([12])$/) {
          my $pe_id1 = $seq->display_id."/1";
          my $pe_id2 = $seq->display_id."/2";
          if (defined($self->{reads_pe}{$pe_id1}) && defined ($self->{reads_pe}{$pe_id2})) {
            my $from1 = $self->{reads_pe}{$pe_id1}{hmm_from};
            my $to1   = $self->{reads_pe}{$pe_id1}{hmm_to};
            my $from2 = $self->{reads_pe}{$pe_id2}{hmm_from};
            my $to2   = $self->{reads_pe}{$pe_id2}{hmm_to};
            if ($from1 > $from2 && $to1 > $to2) {
              my $temp = $from1; $from1 = $from2; $from2 = $temp;
              $temp = $to1; $to1 = $to2; $to2 = $temp;
              $flip = 'rev_pair';
            } else {$flip = 'fwd_pair';}
            if ($from1 < $from2 && $to1 < $to2 && $to1 < $from2) {
              print STDERR "$flip H $from1:$to1:$from2:$to2\n";
              # they dont overlap in hmmsearch, so we can go ahead an relabel it
              my $pair = $1;
              $pair = 2 if (1 == $flip && 1 == $pair);
              $pair = 1 if (1 == $flip && 2 == $pair);
              my $pe_id = $seq->display_id."/$pair";
              my $inter_distance = $to1 - $from2;
              $self->_store_tagvalue($seq->display_id,$node_id,$flip,$inter_distance);
              # $seq->description($seq->description . " H$from1:$to1:$from2:$to2");
              $seq->display_id($pe_id);
              $pairs->{$pe_id} = $seq;

              next;
            }
          } else {
            my $pe_id = $seq->display_id.".$1";
            $seq->display_id($pe_id);
          }
        }
        my $mapback_stable_id = $seq->display_id; $mapback_stable_id =~ s/\:/\@/g; # mapback
        $seq->display_id($mapback_stable_id);
        $out->write_seq($seq);
        $self->{454} = 1 if ($seq->length > 200);
      }
    }
  }
  # Write the pairs at the end to make sure they are sorted
  foreach my $id (sort keys %$pairs) { my $seq = $pairs->{$id};  $out->write_seq($seq);  }
  $out->close;

  $self->{tree}->store_tag('num_input_reads',$num_input_reads);
  my $runtime = time()*1000-$starttime;
  printf("  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
  $self->{'tree'}->store_tag('fetchreads_runtime_msec', $runtime);

  return $readsfilename;
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

sub parse_alignment {
  my $self = shift;
  print STDERR "parse alignment...\n" if ($self->debug);

  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';
  my $worker_temp_directory = $self->worker_temp_directory;

  my $pachinko_output = $self->{outfile};
  return undef unless (-e $self->{xmlfile} && -e $self->{outfile} &&
                  !(-z $self->{xmlfile}) && !(-z $self->{outfile}));

  # $self->get_alignment_scores; #FIXME

  $self->parse_tree_in_xml($self->{xmlfile});

  $DB::deep = 999988; # This is only useful for debugging purposes, causes no harm
  my $nhx_tree_string = $self->_slurp($self->{nhx_tree_file});
  my $nhx_tree;
  eval {
    $nhx_tree = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($nhx_tree_string);
    my @subnodes = $nhx_tree->get_all_subnodes;
  };

  if($@) {
    throw("Error fetching nhx_tree from output");
  } else {
    $self->{mappedtree} = $nhx_tree;
  }

#   $self->parse_tree_in_nhx_tree($self->{nhx_tree_file});

  my $format = 'fasta';
  my $alignio = Bio::AlignIO->new
    (-file => "$pachinko_output",
     -format => "$format");
  my $aln = $alignio->next_aln();

  $self->{pachinko_aln} = $aln;

  $self->{align_hash} = undef;
  my $count = 0;
  foreach my $seq ($aln->each_seq) {
    my $id = $seq->display_id;
    my $sequence = $seq->seq;
    throw("Error fetching sequence from output alignment") unless(defined($sequence));
    $self->{align_hash}{$id} = $sequence;
    print STDERR "[$count]\n" if ($self->debug && $count++ % 100);
  }

  return 1;
}

sub basic_stats {
  my $self = shift;
  print STDERR "basic stats...\n" if ($self->debug);

  my $aln = $self->{pachinko_aln};

  # Basic stats
  my $pachinko_length       = $aln->length;
  my $pachinko_no_residues  = $aln->no_residues;
  my $pachinko_no_sequences = $aln->no_sequences;
  my $original_length       = $self->{gtsa}->length;
  my $original_no_residues  = $self->{gtsa}->no_residues;
  $self->{original_no_sequences} = $self->{gtsa}->no_sequences;
  $self->{coef_no_residues} = $pachinko_no_residues/$original_no_residues;
  my $coef_no_sequences = $pachinko_no_sequences/$self->{original_no_sequences};
  my $coef_length = $pachinko_length/$original_length;
  $self->{tree}->store_tag('coef_no_residues',$self->{coef_no_residues});
  $self->{tree}->store_tag('coef_no_sequences',$coef_no_sequences);
  $self->{tree}->store_tag('coef_length',$coef_length);

  return 1;
}
sub model_coverage_stats {
  my $self = shift;
  return 1 unless (defined($self->{merge_hash}));

  my $merge = $self->{merge_hash};

  foreach my $prefix (keys %{$merge->{_prefixes}}) {
    # do something
  }
}

sub indel_stats {
  my $self = shift;
  my $tag = $self->{tag} || '\wRR';
  print STDERR "indel stats...\n" if ($self->debug);

  my $node_id = $self->{tree}->node_id;
  my $tag = $self->{tag} || '\wRR';

  return 1 unless (grep(/$tag/, keys %{$self->{align_hash}}));

  # Indel stats
  my $num_aligned_reads = 0;
  foreach my $stable_id (keys %{$self->{align_hash}}) {
    next unless ($stable_id =~ /$tag/);
    $num_aligned_reads++;
    my $string = $self->{align_hash}{$stable_id};
    my $tagged = 0;
    my @g1 = $string =~ /[ATCG]{6}\-{1}[ATCG]{6}/g;
    if (0 < scalar(@g1)) { $self->_store_tagvalue($stable_id, $node_id,'g1',scalar(@g1)); $tagged = 1; }
    my @g2 = $string =~ /[ATCG]{6}\-{2}[ATCG]{6}/g;
    if (0 < scalar(@g2)) { $self->_store_tagvalue($stable_id, $node_id,'g2',scalar(@g2)); $tagged = 1; }
    my @g3 = $string =~ /[ATCG]{6}\-{3}[ATCG]{6}/g;
    if (0 < scalar(@g3)) { $self->_store_tagvalue($stable_id, $node_id,'g3',scalar(@g3)); $tagged = 1; }
    my @g4 = $string =~ /[ATCG]{6}\-{4}[ATCG]{6}/g;
    if (0 < scalar(@g4)) { $self->_store_tagvalue($stable_id, $node_id,'g4',scalar(@g4)); $tagged = 1; }
    my @g5 = $string =~ /[ATCG]{6}\-{5}[ATCG]{6}/g;
    if (0 < scalar(@g5)) { $self->_store_tagvalue($stable_id, $node_id,'g5',scalar(@g5)); $tagged = 1; }
    my @g6 = $string =~ /[ATCG]{6,}\-{6,}[ATCG]{6,}/g;
    if (0 < scalar(@g6)) { $self->_store_tagvalue($stable_id, $node_id,'g6p',scalar(@g6)); $tagged = 1; }

    if ($tagged) {
      my $coef_score = $self->{coef_score}{$stable_id};
      $self->_store_tagvalue($stable_id, $node_id,'coef_score',$coef_score);
    }
  }

  $self->{tree}->store_tag('num_aligned_reads',$num_aligned_reads);

#   # Only for manual inspection
#   if (($self->{coef_no_residues} > 1.5 && $self->{original_no_sequences} > 10 && $self->{original_no_sequences} < 55) || ($self->{coef_no_residues} > 1.0 && $self->{original_no_sequences} > 50 && $self->{original_no_sequences} < 100) || $self->{runtime_msec} > 6000) {
#     my $number = ($node_id-(int($node_id/10)*10));
#     my $subdir = $self->{fastq_dir} . "/" . $number;
#     print STDERR "copying into $subdir...\n" if ($self->debug);
#     my $worker_temp_directory = $self->worker_temp_directory;
#     my $ret1 = `rm -f $subdir/$node_id.pagan.fas`;
#     my $ret2 = `cp $worker_temp_directory/$node_id.pagan.fas $subdir/`;
#     $self->{store_this} = $subdir;
#   }

  return 1;
}

sub _to_cigar_line {
    my $self = shift;
    my $alignment_string = shift;

    $alignment_string =~ s/\-([A-Z])/\- $1/g;
    $alignment_string =~ s/([A-Z])\-/$1 \-/g;
    my @cigar_segments = split " ",$alignment_string;
    my $cigar_line = "";
    foreach my $segment (@cigar_segments) {
      my $seglength = length($segment);
      $seglength = "" if ($seglength == 1);
      if ($segment =~ /^\-+$/) {
        $cigar_line .= $seglength . "D";
      } else {
        $cigar_line .= $seglength . "M";
      }
    }
    return $cigar_line;
}

sub parse_tree_in_xml {
  my $self = shift;
  my $xmlfile = shift;

  $DB::deep = 999988; # This is only useful for debugging purposes, causes no harm
  open XML,"$xmlfile" or die $!;
  while (<XML>) {
    if ($_ =~ /\<newick\>(\S+)\<\/newick\>/ ) {
      $self->{outnewick_string} = $1;
#      $self->{outtree} = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($self->{outnewick_string});
      $self->{xmltree} = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($self->{outnewick_string});
    } elsif ($_ =~ /\<leaf id\=\"(\S+)\" name\=\"(\S+)\"\>/) {
      my $id = $1; my $name = $2;
      my $mappednode = $self->{xmltree}->find_node_by_name($id);
      $mappednode->name($name);
    }
  }

 # We map the TIDs in the nhx tree onto the new xmltree
  foreach my $node ($self->{tree}->get_all_subnodes) {
    my $tid = $node->get_tagvalue('TID');
    next unless ($tid > 0);
  }
  return 1;
}

sub load_species_tree_from_file {
  my ($self) = @_;

  my $taxonDBA  = $self->{'comparaDBA'}->get_NCBITaxonAdaptor();
  my $genomeDBA = $self->{'comparaDBA'}->get_GenomeDBAdaptor();

  if($self->debug()) {
  	print STDERR "load_species_tree_from_file\n";
  	print STDERR "Using GenomeDB based tree\n" if $self->{use_genomedb_id};
  }

  my $newick = $self->_slurp($self->{species_tree_file});
  my $tree = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($newick);

  my %used_ids;

  foreach my $node (@{$tree->all_nodes_in_graph()}) {

    #Split based on - to remove comments & sub * for internal nodes.
    #The ID assigned by NewickParser is not the real ID therefore we need to subsitute this in
    my ($id) = split('-',$node->name);
    $id =~ s/\*//;

    #If it looks like a number then assume we are working with an ID (Taxon or GenomeDB)
    if (looks_like_number($id)) {
      $node->node_id($id);

      if($self->{use_genomedb_id}) {
      	my $gdb = $genomeDBA->fetch_by_dbID($id);
      	throw("Cannot find a GenomeDB for the ID ${id}. Ensure your tree is correct and you are using use_genomedb_id correctly") if !defined $gdb;
      	$node->name($gdb->name());
      	$used_ids{$id} = 1;
      	$node->add_tag('_found_genomedb', 1);
      }
      else {
      	my $ncbi_node = $taxonDBA->fetch_node_by_taxon_id($id);
      	$node->name($ncbi_node->name) if (defined $ncbi_node);
      }
    } 
    else {
      $node->name($id);
    }
    $node->add_tag('taxon_id', $id);
  }
  print STDERR "load_species_tree_from_file loaded.\n" if ($self->debug);

  if($self->{use_genomedb_id}) {
  	print "Searching for overlapping identifiers\n" if $self->debug();
  	my $max_id = max(keys(%used_ids));
  	foreach my $node (@{$tree->all_nodes_in_graph()}) {
  		if($used_ids{$node->node_id()} && ! $node->get_tagvalue('_found_genomedb')) {
  			$max_id++;
  			$node->node_id($max_id);
  		}
  	}
  }

  return $tree;
}

sub _slurp {
  my ($self, $file_name) = @_;
  my $slurped;
  {
    local $/ = undef;
    open(my $fh, '<', $file_name);
    $slurped = <$fh>;
    close($fh);
  }
  return $slurped;
}

sub run_prankancestor_pre {
  my $self = shift;
  my $type = shift;

  my $node_id = $self->{tree}->node_id;
  my $tree = $self->{tree};

  # FETCHING THE REFERENCE CDS ALIGNMENT, it may have been re-run so
  # we need to clean up existing attached nodes
  $self->pre_clean_alignment;

  # FETCHING AA ALN
  print STDERR "Fetching $type simplealign...\n" if ($self->debug);
  my $cdna_type = 0; $cdna_type = 1 if ('cds' eq $type);
  my $aln = $self->{tree}->get_SimpleAlign
    (
     -id_type => 'MEMBER',
     -stop2x => 1,
     -cdna => $cdna_type,
     -append_taxon_id => 1,
    );
  $aln->set_displayname_flat(1);
  $aln = $aln->remove_gaps(undef,1);
  $aln->set_displayname_flat(1);

  my $aln_name = $self->worker_temp_directory . "$node_id.prank_pre.$type";

  my $alignIO = Bio::AlignIO->new
    (-file => ">$aln_name",
     -format => 'fasta'
    );
  $alignIO->write_aln($aln);
  $alignIO->close;
  ########################################

  # REFERENCE TREE
  my $treefile = $self->worker_temp_directory . "$node_id.prank_pre.$type.nh";
  foreach my $subnode ($self->{tree}->get_all_subnodes) {
    $subnode->distance_to_parent(0.0001) if (0 == $subnode->distance_to_parent);
    $subnode->distance_to_parent(2) if (1000 < $subnode->distance_to_parent);
  }
  my $newick_string = $self->{tree}->newick_format('member_id_taxon_id');
#  my $newick_string = $self->{tree}->newick_format;
  $newick_string =~ s/\#\d+\#//g; # get rid of the internal node ids

  open NEWICK, ">$treefile" or die "$!";
  print NEWICK "$newick_string\n";
  close NEWICK;
  ########################################

  my $prankancestor_executable = $self->{prankancestor_exe} || "/software/ensembl/compara/prank/latest/src/prank";
  my $outprefix = $self->worker_temp_directory . "$node_id" . ".prank_pre.$type";

  # /software/ensembl/compara/prank/latest/src/prank -d=<input_alignment> -t=<tree> -shortnames -e -notree -nopost -f=11 -writeanc -o=<output_name>
  my $cmd = "$prankancestor_executable ";
  $cmd .= " -d=$aln_name ";
  $cmd .= " -t=$treefile ";
  $cmd .= " -shortnames -e -notree -nopost -f=11 -writeanc ";
  $cmd .= " -o=$outprefix";
  $cmd .= " 2>&1 > /dev/null" unless($self->debug && defined($self->{verbose}));

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  my $starttime = time();
  print STDERR "prankancestor...\n" if ($self->debug);
  my $worker_temp_directory = $self->worker_temp_directory;

  unless(system("cd $worker_temp_directory; $cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running prankancestor_pre $type $!\n");
  }

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

  my $ancfile = "$worker_temp_directory"."$node_id.prank_pre.$type.0.ancseq";
  throw("Output file not found") unless (-e $ancfile);

  my $runtime = time()*1000-$starttime;
  printf("  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
  my $this_tag = 'prankancestor_' . $type . '_runtime_msec';
  $self->{'tree'}->store_tag($this_tag, $runtime);

  my $result = $self->parse_prankancestor($ancfile,$type);
  $result .= ';';
  my ($comment,$rectree) = split(" ",$result);

  return $rectree;
}

sub run_prankancestor {
  my $self = shift;

  my $node_id = $self->{tree}->node_id;
  my $tag  = $self->{tag} || '\wRR';
  my $tree = $self->{tree};

  print STDERR "Tree pruning...\n" if ($self->debug);

  my $alnfile = $self->{outfile};

#   my $outalnio = Bio::AlignIO->new
#     (-file => $alnfile,
#      -format => 'fasta');
#   my $outaln = $outalnio->next_aln;
#   foreach my $seq ($outaln->each_seq) {
#     if ($seq->display_id =~ /$tag/) {
#       my $all_gaps = '-' x ($seq->length-1); $all_gaps .= 'N';
#       $seq->seq($all_gaps);
#     }
#   }
#   $outaln->set_displayname_flat(1);
#   $outaln = $outaln->remove_gaps(undef,1);
#   $outaln->set_displayname_flat(1);
#   my $collapsedfile = "$alnfile.collapsed";
#   my $collapsedalnio = Bio::AlignIO->new
#     (-file => ">$collapsedfile",
#      -format => 'fasta');
#   $collapsedalnio->write_aln($outaln);
#   $collapsedalnio->close;
#   $outalnio->close;

  my $treefile = $self->worker_temp_directory . "$node_id.nh";
  foreach my $subnode ($self->{mappedtree}->get_all_subnodes) {
    $subnode->distance_to_parent(0.0001) if (0 == $subnode->distance_to_parent);
  }
  my $newick_string = $self->{mappedtree}->newick_format;
  $newick_string =~ s/\#\d+\#//g; # get rid of the internal node ids

  open NEWICK, ">$treefile" or die "$!";
  print NEWICK "$newick_string\n";
  close NEWICK;

  my $prankancestor_executable = $self->{prankancestor_exe} || "/software/ensembl/compara/prank/latest/src/prank";
#    unless (-e $pachinko_executable);

  my $outprefix = $self->worker_temp_directory . "$node_id";

  # /software/ensembl/compara/prank/latest/src/prank -d=<input_alignment> -t=<tree> -shortnames -e -notree -nopost -f=11 -writeanc -o=<output_name>
  my $cmd = "$prankancestor_executable ";
  $cmd .= " -d=$alnfile ";
  $cmd .= " -t=$treefile ";
  $cmd .= " -shortnames -e -notree -nopost -f=11 -writeanc ";
  $cmd .= " -o=$outprefix";
  $cmd .= " 1>/dev/null 2>/dev/null" unless($self->debug && defined($self->{verbose}));

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  my $starttime = time();
  print STDERR "prankancestor...\n" if ($self->debug);
  my $worker_temp_directory = $self->worker_temp_directory;

  unless(system("cd $worker_temp_directory; $cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running prankancestor, $!\n");
  }

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

  $self->{ancfile} = "$worker_temp_directory"."$node_id.0.ancseq";
  throw("Output file not found") unless (-e $self->{ancfile});

  my $runtime = time()*1000-$starttime;
  printf("  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
  $self->{'tree'}->store_tag('prankancestor_runtime_msec', $runtime);

  return 1;
}

sub parse_prankancestor {
  my $self    = shift;
  my $ancfile = shift;
  my $type    = shift;
  print STDERR "parse_prankancestor...\n" if ($self->debug);

  my ($Ptreeref, $PRKnodref);
  my $node_id = $self->{tree}->node_id;
  my $tag  = $self->{tag} || '\wRR';
  my $tree = $self->{tree};

  my @names = ();
  my ($name, $seq) = ('','');
  my $length = -1;
  my $different = 0;

  open(FILE, $ancfile) or die "error opening the prankancestor file: $ancfile\n";
  while (<FILE>) {
    chomp();
    next if /^\s+$/;
    if (/^#\s+\(/) {
      $Ptreeref = $_;
      while (/(#\d+#)/g) { push @{$PRKnodref}, $1; }
      next;
    }
    if (/^\>/) {
      s/\>//;
      s/\s+#.*//;			# Get rid of possible comments (?)
      s/\s+//g;
      push @names, $_;
      unless ($seq =~ /^$/) {
        $seq =~ tr/[a-z]/[A-Z]/;
        if ($length < 0) { $length = length($seq); }
        if ($length != length($seq)) {
          warn "$name is of different length: $length and ".length($seq).".\n";
          $different++;
        }
        $seq =~ s/\-//g;
        $self->{aa_ancestral_seqs}{$name} = $seq    if (   'aa' eq $type);
        $self->{cds_ancestral_seqs}{$name} = $seq   if (  'cds' eq $type);
        $self->{final_ancestral_seqs}{$name} = $seq if ('final' eq $type);
      }
      $name = $_;
      $seq = "";
      next;
    }
    $seq .= $_;
  }

  return $Ptreeref;
}

sub run_this_columbus {
  my $self = shift;
  my $this_tid = shift;
  my $id1 = shift;
  print STDERR "run_this_columbus...\n" if ($self->debug);

  $self->{subset_aa_ancestral_seqs} = $self->produce_list_of_ancestors($id1, $self->{ancestor_depth}, $self->{rectree_aa}, $self->{aa_ancestral_seqs});
  my $samfile = $self->generate_sam_and_ref_file($this_tid,$id1);
  my $velvetfile = $self->columbus_merge_fragments($this_tid,$id1) if (defined($samfile));
  #       $self->parse_contigs_file($this_tid,$id1) if (defined($velvetfile));
  #      $self->run_maq($this_tid,$id1) if (defined($velvetfile) && 1 != $self->{454});
  #      $self->run_bwasw($this_tid,$id1) if (defined($velvetfile) && 1 == $self->{454});
  $self->run_bwasw($this_tid,$id1,undef,$velvetfile) if (defined($velvetfile));

  $self->run_exonerate($this_tid,$id1) if (defined($velvetfile));
}


sub generate_sam_and_ref_file {
  my $self = shift;
  my $this_tid = shift;
  my $id1 = shift;
  my $node_id = $self->{tree}->node_id;
  my $tag  = $self->{tag} || '\wRR';

  print STDERR sprintf("%10d %05d",$this_tid,$id1), " generate_sam_and_ref_file...\n" if ($self->debug);

  my $samfile    = $self->worker_temp_directory . "$this_tid.$id1.$node_id.sam";
  my $reffile    = $self->worker_temp_directory . "$this_tid.$id1.$node_id.ref.fasta";
  my $leavesfile = $self->worker_temp_directory . "$this_tid.$id1.$node_id.leaves.fasta";
  my $tleavesfile = $self->worker_temp_directory . "$this_tid.$id1.$node_id.tleaves.fasta";

  my $outfile = $self->{ancfile};
  throw("no output file") unless (!-z $outfile);

  my $alnio    = Bio::AlignIO->new
    (-file => $outfile,
     -format => 'fasta');
  my $aln = $alnio->next_aln;

  my $refio    = Bio::SeqIO->new
    (-file => ">$reffile",
     -format => 'fasta');

  my $leavesio = Bio::SeqIO->new
    (-file => ">$leavesfile",
     -format => 'fasta');
  my $tleavesio = Bio::SeqIO->new
    (-file => ">$tleavesfile",
     -format => 'fasta');

  foreach my $id (keys %{$self->{subset_aa_ancestral_seqs}}) {
    my $seq = $self->{subset_aa_ancestral_seqs}{$id};
    my $outseq = Bio::Seq->new
      (-display_id => $id,
       -seq => $seq);
    $self->{id_subset_cds_ancestral_seq}{$self->{cds_ancestral_seqs}{$id}}{$id} = 1;
    $tleavesio->write_seq($outseq);
  }
  $tleavesio->close;

  foreach my $seq ($aln->each_seq) {
    # is not a read or fragment
    unless ($seq->display_id =~ /$tag/ || $seq->display_id =~ /^NODE/ || $seq->display_id =~ /^SEQUENCE/) {
      # sequence is a reference
      my $ungapped_seq = $seq->seq; $ungapped_seq =~ s/\-//g;

      #       my $is_subset = 0;
      #       foreach my $id (keys %{$self->{subset_aa_ancestral_seqs}}) {

      #         if (defined ($self->{id_subset_cds_ancestral_seq}{$ungapped_seq}{$id})) {
      #           $is_subset = 1;
      #           # I am thinking that instead of comparing sequences, we
      #           # should really be doing a sorted group_concat of the leaves
      #           # on both the final rectree and the pre rectree, and use the
      #           # cross the ancestral nodes from subset_ to the $aln->each_seq

      #         }
      #       }
      #       next unless (1 == $is_subset);

      my $outseq = Bio::Seq->new
        (-display_id => $seq->display_id,
         -description => $seq->description,
         -seq => $ungapped_seq,
         -alphabet => 'dna');
      my $tseq = $outseq->translate;
      $self->{exonerate_seqs}{$tseq->display_id} = $tseq->seq;
      $refio->write_seq($outseq);
      $leavesio->write_seq($outseq);
      if (defined ($self->{subset_aa_ancestral_seqs}{$seq->display_id})) {

        $self->{ref_ids}{$seq->display_id} = 1;
      }
    }
  }


  $leavesio->close;

  my @sam_entries;
  foreach my $ref_id (keys %{$self->{ref_ids}}) {
    my $this_ref_aln = Bio::SimpleAlign->new();
    my $ref_seq = $aln->get_seq_by_id($ref_id);
    $this_ref_aln->add_seq($ref_seq);
    if (!defined($ref_seq)) {
      throw("refseq not defined");
    }
    foreach my $seq ($aln->each_seq) {
      # This adds a fragment
      if ($seq->display_id =~ /$tag/ || $seq->display_id =~ /^NODE/ || $seq->display_id =~ /^SEQUENCE/) {
        $self->{mapped_reads}{$seq->display_id} = 1;
        $this_ref_aln->add_seq($seq);
      }
    }
    $this_ref_aln = $this_ref_aln->remove_gaps(undef,1);
    my $this_ref_seq = $this_ref_aln->get_seq_by_id($ref_id);
    foreach my $seq ($this_ref_aln->each_seq) {
      # This add a reference
      if ($seq->display_id =~ /$tag/ || $seq->display_id =~ /^NODE/ || $seq->display_id =~ /^SEQUENCE/) {
        my $this_pair = Bio::SimpleAlign->new();
        $this_pair->add_seq($this_ref_seq);
        $this_pair->add_seq($seq);
        $this_pair = $this_pair->remove_gaps(undef,1);
        $self->{this_pair} = $this_pair;
        my ($seq1,$seq2) = $this_pair->each_seq;
        if (!defined $seq1 || !defined $seq2) {
            throw("Undefined sequences in sam entries pair");
        }
        my $cigar_line = $self->to_sam_cigar_line($seq1->seq,$seq2->seq);
        my $sam_entry = $self->sam_entry($seq1,$seq2,$cigar_line);
        push @sam_entries, "$sam_entry\n";
      }
    }
  }

  # Add extra unaligned or discarded reads to sam file
  my $readsfilename = $self->worker_temp_directory . $node_id . ".reads";
  my $format = 'fastq';
  my $readsio = Bio::SeqIO->new(-format => $format,
                                -file => "$readsfilename");
  while (my $seq1 = $readsio->next_seq) {
    next if (defined($self->{mapped_reads}{$seq1->display_id}));
    my $unmapped_sam_entry = $self->sam_entry(undef,$seq1,'*');
    push @sam_entries, "$unmapped_sam_entry\n";
  }
  $readsio->close;

  open SAM,">$samfile" or die $!;
  print SAM join("", sort @sam_entries);
  close SAM;

  $self->{samfile}    = $samfile;
  $self->{reffile}    = $reffile;
  $self->{leavesfile} = $leavesfile;

  if (-z $samfile || -z $reffile) {
    print STDERR "no samfile or reffile produced\n";
    return undef;
  } else {
    return $self->{samfile};
  }

  return;
}

sub produce_list_of_ancestors {
  my $self    = shift;
  my $this_node_id = shift;
  my $depth   = shift;
  my $rectree = shift;
  my $aa_ancestral_seqs = shift;

  my $tree = $self->parse_tree_string($rectree);
  my $this_node = $tree->find_node_by_name($this_node_id);
  my $subset_aa_ancestral_seqs;
  # (1) itself
  $subset_aa_ancestral_seqs->{$this_node->name} = $aa_ancestral_seqs->{$this_node->name};
  # (2) deeper nodes as much as $depth
  while (0 < $depth) {
    $depth--;
    my $try_parent;
    eval { $try_parent = $this_node->parent; };
    if (defined $try_parent) {
      # The following 'next if' is to make sure we dont mix different
      # placement nodes when going back to the parents
      next if (defined ($self->{read_placement_sets}{$this_node->parent->name}));
      next if (1 > length($aa_ancestral_seqs->{$this_node->parent->name}));
      $subset_aa_ancestral_seqs->{$this_node->parent->name} = $aa_ancestral_seqs->{$this_node->parent->name};
    }
    eval { $this_node = $this_node->parent; } ;
  }

  return $subset_aa_ancestral_seqs;
}

sub parse_tree_string {
  my $self = shift;
  my $treestring = shift;
  my $tree;
  eval {
    $tree = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($treestring);
    my @leaves = @{$tree->get_all_leaves};
  };

  throw("Error parsing tree $treestring: $@") if ($@);
  return $tree;
}

sub sam_entry {
  my $self = shift;
  my $ref = shift;
  my $query = shift;
  my $cigar_line = shift;

  my $query_ungapped_seq = $query->seq; $query_ungapped_seq =~ s/\-//g;

  my $qname = $query->display_id;
  my $flag  = 0;   $flag  = 8                if (!defined($ref)); # 0x008 = 8
  my $rname = '*'; $rname = $ref->display_id if (defined($ref));
  my $pos   = 1;
  my $mapq  = 0;
  my $cigar = $cigar_line;
  my $mrnm  = '*';
  my $mpos  = 0;
  my $isize = 0;
  my $seq   = $query_ungapped_seq;
  my $qual  = '';

  my $alignment_section = join ("\t",$qname,$flag ,$rname,$pos  ,$mapq ,$cigar,$mrnm ,$mpos ,$isize,$seq  ,$qual);

  return "$alignment_section";
}

sub to_sam_cigar_line {
  my $self = shift;
  my $ref_seq = shift;
  my $query_seq = shift;
  my $cigar_line = undef;

  $ref_seq =~ s/[^-]/1/g;    $ref_seq =~ s/-/0/g;
  $query_seq =~ s/[^-]/4/g;  $query_seq =~ s/-/2/g;
  my $intersection = $ref_seq | $query_seq;

  while (0 < length($intersection)) {
    if ($intersection =~ /^(4+)/) {
      my $segment = substr($intersection,0,length($1),'');
      $cigar_line .= length($segment) if (1 < length($segment)); $cigar_line .= "I";
    } elsif ($intersection =~ /^(3+)/) {
      my $segment = substr($intersection,0,length($1),'');
      $cigar_line .= length($segment) if (1 < length($segment)); $cigar_line .= "D";
    } elsif ($intersection =~ /^(5+)/) {
      my $segment = substr($intersection,0,length($1),'');
      $cigar_line .= length($segment) if (1 < length($segment)); $cigar_line .= "M";
    } else {
      warn "problem with cigar_line: $intersection\n"; die;
    }
  }

  return $cigar_line;
}

sub columbus_merge_fragments {
  my $self = shift;
  my $this_tid = shift;
  my $id1  = shift;
  my $node_id = $self->{tree}->node_id;
  my $tag  = $self->{tag} || '\wRR';

  print STDERR sprintf("%10d %05d",$this_tid,$id1), " columbus_merge_fragments...\n" if ($self->debug);

  my $samfile = $self->{samfile};
  my $reffile = $self->{reffile};

  my $velvet_dir = $self->worker_temp_directory . "$id1.columbus";
  my $readtype = '-short'; $readtype = '-long' if (1 == $self->{454});
  my $velvet_executable = $self->{velvet_exe} || "/nfs/users/nfs_a/avilella/src/velvet/latest/velvet/velveth";
  # my $cmd = "$velvet_executable $velvet_dir 21 -reference $reffile $readtype -sam $samfile";
  my $cmd = "$velvet_executable $velvet_dir 31 -reference $reffile $readtype -sam $samfile";
  print("$cmd\n");
  unless(system("$cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running velveth, $!\n");
  }

  $velvet_executable =~ s/velveth/velvetg/;
  $cmd = "$velvet_executable $velvet_dir -exp_cov 1";
  print("$cmd\n");
  unless(system("$cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running velvetg, $!\n");
  }

  my $velvetfile = $velvet_dir . "/" . "contigs.fa";
  if (-z $velvetfile) {
    print STDERR "no contigs produced\n";
      return undef;
  } else {
    $self->{velvetfile} = $velvetfile;
    return $self->{velvetfile};
  }
}

# sub parse_contigs_file {
#   my $self = shift;
#   my $this_tid = shift;
#   my $id1  = shift;

#   print STDERR sprintf("%10d %05d",$this_tid,$id1), " parse_contigs_file...\n" if ($self->debug);

#   my $node_id = $self->{tree}->node_id;
#   my $tag  = $self->{tag} || '\wRR';
#   my $velvetfile = $self->{velvetfile};

#   my $velvetio = Bio::SeqIO->new
#     (-format => 'fasta',
#      -file => $velvetfile);

# #   my $table_name = 'pachinko_tree_member';
# #   my $sth = $self->{treeDBA}->prepare
# #     ("INSERT ignore INTO $table_name
# #      (node_id,root_id,member_id,method_link_species_set_id,cigar_line)  VALUES (?,?,?,?,?)");

#   while (my $seq  = $velvetio->next_seq) {
#     my $stable_id_name = $node_id . "_i" . $id1 . "_t" . $this_tid . "_" . $seq->display_id;
#     $stable_id_name =~ s/NODE_//; $stable_id_name =~ s/length_/l_/;
#     $stable_id_name =~ s/\#//g;

#     my $this_qtaxon_id = $self->{qtaxon_id};
#     my $fragment_member = Bio::EnsEMBL::Compara::Member->new
#       (
#        -stable_id    => $stable_id_name,
#        -source_name  => 'EXTERNALCDS',
#        -taxon_id     => $this_qtaxon_id,
#        -genome_db_id => 999999,
#        -description  => $this_tid,
#       );
#     # FIXME -- query taxon_id and genome_db_id
#     $fragment_member->sequence($seq->seq);
#     $self->{memberDBA}->store($fragment_member);

#     # We store it as a child node to the root
#     my $node = new Bio::EnsEMBL::Compara::NestedSet;
#     $self->{tree}->add_child($node);  #FIXME add to the right subnode
#     #leaves are NestedSet objects, bless to make into AlignedMember objects
#     bless $node, "Bio::EnsEMBL::Compara::AlignedMember";
#     $node->member_id($fragment_member->member_id);
#     $node->method_link_species_set_id($self->{mlssid});
#     $self->{treeDBA}->store_node($node);
#     $fragment_member->node_id($node->node_id);
#     print STDERR "New node " . $node->node_id . " for member " . $fragment_member->stable_id . "\n" if ($self->debug);

# #     # We store the cigar_line
# #     my $cigar_line = $self->_to_cigar_line(uc($alignment_string));
# #     my @cigar_match_lengths = map { if ($_ eq '') {$_ = 1} else {$_ = $_;} } map { $_ =~ /^(\d*)/ } ( $cigar_line =~ /(\d*[M])/g );
# #     # Sum up the M lengths
# #     my $seq_cigar_length; map { $seq_cigar_length += $_ } @cigar_match_lengths;
# #     if ($seq_cigar_length != length($fragment_cds)) {
# #       print $fragment_cds."\n".$cigar_line."\n" if ($self->debug);
# #       throw("While storing the cigar line, the returned cigar length did not match the sequence length\n");
# #     }

#     # my $mlss_id = $self->{mlssid};
#     # $sth->execute($node->node_id,$root_id,$fragment_member->member_id,$mlss_id,$cigar_line);

#   }
# #  $sth->finish;
# }

sub store_externalcds {
  my $self = shift;
  my $seq  = shift;
  my $this_tid = shift;

  my $this_qtaxon_id = $self->{qtaxon_id};
  my $fragment_member = Bio::EnsEMBL::Compara::Member->new
    (
     -stable_id    => $seq->display_id,
     -source_name  => 'EXTERNALCDS',
     -taxon_id     => $this_qtaxon_id,
     -genome_db_id => 999999,
     -description  => $this_tid,
    );
  # FIXME -- query taxon_id and genome_db_id
  $fragment_member->sequence($seq->seq);
  $self->{memberDBA}->store($fragment_member);

  # We store it as a child node to the root
  my $node = new Bio::EnsEMBL::Compara::NestedSet;
  $self->{tree}->add_child($node);  #FIXME add to the right subnode
  #leaves are NestedSet objects, bless to make into AlignedMember objects
  bless $node, "Bio::EnsEMBL::Compara::AlignedMember";
  $node->member_id($fragment_member->member_id);
  $node->method_link_species_set_id($self->{mlssid});
  $self->{treeDBA}->store_node($node);
  $fragment_member->node_id($node->node_id);
  print STDERR "New node " . $node->node_id . " for member " . $fragment_member->stable_id . "\n" if ($self->debug);

  return;
}

sub get_alignment_scores {
  #   # Get scores
  #   my $cmd = "/nfs/users/nfs_a/avilella/src/tcoffee/latest/T-COFFEE_distribution_Version_8.27/t_coffee_source/t_coffee ";
  #   $cmd   .= "-cpu=1 ";
  #   $cmd   .= "-infile=$pachinko_output ";
  #   $cmd   .= "-special_mode=evaluate";
  #   $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  #   unless(system("cd $worker_temp_directory; $cmd") == 0) {
  #     print("$cmd\n");
  #     $self->throw("error running t_coffee evaluate, $!\n");
  #   }
  #   $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);
}

sub run_bwasw {
  my $self = shift;
  my $this_tid = shift;
  my $id1  = shift;
  my $contigs_dir = shift;
  my $contigsfile = shift;
  my $node_id = $self->{tree}->node_id;
  my $tag  = $self->{tag} || '\wRR';

  print STDERR sprintf("%10d %05d",$this_tid,$id1), " run_bwasw...\n" if ($self->debug);
  my $readsfilename = $self->worker_temp_directory . $node_id . ".reads";

  $contigs_dir = $self->worker_temp_directory . "$id1.columbus" if (!defined $contigs_dir);
  my $ret = system("rm -f $contigs_dir/cns.filter.snp");

  my $bwa_executable = $self->{bwa_exe} || "/nfs/users/nfs_a/avilella/src/bwa/latest/bwa-0.5.8a/bwa";
  my $samtools_executable = $self->{samtools_exe} || "/nfs/users/nfs_a/avilella/src/samtools/latest/samtools/samtools";

  my $cmd;
  # Index
  $cmd = "$bwa_executable index $contigsfile";
  print STDERR "$cmd && " if ($self->debug);
  unless(system("$cmd") == 0) {    print("$cmd\n");    $self->throw("error running bwa index: $!\n");  }

  # bwasw
  $cmd = "$bwa_executable bwasw -f $contigsfile.bwasw $contigsfile $readsfilename";
  print STDERR "$cmd && " if ($self->debug);
  unless(system("$cmd") == 0) {    print("$cmd\n");    $self->throw("error running bwa bwasw: $!\n");  }

  # view
  $cmd = "$samtools_executable view -bt $contigsfile -o $contigsfile.bam $contigsfile.bwasw";
  print STDERR "$cmd && " if ($self->debug);
  unless(system("$cmd") == 0) {    print("$cmd\n");    $self->throw("error running samtools view: $!\n");  }

  # sort
  $cmd = "$samtools_executable sort $contigsfile.bam $contigsfile.bam.srt";
  print STDERR "$cmd && " if ($self->debug);
  unless(system("$cmd") == 0) {    print("$cmd\n");    $self->throw("error running samtools sort: $!\n");  }

  my $samtools_dir = $samtools_executable; $samtools_dir =~ s/\/samtools$//;
  # pileup & consensus
  # FIXME: -D1 is depth 1, which is not conservative
  my $consensus_fastq_filename = "$contigs_dir/cns.fastq";
  $cmd = "$samtools_executable pileup -cf $contigsfile $contigsfile.bam.srt.bam | $samtools_dir/misc/samtools.pl pileup2fq -D1 > $consensus_fastq_filename";
  print STDERR "$cmd && " if ($self->debug);
  $DB::single=1;1;#??
  unless(system("$cmd") == 0) {    print("$cmd\n");    $self->throw("error running samtools pileup and consensus: $!\n");  }

  my $consio = Bio::SeqIO->new
    (-format => 'fastq',
     -file => $consensus_fastq_filename);

  my $consfile = $self->worker_temp_directory . "$this_tid.$id1.cns.fastq";
  my $onelineconsout = Bio::SeqIO->new
    (-format => 'fastq',
     -file => ">$consfile");

  while (my $seq = $consio->next_seq) {
    next if ($seq->seq =~ /^n+$/);

    $seq->seq(uc($seq->seq));
    $onelineconsout->write_seq($seq);
  }
  $consio->close;
  $onelineconsout->close;

  $self->{consfile} = $consfile;

  return;
}


sub run_maq {
  my $self = shift;
  my $this_tid = shift;
  my $id1  = shift;
  my $node_id = $self->{tree}->node_id;
  my $tag  = $self->{tag} || '\wRR';

  print STDERR sprintf("%10d %05d",$this_tid,$id1), " run_maq...\n" if ($self->debug);

  my $readsfilename = $self->worker_temp_directory . $node_id . ".reads";
  my $velvetfile = $self->{velvetfile};

  my $velvet_dir = $self->worker_temp_directory . "$id1.columbus";
  my $ret = system("rm -f $velvet_dir/cns.filter.snp");

  my $cmd = "perl /nfs/users/nfs_a/avilella/src/maq/maq-0.7.1/scripts/maq.pl easyrun -d $velvet_dir $velvetfile $readsfilename";
  print("$cmd\n");
  unless(system("$cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running maq, $!\n");
  }

  my $consensus_fastq_filename = $velvet_dir . "/" . "cns.fq";
  my $consio = Bio::SeqIO->new
    (-format => 'fastq',
     -file => $consensus_fastq_filename);

  my $consfile = $self->worker_temp_directory . "$this_tid.$id1.cns.fasta";
  my $consout = Bio::SeqIO->new
    (-format => 'fasta',
     -file => ">$consfile");

  while (my $seq = $consio->next_seq) {
    next if ($seq->seq =~ /^n+$/);
    my $stable_id_name = $node_id . "_i" . $id1 . "_t" . $this_tid . "_" . $seq->display_id;
    $stable_id_name =~ s/NODE_//; $stable_id_name =~ s/length_/l_/;
    $stable_id_name =~ s/\#//g;
    my $filtseq = lc($seq->seq); $filtseq =~ s/n//g; my $filtlen = length($filtseq);
    $self->_store_tagvalue($stable_id_name,$node_id,'filtlen',$filtlen);

    my $this_qtaxon_id = $self->{qtaxon_id};
    my $fragment_member = Bio::EnsEMBL::Compara::Member->new
      (
       -stable_id    => $stable_id_name,
       -source_name  => 'EXTERNALCDS',
       -taxon_id     => $this_qtaxon_id,
       -genome_db_id => 999999,
       -description  => $this_tid,
      );
    # FIXME -- query taxon_id and genome_db_id
    $fragment_member->sequence($seq->seq);
    $self->{memberDBA}->store($fragment_member);

    # We store it for transcript translation
    my $consseq = Bio::Seq->new
      (-display_id => $stable_id_name,
       -seq => $seq->seq,
       -alphabet => 'dna');
    $consout->write_seq($consseq);

    # We store it as a child node to the root
    my $node = new Bio::EnsEMBL::Compara::NestedSet;
    $self->{tree}->add_child($node);  #FIXME add to the right subnode
    #leaves are NestedSet objects, bless to make into AlignedMember objects
    bless $node, "Bio::EnsEMBL::Compara::AlignedMember";
    $node->member_id($fragment_member->member_id);
    $node->method_link_species_set_id($self->{mlssid});
    $self->{treeDBA}->store_node($node);
    $fragment_member->node_id($node->node_id);
    print STDERR "New node " . $node->node_id . " for member " . $fragment_member->stable_id . "\n" if ($self->debug);
  }
  $consio->close;
  $consout->close;

  $self->{consfile} = $consfile;
  $self->tag_amount_sequence($consfile,'fasta','maq');

  return;
}

sub tag_amount_sequence {
  my $self   = shift;
  my $file   = shift;
  my $format = shift;
  my $tag    = shift;

  my $seqio = Bio::SeqIO->new
    (-file => $file,
     -format => $format);
  my $seqlen = 0; my $numseq = 0;
  while (my $seq = $seqio->next_seq) {
    $seqlen += $seq->length;
    $numseq++;
  }
  $seqio->close;

  $self->{tree}->store_tag('stats_seqlen_'.$tag,$seqlen);
  $self->{tree}->store_tag('stats_numseq_'.$tag,$numseq);
  return 0;
}

sub run_exonerate {
  my $self = shift;
  my $this_tid = shift;
  my $id1  = shift;


  print STDERR sprintf("%10d %05d",$this_tid,$id1), " run_exonerate...\n" if ($self->debug);

  my $node_id = $self->{tree}->node_id;
  my $tag  = $self->{tag} || '\wRR';
  my $velvetfile = $self->{velvetfile};
  my $consfile = $self->{consfile};
  if (-z $consfile) {
    throw ("Consensus file is empty: $consfile\n");
  }
  my $leavesfile = $self->worker_temp_directory . "$this_tid.$id1.$node_id.leaves.fasta";
  my $tleavesfile = $self->worker_temp_directory . "$this_tid.$id1.$node_id.tleaves.fasta";

  my $inputio = Bio::SeqIO->new
    (-format => 'fastq',
     -file => $consfile);
  my $outio = Bio::SeqIO->new
    (-format => 'fasta',
     -file => ">$consfile.fasta");

  while (my $seq  = $inputio->next_seq) {
    $self->{exonerate_seqs}{$seq->display_id} = $seq->seq;
    $self->store_externalcds($seq,$this_tid);
    $outio->write_seq($seq);
  }
  $inputio->close;
  $outio->close;

  my $exonerate_executable = $self->{exonerate_exe} || "/nfs/users/nfs_a/avilella/src/exonerate/latest/exonerate-2.2.0-x86_64/bin/exonerate";
#    unless (-e $exonerate_executable);

  # my $model = 'cdna2genome';
  # $self->{query_type} = 'dna';
  $self->{model} = 'protein2genome'; my $model = $self->{model};
  $self->{query_type} = 'protein';
  my $worker_temp_directory = $self->worker_temp_directory;
  my $resultsfile = $worker_temp_directory . $node_id . ".$model";

  my $cmd = "$exonerate_executable ";
  $cmd .= " --model $model";
  #  $cmd .= " --useaatla FALSE --showalignment" if ($self->debug);
  $cmd .= " --frameshift -5" if (1 == $self->{454});
#   $cmd .= " --forwardcoordinates FALSE";
   $cmd .= " --refine FULL";
#   $cmd .= " --exhaustive FALSE";
#   $cmd .= " --score 50 --dnahspthreshold 60 --dnawordlen 12";
  $cmd .= " --alignmentwidth 200 --showalignment" if ($self->debug);
  $cmd .= " --query $leavesfile" if ($model eq 'cdna2genome');
  $cmd .= " --query $tleavesfile" if ($model eq 'protein2genome');
  $cmd .= " --target $consfile.fasta ";
  # $cmd .= " > $resultsfile ";
  print "# $cmd \n" if ($self->debug);

  my $fh;
  eval { open($fh, "$cmd |") || die $!; };
  if ($@) {
    warn("problem with exonerate $@ $!");
    return;
  }

  my $starttime = time();
  print STDERR "exonerate...\n" if ($self->debug);


  my $transcripts = $self->parse_exonerate_results($fh);
  $self->{tree}->store_tag('num_exonerate_tr',scalar @$transcripts);
  my $runtime = time()*1000-$starttime;
  printf("  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
  $self->{'tree'}->store_tag('exonerate_runtime_msec', $runtime);

  my $this_id;
  my $worker_temp_directory = $self->worker_temp_directory;
  my $transcriptsfile = $worker_temp_directory . $node_id . ".transcripts";
  open TRANSCRIPTSFILE,">$transcriptsfile" or die $!;
  foreach my $transcript (@$transcripts) {
    # Only the best scoring transcript for each input seq is stored
    next if (defined($this_id->{$transcript->seqname}));
    my $trseq = $transcript->seq;
    my $this_seq;
    $this_seq = $trseq     unless ($transcript->strand eq -1);
    $this_seq = $trseq->revcom if ($transcript->strand eq -1);

    # Store transcript
    print STDERR "Storing transcript\n" if ($self->debug);
    my $stable_id_name = $node_id . "_i" . $id1 . "_t" . $this_tid . "_" . $transcript->seqname;
    $stable_id_name =~ s/\#//g;

    my $this_qtaxon_id = $self->{qtaxon_id};
    my $fragment_member = Bio::EnsEMBL::Compara::Member->new
      (
       -stable_id    => $stable_id_name,
       -source_name  => 'ENSEMBLTRANS',
       -taxon_id     => $this_qtaxon_id,
       -genome_db_id => 999999,
       -description  => $this_tid,
      );
    # FIXME -- query taxon_id and genome_db_id
    $fragment_member->sequence($this_seq->seq);
    $self->{memberDBA}->store($fragment_member);

    # We store it as a child node to the root
    my $node = new Bio::EnsEMBL::Compara::NestedSet;
    $self->{tree}->add_child($node);  #FIXME add to the right subnode
    #leaves are NestedSet objects, bless to make into AlignedMember objects
    bless $node, "Bio::EnsEMBL::Compara::AlignedMember";
    $node->member_id($fragment_member->member_id);
    $node->method_link_species_set_id($self->{mlssid});
    $self->{treeDBA}->store_node($node);
    $fragment_member->node_id($node->node_id);
    print STDERR "New node " . $node->node_id . " for member " . $fragment_member->stable_id . "\n" if ($self->debug);

    my $sql = "UPDATE protein_tree_member SET";
    $sql .= " root_id=" . $node_id;
    $sql .= " WHERE member_id=". $fragment_member->member_id;
    $self->dbc->do($sql);

    print STDERR "Storing peptide\n" if ($self->debug);
    my $translated_this_seq = $this_seq->translate;
    my $translated_fragment_member = Bio::EnsEMBL::Compara::Member->new
      (
       -stable_id    => $stable_id_name,
       -source_name  => 'ENSEMBLPEP',
       -taxon_id     => $self->{qtaxon_id},
       -description  => $self->{qtaxon_id},
       -genome_db_id => 999999,
       -description  => $this_tid,
      );
    # FIXME -- query taxon_id and genome_db_id
    $translated_fragment_member->sequence($translated_this_seq->seq);
    $self->{memberDBA}->store($translated_fragment_member);

    # We store it as a child node to the root
    my $node = new Bio::EnsEMBL::Compara::NestedSet;
    $self->{tree}->add_child($node);  #FIXME add to the right subnode
    #leaves are NestedSet objects, bless to make into AlignedMember objects
    bless $node, "Bio::EnsEMBL::Compara::AlignedMember";
    $node->member_id($translated_fragment_member->member_id);
    $node->method_link_species_set_id($self->{mlssid});
    $self->{treeDBA}->store_node($node);
    $translated_fragment_member->node_id($node->node_id);
    print STDERR "New node " . $node->node_id . " for member " . $translated_fragment_member->stable_id . "\n" if ($self->debug);

    my $sql = "UPDATE protein_tree_member SET";
    $sql .= " root_id=" . $node_id;
    $sql .= " WHERE member_id=". $translated_fragment_member->member_id;
    $self->dbc->do($sql);

    print TRANSCRIPTSFILE ">" . $transcript->seqname , "\n", $translated_this_seq->seq, "\n";
    $this_id->{$transcript->seqname} = 1;
    # my $trseq = $transcript->translation->end_Exon->seq->translate->seq;
  }

  close TRANSCRIPTSFILE;
  return 1;
}

1;


sub parse_exonerate_results {
	my ($self, $fh) = @_;
	my %strand_lookup = ( '+' => 1, '-' => -1, '.' => 1 );

        my @transcripts;
        
        # Parse output looking for lines beginning with 'RESULT:'.
        # Each line represents a distinct match to one sequence
        # containing multiple 'exons'.

	while (<$fh>){
          # next unless /^RESULT:/;
		next unless /^vulgar\:/;
                # print STDERR "$_\n" if ($self->debug);
                chomp;

                my ($perc_id, $q_length, $t_length, $gene_orientation);
                # my ($tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score, $perc_id, $q_length, $t_length, $gene_orientation,
                my ($tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score,
                    @align_components) = split;
                $q_length = length($self->{exonerate_seqs}{$q_id});
                $t_length = length($self->{exonerate_seqs}{$t_id});

                $t_strand = $strand_lookup{$t_strand};
                $q_strand = $strand_lookup{$q_strand};
                $gene_orientation = $strand_lookup{$gene_orientation};

                # Read vulgar information and extract exon regions.
                my $exons = $self->_parse_vulgar_block($t_start,
                                                       $t_end,
                                                       $t_strand,
                                                       $t_length,
                                                       $q_start, 
                                                       $q_end,
                                                       $q_strand,
                                                       $q_length,
                                                       \@align_components);

                # now we have extracted the exons and the coordinates are with 
                # reference to the forward strand of the query and target, we can 
                # use the gene_orienation to flip the strands if necessary
                if ($gene_orientation == -1 and $t_strand == 1) {
                  $t_strand *= -1;
                  $q_strand *= -1;
                }

                my $covered_count = 0;
                if ($self->{coverage_as_proportion_of_aligned_residues}) {
                  foreach my $exon (@$exons) {
                    foreach my $sf (@{$exon->{sf}}) {
                      $covered_count += $sf->{query_end} - $sf->{query_start} + 1;
                    }
                  }
                } else {
                  $covered_count = abs($q_end - $q_start);
                }

                my $coverage = sprintf("%.2f", 100 * $covered_count / $q_length);

                # Build FeaturePairs for each region of query aligned to a single
                # Exon.  Create a DnaDnaAlignFeature from these FeaturePairs and then
                # attach this to our Exon.
                my $transcript = Bio::EnsEMBL::Transcript->new();

                my $slice_sequence;
                my $cs; my $slice;
                if ($self->{model} eq 'protein2genome') {
                  $slice_sequence = $self->{exonerate_seqs}{$t_id};
                  $cs = Bio::EnsEMBL::CoordSystem->new(-NAME    => $t_id, -RANK    => 1);
                  $slice = Bio::EnsEMBL::Slice->new
                    (-seq_region_name  => 'misc',
                     -name             => $t_id,
                     -seq_region_length => $t_length,
                     -start            => 1,
                     -end              => $t_length,
                     -strand           => $t_strand,
                     -coord_system     => $cs,
                     -seq              => $slice_sequence,
                    );
                } else {
                  throw("only protein2genome works right now");
#                   $slice_sequence = $self->{exonerate_seqs}{$q_id};
#                   $cs = Bio::EnsEMBL::CoordSystem->new(-NAME    => $q_id, -RANK    => 1);
#                   $slice = Bio::EnsEMBL::Slice->new
#                     (-seq_region_name  => 'misc',
#                      -name             => $q_id,
#                      -seq_region_length => $q_length,
#                      -start            => 1,
#                      -end              => $q_length,
#                      -strand           => $q_strand,
#                      -coord_system     => $cs,
#                      -seq              => $slice_sequence,
#                     );
                }
                my (@tran_feature_pairs, 
                    $cds_start_exon, $cds_start,
                    $cds_end_exon, $cds_end);

                foreach my $proto_exon (@$exons){

                  # Build our exon and set its key values.
                  my $exon = Bio::EnsEMBL::Exon->new();

                  $exon->seqname($t_id);
                  $exon->start($proto_exon->{exon_start});
                  $exon->end($proto_exon->{exon_end});
                  $exon->phase($proto_exon->{phase});
                  $exon->end_phase($proto_exon->{end_phase});
                  $exon->strand($t_strand);

                  my @exon_feature_pairs;
                  foreach my $sf (@{$proto_exon->{sf}}){
                    my $feature_pair = Bio::EnsEMBL::FeaturePair->new
                      (-seqname    => $t_id,
                       -start      => $sf->{target_start},
                       -end        => $sf->{target_end},
                       -strand     => $t_strand,
                       -hseqname   => $q_id,
                       -hstart     => $sf->{query_start},
                       -hend       => $sf->{query_end},
                       -hstrand    => $q_strand,
                       -hcoverage  => $coverage,
                       -score      => $coverage,
                       -percent_id => $perc_id);


                    push @exon_feature_pairs, $feature_pair;
                    push @tran_feature_pairs, $feature_pair;
                  }

                  if (@exon_feature_pairs) {
                    # Use our feature pairs for this exon to create a single 
                    # supporting feature (with cigar line).
                    my $supp_feature;

                    eval{
                      if ($self->{query_type} eq 'protein') {
                        $supp_feature =
                          Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@exon_feature_pairs);
                      } else {
                        $supp_feature = 
                          Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@exon_feature_pairs);
                      }
                    };
                    if ($@) {
                      warning($@);
                      next TRANSCRIPT;
                    }      
                    $exon->add_supporting_features($supp_feature);
                  }      

                  if (exists $proto_exon->{coding_start}) {
                    if (not defined $cds_start_exon) {
                      $cds_start_exon = $exon;
                      $cds_start = $proto_exon->{coding_start};
                    }
                    $cds_end_exon = $exon;
                    $cds_end = $proto_exon->{coding_end};
                  }


                  $transcript->add_Exon($exon);
                }

                # Create a single supporting feature for the whole transcript
                my $t_supp_feat;
                eval{
                  if ($self->{query_type} eq 'protein') {
                    $t_supp_feat =
                      Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@tran_feature_pairs);
                  } else {
                    $t_supp_feat = 
                      Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@tran_feature_pairs);
                  }
                };
                if ($@) {
                  warning("Could not create Transcript supporting feature");
                } else {
                  $transcript->add_supporting_features($t_supp_feat);
                }

                my @exons = @{$transcript->get_all_Exons};

                if (scalar(@exons)) {
                  if (defined $cds_start_exon) {
                    my $translation = Bio::EnsEMBL::Translation->new();
                    $translation->start_Exon($cds_start_exon);
                    $translation->start($cds_start);      
                    $translation->end_Exon($cds_end_exon);
                    $translation->end($cds_end);

                    $transcript->translation($translation);
                  }

                  calculate_exon_phases($transcript, 0);

                  $transcript->slice($slice);
                  foreach my $sf(@{$transcript->get_all_supporting_features}){
                    $sf->slice($slice);
                  }
                  foreach my $exon(@{$transcript->get_all_Exons}){
                    $exon->slice($slice);
                    foreach my $sf(@{$exon->get_all_supporting_features}){
                      $sf->slice($slice);
                    }
                  }
#                   $transcript = $transcript->transfer($t_slice);

                  push @transcripts, [$score, $transcript];
                }
              }

        @transcripts = sort { $b->[0] <=> $a->[0] } @transcripts;
        @transcripts = map { $_->[1] } @transcripts;

        return \@transcripts;
}


sub _parse_vulgar_block {
  my ($self, 
      $target_start, $target_end, $target_strand, $target_length,
      $query_start, $query_end,  $query_strand, $query_length,
      $vulgar_components) = @_;

  # This method works along the length of a vulgar line 
  # exon-by-exon.  Matches that comprise an exon are 
  # grouped and an array of 'proto-exons' is returned.
  # Coordinates from the vulgar line are extrapolated 
  # to actual genomic/query coordinates.

  my @exons;
  my $exon_number = 0;

  # We sometimes need to increment all our start coordinates. Exonerate 
  # has a coordinate scheme that counts _between_ nucleotides at the start.
  # However, for reverse strand matches 
  
  my ($query_in_forward_coords, $target_in_forward_coords);
  my ($cumulative_query_coord, $cumulative_target_coord);

  if ($target_start > $target_end) {
    warn("For target, start and end are in thew wrong order for a reverse strand match")
        if $target_strand != -1;
    $cumulative_target_coord = $target_start;
    $target_in_forward_coords = 1;
  } else {
    $cumulative_target_coord = $target_start + 1;
    $target_in_forward_coords = 0;
  }
  if ($query_start > $query_end) {
    warn("For query, start and end are in thew wrong order for a reverse strand match")
        if $query_strand != -1;
    $cumulative_query_coord = $query_start;
    $query_in_forward_coords = 1;
  } else {
    $cumulative_query_coord = $query_start + 1;
    $query_in_forward_coords = 0;
  }


  while (@$vulgar_components){
    throw("Something funny has happened to the input vulgar string." .
		 "  Expecting components in multiples of three, but only have [" .
		 scalar @$vulgar_components . "] items left to process.")
      unless scalar @$vulgar_components >= 3;

    my $type                = shift @$vulgar_components;
    my $query_match_length  = shift @$vulgar_components;
    my $target_match_length = shift @$vulgar_components;

    throw("Vulgar string does not start with a match.  Was not " . 
		 "expecting this. (Have type $type)")
      if (scalar @exons == 0) && $type ne 'M' && $type ne 'S' && $type ne 'C';

    if ($type eq 'M' or $type eq 'S' or $type eq 'C'){
      my %hash = (type => $type);

      if ($target_strand == -1) {
        if ($target_in_forward_coords) {
          $hash{target_start} = $cumulative_target_coord - ($target_match_length - 1);
          $hash{target_end}   = $cumulative_target_coord;
        } else {
          $hash{target_end}   = $target_length - ($cumulative_target_coord - 1);
          $hash{target_start} = $hash{target_end} - ($target_match_length - 1);
        }
      } else {
        $hash{target_start} = $cumulative_target_coord;
        $hash{target_end}   = $cumulative_target_coord + ($target_match_length - 1);
      }

      if ($query_strand == -1) {
        if ($query_in_forward_coords) {
          $hash{query_start} = $cumulative_query_coord - ($query_match_length - 1);
          $hash{query_end}   = $cumulative_query_coord;
        } else {
          $hash{query_end}   = $query_length - ($cumulative_query_coord - 1);
          $hash{query_start} = $hash{query_end} - ($query_match_length - 1);
        }
      } else {
        $hash{query_start} = $cumulative_query_coord;
        $hash{query_end}   = $cumulative_query_coord + ($query_match_length - 1);
      }

      # there is nothing to add if this is the last state of the exon
      $exons[$exon_number]->{gap_end}   = 0;
      push @{$exons[$exon_number]->{sf}}, \%hash;
    }
    elsif ($type eq "G") {
      if (exists($exons[$exon_number]->{sf})) {
        # this is the gap in the middle of an exon, or at the end. Assume it is 
        # at the end, and then reset if we see another match state in this exon
        $exons[$exon_number]->{gap_end}   = $target_match_length;
      } else {
        # this is a gap at the start of an exon; 
        $exons[$exon_number]->{gap_start} = $target_match_length;
      }
    }
    elsif ($type eq "I" or
           $type eq "F") {

      # in protein mode, any insertion on the genomic side should be treated as 
      # an intron to ensure that the result translates. However, we allow for
      # codon insertions in the genomic sequence with respect to the protein. 
      # This introduces the possibility of in-frame stops, but I don't
      # think "introning over" these insertions is appropriate here. 

      # if we see a gap/intron immediately after an intron, the current exon is "empty"
      if ($exons[$exon_number]) {
        $exon_number++;
      }
    }

    if ($target_in_forward_coords and $target_strand == -1) {
      $cumulative_target_coord -= $target_match_length;
    } else {
      $cumulative_target_coord += $target_match_length;
    }
    if ($query_in_forward_coords and $query_strand == -1) {
      $cumulative_query_coord  -= $query_match_length;
    }
    else {
      $cumulative_query_coord  += $query_match_length;
    }

  }

  for(my $i = 0; $i < @exons; $i++) {
    my $ex = $exons[$i];
    my @ex_sf = @{$ex->{sf}};

    if ($target_strand == -1) {
      $ex->{exon_start} = $ex_sf[-1]->{target_start};
      $ex->{exon_end}   = $ex_sf[0]->{target_end};

      if (exists $ex->{gap_start}) {
        $ex->{exon_end} += $ex->{gap_start};
      }
      if (exists $ex->{gap_end}) {
        $ex->{exon_start} -= $ex->{gap_end};
      }

    } else {
      $ex->{exon_start} = $ex_sf[0]->{target_start};
      $ex->{exon_end}   = $ex_sf[-1]->{target_end};

      if (exists $ex->{gap_start}) {
        $ex->{exon_start} -= $ex->{gap_start};
      }
      if (exists $ex->{gap_end}) {
        $ex->{exon_end} += $ex->{gap_end};
      }
    }

    # split codons are a pain. If the query is dna, we must be in the 
    # cdna2genome model so they must be part of the supporting feature. 
    # If query is protein, they need to be removed from the supporting feature

    if ($self->{query_type} eq 'dna') {
      map { $_->{type} = 'C' if $_->{type} eq 'S' } @ex_sf;
      if (my @cod = grep { $_->{type} eq 'C' } @ex_sf) {
        # at least part of this exon is coding

        @cod = sort { $a->{target_start} <=> $b->{target_start} } @cod;
        my $cod_start = $cod[0]->{target_start};
        my $cod_end   = $cod[-1]->{target_end};
        my $cod_len   = $cod_end - $cod_start + 1;

        if ($target_strand == -1) {
          $ex->{coding_start} = $ex->{exon_end} - $cod_end + 1;
        } else {
          $ex->{coding_start} = $cod_start - $ex->{exon_start} + 1;
        }

        $ex->{coding_end}   = $ex->{coding_start} + $cod_len - 1;
      }
      # merge together abutting ungapped features
      my @nr_sf;
      foreach my $sf (@ex_sf) {
        my $merged = 0;
        if (@nr_sf and
            $sf->{type} eq 'C' and
            $nr_sf[-1]->{type} eq 'C' and
            $sf->{query_start} == $nr_sf[-1]->{query_end} + 1) {
          if ($target_strand == -1) {
            if ($sf->{target_end} == $nr_sf[-1]->{target_start} - 1) {
              $nr_sf[-1]->{target_start} = $sf->{target_start};
              $nr_sf[-1]->{query_end}  = $sf->{query_end};
              $merged = 1;
            }
          } else {
            if ($sf->{target_start} == $nr_sf[-1]->{target_end} + 1) {
              $nr_sf[-1]->{target_end} = $sf->{target_end};
              $nr_sf[-1]->{query_end}  = $sf->{query_end};
              $merged = 1;
            }
          }         
        }
        if (not $merged) {
          push @nr_sf, $sf;
        }
      }
      @ex_sf = @nr_sf;
    } else {
      # query type is protein
      @ex_sf = grep { $_->{type} ne 'S' } @ex_sf;
      $ex->{coding_start} = 1;
      $ex->{coding_end} = $ex->{exon_end} - $ex->{exon_start} + 1;
    }

    $ex->{sf} = \@ex_sf;
  }

  return \@exons;
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

sub illumina_same_set {
  my $self = shift;
  my $assembly = shift;
  my $stable_id = shift;

  my $discard = 0;
  my $read_lane_id;
  my $assembly_lane_id;
  if ($assembly =~ /\S+\_([12])\_([12])\_pf/) {
    $assembly_lane_id = $1 . "_" . $2;
  }
  if ($stable_id =~ /HWUSI-EAS\S+:(\d+):\d+:\d+:\d+#\d+\/([12])/) {
    $read_lane_id = int($1) . "_" . $2;
  }
  if ($read_lane_id ne $assembly_lane_id) {
    $discard = 1;
  }

  return $discard;
}

sub pre_clean_alignment {
  my $self = shift;
  my $tag  = $self->{tag} || '\wRR';
  my @preclean;
  foreach my $leaf (@{$self->{tree}->get_all_leaves}) {
#     if ($leaf->stable_id =~ /$tag/ || $leaf->stable_id eq '' || $leaf->source_name eq 'EXTERNALCDS' || !defined($leaf->genome_db)) {
    my $this_taxon_id;
    eval { $this_taxon_id = $leaf->taxon_id; };
    if ($@) {
      throw("couldnt find taxon for member\n");
    }
    if (!defined $self->{species_set_hash}{$leaf->taxon_id}) {
      push @preclean, $leaf;
      # self_collapse option is only for testing purposes
    } elsif (1 == $self->{self_collapse} && defined($self->{qtaxon_list}{$leaf->taxon_id})) {
      push @preclean, $leaf;
    }
  }

  $self->{tree} = $self->{tree}->remove_nodes(\@preclean);
  $self->{tree} = $self->{tree}->minimize_tree;
  foreach my $node (@preclean) {
    if ($node->taxon_id == $self->{qtaxon_id} || 0 == $node->taxon_id || $node->description == $self->{qtaxon_id}) {

      # Only delete the nodes for the same qtaxon_id as in this job
      $self->delete_old_pachinko_member($node);
    }
  }

  return 0;
}

sub delete_old_pachinko_member {
  my $self = shift;
  my $node = shift;

  my $this_old_member_id = $node->member_id;
  $self->{treeDBA}->delete_flattened_leaf($node);
  $self->{comparaDBA}->dbc->do("DELETE from protein_tree_member WHERE member_id = $this_old_member_id");
  $self->{comparaDBA}->dbc->do("DELETE from member WHERE member_id = $this_old_member_id");

  return;
}

sub dataflow_pachinko_tree {
  my $self = shift;
  my $node_id = $self->{tree}->node_id;
  my $clusterset_id = $self->{tree}->parent->node_id;

  # Dataflow clusters
  # This will create a new jobs for each of the newly extended gene family
  my $output_id = sprintf("{'protein_tree_id'=>%d, 'clusterset_id'=>%d}", 
                          $node_id, $clusterset_id);
  print STDERR "Dataflow # $output_id\n" if ($self->debug);
  $self->dataflow_output_id($output_id, 0);
}


# sub merge_fragments {
#   my $self = shift;

#   my $alnio = $self->{pachinko_aln};
#   my $tag = shift || '\wRR';

#   my $merge;
#   my $reference;

#   print STDERR "Sorting references and non-references\n" if ($self->debug);
#   my ($this_prefix,$dummy) = ('reference',undef);

#   # Label the reads from the references in the aln -- we only want to
#   # merge reads

#   foreach my $this_seq ($alnio->each_seq) {
#     if ($this_seq->display_id =~ /$tag/) {
#       ($this_prefix,$dummy) = split('\.',$this_seq->display_id) if ($tag =~ /RR/);
#       if ($tag !~ /RR/) {
#         my $display_id = $this_seq->display_id; $display_id =~ /($tag.+)\|/;
#         $this_prefix = $1;
#       }
#       $merge->{_prefixes}{$this_prefix}{$this_seq->display_id} = 1;
#       $merge->{_seqs}{$this_seq->display_id}{$this_prefix} = $this_seq->seq;
#     } else {
#       $merge->{_refseqs}{$this_seq->display_id} = $this_seq->seq;
#     }
#   }
#   return 1 unless (0 < scalar keys %{$merge->{_seqs}});

#   my $merge_stats;
#   print STDERR "Merging...\n" if ($self->debug);
#   foreach my $this_seq ( $alnio->each_seq ) {
#     next unless (defined($merge->{_seqs}{$this_seq->display_id}));
#     my @temp = keys %{$merge->{_seqs}{$this_seq->display_id}}; my $read_id = $temp[0];
#     next if (2 > scalar keys %{$merge->{_prefixes}{$read_id}});
#     # print STDERR "." if ($self->debug);
#     my ($this_prefix,$dummy) = split('\.',$read_id) if ($tag =~ /RR/);
#     if ($tag !~ /RR/) {
#       my $display_id = $this_seq->display_id; $display_id =~ /($tag.+)\|/;
#       $this_prefix = $1;
#     }

#     my @read_ids = keys %{$merge->{_prefixes}{$read_id}};
#     while (my $read_id1 = shift (@read_ids)) {
#       foreach my $read_id2 (@read_ids) {
#         my $read1 = $merge->{_seqs}{$read_id1}{$this_prefix};
#         my $read2 = $merge->{_seqs}{$read_id2}{$this_prefix};
#         my $pre1 = $read1; $pre1 =~ s/\-//g; my $pre2 = $read2; $pre2 =~ s/\-//g;

#         my $lcread1 = lc($read1); my $lcread2 = lc($read2);
#         $read1 =~ s/-/substr($lcread2, pos($read1), 1)/eg;
#         $read2 =~ s/-/substr($lcread1, pos($read2), 1)/eg;

#         my $post1 = $read1; $post1 =~ s/\-//g; my $post2 = $read2; $post2 =~ s/\-//g;
#         if (
#             length($post1)>$merge->{_seqs}{$read_id1}{$this_prefix}
#             ||
#             length($post2)>$merge->{_seqs}{$read_id2}{$this_prefix}
#            ) {
#           # $merge_stats->{extension_events}{$this_prefix}++;
#           my $pre_blocks  = scalar grep {length $_ > 2} $pre1 =~ /([ATCGN]+|[atcgn]+)/g;
#           my $post_blocks = scalar grep {length $_ > 2} $post1 =~ /([ATCGN]+|[atcgn]+)/g;
#           my $block_increase = $post_blocks - $pre_blocks;
#           if (1 == $block_increase) {
#             $merge_stats->{bridging_events}{$this_prefix}{1}++;
#           } elsif (1 < $block_increase) {
#             $merge_stats->{bridging_events}{$this_prefix}{$block_increase}++;
#           }
#         }
#         $merge->{_seqs}{$read_id1}{$this_prefix} = $read1;
#         $merge->{_seqs}{$read_id2}{$this_prefix} = $read2;
#       }
#     }
#   }

#   my $node_id = $self->{tree}->node_id;
#   my $tag  = $self->{tag} || '\wRR';
#   foreach my $prefix (keys %{$merge_stats->{bridging_events}}) {
#     next unless (defined($merge_stats->{bridging_events}{$prefix}));
#     foreach my $event_increase (keys %{$merge_stats->{bridging_events}{$prefix}}) {
#       my $tag = 'bridging_' . $event_increase;
#       my $value = $merge_stats->{bridging_events}{$prefix}{$event_increase};
#       $self->_store_tagvalue($prefix, $node_id,$tag,$value);
#     }
#   }

#   my $nodups_merge_num = 0;
#   unless (defined($self->{allow_dups})) {
#     my @ids = keys %{$merge->{_seqs}};
#     while (my $read_id1 = shift (@ids)) {
#       my ($this_prefix1,$dummy) = split('\.',$read_id1);
#       foreach my $read_id2 (@ids) {
#         my ($this_prefix2,$dummy) = split('\.',$read_id2);
#         next unless ($this_prefix1 eq $this_prefix2);
#         next unless (defined($merge->{_seqs}{$read_id1}{$this_prefix1}));
#         next unless (defined($merge->{_seqs}{$read_id2}{$this_prefix2}));
#         my $read1 = $merge->{_seqs}{$read_id1}{$this_prefix1};
#         my $read2 = $merge->{_seqs}{$read_id2}{$this_prefix2};
#         $nodups_merge_num++;
#         if (uc($read1) eq uc($read2)) {
#           delete $merge->{_seqs}{$read_id2};
#           delete $merge->{_prefixes}{$this_prefix1}{$read_id2};
#           $nodups_merge_num--;
#           $merge->{_prefixes}{$this_prefix1}{$read_id1}++;
#         }
#       }
#     }
#   }
#   my $num_merged_fragments = scalar keys %{$merge->{_seqs}};
#   $self->{tree}->store_tag('num_merged_fragments',$num_merged_fragments);
#   $self->{merge_hash} = $merge;


#   return 1;
# }


1;
