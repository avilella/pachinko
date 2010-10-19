#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::Simremoval

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $simremoval = Bio::EnsEMBL::Compara::RunnableDB::Simremoval->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$simremoval->fetch_input(); #reads from DB
$simremoval->run();
$simremoval->output();
$simremoval->write_output(); #writes to DB

=cut


=head1 DESCRIPTION

Simulate alignments and reads

=cut


=head1 CONTACT

  Contact Albert Vilella on module implementation/design detail: avilella@ebi.ac.uk
  Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=cut


=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Compara::RunnableDB::Simremoval;

use strict;
use Getopt::Long;
use Time::HiRes qw(time gettimeofday tv_interval);

use IO::File;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use List::Util qw(max);
use Statistics::Descriptive;
use Algorithm::Diff qw(sdiff);

use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SimpleAlign;
use Bio::EnsEMBL::Compara::Graph::NewickParser;
use Bio::EnsEMBL::Compara::SeqGrinder;
#use Bio::EnsEMBL::Compara::RunnableDB::Pachinkoalign;
#use Bio::EnsEMBL::Compara::RunnableDB::TreeBeSTMagicMirror;

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);

our @ISA = qw(Bio::EnsEMBL::Hive::Process);
#our @ISA = qw(Bio::EnsEMBL::Hive::Process Bio::EnsEMBL::Compara::RunnableDB::Pachinkoalign);


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
  $self->{huge_number} = 400000000;
  $self->{insert_size} = 500;
  $self->{tag} = '\|';

  #create a Compara::DBAdaptor which shares the same DBI handle
  #with the Pipeline::DBAdaptor that is based into this runnable
  $self->{'comparaDBA'} = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
    (
     -DBCONN=>$self->db->dbc
    );

  $self->{'clusterset_id'} = 1;

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

  $self->{mlssid} = @{$self->{mlssDBA}->fetch_all_by_method_link_type('PROTEIN_TREES')}->[0]->dbID;
  throw("no mlssid defined for PROTEIN_TREES") unless (0 < $self->{mlssid});

  $self->get_params($self->parameters);
  $self->get_params($self->input_id);

# # For long parameters, look at analysis_data
#   if($self->{analysis_data_id}) {
#     my $analysis_data_id = $self->{analysis_data_id};
#     my $analysis_data_params = $self->db->get_AnalysisDataAdaptor->fetch_by_dbID($analysis_data_id);
#     $self->get_params($analysis_data_params);
#   }

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

  foreach my $key (qw[species_set gdb fasta_dir protein_tree_id clusterset_id type analysis_data_id]) {
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
  $self->run_simremoval;

  return 0;
}

sub run_simremoval {
  my $self = shift;

  $self->{tree} = 
    $self->{treeDBA}->fetch_node_by_node_id($self->{protein_tree_id});
  my @seqs = $self->pre_clean_alignment;
  # STEP1
  $self->{rectree_aa} = $self->run_prankancestor_pre('aa');

  foreach my $this_seq (@seqs) {
    # STEP2
    my $file_tosimulate = $self->separate_seq_in_genetree($this_seq);
    next if (!defined $file_tosimulate);

    ##  foreach my $q (sort {$b <=> $a} 1..10) {
    foreach my $coverage (sort {$a <=> $b} 1..5) {
#       $self->{qscore} = $q * 4;
#       $self->{qscore} = 1 if (0 == $q);
      $self->{coverage} = $coverage;
      my $q = $coverage;
      print STDERR "# QSCORE " ,$self->{qscore},"\n";

      # STEP3
      my $grinder_file = $self->single_simulate_from_newsamples($coverage,$file_tosimulate);

      # STEP3b
      $self->{rectree_aa} = $self->run_prankancestor_pre('aa');
      $self->{rectree_cds} = $self->run_prankancestor_pre('cds');

      # STEP4
      $self->align_removal;

      next unless defined($self->{outfile});
      return 1 unless($self->parse_alignment);
      $self->define_read_placement_sets;
      #    $self->call_parse_alignment;
      next unless defined($self->{pachinko_aln});
      $self->basic_stats;
      # STEP4
      $self->run_treebest;
      print STDERR "parse...\n" if ($self->debug);
      $self->parse_treebest;
      print STDERR "generate_sam_and_ref_file_placement ...\n" if ($self->debug);
      $self->{subset_aa_ancestral_seqs} = $self->produce_list_of_ancestors_sim('sim', 2, $self->{rectree_aa}, $self->{aa_ancestral_seqs});
      # STEP5
      $self->generate_sam_and_ref_file_placement;

      # STEP6
      print STDERR "columbus_merge_fragments...\n" if ($self->debug);
      my $samfile = $self->{samfile};
      my $velvetfile = $self->columbus_merge_fragments('sim',$q) if (defined($samfile));

      # STEP7
      $self->run_bwasw('sim',$q,undef,$velvetfile) if (defined($velvetfile));

      # STEP8
      $self->run_exonerate('sim',$q) if (defined($velvetfile));

      # STEP9
      $self->run_jackhmmer($self->{transcriptsfile},$self->{pepfile});
      delete $self->{read_placement_sets};
    }
  }
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

  1;
}


##########################################
#
# internal methods
#
##########################################

sub pre_clean_alignment {
  my $self = shift;
  my $tag  = $self->{tag} || '\wRR';
  my @preclean;
  foreach my $leaf (@{$self->{tree}->get_all_leaves}) {
    if ($leaf->stable_id =~ /$tag/ || $leaf->stable_id eq '' || $leaf->source_name eq 'EXTERNALCDS' || !defined($leaf->genome_db)) {

      push @preclean, $leaf;
      # self_collapse option is only for testing purposes
    } elsif (1 == $self->{self_collapse} && defined($self->{qtaxon_list}{$leaf->taxon_id})) {
      push @preclean, $leaf;
    }
  }

  $self->{tree} = $self->{tree}->remove_nodes(\@preclean);
  $self->{tree} = $self->{tree}->minimize_tree;
  foreach my $node (@preclean) {
    $self->{treeDBA}->delete_node($node);
  }

  return 0;
}

sub separate_seq_in_genetree {
  my $self     = shift;
  my $seq = shift;

  my $node_id = $self->{tree}->node_id;

  $self->{this_tid} = undef;
  eval { 
    # We query the parent->parent of the sequence we want to delete,
    # since it's going to be the next MRCA
    $self->{this_tid} = $seq->parent->parent->get_tagvalue('taxon_id');
  };

  my @clean;
  push @clean, $seq;

  $self->{tree} = $self->{tree}->remove_nodes(\@clean);
  return undef if (!defined $self->{tree});
  $self->{tree} = $self->{tree}->minimize_tree;

  if ($self->debug && defined($self->{this_tid})) {
      print STDERR "TID = ", $self->{this_tid}, "\n";
      foreach my $node ($self->{tree}->get_all_subnodes) {
        my $tid = $node->get_tagvalue('taxon_id');
        $node->add_tag('TID',$self->{this_tid}) if ($tid eq $self->{this_tid});
      }
  }

  my $cds_seq = $seq->sequence_cds;

  my $file_tosimulate = $self->worker_temp_directory . $node_id . ".tosimulate";
  my $outseq = Bio::Seq->new
    (-display_id => $seq->member_id . '_' . $seq->taxon_id,
     -seq => $cds_seq,
     -alphabet => 'dna');
  $self->{minlength} = length($cds_seq);
  my $outio = Bio::SeqIO->new
    (-file => ">$file_tosimulate",
    -format => 'fasta');
  $outio->write_seq($outseq);
  $outio->close;

  my $aa_seq = $seq->sequence;
  my $pepfile = $self->worker_temp_directory . $node_id . ".pepfile";
  my $pepseq = Bio::Seq->new
    (-display_id => $seq->member_id . '_' . $seq->taxon_id,
     -seq => $aa_seq,
     -alphabet => 'protein');
  my $pepio = Bio::SeqIO->new
    (-file => ">$pepfile",
    -format => 'fasta');
  $pepio->write_seq($pepseq);
  $pepio->close;
  $self->{pepfile} = $pepfile;

  return $file_tosimulate;
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
  print STDERR "loaded\n" if ($self->debug);
  return $tree;
}

sub reference_alignment {
  my $self = shift;
  my $aln_f = shift;
  my $cds_f = shift;

  my $tree_id = $self->{tree}->node_id;
  my $worker_temp_directory = $self->worker_temp_directory;
  my $alignio = Bio::AlignIO->new
    (-file => "$aln_f",
     -format => "fasta");
  my $aln = $alignio->next_aln();

  my $refalnfile = $worker_temp_directory."$tree_id.cds";
  my $refaln =  Bio::AlignIO->new
    (-file => ">$refalnfile",
     -format => "fasta");

  my $sa = Bio::SimpleAlign->new();
  foreach my $seq ($aln->each_seq) {
    next if ($seq->display_id >= $self->{huge_number});
    $sa->add_seq($seq);
  }
  $sa = $sa->remove_gaps(undef,1);
  $sa->set_displayname_flat(1);
  $refaln->write_aln($sa);
  $refaln->close;

  $self->{gtsa} = $sa;

  my $seqio = Bio::SeqIO->new
    (-file => "$cds_f",
     -format => "fasta");

  my $newsamplesfile = $worker_temp_directory."$tree_id.newsamples";
  my $newsamplesio =  Bio::SeqIO->new
    (-file => ">$newsamplesfile",
     -format => "fasta");

  my $minlength = 99999999;
  while (my $seq = $seqio->next_seq) {
    next unless ($seq->display_id >= $self->{huge_number});
    $minlength = length($seq->seq) if ($minlength > length($seq->seq));

    $self->{seqs}{$seq->display_id} = $seq->seq;
    $newsamplesio->write_seq($seq);
   }
  $newsamplesio->close;

  $self->{refaln} = $refalnfile;
  $self->{newsamplesfile} = $newsamplesfile;
  $self->{minlength} = $minlength;

  return 1;
}


sub single_simulate_from_newsamples {
  my $self = shift;
  my $coverage = shift || 1;
  my $newsamplesfile = shift;

  my $seq_length = 75;
  my $insert_size = $self->{insert_size};

  if ($insert_size > $self->{minlength}) {
      $insert_size = $self->{minlength};
      print STDERR "insert_size $insert_size\n";
  }
  my $mate_epsil = 50;

  my $node_id = $self->{tree}->node_id;
  my $worker_temp_directory = $self->worker_temp_directory;

  print STDERR "Run Grinder ...\n" if ($self->debug);
  my $ret = Bio::EnsEMBL::Compara::SeqGrinder->Grinder
    ($newsamplesfile, #fastafile
     1,               #indep
     undef,           #nof_seqs
     $coverage,       #coverage
     $seq_length,     #seq_length
     0,               #seq_delta
     $insert_size,    #mate pair length (insert size)
     $mate_epsil,     #mate epsil
     'powerlaw',      #distrib
     1, 
     1, 
     0, 
     0, 
     $node_id, 
     $node_id, 
     $worker_temp_directory
    );

  return 1 if (2 == $ret); # Didnt work, maybe too short
  print STDERR "Done.\n" if ($self->debug);

#  my $simfastqfile = $worker_temp_directory . "$node_id.sim.fastq";
  my $simfastqfile = $worker_temp_directory . "$node_id.reads";

  my $seqgrindersamplesfile = $worker_temp_directory . $node_id . "-shotgun.fa";

  # Concatenate the reads and give them to simNGS to process with the -p option
  my $seqio = Bio::SeqIO->new
    (-format => 'fasta',
     -file => $seqgrindersamplesfile);
  my $pairs;
  while (my $seq = $seqio->next_seq) {
    my $root_display_id = $seq->display_id; $root_display_id =~ s/\/[12]$//g;
    push @{$pairs->{$root_display_id}}, $seq;
  }
  my $concatsamplesfile = $worker_temp_directory . $node_id . ".concat.fa";
  my $concatio = Bio::SeqIO->new
    (-format => 'fasta',
     -file => ">$concatsamplesfile");
  foreach my $root_display_id (keys %$pairs) {
    my ($p1seq,$p2seq) = @{$pairs->{$root_display_id}};

    my $concatstring = $p1seq->seq . $p2seq->seq;
    my $concatseq = Bio::Seq->new
      (-display_id => $root_display_id,
       -description => $p2seq->description,
       -seq => $concatstring,
       -alphabet => 'dna');
    $concatio->write_seq($concatseq);
  }
  $concatio->close;

  my $cmd = "/nfs/users/nfs_a/avilella/src/simNGS/latest/simNGS/bin/simNGS ";
  $cmd   .= " -o fastq";
  $cmd   .= " -p"; #PE
  $cmd   .= " -n $seq_length"; # number of cycles
  $cmd   .= " /nfs/users/nfs_a/avilella/src/simNGS/latest/s_4_0033.runfile";
  $cmd   .= " > $simfastqfile";
  $cmd    = " cat $concatsamplesfile | $cmd";
  print STDERR "# $cmd\n" if ($self->debug);

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  unless(system("cd $worker_temp_directory; $cmd") == 0) {
    print("## $cmd\n"); throw("error running simngs, $!\n");
  }
  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

  throw ("No simfastqfile [$simfastqfile]") if (-z $simfastqfile);


  # FIXME do properly
#   my @tids = keys %{$self->{all_tids}};
#   my $this_tid = $tids[0];
  my $this_tid = $self->{this_tid} || "0";

  my $format = 'fastq';
  my $in = Bio::SeqIO->new(-format => $format,
                            -file => "$simfastqfile");
  my $readsfilename = $self->worker_temp_directory . $node_id . ".simreads";
  my $out = Bio::SeqIO->new(-format => $format,
                            -file => ">$readsfilename");
  my @ranked_list;
  # Also store sum_qscores just for QC
  while (my $this_read = $in->next_seq) {

    my $sum_qscores; map {$sum_qscores += $_} @{$this_read->qual};
    $self->{sum_qscores}{$this_read->display_id} = $sum_qscores;
    push @ranked_list, [$sum_qscores, $this_read];
  }
  @ranked_list = sort { $b->[0] <=> $a->[0] } @ranked_list;
  @ranked_list = map { $_->[1] } @ranked_list;

  my $stats = new Statistics::Descriptive::Full;
  my @values = values %{$self->{sum_qscores}};
  $stats->add_data(\@values);
  my $median_sum_qscores = $stats->median;
  $self->{median_sum_qscores} = $median_sum_qscores;

  foreach my $read (@ranked_list) {
    my @tempseq1 = split("",$read->seq); my $p1string = join('',splice(@tempseq1,                   0,length($read->seq)/2));
    my @tempseq2 = split("",$read->seq); my $p2string = join('',splice(@tempseq2,length($read->seq)/2,length($read->seq)  ));
    my @tempqual1 = @{$read->qual}; my @p1qual = splice(@tempqual1,                   0,length($read->seq)/2);
    my @tempqual2 = @{$read->qual}; my @p2qual = splice(@tempqual2,        length($read->seq)/2,length($read->seq)  );
    my $revcomp2string = Bio::Seq->new(-seq => $p2string)->revcom->seq;
    my $p1qual_text = join(' ',@p1qual);
    my $revcomp2qual_text = join(' ',reverse(@p2qual));
    my $p1_desc = $read->description; $p1_desc =~ s/\ P2\=\d+\-\d+//;
    my $p1seq = Bio::Seq::Quality->new
      (-display_id => $read->display_id . "/1",
       -description => $p1_desc,
       -qual => $p1qual_text,
       -seq => $p1string,
       -alphabet => 'dna');
    $p1seq->description(" TID=$this_tid " . $p1seq->description);
    $out->write_seq($p1seq);

    my $p2seq = Bio::Seq::Quality->new
      (-display_id => $read->display_id . "/2",
       -description => $read->description,
       -qual => $revcomp2qual_text,
       -seq => $revcomp2string,
       -alphabet => 'dna');
    $p2seq->description(" TID=$this_tid " . $p2seq->description);
    $out->write_seq($p2seq);
  }
  $out->close;

  $self->{simfastqfile} = $readsfilename;

  return 1;
}

sub align_removal {
  my $self = shift;

  my $node_id = $self->{tree}->node_id;
  return 1 if (!defined($self->{simfastqfile}));

  my $readsfilename = $self->{simfastqfile};
  my $cds_seqfile = $self->{'prank_pre.cds'};

  my $nhx_string = $self->{tree}->nhx_format('member_id_taxon_id');
  my $nhx_treefile = $self->worker_temp_directory . "$node_id.nhx";
  open NHX, ">$nhx_treefile" or die "$!";
  print NHX "$nhx_string\n";
  close NHX;
#  my $nhx_treefile = $self->{'prank_pre.cds.nh'};

#  my $pachinko_executable = "/nfs/users/nfs_a/avilella/src/pagan/latest/pagan";
  my $pachinko_executable = "/nfs/users/nfs_a/avilella/src/pagan/latest/pagan/pagan";

  my $outfile_prefix = $self->worker_temp_directory . $node_id . ".pagan";
  my $cmd = "$pachinko_executable ";
  # ~/src/pagan/latest/test_git/pagan --cds-seqfile $node_id.cds.mfa --cds-treefile $node_id.nh --readsfile $node_id.reads --outfile $node_id
#  $cmd .= " --min-reads-overlap 0.75 "; # default 0.1
#  $cmd .= " --reads-distance 0.75 ";  # default 0.01
#  $cmd .= " --454 "; # only for homopolymer issues
  $cmd .= " --pair-end ";  # PEs
  $cmd .= " --output-nhx-tree";
  my $this_qscore = $self->{qscore} || 20;
  $cmd .= " --qscore-minimum $this_qscore";
  $cmd .= " --trim-read-ends --rank-reads-for-nodes";
#  $cmd .= " --discard-overlapping-identical-reads";
#  $cmd .= " --trim-mean-qscore arg 15";       # sliding window average Q-score to be
#  $cmd .= " --trim-window-width arg 5";       # sliding window width for trimming
#  $cmd .= " --minimum-trimmed-length 20";     # minimum trimmed read length
#  $cmd .= " --allow-skip-low-qscore --align-bad-reads-at-root";
  $cmd .= " --cds-seqfile $cds_seqfile";
  $cmd .= " --cds-treefile $nhx_treefile";
  $cmd .= " --readsfile $readsfilename";
  $cmd .= " --outfile $outfile_prefix";
  $cmd .= " --placement-file $outfile_prefix.placement";
  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  my $starttime = time();
  print STDERR "pachinkoalign ...\n" if ($self->debug);
  my $worker_temp_directory = $self->worker_temp_directory;
  #  $DB::single=1;1;
  print("$cmd\n") if ($self->debug);

  unless(system("cd $worker_temp_directory; $cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running pachinkoalign $!\n");
  }
  $self->{runtime_msec} = time() - $starttime;
  print STDERR $self->worker_temp_directory,"\n";
  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

  $self->{outfile} = "$worker_temp_directory"."$node_id.pagan.fas";
  $self->{xmlfile} = "$worker_temp_directory"."$node_id.pagan.xml";
  $self->{nhx_tree_file} = "$worker_temp_directory"."$node_id.pagan.nhx_tree";
  $self->{placement_file} = "$outfile_prefix.placement.tsv";

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

  $DB::deep = 999999; # This is only useful for debugging purposes, causes no harm
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

  $self->{gtsa} = $aln;
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

sub parse_tree_in_xml {
  my $self = shift;
  my $xmlfile = shift;

  $DB::deep = 999999; # This is only useful for debugging purposes, causes no harm
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

sub run_prankancestor_pre_sim {
  my $self = shift;
  my $type = shift;

  my $node_id = $self->{tree}->node_id;


  # PRE SIM: USING THE SIMULATED ALIGNMENT

  my $aln_name = $self->{aln_f};

  # PRE EXPANDED NEWICK TREE
  my $treefile = $self->{pre_treefile};
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

sub produce_list_of_ancestors_sim {
  my $self    = shift;
  my $this_node_id = shift;
  my $depth   = shift;
  my $rectree = shift;
  my $aa_ancestral_seqs = shift;

    my $subset_aa_ancestral_seqs;
  my $tree = $self->parse_tree_string($rectree);
  foreach my $this_node_id (keys %{$self->{read_placement_sets}}) {
    my $this_node = $tree->find_node_by_name($this_node_id);
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
    my $fragment_member = Bio::EnsEMBL::Compara::Member->new
      (
       -stable_id    => $stable_id_name,
       -source_name  => 'EXTERNALCDS',
       -taxon_id     => 0,
       -genome_db_id => 999999,
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
  my $tree_key = "prank_pre.$type.nh";
  $self->{$tree_key} = $treefile;
  my $cds_key = "prank_pre.$type";
  $self->{$cds_key} = $aln_name;

  foreach my $subnode ($self->{tree}->get_all_subnodes) {
    $subnode->distance_to_parent(0.0001) if (0 == $subnode->distance_to_parent);
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
  my $tleavesfile = $self->worker_temp_directory . "sim.1.$node_id.tleaves.fasta";

  my $inputio = Bio::SeqIO->new
    (-format => 'fastq',
#     -file => $velvetfile);
     -file => $consfile);
  my $outio = Bio::SeqIO->new
    (-format => 'fasta',
#     -file => $velvetfile);
     -file => ">$consfile.fasta");

  while (my $seq  = $inputio->next_seq) {
    $self->{exonerate_seqs}{$seq->display_id} = $seq->seq;
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
  # $DB::single=1;1;
  my $fh;
  eval { open($fh, "$cmd |") || die $!; };
  if ($@) {
    warn("problem with exonerate $@ $!");
    return;
  }

  my $starttime = time();
  print STDERR "exonerate...\n" if ($self->debug);

#  $self->create_template_objects;
  # $DB::single=1;1;
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
    $this_seq = $trseq unless ($transcript->strand eq -1);
    $this_seq = $trseq->revcom if ($transcript->strand eq -1);
    # Store transcript
    # print STDERR "Storing transcript\n" if ($self->debug);
    my $stable_id_name = $node_id . "_i" . $id1 . "_t" . $this_tid . "_" . $transcript->seqname;
    $stable_id_name =~ s/\#//g;

    my $translated_this_seq = $this_seq->translate;
    print TRANSCRIPTSFILE ">" . $transcript->seqname , "\n", $translated_this_seq->seq, "\n";
    $this_id->{$transcript->seqname} = 1;
    # my $trseq = $transcript->translation->end_Exon->seq->translate->seq;
  }

  close TRANSCRIPTSFILE;
  $self->{transcriptsfile} = $transcriptsfile;
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
                unless ($q_length > 0 && $t_length > 0) {
                  print STDERR "Attempting exonerate $q_id\t$t_id\n" if ($self->debug);
                  next;
                }

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

                my $coverage = 0;
                $coverage = sprintf("%.2f", 100 * $covered_count / $q_length) if ($covered_count > 0 && $q_length > 0);

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


sub prepare_fragments_for_assessment {
  my $self = shift;
  my $hash = $self->{merge_hash};

  my $node_id = $self->{tree}->node_id;

  my $pepfile = $self->worker_temp_directory . "$node_id.pep";
  my $pepio = Bio::SeqIO->new
    (-file => ">$pepfile",
     -format => 'fasta',
    );
  my $seq_ids;
  my $count = 1;
  foreach my $stable_id (keys %{$hash->{_seqs}}) {
    my @prefixes = keys %{$hash->{_seqs}{$stable_id}};
    my $prefix = $prefixes[0];
    next unless (defined($prefix));
    my $alignment_string = $hash->{_seqs}{$stable_id}{$prefix};
    my $fragment_cds = $alignment_string; $fragment_cds =~ s/\-//g;
#    my $fragment_cds = $alignment_string; $fragment_cds =~ s/\-//g;
#     my $taxon_id;
#     if ($prefix =~ /(\d+)\_(\d+)/) {
#       $taxon_id = $2;
#     }
    my $is_merged = 0;
    $is_merged = 1 if ($alignment_string =~ /[atcg]/);
    my $stable_id_name = $stable_id;
    my $stable_id_name = $prefix . "." . $count++;
    $stable_id_name .= ".M" if ($is_merged);

    my $seq_id = $stable_id_name;
    my $this_seq = $fragment_cds;
    my $cseq = Bio::Seq->new
      (-display_id => $seq_id,
       -seq => $this_seq,
       -alphabet => 'dna');
    # last if ($count++ > 3 && $self->debug); # FIXME only debug

    $self->{seqs}{$cseq->display_id} = $this_seq;
    $pepio->write_seq($cseq);

  }
  $pepio->close;

  $self->{pepfile} = $pepfile;
  return 1;
}


sub run_magicmirror {
  my $self = shift;
  my $pepfile       = shift || $self->{pepfile};
  my $referencefile = shift || $self->{referencefile};
  my $node_id = $self->{tree}->node_id;

  my $magicmirror_executable = $self->analysis->program_file;
  $magicmirror_executable = "/nfs/users/nfs_a/avilella/src/exonerate/latest/exonerate-2.2.0-x86_64/bin/exonerate";

  my $model = 'cdna2genome';
  $self->{query_type} = 'dna';
  my $worker_temp_directory = $self->worker_temp_directory;
  my $resultsfile = $worker_temp_directory . $node_id . ".$model";

  my $cmd = "$magicmirror_executable ";
  $cmd .= " --model $model";
  $cmd .= " --useaatla FALSE --alignmentwidth 200 --showalignment" if ($self->debug);
  # $cmd .= " --frameshift -1";
  $cmd .= " --forwardcoordinates FALSE";
  $cmd .= " --exhaustive FALSE";
  $cmd .= " --score 50 --dnahspthreshold 60 --dnawordlen 12";
  $cmd .= " --query $pepfile";
  $cmd .= " --target $referencefile ";
  # $cmd .= " > $resultsfile ";
  # print "# $cmd \n" if ($self->debug);

  my $fh;
  eval { open($fh, "$cmd |") || die $!; };
  if ($@) {
    warn("problem with magicmirror $@ $!");
    return;
  }

  my $starttime = time();
  # print STDERR "magicmirror...\n" if ($self->debug);

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  my $transcripts = $self->parse_results($fh);
  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

  my $runtime = time()*1000-$starttime;
  # printf("  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
  $self->{'tree'}->store_tag('magicmirror_runtime_msec', $runtime);

  my $this_id;
  my $worker_temp_directory = $self->worker_temp_directory;
  my $transcriptsfile = $worker_temp_directory . $node_id;
  $transcriptsfile .= ".c" . $self->{coverage} if (defined($self->{coverage}));
  $transcriptsfile .= ".q" . $self->{qscore} if (defined($self->{qscore}));
  $transcriptsfile .= ".transcripts";

  open TRANSCRIPTSFILE,">$transcriptsfile" or die $!;
  foreach my $transcript (@$transcripts) {
    my $seq = $transcript->seq;
    my $this_seq;
    $this_seq = $seq->seq unless ($transcript->strand eq -1);
    $this_seq = $seq->revcom->seq if ($transcript->strand eq -1);
    my $bioperl_seq = Bio::Seq->new
      (-display_id => $transcript->seqname,
       -seq => $this_seq,
       -alphabet => 'dna');
    my $tseq = $bioperl_seq->translate;
    next if (defined($this_id->{$bioperl_seq->display_id}));
    print TRANSCRIPTSFILE ">" . $transcript->seqname , "\n", $this_seq, "\n";
    #    print ">" . $bioperl_seq->display_id , "\n", $bioperl_seq->translate->seq, "\n" if ($self->debug);
    $this_id->{$bioperl_seq->display_id} = 1;
    # my $trseq = $transcript->translation->end_Exon->seq->translate->seq;
  }

  close TRANSCRIPTSFILE;

  return $transcriptsfile;
}

# sub parse_results {
# 	my ($self, $fh) = @_;
# 	my %strand_lookup = ( '+' => 1, '-' => -1, '.' => 1 );

#         my @transcripts;
        
#         # Parse output looking for lines beginning with 'RESULT:'.
#         # Each line represents a distinct match to one sequence
#         # containing multiple 'exons'.

# 	while (<$fh>){
#           # next unless /^RESULT:/;
# 		next unless /^vulgar\:/;
#                 # print STDERR "$_\n" if ($self->debug);
#                 chomp;

#                 my ($perc_id, $q_length, $t_length, $gene_orientation);
#                 # my ($tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score, $perc_id, $q_length, $t_length, $gene_orientation,
#                 my ($tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score,
#                     @align_components) = split;
#                 $q_length = length($self->{seqs}{$q_id});
#                 $t_length = length($self->{seqs}{$t_id});

#                 $t_strand = $strand_lookup{$t_strand};
#                 $q_strand = $strand_lookup{$q_strand};
#                 $gene_orientation = $strand_lookup{$gene_orientation};

#                 # Read vulgar information and extract exon regions.
#                 my $exons = $self->_parse_vulgar_block($t_start,
#                                                        $t_end,
#                                                        $t_strand,
#                                                        $t_length,
#                                                        $q_start, 
#                                                        $q_end,
#                                                        $q_strand,
#                                                        $q_length,
#                                                        \@align_components);

#                 # now we have extracted the exons and the coordinates are with 
#                 # reference to the forward strand of the query and target, we can 
#                 # use the gene_orienation to flip the strands if necessary
#                 if ($gene_orientation == -1 and $t_strand == 1) {
#                   $t_strand *= -1;
#                   $q_strand *= -1;
#                 }

#                 my $covered_count = 0;
#                 if ($self->{coverage_as_proportion_of_aligned_residues}) {
#                   foreach my $exon (@$exons) {
#                     foreach my $sf (@{$exon->{sf}}) {
#                       $covered_count += $sf->{query_end} - $sf->{query_start} + 1;
#                     }
#                   }
#                 } else {
#                   $covered_count = abs($q_end - $q_start);
#                 }

#                 my $coverage = sprintf("%.2f", 100 * $covered_count / $q_length);

#                 # Build FeaturePairs for each region of query aligned to a single
#                 # Exon.  Create a DnaDnaAlignFeature from these FeaturePairs and then
#                 # attach this to our Exon.
#                 my $transcript = Bio::EnsEMBL::Transcript->new();

#                 my $query_sequence = $self->{seqs}{$q_id};
#                 my $cs = Bio::EnsEMBL::CoordSystem->new(-NAME    => $q_id, -RANK    => 1);
#                 my $slice = Bio::EnsEMBL::Slice->new
#                   (-seq_region_name  => 'misc',
#                    -name             => $q_id,
#                    -seq_region_length => $q_length,
#                    -start            => 1,
#                    -end              => $q_length,
#                    -strand           => $q_strand,
#                    -coord_system     => $cs,
#                    -seq              => $query_sequence,
#                   );

#                 my (@tran_feature_pairs, 
#                     $cds_start_exon, $cds_start,
#                     $cds_end_exon, $cds_end);

#                 foreach my $proto_exon (@$exons){

#                   # Build our exon and set its key values.
#                   my $exon = Bio::EnsEMBL::Exon->new();

#                   $exon->seqname($t_id);
#                   $exon->start($proto_exon->{exon_start});
#                   $exon->end($proto_exon->{exon_end});
#                   $exon->phase($proto_exon->{phase});
#                   $exon->end_phase($proto_exon->{end_phase});
#                   $exon->strand($t_strand);

#                   my @exon_feature_pairs;
#                   foreach my $sf (@{$proto_exon->{sf}}){
#                     my $feature_pair = Bio::EnsEMBL::FeaturePair->new
#                       (-seqname    => $t_id,
#                        -start      => $sf->{target_start},
#                        -end        => $sf->{target_end},
#                        -strand     => $t_strand,
#                        -hseqname   => $q_id,
#                        -hstart     => $sf->{query_start},
#                        -hend       => $sf->{query_end},
#                        -hstrand    => $q_strand,
#                        -hcoverage  => $coverage,
#                        -score      => $coverage,
#                        -percent_id => $perc_id);


#                     push @exon_feature_pairs, $feature_pair;
#                     push @tran_feature_pairs, $feature_pair;
#                   }

#                   if (@exon_feature_pairs) {
#                     # Use our feature pairs for this exon to create a single 
#                     # supporting feature (with cigar line).
#                     my $supp_feature;

#                     eval{
#                       if ($self->{query_type} eq 'protein') {
#                         $supp_feature =
#                           Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@exon_feature_pairs);
#                       } else {
#                         $supp_feature = 
#                           Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@exon_feature_pairs);
#                       }
#                     };
#                     if ($@) {
#                       warn($@);
#                       next TRANSCRIPT;
#                     }      
#                     $exon->add_supporting_features($supp_feature);
#                   }      

#                   if (exists $proto_exon->{coding_start}) {
#                     if (not defined $cds_start_exon) {
#                       $cds_start_exon = $exon;
#                       $cds_start = $proto_exon->{coding_start};
#                     }
#                     $cds_end_exon = $exon;
#                     $cds_end = $proto_exon->{coding_end};
#                   }


#                   $transcript->add_Exon($exon);
#                 }

#                 # Create a single supporting feature for the whole transcript
#                 my $t_supp_feat;
#                 eval{
#                   if ($self->{query_type} eq 'protein') {
#                     $t_supp_feat =
#                       Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@tran_feature_pairs);
#                   } else {
#                     $t_supp_feat = 
#                       Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@tran_feature_pairs);
#                   }
#                 };
#                 if ($@) {
#                   warn("Could not create Transcript supporting feature");
#                 } else {
#                   $transcript->add_supporting_features($t_supp_feat);
#                 }

#                 my @exons = @{$transcript->get_all_Exons};

#                 if (scalar(@exons)) {
#                   if (defined $cds_start_exon) {
#                     my $translation = Bio::EnsEMBL::Translation->new();
#                     $translation->start_Exon($cds_start_exon);
#                     $translation->start($cds_start);
#                     $translation->end_Exon($cds_end_exon);
#                     $translation->end($cds_end);

#                     $transcript->translation($translation);
#                   }

#                   calculate_exon_phases($transcript, 0);
#                   $transcript->slice($slice);
#                   foreach my $sf(@{$transcript->get_all_supporting_features}){
#                     $sf->slice($slice);
#                   }
#                   foreach my $exon(@{$transcript->get_all_Exons}){
#                     $exon->slice($slice);
#                     foreach my $sf(@{$exon->get_all_supporting_features}){
#                       $sf->slice($slice);
#                     }
#                   }
# #                   $transcript = $transcript->transfer($t_slice);

#                   push @transcripts, [$score, $transcript];
# #                  last if (5 > scalar @transcripts); #FIXME debug
# #                   my $seq = $transcript->seq;
# #                   my $this_seq;
# #                   $this_seq = $seq->seq unless ($transcript->strand eq -1);
# #                   $this_seq = $seq->revcom->seq if ($transcript->strand eq -1);
# #                   my $bioperl_seq = Bio::Seq->new
# #                     (-display_id => $transcript->seqname,
# #                      -seq => $this_seq,
# #                      -alphabet => 'dna');
# #                   my $tseq = $bioperl_seq->translate;
#                   # next if (defined($self->{transcript_ids}->{$bioperl_seq->display_id}));
# #                  $DB::single=1;1;
# #                  eval {require Bio::EnsEMBL::Compara::Member;};
# #                   #if ($score > 500 && $tseq->seq != /\*/) {
# #                     my $transcript_member = Bio::EnsEMBL::Compara::Member->new
# #                       (
# #                        -stable_id    => $transcript->seqname,
# #                        -source_name  => 'EXTERNALTRANS',
# #                        -taxon_id     => $taxon_id,
# #                        -genome_db_id => $genome_db_id,
# #                       );
# #                     $fragment_member->sequence($fragment_cds);
# #                     $self->{memberDBA}->store($fragment_member);

# #                     print STDERR ">" . $bioperl_seq->display_id , "\n", $tseq->seq, "\n";
# #                     $self->{transcript_ids}{$bioperl_seq->display_id} = 1;
#                   #}
#                 }
#               }

#         @transcripts = sort { $b->[0] <=> $a->[0] } @transcripts;
#         @transcripts = map { $_->[1] } @transcripts;

#         return \@transcripts;
# }


# sub _parse_vulgar_block {
#   my ($self, 
#       $target_start, $target_end, $target_strand, $target_length,
#       $query_start, $query_end,  $query_strand, $query_length,
#       $vulgar_components) = @_;

#   # This method works along the length of a vulgar line 
#   # exon-by-exon.  Matches that comprise an exon are 
#   # grouped and an array of 'proto-exons' is returned.
#   # Coordinates from the vulgar line are extrapolated 
#   # to actual genomic/query coordinates.

#   my @exons;
#   my $exon_number = 0;

#   # We sometimes need to increment all our start coordinates. Exonerate 
#   # has a coordinate scheme that counts _between_ nucleotides at the start.
#   # However, for reverse strand matches 
  
#   my ($query_in_forward_coords, $target_in_forward_coords);
#   my ($cumulative_query_coord, $cumulative_target_coord);

#   if ($target_start > $target_end) {
#     warn("For target, start and end are in thew wrong order for a reverse strand match")
#         if $target_strand != -1;
#     $cumulative_target_coord = $target_start;
#     $target_in_forward_coords = 1;
#   } else {
#     $cumulative_target_coord = $target_start + 1;
#     $target_in_forward_coords = 0;
#   }
#   if ($query_start > $query_end) {
#     warn("For query, start and end are in thew wrong order for a reverse strand match")
#         if $query_strand != -1;
#     $cumulative_query_coord = $query_start;
#     $query_in_forward_coords = 1;
#   } else {
#     $cumulative_query_coord = $query_start + 1;
#     $query_in_forward_coords = 0;
#   }


#   while (@$vulgar_components){
#     throw("Something funny has happened to the input vulgar string." .
# 		 "  Expecting components in multiples of three, but only have [" .
# 		 scalar @$vulgar_components . "] items left to process.")
#       unless scalar @$vulgar_components >= 3;

#     my $type                = shift @$vulgar_components;
#     my $query_match_length  = shift @$vulgar_components;
#     my $target_match_length = shift @$vulgar_components;

#     throw("Vulgar string does not start with a match.  Was not " . 
# 		 "expecting this. (Have type $type)")
#       if (scalar @exons == 0) && $type ne 'M' && $type ne 'S' && $type ne 'C';

#     if ($type eq 'M' or $type eq 'S' or $type eq 'C'){
#       my %hash = (type => $type);

#       if ($target_strand == -1) {
#         if ($target_in_forward_coords) {
#           $hash{target_start} = $cumulative_target_coord - ($target_match_length - 1);
#           $hash{target_end}   = $cumulative_target_coord;
#         } else {
#           $hash{target_end}   = $target_length - ($cumulative_target_coord - 1);
#           $hash{target_start} = $hash{target_end} - ($target_match_length - 1);
#         }
#       } else {
#         $hash{target_start} = $cumulative_target_coord;
#         $hash{target_end}   = $cumulative_target_coord + ($target_match_length - 1);
#       }

#       if ($query_strand == -1) {
#         if ($query_in_forward_coords) {
#           $hash{query_start} = $cumulative_query_coord - ($query_match_length - 1);
#           $hash{query_end}   = $cumulative_query_coord;
#         } else {
#           $hash{query_end}   = $query_length - ($cumulative_query_coord - 1);
#           $hash{query_start} = $hash{query_end} - ($query_match_length - 1);
#         }
#       } else {
#         $hash{query_start} = $cumulative_query_coord;
#         $hash{query_end}   = $cumulative_query_coord + ($query_match_length - 1);
#       }

#       # there is nothing to add if this is the last state of the exon
#       $exons[$exon_number]->{gap_end}   = 0;
#       push @{$exons[$exon_number]->{sf}}, \%hash;
#     }
#     elsif ($type eq "G") {
#       if (exists($exons[$exon_number]->{sf})) {
#         # this is the gap in the middle of an exon, or at the end. Assume it is 
#         # at the end, and then reset if we see another match state in this exon
#         $exons[$exon_number]->{gap_end}   = $target_match_length;
#       } else {
#         # this is a gap at the start of an exon; 
#         $exons[$exon_number]->{gap_start} = $target_match_length;
#       }
#     }
#     elsif ($type eq "I" or
#            $type eq "F") {

#       # in protein mode, any insertion on the genomic side should be treated as 
#       # an intron to ensure that the result translates. However, we allow for
#       # codon insertions in the genomic sequence with respect to the protein. 
#       # This introduces the possibility of in-frame stops, but I don't
#       # think "introning over" these insertions is appropriate here. 

#       # if we see a gap/intron immediately after an intron, the current exon is "empty"
#       if ($exons[$exon_number]) {
#         $exon_number++;
#       }
#     }

#     if ($target_in_forward_coords and $target_strand == -1) {
#       $cumulative_target_coord -= $target_match_length;
#     } else {
#       $cumulative_target_coord += $target_match_length;
#     }
#     if ($query_in_forward_coords and $query_strand == -1) {
#       $cumulative_query_coord  -= $query_match_length;
#     }
#     else {
#       $cumulative_query_coord  += $query_match_length;
#     }

#   }

#   for(my $i = 0; $i < @exons; $i++) {
#     my $ex = $exons[$i];
#     my @ex_sf = @{$ex->{sf}};

#     if ($target_strand == -1) {
#       $ex->{exon_start} = $ex_sf[-1]->{target_start};
#       $ex->{exon_end}   = $ex_sf[0]->{target_end};

#       if (exists $ex->{gap_start}) {
#         $ex->{exon_end} += $ex->{gap_start};
#       }
#       if (exists $ex->{gap_end}) {
#         $ex->{exon_start} -= $ex->{gap_end};
#       }

#     } else {
#       $ex->{exon_start} = $ex_sf[0]->{target_start};
#       $ex->{exon_end}   = $ex_sf[-1]->{target_end};

#       if (exists $ex->{gap_start}) {
#         $ex->{exon_start} -= $ex->{gap_start};
#       }
#       if (exists $ex->{gap_end}) {
#         $ex->{exon_end} += $ex->{gap_end};
#       }
#     }

#     # split codons are a pain. If the query is dna, we must be in the 
#     # cdna2genome model so they must be part of the supporting feature. 
#     # If query is protein, they need to be removed from the supporting feature

#     if ($self->{query_type} eq 'dna') {
#       map { $_->{type} = 'C' if $_->{type} eq 'S' } @ex_sf;
#       if (my @cod = grep { $_->{type} eq 'C' } @ex_sf) {
#         # at least part of this exon is coding

#         @cod = sort { $a->{target_start} <=> $b->{target_start} } @cod;
#         my $cod_start = $cod[0]->{target_start};
#         my $cod_end   = $cod[-1]->{target_end};
#         my $cod_len   = $cod_end - $cod_start + 1;

#         if ($target_strand == -1) {
#           $ex->{coding_start} = $ex->{exon_end} - $cod_end + 1;
#         } else {
#           $ex->{coding_start} = $cod_start - $ex->{exon_start} + 1;
#         }

#         $ex->{coding_end}   = $ex->{coding_start} + $cod_len - 1;
#       }
#       # merge together abutting ungapped features
#       my @nr_sf;
#       foreach my $sf (@ex_sf) {
#         my $merged = 0;
#         if (@nr_sf and
#             $sf->{type} eq 'C' and
#             $nr_sf[-1]->{type} eq 'C' and
#             $sf->{query_start} == $nr_sf[-1]->{query_end} + 1) {
#           if ($target_strand == -1) {
#             if ($sf->{target_end} == $nr_sf[-1]->{target_start} - 1) {
#               $nr_sf[-1]->{target_start} = $sf->{target_start};
#               $nr_sf[-1]->{query_end}  = $sf->{query_end};
#               $merged = 1;
#             }
#           } else {
#             if ($sf->{target_start} == $nr_sf[-1]->{target_end} + 1) {
#               $nr_sf[-1]->{target_end} = $sf->{target_end};
#               $nr_sf[-1]->{query_end}  = $sf->{query_end};
#               $merged = 1;
#             }
#           }         
#         }
#         if (not $merged) {
#           push @nr_sf, $sf;
#         }
#       }
#       @ex_sf = @nr_sf;
#     } else {
#       # query type is protein
#       @ex_sf = grep { $_->{type} ne 'S' } @ex_sf;
#       $ex->{coding_start} = 1;
#       $ex->{coding_end} = $ex->{exon_end} - $ex->{exon_start} + 1;
#     }

#     $ex->{sf} = \@ex_sf;
#   }

#   return \@exons;
# }

# sub calculate_exon_phases {
#   my ($transcript, $start_phase) = @_;

#   #$transcript->sort; Commented out as it has been deprecated. Transcript are now sorted by default.
#   foreach my $e (@{$transcript->get_all_Exons}) {
#     $e->phase(-1);
#     $e->end_phase(-1);
#   }

#   if ($transcript->translation) {
#     my $tr = $transcript->translation;

#     my @exons = @{$transcript->get_all_Exons};

#     while($exons[0] != $tr->start_Exon) {
#       shift @exons;
#     }
#     while($exons[-1] != $tr->end_Exon) {
#       pop @exons;
#     }
    
#     # set phase of for first coding exon
#     my $cds_len = $exons[0]->length;
#     if ($tr->start == 1) {
#       $exons[0]->phase($start_phase);
#       $cds_len += $start_phase;
#     } else {
#       $cds_len -= ($tr->start - 1);
#     }
#     $exons[0]->end_phase($cds_len % 3);

#     # set phase for internal coding exons      
#     for(my $i=1; $i < @exons; $i++) {
#       $exons[$i]->phase($exons[$i-1]->end_phase);
#       $exons[$i]->end_phase(($exons[$i]->length + $exons[$i]->phase) % 3);
#     }
        
#     # set phase for last coding exon
#     if ($exons[-1]->length > $tr->end) {
#       $exons[-1]->end_phase(-1);
#     }
#   }
# }

sub pre_clean_alignment {
  my $self = shift;
  my $tag  = $self->{tag} || '\wRR';
  my @preclean;
  foreach my $leaf (@{$self->{tree}->get_all_leaves}) {
    if ($leaf->stable_id =~ /$tag/ || $leaf->stable_id eq '' || $leaf->source_name eq 'EXTERNALCDS' || !defined($leaf->genome_db)) {

      push @preclean, $leaf;
      # self_collapse option is only for testing purposes
    } elsif (1 == $self->{self_collapse} && defined($self->{qtaxon_list}{$leaf->taxon_id})) {
      push @preclean, $leaf;
    }
  }

  $self->{tree} = $self->{tree}->remove_nodes(\@preclean);
  $self->{tree} = $self->{tree}->minimize_tree;
  foreach my $node (@preclean) {
    $self->{treeDBA}->delete_node($node);
  }

  return @{$self->{tree}->get_all_leaves};
}

sub run_jackhmmer {
  my $self = shift;
  my $transcriptsfile = shift || $self->{transcriptsfile};
  my $pepfile         = shift || $self->{pepfile};
  my $node_id = $self->{tree}->node_id;

  print STDERR "run_jackhmmer...\n" if ($self->debug);
  my $jackhmmer_executable = $self->analysis->program_file;
  $jackhmmer_executable = "/software/ensembl/compara/hmmer3/hmmer-3.0/src/jackhmmer";
#    unless (-e $jackhmmer_executable);

  my $worker_temp_directory = $self->worker_temp_directory;
  my $resultsfile = $worker_temp_directory . $node_id;
  $resultsfile .= ".c" . $self->{coverage} if (defined($self->{coverage}));
  $resultsfile .= ".q" . $self->{qscore} if (defined($self->{qscore}));
  $resultsfile .= ".jackhmmer";

  my $cmd = "$jackhmmer_executable ";
  $cmd .= " --cpu 1 --noali";
  $cmd .= " --tblout $resultsfile";
  $cmd .= " $transcriptsfile $pepfile";
  $cmd .= " 1>/dev/null 2>/dev/null";
  # $cmd .= " 1>/dev/null 2>/dev/null" unless ($self->{debug});
  # print "# $cmd \n" if ($self->debug);

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  print STDERR "# $cmd\n" if ($self->debug);
  unless(system("cd $worker_temp_directory; $cmd") == 0) {
    print("## $cmd\n"); warn("error running jackhmmer, $!\n");
    return 0;
  }
  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

  open RESULTSFILE,"$resultsfile" or die $!;
  while (<RESULTSFILE>) {
    next if ($_ =~ /^#/);

    if ($_ =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/) {
      my ($query_id, $q_accession, $hit_id, $q_acc, $e_value, $score, $bias) = ($1,$2,$3,$4,$5,$6);
      $self->{jackmmer_hits}{$query_id}{score}{$score}     = $hit_id;
      $self->{jackmmer_hits}{$query_id}{e_value}{$e_value} = $hit_id;
      $self->{jackmmer_hits}{$query_id}{bias}{$bias}       = $hit_id;
    } else {
      warn("Error parsing jackhmmer results [$resultsfile]");
    }
  }
  close RESULTSFILE;

  foreach my $query_id (keys %{$self->{jackmmer_hits}}) {
    foreach my $score (sort { $b <=> $a } keys %{$self->{jackmmer_hits}{$query_id}{score}}) {
      # This is the best score for this query_id, store as a tag and go through the next one
      my $tag = "jackhmmer_score";
      $tag   .= ".c" . $self->{coverage} if (defined($self->{coverage}));
      $tag   .= ".q" . $self->{qscore}   if (defined($self->{qscore}));
      my $value = $score;
      $self->_store_tagvalue($query_id,$node_id,$tag,$value);
      last;
    }
  }

  $self->{jackhmmerfile} = $resultsfile;
  return 1;
}

sub run_jackhmmer_merge_by_merge {
  my $self = shift;
  my $node_id = $self->{tree}->node_id;

  my $format = 'fasta';
  my $referencetranscriptsfile = $self->worker_temp_directory . $self->{tree}->node_id . ".sim.fas";
  my $refseqio = Bio::SeqIO->new
      (-file => "$referencetranscriptsfile",
       -format => "$format");
  while (my $refseq = $refseqio->next_seq) {
    $self->{truesim_transcripts}{$refseq->display_id} = $refseq;
  }

  my $aln  = $self->{pachinko_aln};
#   # FIXME debug
#   my $pachinko_output = $self->worker_temp_directory . $self->{tree}->node_id . ".pagan.fas";
#   my $debug_alignio;
#   if (!defined($self->{pachinko_aln})) {
#     $debug_alignio = Bio::AlignIO->new
#       (-file => "$pachinko_output",
#        -format => "$format");
#     $aln = $debug_alignio->next_aln;
#   }
#   # FIXME debug

  my $tag = shift || '\wRR';

  my $merge;
  my $reference;

  print STDERR "Sorting references and non-references\n" if ($self->debug);
  my ($this_prefix,$dummy) = ('reference',undef);

  # Label the reads from the references in the aln -- we only want to
  # merge reads

  foreach my $this_seq ($aln->each_seq) {
    if ($this_seq->display_id =~ /$tag/) {
      ($this_prefix,$dummy) = split('\.',$this_seq->display_id) if ($tag =~ /RR/);
      if ($tag !~ /RR/) {
        my $display_id = $this_seq->display_id; $display_id =~ /($tag.+)\|/;
        $this_prefix = $1;
      }
      $merge->{_prefixes}{$this_prefix}{$this_seq->display_id} = 1;
      $merge->{_seqs}{$this_seq->display_id}{$this_prefix} = $this_seq->seq;
      $merge->{_seqs_tags}{$this_seq->display_id}{$this_prefix} = $this_seq->description;
    } else {
      $merge->{_refseqs}{$this_seq->display_id} = $this_seq->seq;
    }
  }
  return 1 unless (0 < scalar keys %{$merge->{_seqs}});

  my $merge_stats;
  # print STDERR "Merging...\n" if ($self->debug);
  foreach my $this_seq ( $aln->each_seq ) {
    next unless (defined($merge->{_seqs}{$this_seq->display_id}));
    my @temp = keys %{$merge->{_seqs}{$this_seq->display_id}}; my $read_id = $temp[0];
    next if (2 > scalar keys %{$merge->{_prefixes}{$read_id}});
    # print STDERR "." if ($self->debug);
    my ($this_prefix,$dummy) = split('\.',$read_id) if ($tag =~ /RR/);
    if ($tag !~ /RR/) {
      my $display_id = $this_seq->display_id; $display_id =~ /($tag.+)\|/;
      $this_prefix = $1;
    }

    my @read_ids = keys %{$merge->{_prefixes}{$read_id}};
    while (my $read_id1 = shift (@read_ids)) {
      foreach my $read_id2 (@read_ids) {
        my $read1 = $merge->{_seqs}{$read_id1}{$this_prefix};
        my $read2 = $merge->{_seqs}{$read_id2}{$this_prefix};
        my $pre1 = $read1; $pre1 =~ s/\-//g; my $pre2 = $read2; $pre2 =~ s/\-//g;

        my $lcread1 = lc($read1); my $lcread2 = lc($read2);
        $read1 =~ s/-/substr($lcread2, pos($read1), 1)/eg;
        $read2 =~ s/-/substr($lcread1, pos($read2), 1)/eg;

        my $post1 = $read1; $post1 =~ s/\-//g; my $post2 = $read2; $post2 =~ s/\-//g;
        if (
            length($post1)>$merge->{_seqs}{$read_id1}{$this_prefix}
            ||
            length($post2)>$merge->{_seqs}{$read_id2}{$this_prefix}
           ) {
          # $merge_stats->{extension_events}{$this_prefix}++;

          my $pre_blocks  = scalar grep {length $_ > 2} $pre1  =~ /([ATCGN]+|[atcgn]+)/g;
          my $post_blocks = scalar grep {length $_ > 2} $post1 =~ /([ATCGN]+|[atcgn]+)/g;
          my $block_increase = $post_blocks - $pre_blocks;
          if (1 == $block_increase) {
            $merge_stats->{bridging_events}{$this_prefix}{1}++;
          } elsif (1 < $block_increase) {
            $merge_stats->{bridging_events}{$this_prefix}{$block_increase}++;
          }

          # Traverse to see if this pair of simulated sequences are consistently aligned
          my @a1 = split("",$merge->{_seqs}{$read_id1}{$this_prefix});
          my @a2 = split("",$merge->{_seqs}{$read_id2}{$this_prefix});
          my $tags1 = $merge->{_seqs_tags}{$read_id1}{$this_prefix};
          my $tags2 = $merge->{_seqs_tags}{$read_id2}{$this_prefix};
          $tags1 =~ /P1\=(\d+\-\d+) P2\=(\d+\-\d+)/; 
          my ($start1p1,$end1p1) = split("-",$1); my ($start1p2,$end1p2) = split("-",$2);
          $tags2 =~ /P1\=(\d+\-\d+) P2\=(\d+\-\d+)/; 
          my ($start2p1,$end2p1) = split("-",$1); my ($start2p2,$end2p2) = split("-",$2);
          # 1 coordinates always before 2
          if ($start2p1 < $start1p1) {
            my ($temp_start,$temp_end) = ($start1p1,$end1p1); ($start1p1,$end1p1) = ($start2p1,$end2p1);
            ($start2p1,$end2p1) = ($temp_start,$temp_end);
            my ($temp_start,$temp_end) = ($start1p2,$end1p2); ($start1p2,$end1p2) = ($start2p2,$end2p2);
            ($start2p2,$end2p2) = ($temp_start,$temp_end);
            my @temp = @a1; @a1 = @a2; @a2 = @temp;
            my $rtemp = $read_id1; $read_id1 = $read_id2; $read_id2 = $rtemp;
            my $ttemp = $tags1; $tags1 = $tags2; $tags2 = $ttemp;
          }
          my $status1 = 0;
          my $status2 = 0;
          my $is_incompatible = 0;
          foreach my $i (0..(scalar @a1)-1) {
            my $a1 = $a1[$i]; my $a2 = $a2[$i];
            if ($a1 =~ /[ATCGN]/) {
              $status1++;
            }
            if ($a2 =~ /[ATCGN]/) {
              $status2++;
            }
            if ($status2 > $status1) {
#               print "Inconsistent status position [$i]\n";
#               print $merge->{_seqs}{$read_id1}{$this_prefix}, "\n";
#               print $merge->{_seqs}{$read_id2}{$this_prefix}, "\n";
#               print "X" x ($i);
#               print "0" x (scalar @a1 - $i);print "\n";
#               print "$read_id1 $tags1\n"; print "$read_id2 $tags2\n";
              my $root_id1 = $read_id1; $root_id1 =~ s/\/p12//;
              my $root_id2 = $read_id2; $root_id2 =~ s/\/p12//;
              my $sum_qscores1 = $self->{sum_qscores}{$root_id1};
              my $sum_qscores2 = $self->{sum_qscores}{$root_id2};
              my $lowest_qscore_id; my $lowest_qscore;
              if ($sum_qscores1 <= $sum_qscores2) {
                $lowest_qscore_id = $root_id1; $lowest_qscore = $sum_qscores1;
              } else {
                $lowest_qscore_id = $root_id2;  $lowest_qscore = $sum_qscores2;
              }
              my $lowest_avg_degree_coef = $lowest_qscore . ":" . $self->{median_sum_qscores};
              print "incompatible_read\n";
              my $cleanread1 = uc($merge->{_seqs}{$read_id1}{$this_prefix});
              my $cleanread2 = uc($merge->{_seqs}{$read_id2}{$this_prefix});
              my $degree_of_difference = -1; # degree_of_difference($cleanread1,$cleanread2);
              $cleanread1 =~ s/\-//g; $cleanread1 =~ s/N//g;
              $cleanread2 =~ s/\-//g; $cleanread2 =~ s/N//g;
               my $coef_of_difference = sprintf("%.3f",$degree_of_difference/(length($cleanread1)+length($cleanread2))/2);
               $lowest_avg_degree_coef .= ":" . $degree_of_difference . ":" . $coef_of_difference;
              $self->_store_tagvalue($lowest_qscore_id,$node_id,'incompatible_read.q'.$self->{qscore},$lowest_avg_degree_coef);
              $self->{incompatible_read}{$root_id1}++;
              $self->{incompatible_read}{$root_id2}++;
              my ($dummy,$taxon_id) = split("_",$root_id1);
              $self->{ccounts}{$taxon_id}{incompat}++;
              last; # final traverse
            }
          }
          if (!$is_incompatible) {
            print "  compatible_read\n";
            $read1 =~ s/\-//g; $read1 =~ s/N//ig;
            my $double_pair = $read_id1 . ":" . $read_id2;
            next if (defined($self->{compatible_read_pair}{$double_pair}));
            my $this_best_score = $self->jackhmmer_pair($self->{truesim_transcripts}{$this_prefix},$read_id1,$read1);
            $self->{compatible_read_pair}{$double_pair} = 1;
            #            print "Storing:\n", "$double_pair,$node_id,'jackhmmer_pair_score',$this_best_score\n";
            $self->_store_tagvalue($double_pair,$node_id,'jackhmmer_pair_score.q'.$self->{qscore},$this_best_score);
            # The fragment is embedded in the transcript
            if ($self->{truesim_transcripts}{$this_prefix}->seq =~ /$read1/i) {
              $self->_store_tagvalue($read_id1,$node_id,'embedded_fragment.q'.$self->{qscore},$this_best_score);
            }
            my $root_id1 = $read_id1; $root_id1 =~ s/\/p12//;
            my $root_id2 = $read_id2; $root_id2 =~ s/\/p12//;
            my ($dummy,$taxon_id) = split("_",$root_id1);
            $self->{ccounts}{$taxon_id}{compat}++;
          }
        }
      }
    }
  }
  foreach my $read_id (keys %{$self->{incompatible_read}}) {
    my $value = $self->{incompatible_read}{$read_id};
    $self->_store_tagvalue($read_id,$node_id,'num_incompatibilities.q'.$self->{qscore},$value);
  }

  foreach my $taxon_id (keys %{$self->{ccounts}}) {
    my $compat_count   =   $self->{ccounts}{$taxon_id}{compat};
    my $incompat_count = $self->{ccounts}{$taxon_id}{incompat};
    my $coef = 0;
    if (0 < $compat_count && $incompat_count) {
      $coef = $compat_count/$incompat_count;
    }
    $self->_store_tagvalue($taxon_id,$node_id,'num_incompatibilities.q'.$self->{qscore},$coef);
  }

  return 1;
}

sub generate_sam_and_ref_file_placement {
  my $self = shift;
  my $node_id = $self->{tree}->node_id;

  my $samfile = $self->worker_temp_directory . "$node_id.sam";
  my $reffile = $self->worker_temp_directory . "$node_id.ref.fasta";
  my $tleavesfile = $self->worker_temp_directory . "sim.1.$node_id.tleaves.fasta";

  my $outfile = $self->{ancfile};
  throw("no output file") unless (!-z $outfile);

  my $alnio = Bio::AlignIO->new
    (-file => $outfile,
     -format => 'fasta');
  my $aln = $alnio->next_aln;
  my $refio = Bio::SeqIO->new
    (-file => ">$reffile",
     -format => 'fasta');
  my $tleavesio = Bio::SeqIO->new
    (-file => ">$tleavesfile",
     -format => 'fasta');

  foreach my $id (keys %{$self->{subset_aa_ancestral_seqs}}) {
    my $seq = $self->{subset_aa_ancestral_seqs}{$id};
    my $outseq = Bio::Seq->new
      (-display_id => $id,
       -seq => $seq);
    $self->{id_subset_cds_ancestal_seq}{$self->{cds_ancestral_seqs}{$id}}{$id} = 1;
    $self->{exonerate_seqs}{$outseq->display_id} = $outseq->seq;
    $tleavesio->write_seq($outseq);
  }
  $tleavesio->close;

  foreach my $seq ($aln->each_seq) {
    if ($seq->display_id =~ /\|/) {
      # read
    } else {
      # reference
      my $ungapped_seq = $seq->seq; $ungapped_seq =~ s/\-//g;
      my $outseq = Bio::Seq->new
        (-display_id => $seq->display_id,
         -description => $seq->description,
         -seq => $ungapped_seq,
         -alphabet => 'dna');
      $refio->write_seq($outseq);
      $self->{ref_ids}{$seq->display_id} = 1;
    }
  }

  my @sam_entries;
#   foreach my $ref_id (keys %{$self->{ref_ids}}) {
  foreach my $ref_id (keys %{$self->{read_placement_sets}}) {

    my $this_ref_aln = Bio::SimpleAlign->new();
    my $ref_seq = $aln->get_seq_by_id($ref_id);
    $this_ref_aln->add_seq($ref_seq);
    if (!defined($ref_seq)) {
      $DB::single=1;1;
    }
    foreach my $stable_id (keys %{$self->{read_placement_sets}{$ref_id}}) {
      my $read_seq = $aln->get_seq_by_id($stable_id);
      if (!defined $read_seq) {
        print STDERR "$stable_id not found\n" if ($self->debug);
        next;
        $DB::single=1;1;
      }
      $this_ref_aln->add_seq($read_seq);
    }
    $this_ref_aln = $this_ref_aln->remove_gaps(undef,1);

    my $this_ref_seq = $this_ref_aln->get_seq_by_id($ref_id);
    foreach my $seq ($this_ref_aln->each_seq) {
      next if ($seq->display_id !~ /\|/);
      my $this_pair = Bio::SimpleAlign->new();
      $this_pair->add_seq($this_ref_seq);
      $this_pair->add_seq($seq);
      $this_pair = $this_pair->remove_gaps(undef,1);
      $self->{this_pair} = $this_pair;
      my ($seq1,$seq2) = $this_pair->each_seq;
      my $cigar_line = $self->to_sam_cigar_line($seq1->seq,$seq2->seq);
      my $sam_entry = $self->sam_entry($seq1,$seq2,$cigar_line);
      push @sam_entries, "$sam_entry\n";
    }
  }

  open SAM,">$samfile" or die $!;
  print SAM join("", sort @sam_entries);
  close SAM;

  $self->{samfile} = $samfile;
  $self->{reffile} = $reffile;


  return;
}

sub jackhmmer_pair {
  my $self = shift;
  my $reftruetranscript = shift;
  my $fragment_id = shift;
  my $fragment_seq = shift;
  my $node_id = $self->{tree}->node_id;

  my $jackhmmer_executable = $self->analysis->program_file;
  $jackhmmer_executable = "/software/ensembl/compara/hmmer3/hmmer-3.0/src/jackhmmer";
#    unless (-e $jackhmmer_executable);

  my $worker_temp_directory = $self->worker_temp_directory;
  my $resultsfile = $worker_temp_directory . $node_id;
  $resultsfile .= ".c" . $self->{coverage} if (defined($self->{coverage}));
  $resultsfile .= ".q" . $self->{qscore} if (defined($self->{qscore}));
  $resultsfile .= ".jackhmmer";

  my $truecdnafile = $self->worker_temp_directory . $node_id . ".truecdna.fas";
  open TRUECDNAFILE, ">$truecdnafile" or die $!;
  print TRUECDNAFILE ">", $reftruetranscript->display_id, "\n", $reftruetranscript->seq, "\n";;
  close TRUECDNAFILE;
  my $truepepfile = $self->worker_temp_directory . $node_id . ".truepep.fas";
  open TRUEPEPFILE, ">$truepepfile" or die $!;
  print TRUEPEPFILE ">", $reftruetranscript->display_id, "\n", $reftruetranscript->translate->seq, "\n";;
  close TRUEPEPFILE;

  my $fragmfile = $self->worker_temp_directory . $node_id . ".fragm.fas";
  open FRAGMFILE, ">$fragmfile" or die $!;
  print FRAGMFILE ">", $fragment_id, "\n", $fragment_seq, "\n";
  # print ">", $fragment_id, "\n", $fragment_seq, "\n" if ($self->debug); #FIXME
  close FRAGMFILE;
  $self->{seqs}{$fragment_id} = $fragment_seq;

  my $fragmtranscriptsfile = $self->run_magicmirror($truecdnafile,$fragmfile);
  next unless (!(-z $fragmtranscriptsfile));

  my $cmd = "$jackhmmer_executable ";
  $cmd .= " --cpu 1 --noali";
  $cmd .= " --tblout $resultsfile";
  $cmd .= " $fragmfile $truecdnafile";
  $cmd .= " 1>/dev/null 2>/dev/null";
  # print "# $cmd \n" if ($self->debug);

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  unless(system("cd $worker_temp_directory; $cmd") == 0) {
    print("## $cmd\n"); throw("error running jackhmmer, $!\n");
  }
  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

  open RESULTSFILE,"$resultsfile" or die $!;
  my $best_score = 0;
  while (<RESULTSFILE>) {
    next if ($_ =~ /^#/);

    if ($_ =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/) {
      my ($query_id, $q_accession, $hit_id, $q_acc, $e_value, $score, $bias) = ($1,$2,$3,$4,$5,$6);
      $self->{jackmmer_hits}{$query_id}{score}{$score}     = $hit_id;
      $self->{jackmmer_hits}{$query_id}{e_value}{$e_value} = $hit_id;
      $self->{jackmmer_hits}{$query_id}{bias}{$bias}       = $hit_id;
      $best_score = $score if ($best_score < $score);
    } else {
      warn("Error parsing jackhmmer results [$resultsfile]");
    }
  }
#   foreach my $query_id (keys %{$self->{jackmmer_hits}}) {
#     foreach my $score (sort { $b <=> $a } keys %{$self->{jackmmer_hits}{$query_id}{score}}) {
#       # This is the best score for this query_id, store as a tag and go through the next one
#       my $tag = "jackhmmer_score";
#       $tag   .= ".c" . $self->{coverage} if (defined($self->{coverage}));
#       $tag   .= ".q" . $self->{qscore}   if (defined($self->{qscore}));
#       my $value = $score;
#       $self->_store_tagvalue($query_id,$node_id,$tag,$value);
#       last;
#     }
#   }

  return $best_score;
}

sub degree_of_difference {
  my ($string_x, $string_y) = @_;

  s{\s}''xmsg for $string_x, $string_y;

  # the longest string always comes first:
  if (length($string_x) < length($string_y)) {
    my $temp  = $string_x;
    $string_x = $string_y;
    $string_y = $temp;
  }

  my @chain_x = split m{}xms, $string_x;
  my @chain_y = split m{}xms, $string_y;

  my @sd = sdiff(\@chain_x, \@chain_y);

  my $inserts   = () = grep {$_->[0] eq '+'} @sd;
  my $deletes   = () = grep {$_->[0] eq '-'} @sd;
  my $changes   = () = grep {$_->[0] eq 'c'} @sd;
  my $unchanged = () = grep {$_->[0] eq 'u'} @sd;

  $inserts + $changes;
}


sub call_parse_alignment {
  my $self = shift;

  # We inherit from Pachinkoalign
  $self->parse_alignment;

  return 1;
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

sub run_treebest {
  my $self = shift;

  my $node_id = $self->{tree}->node_id;
  my $tree = $self->{tree};

  print STDERR "Tree pruning...\n" if ($self->debug);

  my $alnfile = $self->{outfile};

  my $treefile = $self->worker_temp_directory . "$node_id.nh";
  foreach my $subnode ($self->{mappedtree}->get_all_subnodes) {
    $subnode->distance_to_parent(0.0001) if (0 == $subnode->distance_to_parent);
  }
  my $newick_string = $self->{mappedtree}->newick_format;
  $newick_string =~ s/\#\d+\#//g; # get rid of the internal node ids

  open NEWICK, ">$treefile" or die "$!";
  print NEWICK "$newick_string\n";
  close NEWICK;

  my $treebestmagicmirror_executable = $self->analysis->program_file;
  $treebestmagicmirror_executable = "/software/ensembl/compara/prank/latest/src/prank";
#    unless (-e $pachinko_executable);

  my $outprefix = $self->worker_temp_directory . "$node_id";

  # /software/ensembl/compara/prank/latest/src/prank -d=<input_alignment> -t=<tree> -shortnames -e -notree -nopost -f=11 -writeanc -o=<output_name>
  my $cmd = "$treebestmagicmirror_executable ";
  $cmd .= " -d=$alnfile ";
  $cmd .= " -t=$treefile ";
  $cmd .= " -shortnames -e -notree -nopost -f=11 -writeanc ";
  $cmd .= " -o=$outprefix";
  $cmd .= " 2>&1 > /dev/null" unless($self->debug && defined($self->{verbose}));

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(1);
  my $starttime = time();
  print STDERR "treebestmagicmirror...\n" if ($self->debug);
  my $worker_temp_directory = $self->worker_temp_directory;
  # $DB::single=1;1;
  unless(system("cd $worker_temp_directory; $cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running treebestmagicmirror), $!\n");
  }

  $self->{'comparaDBA'}->dbc->disconnect_when_inactive(0);

  $self->{ancfile} = "$worker_temp_directory"."$node_id.0.ancseq";
  throw("Output file not found") unless (-e $self->{outfile});

  my $runtime = time()*1000-$starttime;
  printf("  %1.3f secs\n", (time()-$starttime)) if ($self->debug);
  $self->{'tree'}->store_tag('treebestmagicmirror_runtime_msec', $runtime);

  return 1;
}

sub parse_treebest {
  my $self = shift;

  my ($Ptreeref, $PRKnodref);
  my $node_id = $self->{tree}->node_id;
  my $tree = $self->{tree};
  my $ancfile = $self->{ancfile};

  my @names = ();
  my ($name, $seq) = ('','');
  my $length = -1;
  my $different = 0;

  open(FILE, $ancfile) or die "error opening the treebestmagicmirror file: $ancfile\n";
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
        $self->{seqs}{$name} = $seq;
        if ($length < 0) { $length = length($seq); }
        if ($length != length($seq)) {
          warn "$name is of different length: $length and ".length($seq).".\n";
          $different++;
        }
      }
      $name = $_;
      $seq = "";
      next;
    }
    $seq .= $_;
  }

  $self->{rectree} = $Ptreeref;

  return 1;
}


sub sam_entry {
  my $self = shift;
  my $ref = shift;
  my $query = shift;
  my $cigar_line = shift;

  my $query_ungapped_seq = $query->seq; $query_ungapped_seq =~ s/\-//g;

  my $qname = $query->display_id;
  my $flag  = 0;
  my $rname = $ref->display_id;
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

sub columbus_merge_fragments_old {
  my $self = shift;
  my $node_id = $self->{tree}->node_id;

  my $samfile = $self->{samfile};
  my $reffile = $self->{reffile};

  my $velvet_dir = $self->worker_temp_directory . "columbus";

  my $cmd = "/nfs/users/nfs_a/avilella/src/velvet/latest/velvet/velveth $velvet_dir 21 -reference $reffile -short -sam $samfile";
  print("$cmd\n");
  unless(system("$cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running velveth), $!\n");
  }

  $cmd = "/nfs/users/nfs_a/avilella/src/velvet/latest/velvet/velvetg $velvet_dir -exp_cov 1";
  print("$cmd\n");
  unless(system("$cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running velvetg), $!\n");
  }

  my $velvetfile = $velvet_dir . "/" . "contigs.fa";
  throw("no contigs produced") if (-z $velvetfile);

  $self->{velvetfile} = $velvetfile;

  return;
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


# sub label_taxa_genetree {
#   my $self = shift;

#   my $node_id = $self->{tree}->node_id;

#   my $sql1 = "select value from protein_tree_tag where tag='species_tree_string'";
#   my $sth1 = $self->dbc->prepare($sql1);
#   $sth1->execute;
#   my $species_tree_string = $sth1->fetchrow_hashref;
#   $sth1->finish;
#   my $eval_species_tree;
#   eval {
#     $eval_species_tree = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($species_tree_string->{value});
#     my @leaves = @{$eval_species_tree->get_all_leaves};
#   };
#   if($@) {
#     unless(-e $self->{'species_tree_file'}) {
#       throw("can't find species_tree\n");
#     }
#   } else {
#     $self->{species_tree_string} = $species_tree_string->{value};
#     my $spfilename = $self->worker_temp_directory . "spec_tax.nh";
#     open SPECIESTREE, ">$spfilename" or die "$!";
#     print SPECIESTREE $self->{species_tree_string};
#     close SPECIESTREE;
#     $self->{'species_tree_file'} = $spfilename;
#   }
#   my $taxon_tree = $self->load_species_tree_from_file;
#   $self->{'taxon_tree'} = $taxon_tree;
#   $DB::single=1;1;
#   my @leaves = @{$self->{tree}->get_all_leaves};
#   if (scalar @leaves < 2) {return 0;}

#   # Foreach leaf in the extended species tree
#   my $flagged_global;
#   my $num_new_taxa = 0;
#   foreach my $taxon_subnode (@{$self->{taxon_tree}->get_all_leaves}) {
#     print STDERR "Sampling " . $taxon_subnode->name,"\n" if ($self->debug);

#     # If it's a new taxon
#     next if (defined ($self->{species_set_hash}{$taxon_subnode->node_id}));
#     $num_new_taxa++;
#     # Traverse from the leaf to the root of the species tree, then try
#     # to find the taxon level in the reference genetree
#     my $taxon_id = $taxon_subnode->node_id;
#     my $parent = $taxon_subnode;
#     while ($parent->has_parent) {
#       my @ref_subnodes = @{$self->{tree}->get_all_nodes_by_tag_value('taxon_id',$parent->node_id)};
#       $parent = $parent->parent;

#       # We found a node in the reference genetree that is the smallest
#       # taxon level in the species tree for the new taxon
#       next unless (defined(@ref_subnodes) && 0 < scalar @ref_subnodes);
#       my $flagged_nested;
#       # Check if the reference subnodes are overlapping
#       my $num_nonflagged;
#       if (1 < scalar @ref_subnodes) {
#         my $i = 0;
#         # If more than one subnode, check if they are nested, and
#         # only use the smallest
#         # If not nested, use all
#         foreach my $j ($i+1 .. scalar(@ref_subnodes)-1 ) {
#           my $a = $ref_subnodes[$i];
#           my $b = $ref_subnodes[$j];
#           my $ancestor = $a->find_first_shared_ancestor($b);
#           if ($ancestor->node_id eq $a->node_id || $ancestor->node_id eq $b->node_id) {
#             # They are nested, flag out the biggest subtree
#             $flagged_nested->{$a->node_id} = 1 if ($ancestor->node_id eq $a->node_id);
#             $flagged_nested->{$b->node_id} = 1 if ($ancestor->node_id eq $b->node_id);
#           }
#         }
#         $num_nonflagged = scalar @ref_subnodes - scalar keys %{$flagged_nested};
#       }
#       # Then go through each of them and attach a new child
#       foreach my $ref_subnode (@ref_subnodes) {
#         if ($ref_subnode->node_id eq $self->{tree}->node_id) {
#           next;
#         }
#         if (defined($flagged_nested->{$ref_subnode->node_id})) {
#           next;
#         }
#         if (defined($flagged_global->{$ref_subnode->node_id})) {
#           next;
#         }
#         # We are in, so let's tag them with TID first
#         # my $tid = $ref_subnode->get_tagvalue('taxon_id');
#         my $tid = $taxon_id; # tid is from the new species onto this node
#         $ref_subnode->store_tag('TID',$tid) if ('' ne $tid);
#         $self->{all_tids}{$tid} = 1;
#       }
#     }
#   }
#   return 1 if ($num_new_taxa < 1);

#   my $reloaded_ref_tree = $self->{treeDBA}->fetch_node_by_node_id($self->{tree}->node_id);
#   my $nhx_treefile = $self->worker_temp_directory . "$node_id.nhx";
#   my $nhx_string = $reloaded_ref_tree->nhx_format('member_id_taxon_id');
#   open REFNHX, ">$nhx_treefile" or die "$!";
#   print REFNHX "$nhx_string\n";
#   close REFNHX;
#   $self->{nhx_treefile} = $nhx_treefile;

#   return 1;
# }

1;

