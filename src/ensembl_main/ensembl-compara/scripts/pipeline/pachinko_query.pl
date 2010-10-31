#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Hive::URLFactory;

use File::Basename;
use Time::HiRes qw(time gettimeofday tv_interval);

use POSIX; # stemplot_transcripts

# ok this is a hack, but I'm going to pretend I've got an object here
# by creating a blessed hash ref and passing it around like an object
# this is to avoid using global variables in functions, and to consolidate
# the globals into a nice '$self' package
my $self = bless {};

GetOptions('help'                                     => \$self->{help},
           'url=s'                                    => \$self->{url},
           'xyz=s'                                    => \$self->{xyz},
           'n50=s'                                    => \$self->{n50},
           'transcripts|dump_transcripts=s'           => \$self->{dump_transcripts},
           'stemplot_transcripts=s'                   => \$self->{stemplot_transcripts},
           'family_stats=s'                      => \$self->{family_stats},
           'set_numreads=s'                      => \$self->{set_numreads},
);

if ($self->{help}) {
  usage();
}

sub usage {
  print 'pachinko_query.pl -url mysql://username:password@host:port/dbname --method=parameter\n';
  print "  -help                  : print this help\n";
  print "  -url <url>             : connect to compara at url\n";
  print "\n";
  exit(1);
}

usage() unless (defined($self->{url}));
$self->{'comparaDBA'}  = Bio::EnsEMBL::Hive::URLFactory->fetch($self->{url}, 'compara');

$self->xyz if ($self->{xyz});
$self->n50 if ($self->{n50});
$self->dump_transcripts if ($self->{dump_transcripts});
$self->stemplot_transcripts if ($self->{stemplot_transcripts});
$self->set_numreads if ($self->{set_numreads});
$self->family_stats if ($self->{family_stats});

sub xyz {
  my $self = shift;
  my $xyz = $self->{xyz};

  exit(0);
}

sub n50 {
  my $self = shift;
  my $n50 = $self->{n50};

  open FILE,$n50 or die $!;
  my ($len,$total,$ncontigs,$p100,$maxlen)=(0,0,0,0,0);
  my @x;
  while (<FILE>) {
    if (/^[\>\@]/) {
      if ($len>0) {
        $total+=$len;
        push @x,$len;
        $p100++ if ($len > 100);
        $ncontigs++;
        $maxlen = $len if ($len > $maxlen);
      }
      $len=0;
    } else {
      s/\s//g;
      $len+=length($_);
    }
  }
  if ($len>0) {
    $total+=$len;
    push @x,$len;
  }
  @x=sort{$b<=>$a} @x; 
  my ($count,$half)=(0,0);
  for (my $j=0;$j<@x;$j++) {
    $count+=$x[$j];
    if (($count>=$total/2)&&($half==0)) {
      print "N50,$x[$j]\n";
      $half=$x[$j]
    } elsif ($count>=$total*0.9) {
      print "N90,$x[$j]\n";
      last;
    }
  }
  print "Total_length,$total\n";
  print "No_of_contigs,$ncontigs\n";
  print "Contigs>100,$p100\n";
  print "Max_length,$maxlen\n";
  close FILE;

  exit(0);
}

sub dump_transcripts {
  my $self = shift;
  my $dump_transcripts = $self->{dump_transcripts};

  $self->{memberDBA} = $self->{comparaDBA}->get_MemberAdaptor;
  my $source_name = 'EXTERNALCDS';
  $source_name = 'ENSEMBLTRANS' if ($dump_transcripts =~ /ENSEMBLTRANS/);
  $source_name = 'ENSEMBLPEP'   if ($dump_transcripts =~ /ENSEMBLPEP/);
  $source_name = 'ENSEMBLGENE'  if ($dump_transcripts =~ /ENSEMBLGENE/);
  my $taxon_id = 0; my $genome_db_id = 999999; my @transcripts;
  if ($dump_transcripts =~ /taxon=(\d+)/) { $taxon_id = $1; }
  @transcripts = @{$self->{memberDBA}->fetch_all_by_source_genome_db_id($source_name,$genome_db_id)} unless (0 < $taxon_id);
  @transcripts = @{$self->{memberDBA}->fetch_all_by_source_taxon($source_name,$taxon_id)} if (0 < $taxon_id);
  foreach my $transcript (@transcripts) {
    my $sequence;
    $sequence = $transcript->sequence;
    eval { $sequence = $transcript->sequence_cds if ($dump_transcripts =~ /sequence_cds/);};
    next unless (defined $sequence && length($sequence)>3);
    my $transcript_id = $transcript->stable_id;
    my $description = $transcript->description;
    my $taxon_id = $transcript->taxon_id;
    $sequence =~ s/N//g unless ($dump_transcripts =~ /full/ || $dump_transcripts =~ /ENSEMBLPEP/);
    print ">$transcript_id DESCRIPTION=$description TAXON=$taxon_id\n$sequence\n";
  }
  exit(0);
}

sub stemplot_transcripts {
  my $self = shift;
  my $taxon_id;
  if (1 < $self->{stemplot_transcripts}) {
    $taxon_id = $self->{stemplot_transcripts};
  } else {
    $taxon_id = "select taxon_id from genome_db where length(genebuild)<2";
  }
  my @values;
  my $sql1 = "SELECT s.length from member m, sequence s where m.sequence_id=s.sequence_id and m.genome_db_id=999999 and m.taxon_id in ($taxon_id)";
  my $sth1 = $self->{comparaDBA}->dbc->prepare($sql1);
  $sth1->execute;
  my $result1;
  while ($result1 = $sth1->fetchrow) {
    push @values, $result1;
  }
  $sth1->finish;

  stemplot(\@values);
  exit(0);
}

sub set_numreads {
  my $self = shift;
  exit(0);
}
sub family_stats {
  my $self = shift;
  my $tag = $self->{family_stats};
  my $header = "node_id, tag, value";
  print STDOUT $header,"\n";
  my $sql1 = "SELECT $header from protein_tree_tag where tag=\'$tag\'";
  my $sth1 = $self->{comparaDBA}->dbc->prepare($sql1);
  $sth1->execute;
  my $result1;
  while ($result1 = $sth1->fetchrow_arrayref) {
    print STDOUT join(", ",@$result1) ,"\n";
  }
  $sth1->finish;
  exit(0);
}

sub stemplot {

  my $p = shift;
  my @sortedData = sort { $a <=> $b } @$p;

  my $div = 10 ** (floor log10 $sortedData[ int ($#sortedData / 2)]);
  my $range = $sortedData[-1] - $sortedData[0];

  if (($range / $div) < 5.0) {
    $div /= 10;
  }
  elsif (($range / $div) > 35.0) {
    $div *= 10;
  }

  my @dataTo1DecimalPlace = map { sprintf "%.1f", ($_ / $div) } @sortedData;

  my %leaves;

  foreach (@dataTo1DecimalPlace) {
    my $stem = floor $_;
    $leaves{$stem} .= sprintf "%1.0f", 10 * ($_ - $stem);
  }

  my $minStem = floor $dataTo1DecimalPlace[0];
  my $maxStem = floor $dataTo1DecimalPlace[-1];
  print "key: 10 stem = ", 10 * $div, " data\n";

  foreach ($minStem .. $maxStem) {
    printf "%3d: %s\n", $_, defined $leaves{$_} ? $leaves{$_} : "";
  }

  print "\n";

}



1;
