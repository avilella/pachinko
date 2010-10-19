#!/usr/local/ensembl/bin/perl -w

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Hive::URLFactory;

use File::Basename;
use Time::HiRes qw(time gettimeofday tv_interval);

# ok this is a hack, but I'm going to pretend I've got an object here
# by creating a blessed hash ref and passing it around like an object
# this is to avoid using global variables in functions, and to consolidate
# the globals into a nice '$self' package
my $self = bless {};

GetOptions('help'                                     => \$self->{help},
           'url=s'                                    => \$self->{url},
           'xxx=s'                      => \$self->{xxx},
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

$self->xxx if ($self->{xxx});
$self->set_numreads if ($self->{set_numreads});
$self->family_stats if ($self->{family_stats});

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
$DB::single=1;1;

1;
