#!/usr/local/bin/perl
use Getopt::Long;

my ($databasefile,$readsfile,$onlyindex,$tmpdir,$tag,$debug);

my $bwa_executable      = "/nfs/users/nfs_a/avilella/src/bwa/latest/bwa-0.5.9rc1/bwa";
my $samtools_executable = "/nfs/users/nfs_a/avilella/src/samtools/latest/samtools/samtools";
my $indexmethod = '';
my $tmpdir = '/tmp';

GetOptions(
	   'db|input|databasefile:s' => \$databasefile,
           'tmpdir:s'  => \$tmpdir,
	   'r|reads|readsfile:s' => \$readsfile,
	   't|tag:s' => \$tag,
           'd|debug:s' => \$debug,
           'onlyindex:s' => \$onlyindex,
           'indexmethod:s' => \$indexmethod,
	   'bwa_exe:s' => \$bwa_executable,
	   'samtools_exe:s' => \$samtools_executable,
          );

my $self = bless {};
my $tmp_dir = $self->worker_process_temp_directory;

my $cmd;
# bwa index
my $indexmethod = "-a $indexmethod" unless ('' eq $indexmethod);
$cmd = "$bwa_executable index $indexmethod $databasefile";
print STDERR "# $cmd\n" if ($debug);
if (1 != $onlyindex && ((-e "$databasefile.bwt" && !-z "$databasefile.bwt") || (-e "$databasefile.pac" && !-z "$databasefile.pac"))) {
  print STDERR "$databasefile already indexed\n";
} else {
  unless(system("$cmd") == 0) {    print("$cmd\n");    throw("error running bwa index: $!\n");  }
}

exit 0 if ($onlyindex);

# bwa aln
$cmd = "$bwa_executable aln $databasefile $readsfile > $tmp_dir/$$.sai";
print STDERR "$cmd\n" if ($debug);
unless(system("$cmd") == 0) {    print("$cmd\n");    throw("error running bwa aln: $!\n");  }

# bwa samse
$cmd = "$bwa_executable samse $databasefile $tmp_dir/$$.sai $readsfile > $tmp_dir/$$.sam";
print STDERR "$cmd\n" if ($debug);
unless(system("$cmd") == 0) {    print("$cmd\n");    throw("error running bwa samse: $!\n");  }

# samtools view
$cmd = "$samtools_executable view -S -b -o $tmp_dir/$$.bam $tmp_dir/$$.sam";
print STDERR "$cmd\n" if ($debug);
unless(system("$cmd") == 0) {    print("$cmd\n");    throw("error running samtools view: $!\n");  }

# samtools sort
$cmd = "$samtools_executable sort $tmp_dir/$$.bam $tag";
print STDERR "$cmd\n" if ($debug);
unless(system("$cmd") == 0) {    print("$cmd\n");    throw("error running samtools sort: $!\n");  }

# samtools index
$cmd = "$samtools_executable index $tag.bam";
print STDERR "$cmd\n" if ($debug);
unless(system("$cmd") == 0) {    print("$cmd\n");    throw("error running samtools index: $!\n");  }

print "$tag.bam\n";
print "$tag.bam.bai\n";

$self->cleanup_worker_process_temp_directory;

########################################
#### METHODS
########################################

sub DESTROY {
  my $self = shift;
  $self->cleanup_worker_process_temp_directory;
}

sub worker_process_temp_directory {
  my $self = shift;
  
  unless(defined($self->{'_tmp_dir'}) and (-e $self->{'_tmp_dir'})) {
    #create temp directory to hold fasta databases
    $self->{'_tmp_dir'} = $tmpdir . "/worker.$$/";
    mkdir($self->{'_tmp_dir'}, 0777);
    throw("unable to create ".$self->{'_tmp_dir'}) unless(-e $self->{'_tmp_dir'});
  }
  return $self->{'_tmp_dir'};
}


sub cleanup_worker_process_temp_directory {
  my $self = shift;
  if($self->{'_tmp_dir'}) {
    my $cmd = "rm -r ". $self->{'_tmp_dir'};
    system($cmd) if (-e $self->{'_tmp_dir'});
  }
}
