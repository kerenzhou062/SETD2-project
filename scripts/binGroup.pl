#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for grouping the bin files.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $file;
my $group;
my $skip = 0;
my $output = './groupResult.txt';

GetOptions (
  "f|file=s{1,1}"                    =>\$file,
  "g|group=s{1,1}"                   =>\$group,
  "o|output=s{1,1}"                  =>\$output,
  "skip=s"                           =>\$skip,
  "h|help"                           =>\$help,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

my @files = split (",", $file);
my @groups = split (",", $group);
if (scalar(@files) != scalar(@groups)) {
  print STDERR "The groups do not match the number of files!\n";
  exit;
}
## check the groups
my %hashGroup;
foreach my $eachGroup(@groups) {
  if (exists($hashGroup{$eachGroup})) {
    print STDERR "Reduplicative groups are not allowed!\n";
    exit;
  }else{
    $hashGroup{$eachGroup}++;
  }
}

if ($skip >= 0) {
  print STDERR "The script will skip the fisrt $skip lines from each files.\n";
}else{
  print STDERR "The -skip do not allow negative number!\n";
  exit;
}

open (OUT, ">$output") or die "cannot open file($output), $!\n";
print OUT "region\tbin\tpercentage\tgroup\n";
for (my $i = 0; $i < scalar(@groups); $i++) {
  my $eachFile = $files[$i];
  my $eachGroup = $groups[$i];
  open (IN, "<$eachFile") or die "cannot open file($eachFile), $!\n";
  my $count = 0;
  while (my $line=<IN>) {
    $count++;
    next if ($count <= $skip);
    chomp($line);
    print OUT join ("\t", ($line, $eachGroup)), "\n";
  }
}
close(OUT);

################# Abbreviations for this script #################
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for grouping the bin files.

=head1 CAUTIONS

If "-skip" are not defined, it will be set to "0".

=head1 SYNOPSIS


binGroup.pl [options] -f [file,file] -g [name1,name2] -o [file]

 Options:
    -f | --file        The input files (comma separated)
    -g | --group       The group names (comma separated)
    -o | --output      The output file
    -skip              Skip the first n lines from files;
    -h | --help        Brief help message
    -man               Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

binGroup.pl [options] -f [file,file] -g [name1,name2] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will group the bin files.

=cut
