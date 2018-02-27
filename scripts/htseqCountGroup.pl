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
my @files;
my $group;
my $skip = 0;
my $output = './groupResult.txt';

GetOptions (
  "f|file=s{1,}"                     =>\@files,
  "o|output=s{1,1}"                  =>\$output,
  "h|help"                           =>\$help,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;


my %hashGene;

for (my $i = 0; $i < scalar(@files); $i++) {
  my $file = $files[$i];
  open (IN, "<$file") or die "cannot open file:$file, $!\n";
  while (my $line=<IN>) {
    next if ($line =~ /^__/);
    chomp($line);
    my ($gene, $count) = split ("\t", $line);
    next if ($gene eq "__no_feature");
    next if ($gene eq "__ambiguous");
    next if ($gene eq "__too_low_aQual");
    next if ($gene eq "__not_aligned");
    next if ($gene eq "__alignment_not_unique");
    $hashGene{$gene}->{$i} = $count;
  }
}

open (OUT, ">$output") or die "cannot open file($output), $!\n";
my @genes = sort keys (%hashGene);
foreach my $gene (@genes) {
  print OUT $gene;
  for (my $i = 0; $i < scalar(@files); $i++) {
    my $count = $hashGene{$gene}->{$i};
    print OUT "\t", $count;
  }
  print OUT "\n";
}

################# Abbreviations for this script #################
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for grouping the results of htseq-count (version 0.6.0).

=head1 CAUTIONS

=head1 SYNOPSIS


htseqCountGroup.pl [options] -f [file1 file2] -o [file]

 Options:
    -f | --file        The input files (htseq-count1 htseq-count2...)
    -o | --output      The output file
    -h | --help        Brief help message
    -man               Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

htseqCountGroup.pl [options] -f [file1 file2] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will group the results of htseq-count (version 0.6.0).

=cut
