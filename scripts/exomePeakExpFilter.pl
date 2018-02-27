#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
###develop by K.R.Chow, designed for formatting exomePeak (bed12) to summit (bed6).

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $output;
my $expFile;
my $cutoff = 1;
my $peakXlsFile;
my $sample;
my $verbose;

GetOptions (
  "o|output=s{1,1}"                       =>\$output,
  "expF=s{1,1}"                           =>\$expFile,
  "peakF=s{1,1}"                          =>\$peakXlsFile,
  "cutoff=s{1,1}"                         =>\$cutoff,
  "sample=s{1,1}"                         =>\$sample,
  "verbose"                               =>\$verbose,
  "h|help"                                =>\$help,
  "man"                                   =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

my %hashGeneExp;
open (EXPFILE, "<$expFile") or die "cannot open file:$!";
<EXPFILE>;#GeneID geneName  rep1  rep2  rep3  MeanFPKM  log2MeanFPKM
while (my $line=<EXPFILE>) {
  chomp($line);
  my @lineContents = split("\t", $line);
  my $geneID = $lineContents[0];
  if ($sample eq 'rep1'){
    if ($lineContents[2] >= $cutoff){
      $hashGeneExp{$geneID} = $lineContents[2];
    }
  }elsif ($sample eq 'rep2'){
    if ($lineContents[3] >= $cutoff){
      $hashGeneExp{$geneID} = $lineContents[3];
    }
  }elsif ($sample eq 'rep3'){
    if ($lineContents[4] >= $cutoff){
      $hashGeneExp{$geneID} = $lineContents[4];
    }
  }elsif ($sample eq 'mean'){
    if ($lineContents[5] >= $cutoff){
      $hashGeneExp{$geneID} = $lineContents[5];
    }
  }
}
close EXPFILE;

open (PEAK, "<$peakXlsFile") or die "cannot open file:$!";
<PEAK>;
open (OUT, ">$output") or die "cannot open file:$!";
while (my $line=<PEAK>) {
  my $geneID = (split("\t", $line))[3];
  if (exists($hashGeneExp{$geneID})) {
    print OUT $line;
  }
}

close PEAK;
close OUT;

sub verbose{
  my $warningText = shift;
  if (defined($verbose)) {
    print STDERR $warningText, "\n";
  }
}


################# Abbreviations for this script #################
#
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for filtering the peak with expression file

=head1 CAUTIONS

Defaults: -cutoff=(1)

=head1 SYNOPSIS

exomePeakExpFilter.pl [options] -x [file1 file2] -o [file]

 Options:
    -o | --output             The output file name
    -expF                     The expression file
    -peakF                    The exomePeak xls file (result from exomePeakToSummit.pl -original)
    -cutoff                   The cutoff value for expressed gene
    -sample                   The sample used for estimate gene expression (rep1, rep2, rep3 and mean)
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

exomePeakExpFilter.pl [options] -x [file1 file2] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will filter the peak with expression file

=cut
