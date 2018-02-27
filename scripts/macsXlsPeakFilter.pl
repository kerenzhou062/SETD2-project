#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
###develop by K.R.Chow, designed for filtering and formatting the macs xls result to bed6

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $xlsFile;
my $output;
my $foldEnrich = 0;
my $name4thCol   = 'peak';
my $fdr        = 0.05;
my $pvalue     = 0.05;
my $sort       = "V";
my $macs2Flag;
my $original;
my $verbose;

GetOptions (
  "f|fold=s{1,1}"                         =>\$foldEnrich,
  "o|output=s{1,1}"                       =>\$output,
  "x|xls=s{1,1}"                          =>\$xlsFile,
  "fdr=s{1,1}"                            =>\$fdr,
  "name=s{1,1}"                           =>\$name4thCol,
  "pval=s{1,1}"                           =>\$pvalue,
  "sort=s{1,1}"                           =>\$sort,
  "original"                              =>\$original,
  "macs2"                                 =>\$macs2Flag,
  "verbose"                               =>\$verbose,
  "h|help"                                =>\$help,
  "man"                                   =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

if ($pvalue <= 0) {
  &verbose("Invalid -pval: $pvalue! Will set -pval to 0.05");
  $pvalue = 0.05;
}

if ($fdr <= 0) {
  &verbose("Invalid -pval: $fdr! Will set -fdr to 0.05");
  $fdr = 0.05;
}

if (defined($macs2Flag)) {
  $pvalue = -1 * log($pvalue) / log(10);
  $fdr = -1 * log($fdr);
}else{
  $pvalue = -10 * log($pvalue) / log(10);
  $fdr = $fdr * 100;
}

open (IN, "<$xlsFile") or die "Cannot open file: $xlsFile, $!\n";
open (OUT, ">$output") or die "Cannot open file: $output, $!\n";
my $count = 1;
while (my $line=<IN>) {
  next if ($line =~ /^("?#.+|\s+)/);
  chomp($line);
  my @lineContents = split("\t", $line);
  my $chr          = $lineContents[0];
  my $start        = $lineContents[1];
  my $end          = $lineContents[2];

  next if ($chr eq 'chr' and $start eq 'start' and $end eq 'end');

  my $logPVal;
  my $foldEnrichVal;
  my $fdrVal;
  if (defined($macs2Flag)) {
    $logPVal       = $lineContents[5];# -1*LOG10(pvalue)
    $foldEnrichVal = $lineContents[6];
    $fdrVal        = $lineContents[7];# -1*LOG10(qvalue)
  }else{
    $logPVal       = $lineContents[6];# -1*LOG10(pvalue)
    $foldEnrichVal = $lineContents[7];
    $fdrVal        = $lineContents[8];# FDR(%)
  }
  next if (!&looks_like_number($logPVal) or !&looks_like_number($foldEnrichVal) or !&looks_like_number($fdrVal));
  next if ($logPVal < $pvalue);
  next if ($foldEnrichVal < $foldEnrich);
  if (defined($macs2Flag)) {
    next if ($fdrVal < $fdr);
  }else{
    next if ($fdrVal > $fdr);
  }

  if (defined($original)) {
    print OUT $line, "\n";
  }else{
    $start = $start - 1;
    my $name;
    if (defined($macs2Flag)) {
      my $pileup = $lineContents[4];
      $name      = "$name4thCol=$foldEnrichVal=$pileup=$count";
    }else{
      my $summit = $lineContents[4];
      my $tagNum = $lineContents[5];
      $name      = "$name4thCol=$foldEnrichVal=$tagNum=$summit=$count";
    }
    print OUT join ("\t", ($chr, $start, $end, $name, $foldEnrichVal, '+')), "\n";
    $count++;
  }
}
close(IN);
close(OUT);

&sortBedFile($output, $sort);

sub sortBedFile{
  my ($file, $type) = @_;
  if ($type eq 'n') {
    `sort -t \$'\t' -k 1,1 -k 2,2n $file -o $file`;
  }elsif($type eq 'V') {
    `sort -t \$'\t' -k 1,1V -k 2,2n $file -o $file`;
  }else{
    `sort -t \$'\t' -k 1,1V -k 2,2n $file -o $file`;
  }
}

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

This script is designed for filtering and formatting the macs xls result to bed6

=head1 CAUTIONS

Defaults: -name=(peak), -f=(0), -pval=(0.05), -fdr=(0), -sort=(V)

=head1 SYNOPSIS

macsXlsPeakFilter.pl [options] --name string -x [file] -o [file]

 Options:
    -f | --fold               The cutoff of fold_enrichment
    -o | --output             The output file name
    -x | --xls                The input xls file
    -h | --help               Brief help message
    -fdr                      The cutoff FDR
    -name                     The name for 4th column
    -original                 Output the original record
    -pval                     The cutoff p-value
    -sort                     Sort the ouput ('n' or 'V')
    -macs2                    MACS2 program result input
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

macsXlsPeakFilter.pl [options] --name string -x [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will filter and format the macs xls result to bed6

=cut
