#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
###develop by K.R.Chow, designed for formating MACS xls file into circos format (fold_enrichment)

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $genome = "human";
my $genomeSizeFile;
my $xlsFile;
my $output;
my $foldEnrich = 0;
my $name4thCol   = 'peak';
my $fdr        = 0.05;
my $pvalue     = 0.05;
my $sort       = "V";
my $span;
my $macs2Flag;

my $original;
my $verbose;

GetOptions (
  "f|fold=s{1,1}"                         =>\$foldEnrich,
  "g|genome=s{1,1}"                       =>\$genome,
  "gs|genome_size=s{1,1}"                 =>\$genomeSizeFile,
  "o|output=s{1,1}"                       =>\$output,
  "x|xls=s{1,1}"                          =>\$xlsFile,
  "fdr=s{1,1}"                            =>\$fdr,
  "pval=s{1,1}"                           =>\$pvalue,
  "sort=s{1,1}"                           =>\$sort,
  "span=s{1,1}"                           =>\$span,
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

if (defined($span)) {
  $span = int($span / 1000);
  $span = $span * 1000;
}

my %hashGenomeCircos;
$hashGenomeCircos{'human'} = 'hs';
$hashGenomeCircos{'mouse'} = 'mm';

if (!exists($hashGenomeCircos{$genome})) {
  &verbose("Unkown -g: $genome, 'human' is set.");
  $genome = 'human';
}

$genome = $hashGenomeCircos{$genome};

my %hashChrSplit;
if (defined($span)) {
  open (IN, "$genomeSizeFile") or die "Cannot open file: $genomeSizeFile, $!\n";
  while (my $line=<IN>) {
    chomp($line);
    my ($chr, $size) = split ("\t", $line);
    $chr =~ s/chr/$genome/i;
    my $binSplitMax = int($size / $span);
    $hashChrSplit{$chr}->{'size'} = $size;
    $hashChrSplit{$chr}->{'block'} = $binSplitMax;
  }
  close (IN);
}


my %hashChrBinStats;
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

  $chr =~ s/chr/$genome/i;
  $start = $start - 1;
  if (defined($span)) {
    my $center = int( ($start + $end) / 2 );
    my $bin = int($center / $span);
    $hashChrBinStats{$chr}->{$bin}->{'sum'} += $foldEnrichVal;
  }else{
    print OUT join ("\t", ($chr, $start, $end, $foldEnrichVal)), "\n";
  }
}

close(IN);

foreach my $chr (sort keys %hashChrBinStats) {
  my $totalBin = $hashChrSplit{$chr}->{'block'};
  for (my $i = 0; $i < $totalBin; $i++) {
    my $binVal;
    my $binStart = $i * $span;
    my $binEnd = $binStart + $span - 1;
    if ( $i == ($totalBin - 1) ){
      $binEnd = $hashChrSplit{$chr}->{'size'} - 1;
    }
    if (exists($hashChrBinStats{$chr}->{$i})) {
      my $sum = $hashChrBinStats{$chr}->{$i}->{'sum'};
      $binVal = $sum;
    }else{
      $binVal = 0;
    }
    print OUT join ("\t", ($chr, $binStart, $binEnd, $binVal)), "\n";
  }
}

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

This script is designed for formating MACS xls file into circos format (fold_enrichment)

=head1 CAUTIONS

Defaults: -f=(0), -pval=(0.05), -fdr=(0), -sort=(V)

=head1 SYNOPSIS

macsXlsToCircos.pl [options] -x [file] -o [file]

 Options:
    -f | --fold               The cutoff of fold_enrichment
    -g | --genome             "human" or "mouse"
    -gs| --genome_size        Genome size file
    -o | --output             The output file name
    -x | --xls                The input xls file
    -h | --help               Brief help message
    -fdr                      The cutoff FDR
    -pval                     The cutoff p-value
    -sort                     Sort the ouput ('n' or 'V')
    -span                     Span size (int)
    -macs2                    MACS2 program result input
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

macsXlsToCircos.pl [options] -x [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will format MACS xls file into circos format (fold_enrichment)

=cut
