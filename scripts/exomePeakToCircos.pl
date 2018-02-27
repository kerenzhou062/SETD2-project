#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
###develop by K.R.Chow, designed for formatting exomePeak xls file into circos format (fold_enrichment)

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $genome = 'human';
my $genomeSizeFile;
my $output;
my $pvalue = 0.05;
my $sort = "V";
my $diffFC = 1;
my $fdr = 0.05;
my $foldEnrich = 1;
my $xlsFile;
my $span;
my $verbose;

GetOptions (
  "g|genome=s{1,1}"                       =>\$genome,
  "gs|genome_size=s{1,1}"                 =>\$genomeSizeFile,
  "o|output=s{1,1}"                       =>\$output,
  "p|pvalue=s{1,1}"                       =>\$pvalue,
  "s|sort=s{1,1}"                         =>\$sort,
  "x|xls=s{1,1}"                          =>\$xlsFile,
  "fdr=s{1,1}"                            =>\$fdr,
  "fold=s{1,1}"                           =>\$foldEnrich,
  "span=s{1,1}"                           =>\$span,
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
$fdr = log($fdr);

my %hashGenomeCircos;
$hashGenomeCircos{'human'} = 'hs';
$hashGenomeCircos{'mouse'} = 'mm';

my $diffFlag;
if ($diffFC < 0) {
  $diffFlag = -1;
  $diffFC = - log((-$diffFC)) / log(2);
}elsif ($diffFC > 0) {
  $diffFlag = 1;
  $diffFC = log($diffFC) / log(2);
}else{
  print STDERR "Error in -diff.\n";
  exit;
}

if (defined($span)) {
  $span = int($span / 1000);
  $span = $span * 1000;
}

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
if (defined($xlsFile)) {
  <IN>;#skip the header of xls file
}
while (my $line=<IN>) {
  next if ($line =~ /^#/);
  chomp($line);
  my @lineContents    = split("\t", $line);
  my $chr             = $lineContents[0];
  my $start           = $lineContents[1];
  my $end             = $lineContents[2];
  #my $name            = $lineContents[3];
  my $pvalueVal       = $lineContents[4];
  my $strand          = $lineContents[5];
  #my $cdsStart        = $lineContents[6];
  #my $cdsEnd          = $lineContents[7];
  #my $blockCount     = $lineContents[9];
  my $fdrVal;
  my $foldEnrichVal;
  my $diffFCVal;

  if (scalar(@lineContents) == 15) {
    $fdrVal = $lineContents[-2];
    $foldEnrichVal = $lineContents[-1];
    next if (!&looks_like_number($fdrVal) or !&looks_like_number($foldEnrichVal));
    next if ($fdrVal > $fdr);
    next if ($foldEnrichVal < $foldEnrich);
  }elsif (scalar(@lineContents) == 18) {
    if ($count == 1) {
      $pvalue = log ($pvalue);
    }
    $foldEnrichVal = $lineContents[-4];
    $fdrVal        = $lineContents[-3];
    $pvalueVal     = $lineContents[-2];
    $diffFCVal     = $lineContents[-1];
    next if (!&looks_like_number($pvalueVal) or !&looks_like_number($fdrVal) or !&looks_like_number($foldEnrichVal) or !&looks_like_number($diffFCVal));
    next if ($fdrVal > $fdr);
    next if ($foldEnrichVal < $foldEnrich);
    if ($diffFlag < 0) {
      next if ($diffFCVal > $diffFC);
    }elsif ($diffFlag > 0) {
      next if ($diffFCVal < $diffFC);
    }
  }

  next if ($pvalueVal > $pvalue);
  $chr =~ s/chr/$genome/i;
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

sub verbose{
  my $warningText = shift;
  if (defined($verbose)) {
    print STDERR $warningText, "\n";
  }
}

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


################# Abbreviations for this script #################
#
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for formatting exomePeak xls file into circos format (fold_enrichment)

=head1 CAUTIONS

Defaults: --name=(peak), -p=(0.05), -d=(0), -u=(0), -fdr=(0.05), -fold=(1), -diff=(1)
-format=(bed)
For -diff, negative for hypermethylated while positive for hypomethylated

=head1 SYNOPSIS

exomePeakToCircos.pl [options] -x [file] -o [file]

 Options:
    -g | --genome             'human', 'mouse'
    -gs| --genome_size        Genome size file
    -o | --output             The output file name
    -u | --upstream           Elongate the summit towards upstream
    -p | --pvalue             The cutoff p-value
    -s | --sort               Sort the ouput ('n' or 'V')
    -x | --xls                The input xls file
    -h | --help               Brief help message
    -diff                     The cutoff of diff peak
    -fdr                      The cutoff of FDR
    -fold                     The cutoff of fold_enrichment
    -span                     Span size (int)
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

exomePeakToCircos.pl [options] -x [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will format exomePeak xls file into circos format (fold_enrichment)

=cut
