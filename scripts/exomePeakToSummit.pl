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
my $bed12File;
my $downstream = 0;
my $name4thCol = "peak";
my $output;
my $upstream = 0;
my $pvalue = 0.05;
my $sort = "V";
my $diffFC = 1;
my $fdr = 0.05;
my $foldEnrich = 1;
my $xlsFile;
my $original;
my $format = 'bed';
my $verbose;

GetOptions (
  "b|bed=s{1,1}"                          =>\$bed12File,
  "d|downstream=s{1,1}"                   =>\$downstream,
  "n|name=s{1,1}"                         =>\$name4thCol,
  "o|output=s{1,1}"                       =>\$output,
  "p|pvalue=s{1,1}"                       =>\$pvalue,
  "u|upstream=s{1,1}"                     =>\$upstream,
  "s|sort=s{1,1}"                         =>\$sort,
  "x|xls=s{1,1}"                          =>\$xlsFile,
  "diff=s{1,1}"                           =>\$diffFC,
  "fdr=s{1,1}"                            =>\$fdr,
  "fold=s{1,1}"                           =>\$foldEnrich,
  "format=s{1,1}"                         =>\$format,
  "original"                              =>\$original,
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
$fdr = log($fdr) / log(10);

my $inputFile;
if (defined($bed12File)) {
  $inputFile = $bed12File;
}else{
  $inputFile = $xlsFile;
}

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

open (IN, "<$inputFile") or die "Cannot open file: $inputFile, $!\n";
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
  my $scoreVal       = $lineContents[4];
  my $strand          = $lineContents[5];
  #my $cdsStart        = $lineContents[6];
  #my $cdsEnd          = $lineContents[7];
  #my $blockCount     = $lineContents[9];
  my $blockLengthCol  = $lineContents[10]; $blockLengthCol =~ s/,$//;
  my $blockStartCol   = $lineContents[11]; $blockStartCol =~ s/,$//;
  my @blockLengthList = split(',', $blockLengthCol);
  my @blockStartList  = split(',', $blockStartCol);
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
      $pvalue = log ($pvalue) / log(10);
    }
    $foldEnrichVal = $lineContents[-4];
    $fdrVal        = $lineContents[-3];
    $scoreVal     = $lineContents[-2];
    $diffFCVal     = $lineContents[-1];
    next if (!&looks_like_number($scoreVal) or !&looks_like_number($fdrVal) or !&looks_like_number($foldEnrichVal) or !&looks_like_number($diffFCVal));
    next if ($fdrVal > $fdr);
    next if ($foldEnrichVal < $foldEnrich);
    if ($diffFlag < 0) {
      next if ($diffFCVal > $diffFC);
    }elsif ($diffFlag > 0) {
      next if ($diffFCVal < $diffFC);
    }
  }

  next if ($scoreVal > $pvalue);

  if (defined($original)) {
    if ($format eq 'xls') {
      print OUT $line, "\n";
    }else{
      print OUT join ("\t", @lineContents[0..11]), "\n";
    }
  }else{
    my $exonTotalLength = 0;
    for (my $i=0; $i < scalar(@blockLengthList); $i++) {
      $exonTotalLength += $blockLengthList[$i];
    }
    my $halfExonLength = ceil($exonTotalLength / 2);
    my $exonAccumLength = 0;
    my $halfBLock = 0;
    for (my $i=0; $i < scalar(@blockLengthList); $i++) {
      if ($exonAccumLength < $halfExonLength) {
        $halfBLock = $i;
        if ( ($exonAccumLength + $blockLengthList[$i]) <= $halfExonLength ) {
          $exonAccumLength += $blockLengthList[$i];
        }
      }else{
        last;
      }
    }
    my $halfBLockStart = $start + $blockStartList[$halfBLock];
    my $summitStart;
    my $summitEnd;
    if ($strand eq '+') {
      $summitEnd = ($halfExonLength - $exonAccumLength) + $halfBLockStart;
      $summitStart = $summitEnd - 1 - $upstream;
      $summitEnd += $downstream;
    }else{
      if ( ($exonTotalLength / 2) == int (($exonTotalLength / 2)) ){
        $summitEnd = ($halfExonLength - $exonAccumLength) + $halfBLockStart + 1;
        $summitStart = $summitEnd - 1;
      }else{
        $summitEnd = ($halfExonLength - $exonAccumLength) + $halfBLockStart;
        $summitStart = $summitEnd - 1;
      }
      $summitEnd += $upstream;
      $summitStart -= $downstream;
    }
    my $summitName;
    if ($name4thCol eq 'original') {
      $summitName = $lineContents[3];
    }else{
      $summitName = $name4thCol . '='. $count;
    }
    print OUT join ( "\t", ($chr, $summitStart, $summitEnd, $summitName, $scoreVal, $strand) ), "\n";
    $count++;
  }
}
close(IN);
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

This script is designed for formatting exomePeak (bed12) to summit (bed6).

=head1 CAUTIONS

Defaults: --name=(peak), -p=(0.05), -d=(0), -u=(0), -fdr=(0.05), -fold=(1), -diff=(1)
-format=(bed)
For -diff, negative for hypermethylated while positive for hypomethylated

=head1 SYNOPSIS

exomePeakToSummit.pl [options] --name string -b [file] -o [file]

 Options:
    -b | --bed                The input bed12 file
    -d | --downstream         Elongate the summit towards downstream
    -n | --name               The name for 4th column (or 'original')
    -o | --output             The output file name
    -u | --upstream           Elongate the summit towards upstream
    -p | --pvalue             The cutoff p-value
    -s | --sort               Sort the ouput ('n' or 'V')
    -x | --xls                The input xls file
    -h | --help               Brief help message
    -diff                     The cutoff of diff peak
    -fdr                      The cutoff of FDR
    -fold                     The cutoff of fold_enrichment
    -format                   The output format (bed or xls)
    -original                 Output the original bed12
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

exomePeakToSummit.pl [options] --name string -b [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will format exomePeak (bed12) to summit (bed6).

=cut
