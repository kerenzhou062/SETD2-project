#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
###develop by K.R.Chow, designed for calculating the foldchange between 2 macsPeak files
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $output;
my $pvalue = 0.05;
my $sort = "V";
my $foldEnrich = 1;
my @xlsFiles;
my $cutoff = 0.001;
my $fractionA;
my $fractionB;
my $foldEnrichFlag;
my $fdr = 0.05;
my $log;
my $keepNonOverlapAFlag;
my $keepNonOverlapBFlag;
my $macs2Flag;
my $verbose;

GetOptions (
  "o|output=s{1,1}"                       =>\$output,
  "p|pvalue=s{1,1}"                       =>\$pvalue,
  "s|sort=s{1,1}"                         =>\$sort,
  "x|xls=s{1,2}"                          =>\@xlsFiles,
  "cutoff=s{1,1}"                         =>\$cutoff,
  "f=s{1,1}"                              =>\$fractionA,
  "F=s{1,1}"                              =>\$fractionB,
  "fdr=s{1,1}"                            =>\$fdr,
  "fold=s{1,1}"                           =>\$foldEnrich,
  "log=s{1,1}"                            =>\$log,
  "FE"                                    =>\$foldEnrichFlag,
  "nonOverlapA"                           =>\$keepNonOverlapAFlag,
  "nonOverlapB"                           =>\$keepNonOverlapBFlag,
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

if (defined($log)) {
  if (!&looks_like_number($log) or $log <= 0) {
    &verbose("Invalid -log: $log! Will set -log to 2");
    $log = 2;
  }
  $log = log($log);
}

if (defined($macs2Flag)) {
  $pvalue = -1 * log($pvalue) / log(10);
  $fdr = -1 * log($fdr);
}else{
  $pvalue = -10 * log($pvalue) / log(10);
  $fdr = $fdr * 100;
}

my $fractionACommand;
if (defined($fractionA)) {
  $fractionACommand = "-f $fractionA";
}else{
  $fractionACommand = '';
}

my $fractionBCommand;
if (defined($fractionB)) {
  $fractionBCommand = "-F $fractionB";
}else{
  $fractionBCommand = '';
}

my @outputPath  = split (/\\|\//, $output); splice (@outputPath, -1, 1);
if (scalar(@outputPath) == 0) {
  unshift (@outputPath , ".");
}
my $outputPath  = join "/", @outputPath;

my %hashPeakInfo;
my %hashTempBed6;
my $xlsFilesNum = scalar(@xlsFiles);
for (my $i = 0; $i < $xlsFilesNum; $i++) {
  my $inFile  = $xlsFiles[$i];
  my @XlsCols = split ("\t", `sed -n '2p' $inFile`);
  if (scalar(@XlsCols) != 9) {
    print STDERR "Unkown file format in -xls!\n";
    exit;
  }
  my $tempBed6      = $outputPath .'/' . rand(1000000) . '.tmp';
  $hashTempBed6{$i} = $tempBed6;
  &tempFile6($inFile, $tempBed6, $i);
}

my %hashPeakIDPair;
my $tempBed6Cont  = $hashTempBed6{0};
my $tempBed6Treat = $hashTempBed6{1};
my $overlapTemp   = $outputPath .'/'. rand(1000000) . '_overlap.tmp';
`bedtools intersect -wao -split -s $fractionACommand $fractionBCommand -a $tempBed6Cont -b $tempBed6Treat > $overlapTemp`;
open (IN, "<$overlapTemp") or die "Cannot open file:$overlapTemp";
while (my $line=<IN>) {
  chomp($line);
  my @lineContents  = split ("\t", $line);
  my $overlapLength = $lineContents[-1];
  my $peakIdCont    = $lineContents[3];
  my $peakIdTreat   = $lineContents[9];
  if ($overlapLength != 0) {
    if (!exists($hashPeakIDPair{'overlap'}->{$peakIdCont})) {
      $hashPeakIDPair{'overlap'}->{$peakIdCont}->{'paired'} = $peakIdTreat;
      $hashPeakIDPair{'overlap'}->{$peakIdCont}->{'OverlapLength'} = $overlapLength;
    }else{
      my $existedOverlagLength = $hashPeakIDPair{'overlap'}->{$peakIdCont}->{'OverlapLength'};
      if ($existedOverlagLength < $overlapLength) {
        $hashPeakIDPair{'overlap'}->{$peakIdCont}->{'paired'} = $peakIdTreat;
        $hashPeakIDPair{'overlap'}->{$peakIdCont}->{'OverlapLength'} = $overlapLength;
      }
    }
  }else{
    if (defined($keepNonOverlapAFlag)) {
      $hashPeakIDPair{'nonOverlapA'}->{$peakIdCont}++;
    }
  }
}
close(IN);

`bedtools intersect -wao -split -s $fractionACommand $fractionBCommand -a $tempBed6Treat -b $tempBed6Cont > $overlapTemp`;
open (IN, "<$overlapTemp") or die "Cannot open file:$overlapTemp";
while (my $line=<IN>) {
  chomp($line);
  my @lineContents  = split ("\t", $line);
  my $overlapLength = $lineContents[-1];
  if (defined($keepNonOverlapBFlag)) {
    if ($overlapLength == 0) {
      my $peakIdTreat = $lineContents[3];
      $hashPeakIDPair{'nonOverlapB'}->{$peakIdTreat}++;
    }
  }
}
`rm -f $tempBed6Cont $tempBed6Treat $overlapTemp`;

open (OUT, ">$output") or die "Cannot open file: $output, $!\n";

my @overlappedIDs = sort keys %{$hashPeakIDPair{'overlap'}};
foreach my $contPeakId (@overlappedIDs) {
  my $bed6            = $hashPeakInfo{0}->{$contPeakId}->{'bed6'};
  my $pairedPeakID    = $hashPeakIDPair{'overlap'}->{$contPeakId}->{'paired'};
  my $peakIdContFold  = $hashPeakInfo{0}->{$contPeakId}->{'foldEnrich'};
  my $peakIdTreatFold = $hashPeakInfo{1}->{$pairedPeakID}->{'foldEnrich'};
  if (defined($foldEnrichFlag)) {
    print OUT join ("\t", ($bed6, $peakIdContFold, $peakIdTreatFold)), "\n";
  }else{
    my $enrichLevel = $peakIdTreatFold / $peakIdContFold;
    if (defined($log)) {
      $enrichLevel = log($enrichLevel) / $log;
    }
    print OUT join ("\t", ($bed6, $enrichLevel)), "\n";
  }
}

if (defined($keepNonOverlapAFlag)) {
  my @nonOverlapAIDs = sort keys %{$hashPeakIDPair{'nonOverlapA'}};
  foreach my $contPeakId (@nonOverlapAIDs) {
    my $bed6           = $hashPeakInfo{0}->{$contPeakId}->{'bed6'};
    my $peakIdContFold  = $hashPeakInfo{0}->{$contPeakId}->{'foldEnrich'};
    my $peakIdTreatFold = 0;
    if (defined($foldEnrichFlag)) {
      print OUT join ("\t", ($bed6, $peakIdContFold, $peakIdTreatFold)), "\n";
    }else{
      my $enrichLevel = ($peakIdTreatFold + $cutoff) / ($peakIdContFold + $cutoff);
      if (defined($log)) {
        $enrichLevel = log($enrichLevel) / $log;
      }
      print OUT join ("\t", ($bed6, $enrichLevel)), "\n";
    }
  }
}

if (defined($keepNonOverlapBFlag) and !defined($foldEnrichFlag)) {
  my @nonOverlapBIDs = sort keys %{$hashPeakIDPair{'nonOverlapB'}};
  foreach my $peakId (@nonOverlapBIDs) {
    my $bed6           = $hashPeakInfo{1}->{$peakId}->{'bed6'};
    my $peakIdTreatFold = $hashPeakInfo{1}->{$peakId}->{'foldEnrich'};
    my $peakIdContFold  = 0;
    my $enrichLevel     = ($peakIdTreatFold + $cutoff) / ($peakIdContFold + $cutoff);
    if (defined($foldEnrichFlag)) {
      print OUT join ("\t", ($bed6, $peakIdContFold, $peakIdTreatFold)), "\n";
    }else{
      my $enrichLevel = ($peakIdTreatFold + $cutoff) / ($peakIdContFold + $cutoff);
      if (defined($log)) {
        $enrichLevel = log($enrichLevel) / $log;
      }
      print OUT join ("\t", ($bed6, $enrichLevel)), "\n";
    }
    print OUT join ("\t", ($bed6, $enrichLevel)), "\n";
  }
}
close(OUT);

&sortBedFile($output, $sort);


##################################### subroutines ########################
sub tempFile6{
  my ($inFile, $tempFile6 , $tag) = @_;
  open (IN, "<$inFile") or die "Cannot open file: $inFile, $!\n";
  open (TEMP, ">$tempFile6") or die "Cannot open file: $tempFile6, $!\n";

  my $count = 1;
  while (my $line=<IN>) {
    next if ($line =~ /^(#|\s+)/);
    chomp($line);
    my @lineContents  = split("\t", $line);
    my $chr           = $lineContents[0];
    my $start         = $lineContents[1];
    my $end           = $lineContents[2];

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
    next if (!&looks_like_number($fdrVal) or !&looks_like_number($foldEnrichVal));
    next if ($logPVal < $pvalue);
    next if ($foldEnrichVal < $foldEnrich);
    if (defined($macs2Flag)) {
      next if ($fdrVal < $fdr);
    }else{
      next if ($fdrVal > $fdr);
    }

    $start = $start - 1;
    my $name4thCol;
    if (defined($macs2Flag)) {
      my $pileup = $lineContents[4];
      $name4thCol      = "peak=$foldEnrichVal=$pileup=$count";
    }else{
      my $summit = $lineContents[4];
      my $tagNum = $lineContents[5];
      $name4thCol      = "peak=$foldEnrichVal=$tagNum=$summit=$count";
    }
    my $peakInfo   = join ("\t", ($chr, $start, $end, $name4thCol, $foldEnrichVal, '+') );
    $hashPeakInfo{$tag}->{$name4thCol}->{'bed6'} = $peakInfo;
    $hashPeakInfo{$tag}->{$name4thCol}->{'foldEnrich'} = $foldEnrichVal;
    print TEMP $peakInfo, "\n";
    $count++;
  }
  close(IN);
  close(TEMP);
  &sortBedFile($tempFile6, $sort);
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

This script is designed for calculating the foldchange between 2 macsPeak files.

=head1 CAUTIONS

Defaults: -p=(0.05), -f=(1), -fdr=(0.05), -log=(false), -FE=(false)

=head1 SYNOPSIS

macs2PeakFC.pl [options] -x [file1 file2] -o [file]

 Options:

    -o | --output             The output file name
    -p | --pvalue             The cutoff p-value
    -s | --sort               Sort the ouput ('n' or 'V')
    -x | --xls                The input xls files(Cont.xls Treat.xls)
    -h | --help               Brief help message
    -cutoff                   The cutoff value for zero fold_enrichment
    -f                        Minimum overlap required as a fraction of A
    -F                        Minimum overlap required as a fraction of B
    -fdr                      The cutoff of FDR
    -fold                     The cutoff of fold_enrichment
    -log                      Log scale (postive number)
    -FE                       Flag for report the original fold_enrichment
    -nonOverlapA              keep the non-overlapped part in Control
    -nonOverlapB              keep the non-overlapped part in Treat
    -macs2                    MACS2 program result input
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

macs2PeakFC.pl [options] -x [file1 file2] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will calculate the fold change of control peak under different conditions.

=cut
