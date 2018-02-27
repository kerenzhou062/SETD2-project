#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
###develop by K.R.Chow, designed for formatting exomePeak (bed6) to summit (bed6).
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
my $fdr = 0.05;
my $log;
my @overlaps;
my $foldEnrichFlag;
my $macs2Flag;
my $verbose;

GetOptions (
  "o|output=s{1,1}"                       =>\$output,
  "p|pvalue=s{1,1}"                       =>\$pvalue,
  "s|sort=s{1,1}"                         =>\$sort,
  "x|xls=s{1,}"                           =>\@xlsFiles,
  "cutoff=s{1,1}"                         =>\$cutoff,
  "f=s{1,1}"                              =>\$fractionA,
  "F=s{1,1}"                              =>\$fractionB,
  "fdr=s{1,1}"                            =>\$fdr,
  "fold=s{1,1}"                           =>\$foldEnrich,
  "log=s{1,1}"                            =>\$log,
  "overlap=s{1,}"                         =>\@overlaps,
  "FE"                                    =>\$foldEnrichFlag,
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
  if ($log <= 0) {
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


if (scalar(@xlsFiles) != scalar(@overlaps) and scalar(@overlaps) != 0) {
  print STDERR ("Invalid -overlap! -overlap should equal to --xls!");
  exit;
}

if (scalar(@overlaps) != 0) {
  $overlaps[0] = 0;
  for (my $i = 0; $i < scalar(@overlaps); $i++) {
    if ($overlaps[$i] != 1 and $overlaps[$i] != 0) {
      print STDERR ("Invalid -overlap! Only 1 and 0 are allow in -overlap!\n");
      exit;
    }
  }
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

my $tempCont = $hashTempBed6{0};
my %hashPeakFoldEnrich;
for (my $i = 1; $i < $xlsFilesNum; $i++) {
  my %hashPeak;
  my $tempBed6Cont  = $tempCont;
  my $tempBed6Treat = $hashTempBed6{$i};
  my $overlapTemp    = $outputPath .'/'. rand(1000000) . '_overlap.tmp';
  `bedtools intersect -wao $fractionACommand $fractionBCommand -a $tempBed6Cont -b $tempBed6Treat > $overlapTemp`;
  open (IN, "<$overlapTemp") or die "Cannot open file: $overlapTemp, $!\n";
  while (my $line=<IN>) {
    chomp($line);
    my @lineContents = split("\t", $line);
    my $peakIdCont  = $lineContents[3];
    my $peakIdTreat = $lineContents[9];
    my $overlapLength = $lineContents[12];
    next if ($overlapLength == 0);
    if (exists($hashPeak{$peakIdCont})) {
      if ($hashPeak{$peakIdCont}) {
        if ($overlapLength > $hashPeak{$peakIdCont}) {
          $hashPeak{$peakIdCont} = $overlapLength;
        }else{
          next;
        }
      }
    }else{
      $hashPeak{$peakIdCont} = $overlapLength;
    }
    my $peakIdContFold                     = $hashPeakInfo{0}->{$peakIdCont}->{'foldEnrich'};
    my $peakIdTreatFold                    = $hashPeakInfo{$i}->{$peakIdTreat}->{'foldEnrich'};
    my $enrichLevel;
    if (defined($foldEnrichFlag)) {
      $enrichLevel = $peakIdTreatFold;
    }else{
      $enrichLevel = ($peakIdTreatFold + $cutoff) / ($peakIdContFold + $cutoff);
    }
    if (defined($log)) {
      $enrichLevel = log($enrichLevel) / $log;
    }
    $hashPeakFoldEnrich{$peakIdCont}->{$i} = $enrichLevel;
  }
  close(IN);
  `rm -f $tempBed6Treat $overlapTemp`;
}
`rm -f $tempCont`;

open (OUT, ">$output") or die "Cannot open file: $output, $!\n";
my @peakIdCont = sort keys %{$hashPeakInfo{0}};

for (my $i = 0; $i < scalar(@peakIdCont); $i++) {
  if (scalar(@overlaps) != 0) {
    my $nextFlag = 0;
    for (my $j = 1; $j < $xlsFilesNum; $j++) {
      if (!exists($hashPeakFoldEnrich{$peakIdCont[$i]}->{$j}) and $overlaps[$j] == 1) {
        $nextFlag = 1;
        last;
      }
    }
    next if ($nextFlag);
  }
  my $bed6 = $hashPeakInfo{0}->{$peakIdCont[$i]}->{'bed6'};
  my $peakIdContFold = $hashPeakInfo{0}->{$peakIdCont[$i]}->{'foldEnrich'};
  print OUT $bed6, "\t". $peakIdContFold;
  for (my $j = 1; $j < $xlsFilesNum; $j++) {
    my $enrichLevel;
    if (!exists($hashPeakFoldEnrich{$peakIdCont[$i]}->{$j})) {
      if (defined($foldEnrichFlag)) {
        $enrichLevel = 0;
      }else{
        $enrichLevel = $cutoff / ($peakIdContFold + $cutoff);
      }
      if (defined($log)) {
        if ($enrichLevel == 0) {
          $enrichLevel = -Inf;
        }else{
          $enrichLevel = log($enrichLevel) / $log;
        }
      }
    }else{
      $enrichLevel = $hashPeakFoldEnrich{$peakIdCont[$i]}->{$j};
    }
    print OUT "\t". $enrichLevel;
  }
  print OUT "\n";
}

&sortBedFile($output, $sort);



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
    if ($tag eq '0') {
      $hashPeakInfo{$tag}->{$name4thCol}->{'bed6'} = $peakInfo;
    }
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

This script is designed for calculating the fold change of control peak under different conditions.

=head1 CAUTIONS

Defaults: -p=(0.05), -f=(1), -fdr=(0.05), -log=(false), -FE=(false)

=head1 SYNOPSIS

macsPeakFC.pl [options] -x [file1 file2] -o [file]

 Options:

    -o | --output             The output file name
    -p | --pvalue             The cutoff p-value
    -s | --sort               Sort the ouput ('n' or 'V')
    -x | --xls                The input xls files(Cont.xls Treat1.xls Treat2.xls...)
    -h | --help               Brief help message
    -cutoff                   The cutoff value for zero fold_enrichment
    -f                        Minimum overlap required as a fraction of A
    -F                        Minimum overlap required as a fraction of B
    -fdr                      The cutoff of FDR
    -fold                     The cutoff of fold_enrichment
    -log                      Log scale (postive number)
    -overlap                  Only report overlapped parts (0 0 1 1...)
    -FE                       Flag for report the original fold_enrichment
    -macs2                    MACS2 program result input
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

macsPeakFC.pl [options] -x [file1 file2] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will calculate the fold change of control peak under different conditions.

=cut
