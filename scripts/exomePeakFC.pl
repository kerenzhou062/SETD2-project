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
my @overlaps;
my $keepFlag;
my $verbose;

GetOptions (
  "o|output=s{1,1}"                       =>\$output,
  "p|pvalue=s{1,1}"                       =>\$pvalue,
  "s|sort=s{1,1}"                         =>\$sort,
  "x|xls=s{1,}"                           =>\@xlsFiles,
  "cutoff=s{1,1}"                         =>\$cutoff,
  "f=s{1,1}"                              =>\$fractionA,
  "F=s{1,1}"                              =>\$fractionB,
  "FE"                                    =>\$foldEnrichFlag,
  "fdr=s{1,1}"                            =>\$fdr,
  "fold=s{1,1}"                           =>\$foldEnrich,
  "log=s{1,1}"                            =>\$log,
  "overlap=s{1,}"                         =>\@overlaps,
  "keep"                                  =>\$keepFlag,
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

$fdr = log($fdr) / log(10);

if (scalar(@xlsFiles) != scalar(@overlaps) and scalar(@overlaps) != 0) {
  print STDERR ("Invalid -overlap! -overlap should equal to --xls!\n");
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
my %hashTempBed12;
my $xlsFilesNum = scalar(@xlsFiles);
for (my $i = 0; $i < $xlsFilesNum; $i++) {
  my $inFile  = $xlsFiles[$i];
  my @XlsCols = split ("\t", `sed -n '2p' $inFile`);
  if (scalar(@XlsCols) != 15) {
    print STDERR "Unkown file format in -xls!\n";
    exit;
  }
  my $tempBed12      = $outputPath .'/' . rand(1000000) . '.tmp';
  $hashTempBed12{$i} = $tempBed12;
  &tempFile12($inFile, $tempBed12, $i);
}

my $tempCont = $hashTempBed12{0};
my %hashPeakFoldEnrich;
for (my $i = 1; $i < $xlsFilesNum; $i++) {
  my %hashPeak;
  my $tempBed12Cont  = $tempCont;
  my $tempBed12Treat = $hashTempBed12{$i};
  my $overlapTemp    = $outputPath .'/'. rand(1000000) . '_overlap.tmp';
  `bedtools intersect -wao -split -s $fractionACommand $fractionBCommand -a $tempBed12Cont -b $tempBed12Treat > $overlapTemp`;
  open (IN, "<$overlapTemp") or die "Cannot open file: $overlapTemp, $!\n";
  while (my $line=<IN>) {
    chomp($line);
    my @lineContents = split("\t", $line);
    my $peakIdCont  = $lineContents[3];
    my $peakIdTreat = $lineContents[15];
    my $overlapLength = $lineContents[24];
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
  `rm -f $tempBed12Treat $overlapTemp`;
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
  my $bed12 = $hashPeakInfo{0}->{$peakIdCont[$i]}->{'bed12'};
  my $peakIdContFold = $hashPeakInfo{0}->{$peakIdCont[$i]}->{'foldEnrich'};
  print OUT $bed12, "\t". $peakIdContFold;
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



sub tempFile12{
  my ($inFile, $tempFile12 , $tag) = @_;
  open (IN, "<$inFile") or die "Cannot open file: $inFile, $!\n";
  open (TEMP, ">$tempFile12") or die "Cannot open file: $tempFile12, $!\n";

  my $count = 1;
  <IN>;# remove the header
  while (my $line=<IN>) {
    next if ($line =~ /^#/);
    chomp($line);
    my @lineContents  = split("\t", $line);
    my $chr           = $lineContents[0];
    my $start         = $lineContents[1];
    my $end           = $lineContents[2];
    my $pvalueVal     = $lineContents[4];
    my $fdrVal        = $lineContents[-2];
    my $foldEnrichVal = $lineContents[-1];
    next if (!&looks_like_number($fdrVal) or !&looks_like_number($foldEnrichVal));
    next if ($fdrVal > $fdr);
    next if ($foldEnrichVal < $foldEnrich);
    next if ($pvalueVal > $pvalue);

    my $name4thCol = 'peak='. $count;
    if ($tag eq '0') {
      $hashPeakInfo{$tag}->{$name4thCol}->{'bed12'} = join ("\t", @lineContents[0..11]);
    }
    $hashPeakInfo{$tag}->{$name4thCol}->{'foldEnrich'} = $foldEnrichVal;
    $lineContents[3]                                   = $name4thCol;
    print TEMP join ("\t", @lineContents[0..11]), "\n";
    $count++;
  }
  close(IN);
  close(TEMP);
  &sortBedFile($tempFile12, $sort);
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

exomePeakFC.pl [options] -x [file1 file2] -o [file]

 Options:

    -o | --output             The output file name
    -p | --pvalue             The cutoff p-value
    -s | --sort               Sort the ouput ('n' or 'V')
    -x | --xls                The input xls files(Cont.xls Treat1.xls Treat2.xls...)
    -h | --help               Brief help message
    -cutoff                   The cutoff value for zero fold_enrichment
    -f                        Minimum overlap required as a fraction of A
    -F                        Minimum overlap required as a fraction of B
    -FE                       Flag for report the original fold_enrichment
    -fdr                      The cutoff of FDR
    -fold                     The cutoff of fold_enrichment
    -log                      Log scale (postive number)
    -overlap                  Only report overlapped parts (0 0 1 1...)
    -keep                     Keep unique peaks
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

exomePeakFC.pl [options] -x [file1 file2] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will calculate the fold change of control peak under different conditions.

=cut
