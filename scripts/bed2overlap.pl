#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for reporting the intervals that bFile overlap with aFlie.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $aBed;
my $bBed;
my $strandedness;
my $help   = 0;
my $man    = 0;
my $name  = 'peak';
my $output = './overlapResult.txt';
my $aOverlap;
my $split;
my $fractionA;
my $fractionB;

GetOptions (
  "a|aFile=s{1,1}"                   =>\$aBed,
  "b|bFile=s{1,1}"                   =>\$bBed,
  "n|Name=s{1,1}"                    =>\$name,
  "o|output=s{1,1}"                  =>\$output,
  "aOver"                            =>\$aOverlap,
  "f=s{1,1}"                         =>\$fractionA,
  "F=s{1,1}"                         =>\$fractionB,
  "s"                                =>\$strandedness,
  "split"                            =>\$split,
  "h|help"                           =>\$help,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

my @outputPath  = split (/\\|\//, $output); splice (@outputPath, -1, 1);
if (scalar(@outputPath) == 0) {
  unshift (@outputPath , ".");
}
my $outputPath  = join "/", @outputPath;

my @aBedCols = split ("\t", `sed -n '2p' $aBed`);
my $aBedColNum = scalar(@aBedCols);

my $overlapTemp = $outputPath .'/'. rand(1000000) .'.tmp';
my $splitCommand;
if (defined($split)) {
  $splitCommand = '-split';
}else{
  $splitCommand = '';
}
my $strandCommand;
if (defined($split)) {
  $strandCommand = '-s';
}else{
  $strandCommand = '';
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

`bedtools intersect -wao $splitCommand $strandCommand $fractionACommand $fractionBCommand -a $aBed -b $bBed > $overlapTemp`;

my $overlapCount = 1;
my %hashPeak;
open (INTERSECT, "<$overlapTemp") or die "cannot open file($overlapTemp), $!\n";
open (OUT, ">$output") or die "can't open $output, $!\n";
while (my $line=<INTERSECT>) {
  chomp($line);
  my @lineContents   = split("\t", $line);
  my $overlapLength  = $lineContents[-1];
  next if ($overlapLength == 0);
  if (defined($aOverlap)) {
    my $hashKey = join ("\t", @lineContents[0..($aBedColNum-1)]);
    if (! exists($hashPeak{$hashKey})){
      print OUT join ("\t", @lineContents[0..($aBedColNum-1)]), "\n";
      $hashPeak{$hashKey}++;
    }
  }else{
    splice (@lineContents, 0, $aBedColNum);
    my $hashKey = join ("\t", (@lineContents[0..2], $lineContents[5]));
    if (! exists($hashPeak{$hashKey})){
      $lineContents[3] = $name ."=". $overlapCount;
      $lineContents[4] = 0;
      print OUT join ("\t", @lineContents[0..5]), "\n";
      $hashPeak{$hashKey}++;
      $overlapCount++;
    }
  }
}

`rm -f $overlapTemp`;
`sort -t \$'\t' -k 1,1V -k 2,2n $output -o $output`;

################# Abbreviations for this script #################
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for reporting the intervals that bFile overlap with aFlie.

=head1 CAUTIONS
Defaults: -s=(false), -bpeak=(false), -n=(peak)
Report the overlapped peak by default.
If -bpeak is set, -n will be ignore;

If "-s" are not defined, overlaps are reported without respect to strand.

=head1 SYNOPSIS

bed2overlap.pl [options] -a [file] -b [file] -o [file]

 Options:
    -a | --aFile       The A bed file
    -n|--name          The prefix name for 4th in output
    -b | --bFile       The B bed file
    -o | --output      The output file name
    -aOver             Report the original record in A file
    -f                 Minimum overlap required as a fraction of A
    -F                 Minimum overlap required as a fraction of B
    -s                 Only report overlap on the same strand
    -split             Reporting overlaps with blocked BED features
    -h | --help        Brief help message
    -man               Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

bed2overlap.pl [options] -a [file] -b [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will report the intervals that bFile overlap in aFlie(bedtools intersect).

=cut
