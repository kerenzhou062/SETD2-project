#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for normalizing the counts of peaks(bed3+)

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $bedFile;
my @bamFiles;
my $fileName = "foldEnrichmentResult";
my $outDir = ".";
my $strandedness;
my $scale;
my $verbose;

GetOptions (
  "i|input=s{1,1}"                   =>\$bedFile,
  "b|bam=s{1,}"                      =>\@bamFiles,
  "n|name=s{1,}"                     =>\$fileName,
  "o|outDir=s{1,1}"                  =>\$outDir,
  "s"                                =>\$strandedness,
  "scale=s{1,1}"                     =>\$scale,
  "h|help"                           =>\$help,
  "verbose"                          =>\$verbose,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

if (defined($scale)) {
  if ($scale ne 'log2' and $scale ne 'log10') {
    $scale = 'log2';
    if (defined($verbose)) {
      print STDERR "Unkown -scale argument! -scale will set as log2!\n";
    }
  }
}

my $bamFiles = join " ", @bamFiles;
my $countFile = $outDir .'/'. $fileName .'.count';
#if (defined($strandedness)) {
#  `bedtools multicov -s -bams $bamFiles -bed $bedFile > $countFile`;
#}else{
#  `bedtools multicov -bams $bamFiles -bed $bedFile > $countFile`;
#}

my @columns         = split ("\t", `sed -n 1p $countFile`);
my $colNumber       = scalar(@columns);
my $bamFileNum      = scalar(@bamFiles);
if (($bamFileNum % 2) != 0){
  print STDERR "bam files must be even!\n";
  exit;
}
my $countStartIndex = $colNumber - $bamFileNum;
my %hashTotalCount;
open (IN, "<$countFile") or die "cannot open file: $countFile, $!\n";
while (my $line=<IN>) {
  my @lineContents = split ("\t", $line);
  for (my $i = $countStartIndex; $i < $colNumber; $i++){
    $hashTotalCount{$i} += $lineContents[$i];
  }
}
close(IN);

my $foldFile = $outDir .'/'. $fileName .'.FE';
open (IN, "<$countFile") or die "cannot open file: $countFile, $!\n";
open (OUT, ">$foldFile") or die "cannot open file: $foldFile, $!\n";
while (my $line=<IN>) {
  chomp($line);
  my @lineContents = split ("\t", $line);
  my $peakLength = $lineContents[2] - $lineContents[1];
  my @outputContents = @lineContents[0..($countStartIndex-1)];

  for (my $i = $countStartIndex; $i < $colNumber; $i++){
    my $mappedCounts      = $lineContents[$i];
    my $totalMappedCounts = $hashTotalCount{$i};
    $lineContents[$i]     = sprintf ("%.3f", ($mappedCounts * 1e9)/ ( $totalMappedCounts * $peakLength) );
  }
  for (my $i = $countStartIndex; $i < $colNumber; $i = $i+2) {
    my $foldEnrich = ($lineContents[$i] + 0.001) / ($lineContents[$i+1] + 0.001);
    if (defined($scale)) {
      if ($scale eq 'log10') {
        $foldEnrich = log($foldEnrich);
      }else{
        $foldEnrich = log($foldEnrich) / log(2);
      }
    }
    push (@outputContents, $foldEnrich);
  }
  print OUT join ("\t", (@outputContents)), "\n";
}
close(IN);
close(OUT);

################# Abbreviations for this script #################
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for normalizing (RPKM or RPM) the counts of peaks(bed3+)

=head1 CAUTIONS

Default: -t (rpkm), -s (false), -scale (false)

=head1 SYNOPSIS

peakCountNor.pl [options] --input [file] -b [bam1 bam2] --name name -o [file path]

 Options:
    -i | --input       The reference bed file
    -b | --bam         The bam files (IP1 Input1 IP2 Input2 ...)
    -n | --name        The output file name
    -o | --outDir      The output directory
    -s                 Only report overlap on the same strand
    -scale             Scale the value (log2 or log10)
    -h | --help        Brief help message
    -verbose           Verbose mode
    -man               Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

peakCountNor.pl [options] --input [file] -b [bam1 bam2] --name name -o [file path]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will normalize the counts of peaks(bed3+)

=cut
