#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for normalizing the counts of peaks(bed3+)

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $aBed;
my $bBed;
my $help   = 0;
my $man    = 0;
my $bedFile;
my @bamFiles;
my $fileName = "normalizedResult";
my $outDir = ".";
my $strandedness;
my $log;
my $type = 'fpkm';
my $verbose;

GetOptions (
  "i|input=s{1,1}"                   =>\$bedFile,
  "b|bam=s{1,}"                      =>\@bamFiles,
  "n|name=s{1,}"                     =>\$fileName,
  "o|outDir=s{1,1}"                  =>\$outDir,
  "s"                                =>\$strandedness,
  "log=s{1,1}"                       =>\$log,
  "t|type=s{1,1}"                    =>\$type,
  "h|help"                           =>\$help,
  "verbose"                          =>\$verbose,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

if ($type ne 'fpkm' and $type ne 'rpkm') {
  $type = 'rpkm';
  if (defined($verbose)) {
    print STDERR "Unkown -t argument! -t will be set as default!\n";
  }
}

if (defined($log)) {
  if (defined($verbose)) {
    if ($log <= 0) {
      print STDERR "Invalid -log: $log! Will set -log to 2";
    }
  }
  $log = 2;
  $log = log($log);
}


my $bamFiles = join " ", @bamFiles;
my $countFile = $outDir .'/'. $fileName .'.count';
if (defined($strandedness)) {
  `bedtools multicov -s -bams $bamFiles -bed $bedFile > $countFile`;
}else{
  `bedtools multicov -bams $bamFiles -bed $bedFile > $countFile`;
}

my @columns         = split ("\t", `sed -n 1p $countFile`);
my $colNumber       = scalar(@columns);
my $bamFileNum      = scalar(@bamFiles);
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

my $rpkmFile;
if ($type eq 'rpkm') {
  $rpkmFile = $outDir .'/'. $fileName .'.rpkm';
}else{
  $rpkmFile = $outDir .'/'. $fileName .'.rpm';
}
open (IN, "<$countFile") or die "cannot open file: $countFile, $!\n";
open (OUT, ">$rpkmFile") or die "cannot open file: $rpkmFile, $!\n";
while (my $line=<IN>) {
  my @lineContents = split ("\t", $line);
  my $peakLength = $lineContents[2] - $lineContents[1];
  for (my $i = $countStartIndex; $i < $colNumber; $i++){
    my $mappedCounts      = $lineContents[$i];
    my $totalMappedCounts = $hashTotalCount{$i};
    if ($type eq 'rpkm') {
      $lineContents[$i]     = sprintf ("%.3f", ($mappedCounts * 1e9)/ ( $totalMappedCounts * $peakLength) );
    }else{
      $lineContents[$i]     = sprintf ("%.3f", ($mappedCounts * 1e6)/ $totalMappedCounts );
    }
    if (defined($log)) {
      $lineContents[$i] = log($lineContents[$i] + 0.001) / $log;
    }
  }
  print OUT join ("\t", @lineContents), "\n";
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

Default: -t (rpkm), -s (false), -log (false)

=head1 SYNOPSIS

peakCountNor.pl [options] --input [file] -b [bam1 bam2] --name name -o [file path]

 Options:
    -i | --input       The reference bed file
    -b | --bam         The bam files (bam1 bam2 ...)
    -n | --name        The output file name
    -o | --outDir      The output directory
    -s                 Only report overlap on the same strand
    -log               Log scale the value (postive number)
    -t | --type        The normalize type (rpkm or rpm)
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
