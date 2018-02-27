#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for excluding the overlap intervals from aFile.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $aBed;
my $bBed;
my $strand;
my $split;
my $help   = 0;
my $man    = 0;
my $output = './excludeResult.txt';
my $fractionA;
my $fractionB;

GetOptions (
  "a|aFile=s{1,1}"                   =>\$aBed,
  "b|bFile=s{1,1}"                   =>\$bBed,
  "o|output=s{1,1}"                  =>\$output,
  "f=s{1,1}"                         =>\$fractionA,
  "F=s{1,1}"                         =>\$fractionB,
  "s"                                =>\$strand,
  "split"                            =>\$split,
  "h|help"                           =>\$help,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

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

`bedtools intersect -wa -wb $splitCommand $strandCommand $fractionACommand $fractionBCommand -a $aBed -b $bBed -v > $output`;

`sort -t \$'\t' -k 1,1V -k 2,2n $output -o $output`;

################# Abbreviations for this script #################
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for excluding the overlap intervals from aFile.

=head1 CAUTIONS

If "-s" are not defined, overlaps are treated without respect to strand.

=head1 SYNOPSIS

bed2exclude.pl [options] -a [file] -b [file] -o [file]

 Options:
    -a | --aFile       The A bed file
    -b | --bFile       The B bed file
    -o | --output      The output file name
    -f                 Minimum overlap required as a fraction of A
    -F                 Minimum overlap required as a fraction of B
    -s                 Only report overlap on the same strand
    -split             Reporting overlaps with blocked BED features
    -h | --help        Brief help message
    -man               Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

bed2exclude.pl [options] -a [file] -b [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will exclude the overlap intervals from aFile.

=cut
