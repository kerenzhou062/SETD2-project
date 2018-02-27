#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for extracting the longest isoform of transcripts.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);


my $help   = 0;
my $man    = 0;
my $bed12File;
my $output;
my $type = 'full';
my $verbose;

GetOptions (
  "b|bed12=s{1,1}"                   =>\$bed12File,
  "t|type=s{1,1}"                    =>\$type,
  "o|output=s{1,1}"                  =>\$output,
  "verbose"                          =>\$verbose,
  "h|help"                           =>\$help,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

my %hashGeneTnfo;

if ($type ne 'full' and $type ne 'exon') {
  $type = 'full';
  if (defined($verbose)) {
    print STDERR "Unkown -t, will set as 'full' mode\n";
  }
}

open (IN, "<$bed12File") or die "cannot open file: $bed12File $!";
while (my $line=<IN>) {
  chomp($line);
  my @lineContents                                 = split ("\t", $line);
  my ($geneID, undef , undef, $txID, undef, undef) = split (":", $lineContents[3]);
  my $start                                        = $lineContents[1];
  my $end                                          = $lineContents[2];
  my @blockLengthList                              = split (',', $lineContents[10]);
  my $txLength                                     = $end - $start;
  my $exonLength;
  for (my $i = 0; $i < scalar(@blockLengthList); $i++) {
    $exonLength += $blockLengthList[$i];
  }
  $hashGeneTnfo{$geneID}->{$txID}->{'full'} = $txLength;
  $hashGeneTnfo{$geneID}->{$txID}->{'exon'} = $exonLength;
}

my %hashGeneRemain;

if ($type eq 'full') {
  my @geneIDs = sort keys %hashGeneTnfo;
  foreach my $geneID (@geneIDs) {
    my $longestFull    = 0;
    my $LogestExonLength = 0;
    my @txIDs = sort keys %{$hashGeneTnfo{$geneID}};
    foreach my $txID (@txIDs) {
      my $txLength   = $hashGeneTnfo{$geneID}->{$txID}->{'full'};
      my $exonLength = $hashGeneTnfo{$geneID}->{$txID}->{'exon'};
      if ($txLength > $longestFull) {
        $longestFull             = $txLength;
        $LogestExonLength        = $exonLength;
        $hashGeneRemain{$geneID} = $txID;
      }elsif ($txLength == $longestFull) {
        if ($exonLength >= $LogestExonLength) {
          $longestFull             = $txLength;
          $LogestExonLength        = $exonLength;
          $hashGeneRemain{$geneID} = $txID;
        }
      }
    }
  }
}else{
  my @geneIDs = sort keys %hashGeneTnfo;
  foreach my $geneID (@geneIDs) {
    my $longestFull    = 0;
    my $LogestExonLength = 0;
    my @txIDs = sort keys %{$hashGeneTnfo{$geneID}};
    foreach my $txID (@txIDs) {
      my $txLength   = $hashGeneTnfo{$geneID}->{$txID}->{'full'};
      my $exonLength = $hashGeneTnfo{$geneID}->{$txID}->{'exon'};
      if ($exonLength > $LogestExonLength) {
        $longestFull             = $txLength;
        $LogestExonLength        = $exonLength;
        $hashGeneRemain{$geneID} = $txID;
      }elsif ($exonLength == $LogestExonLength) {
        if ($txLength >= $longestFull) {
          $longestFull             = $txLength;
          $LogestExonLength        = $exonLength;
          $hashGeneRemain{$geneID} = $txID;
        }
      }
    }
  }
}
close(IN);

open (IN, "<$bed12File") or die "cannot open file: $bed12File $!";
open (OUT, ">$output") or die "cannot open file: $output $!";
while (my $line=<IN>) {
  my $_4thCol = (split ("\t", $line) )[3];
  my ($geneID, undef , undef, $txID, undef, undef) = split (":", $_4thCol);
  if ($hashGeneRemain{$geneID} eq $txID) {
    print OUT "$line";
  }
}
close(IN);
close(OUT);

`sort -t \$'\t' -k 1,1V -k 2,2n $output -o $output`;

################# Abbreviations for this script #################
#
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for extracting the longest isoform of transcripts.

=head1 CAUTIONS

Defaults: -t (full)

=head1 SYNOPSIS

longestIsoform.pl [options] -t full -b [file] -o [file]

 Options:
    -b | --bed         The bed12 file
    -o | --output      The output file name
    -t | --type        The the longest type (full or exon)
    -h | --help        Brief help message
    -verbose           Print more error message
    -man               Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

longestIsoform.pl [options] -t full -b [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will extract the longest isoform of transcripts.

=cut