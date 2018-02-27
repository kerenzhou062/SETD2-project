#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for computing the distribution patterns on transcripts.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);


my $help   = 0;
my $man    = 0;
my $input;
my $bed12File;
my $bed6File;
my $output = './Result.gene';
my $sort;
my $strandedness;
my $verbose;

GetOptions (
  "i|input=s{1,1}"                   =>\$input,
  "o|output=s{1,1}"                  =>\$output,
  "bed12=s{1,1}"                     =>\$bed12File,
  "bed6=s{1,1}"                      =>\$bed6File,
  "sort"                             =>\$sort,
  "strand"                           =>\$strandedness,
  "verbose"                          =>\$verbose,
  "h|help"                           =>\$help,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;


my @outputPath  = split (/\\|\//, $output); splice (@outputPath, -1, 1);
if (scalar(@outputPath) == 0) {
  unshift (@outputPath , ".");
}
my $outputPath  = join "/", @outputPath;

my %hashCategoryType = (
  "protein_coding" => ["protein_coding"],
  "lncRNA"         => ["processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding",
                       "sense_intronic", "sense_overlapping", "Retained_intron", "TEC", "known_ncrna", "macro_lncRNA",
                       "bidirectional_promoter_lncrna"],
  "small ncRNA"    => ["rRNA", "Mt_tRNA", "Mt_rRNA", "misc_RNA", "ribozyme", "sRNA", "scaRNA", "vaultRNA"],
  "miRNA"          => ["miRNA"],
  "snoRNA"         => ["snoRNA"],
  "snRNA"          => ["snRNA"],
  "pseudogene"     => ["Mt_tRNA_pseudogene", "tRNA_pseudogene", "snoRNA_pseudogene", "snRNA_pseudogene", "scRNA_pseudogene",
                       "rRNA_pseudogene", "misc_RNA_pseudogene", "miRNA_pseudogene", "pseudogene", "processed_pseudogene", "polymorphic_pseudogene",
                       "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene",
                       "translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene"]
);

my %hashTypeCategory;
foreach my $category(sort keys %hashCategoryType) {
  my @types = @{$hashCategoryType{$category}};
  foreach my $type(@types) {
    $hashTypeCategory{$type} = $category;
  }
}

if (defined($bed12File) and defined($bed6File)) {
  &verbose("-bed6 will be ignore!");
}

my %hashPeak;

my $totalPeak = 1;
my $tempCenter = $outputPath .'/'. rand(1000000) .'_center.tmp';
open (PEAK, "<$input") or die "Cannot open file: $input, $!\n";
open (TEMP, ">$tempCenter") or die "Cannot open file: $tempCenter, $!\n";
while (my $line=<PEAK>) {
  next if ($line =~ /^#/);
  chomp($line);
  my @lineContents = split("\t", $line);
  my $pcStart      = int( ($lineContents[1] + $lineContents[2]) / 2 );
  my $pcEnd        = $pcStart + 1;
  if (scalar(@lineContents) < 6) {
    $lineContents[5] = '+';
  }
  my $peakKey      = join("\t", @lineContents[0..2]) ."\t". $lineContents[5];
  ### judge if the bed peaks contain any duplicates
  if (exists($hashPeak{$peakKey})) {
    next;
  }else{
    $hashPeak{$peakKey}++;
  }
  $lineContents[3] = 'peak='. $totalPeak;
  print TEMP $lineContents[0], "\t", $pcStart, "\t", $pcEnd, "\t", join("\t", @lineContents[3..5]),"\n";
  $totalPeak++;
}
close(PEAK);
close(TEMP);



if (defined($bed12File)) {
  $bed6File = $outputPath .'/'. rand(1000000) .'_bed6.tmp';
  `bed12ToBed6.pl -b $bed12File -o $bed6File`;
}

my $intersectTemp = $outputPath .'/'. rand(1000000) .'_intersect.tmp';
if (defined($sort)) {
  `sort -t \$'\t' -k 1,1V -k 2,2n $tempCenter -o $tempCenter`;
  `sort -t \$'\t' -k 1,1V -k 2,2n $bed6File -o $bed6File`;
}
if (defined($strandedness)) {
  `bedtools intersect -wb -s -nonamecheck -a $tempCenter -b $bed6File > $intersectTemp`;
}else{
  `bedtools intersect -wb -a -nonamecheck -a $tempCenter -b $bed6File > $intersectTemp`;
}


my %hashRecord;
my %hashStats;
open(INTERSECT,"<$intersectTemp") or die "can't open $intersectTemp, $!\n";
while (my $line=<INTERSECT>) {
  chomp($line);
  my @lineContents = split ("\t", $line);
  my $pcEnd = $lineContents[2];
  my $pcName = $lineContents[3];
  my $featureStart = $lineContents[7];
  my $featureEnd = $lineContents[8];
  my ($featureName, $feature, $order) = split (/\|/, $lineContents[9]);
  my $featureStrand = $lineContents[11];
  my ($geneID, $geneName, $geneType, $txID, $txName, $txType) = split (":", $featureName);

  if (!exists($hashRecord{$pcName})) {
    $hashRecord{$pcName}++;
  }else{
    next;
  }
  if (!exists($hashRecord{$txID})) {
    $hashRecord{$txID}++;
  }else{
    next;
  }
  $hashStats{'type'}->{$txType}++;
}
close(INTERSECT);

if (defined($bed12File)) {
  `rm -f $tempCenter $intersectTemp $bed6File`;
}else{
  `rm -f $tempCenter $intersectTemp`;
}

my @txTypes = sort keys %{$hashStats{'type'}};

foreach my $type (@txTypes) {
  my $txCategory;
  if (exists($hashTypeCategory{$type})) {
    $txCategory = $hashTypeCategory{$type};
  }else{
    $txCategory = $type;
  }
  my $txNum = $hashStats{'type'}->{$type};
  $hashStats{'category'}->{$txCategory} += $txNum;
}

my @txCategories = sort {$hashStats{'category'}->{$b} <=> $hashStats{'category'}->{$a} or $a cmp $b} keys %{$hashStats{'category'}};

open (OUT, ">$output") or die "can't open $output, $!\n";
foreach my $txCategory (@txCategories) {
  print OUT join ("\t", $txCategory, $hashStats{'category'}->{$txCategory}), "\n";
}
close(OUT);


sub verbose{
  my $warningText = shift;
  if (defined($verbose)) {
    print STDERR $warningText, "\n";
  }
}

################# Abbreviations for this script #################
#
# bed6File (4thCol) = $geneID:$geneName:$geneType:$txID:$txName:$txType|cds|1
# pcEnd             = end position of a peak center
# pcShiftRe         = length that a peak center shifts from the bigining of a transcript region (5utr, cds, or 3utr)
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for computing the distribution of gene types of transcripts.

=head1 CAUTIONS

Defaults: -sort=(false), -strand=(false).

=head1 SYNOPSIS

geneDistribution.pl [options] --input [file] -bed12 [file] -o [file]

 Options:
    -bed12                The annotation bed12 file
    -bed6                 The annotation bed6 file
    -sort                 Sort bed before using bedtools
    -strand               Only report overlap on the same strand
    -i | --input          The input bed6 file
    -o | --output         The output file
    -h | --help           Brief help message
    -verbose              Print more warnings
    -man                  Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

geneDistribution.pl [options] --input [file] -bed12 [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will compute the distribution of gene types of transcripts.

=cut
