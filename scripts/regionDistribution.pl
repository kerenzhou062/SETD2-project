#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
###develop by K.R.Chow, designed for computing the distribution patterns on transcripts.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);


my $help   = 0;
my $man    = 0;
my $features = '5utr,cds,stopCodon,3utr';
my $input;
my $bed12File;
my $bed6File;
my $output = './Result.region';
my $sort;
my $cover;
my $size = 200;
my $strandedness;
my $verbose;

GetOptions (
  "f|feature=s{1,1}"                 =>\$features,
  "i|input=s{1,1}"                   =>\$input,
  "o|output=s{1,1}"                  =>\$output,
  "bed12=s{1,1}"                     =>\$bed12File,
  "bed6=s{1,1}"                      =>\$bed6File,
  "cover=s{1,1}"                     =>\$cover,
  "size=s{1,1}"                      =>\$size,
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

my $startCodonSize = $size;
my $stopCodonSize = $size;

if (defined($cover) and looks_like_number($cover)) {
  my %hashBed6;
  open(IN,"<$bed6File") or die "can't open $bed6File, $!\n";
  while (my $line=<IN>) {
    chomp($line);
    my @lineContents = split ("\t", $line);
    my $start = $lineContents[1];
    my $end = $lineContents[2];
    my ($name, $feature, $order) = split (/\|/, $lineContents[3]);
    $hashBed6{$name}->{$feature}->{'order'}->{$order} = $end - $start;
    $hashBed6{$name}->{$feature}->{'size'} += $end - $start;
  }
  close(IN);

  my $utr5AveLen = 0;
  my $utr3AveLen = 0;
  my %hashFeatureSize;
  my @coverFeatures = ('5utr', '3utr');
  foreach my $featureName (keys %hashBed6) {
    foreach my $feature (keys %{$hashBed6{$featureName}}) {
      my $eachSize = $hashBed6{$featureName}->{$feature}->{'size'};
      $hashFeatureSize{$feature}->{'count'}++;
      $hashFeatureSize{$feature}->{'sumSize'} += $eachSize;
    }
  }
  my $fullFeatureLength = 0;
  foreach my $feature (@coverFeatures) {
    my $count = $hashFeatureSize{$feature}->{'count'};
    my $sumSize = $hashFeatureSize{$feature}->{'sumSize'};
    $hashFeatureSize{$feature}->{'average'} = $sumSize / $count;
    $fullFeatureLength += $hashFeatureSize{$feature}->{'average'};
  }
  $startCodonSize = $hashFeatureSize{'5utr'}->{'average'} * $cover;
  $stopCodonSize = $hashFeatureSize{'3utr'}->{'average'} * $cover;
}

my %hashCondon;
my %hashNameFeatureLength;
open(IN,"<$bed6File") or die "can't open $bed6File, $!\n";
while (my $line=<IN>) {
  chomp($line);
  my @lineContents = split ("\t", $line);
  my $start = $lineContents[1];
  my $end = $lineContents[2];
  my ($name, $eachFeature, $order) = split (/\|/, $lineContents[3]);
  if ($eachFeature eq 'cds') {
    if ( exists($hashCondon{$name}->{'stop'})) {
      if ($end > $hashCondon{$name}->{'stop'}) {
        $hashCondon{$name}->{'stop'} = $end;
      }
    }else{
      $hashCondon{$name}->{'stop'} = $end;
    }
    if ($order == 1) {
      $hashCondon{$name}->{'start'} = $start;
    }
  }
  $hashNameFeatureLength{$name}->{$eachFeature} += $end - $start;
}
close(IN);

my $intersectTemp = $outputPath .'/'. rand(1000000) .'_intersect.tmp';
if (defined($sort)) {
  `sort -t \$'\t' -k 1,1V -k 2,2n $tempCenter -o $tempCenter`;
  `sort -t \$'\t' -k 1,1V -k 2,2n $bed6File -o $bed6File`;
}
if (defined($strandedness)) {
  `bedtools intersect -wb -s -nonamecheck -a $tempCenter -b $bed6File > $intersectTemp`;
}else{
  `bedtools intersect -wb -nonamecheck -a $tempCenter -b $bed6File > $intersectTemp`;
}


my @features = split (",", $features);
my %hashFeatures;
foreach my $feature (@features) {
  $hashFeatures{$feature}++;
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
    if (exists($hashTypeCategory{$txType})) {
      if ($hashTypeCategory{$txType} eq 'protein_coding') {
        if (exists($hashFeatures{'stopCodon'}) and exists($hashFeatures{'startCodon'})) {
          if ($featureStrand eq '+') {
            if (abs(($pcEnd - $hashCondon{$featureName}->{'stop'})) <= ($stopCodonSize - 1)) {
              $hashStats{'feature'}->{'stopCodon'}++;
            }elsif (abs(($pcEnd - $hashCondon{$featureName}->{'start'})) <= $startCodonSize) {
              $hashStats{'feature'}->{'startCodon'}++;
            }else{
              $hashStats{'feature'}->{$feature}++;
              $hashStats{'count'}->{$featureName}->{$feature}++;
            }
          }else{
            if (abs(($pcEnd - $hashCondon{$featureName}->{'start'})) <= $stopCodonSize) {
              $hashStats{'feature'}->{'stopCodon'}++;
            }elsif (abs(($pcEnd - $hashCondon{$featureName}->{'stop'})) <= ($startCodonSize - 1)) {
              $hashStats{'feature'}->{'startCodon'}++;
            }else{
              $hashStats{'feature'}->{$feature}++;
              $hashStats{'count'}->{$featureName}->{$feature}++;
            }
          }
        }elsif (exists($hashFeatures{'stopCodon'})) {
          if ($featureStrand eq '+') {
            if (abs(($pcEnd - $hashCondon{$featureName}->{'stop'})) <= ($stopCodonSize - 1)) {
              $hashStats{'feature'}->{'stopCodon'}++;
            }else{
              $hashStats{'feature'}->{$feature}++;
              $hashStats{'count'}->{$featureName}->{$feature}++;
            }
          }else{
            if (abs(($pcEnd - $hashCondon{$featureName}->{'start'})) <= $stopCodonSize) {
              $hashStats{'feature'}->{'stopCodon'}++;
            }else{
              $hashStats{'feature'}->{$feature}++;
              $hashStats{'count'}->{$featureName}->{$feature}++;
            }
          }
        }elsif (exists($hashFeatures{'startCodon'})) {
          if ($featureStrand eq '+') {
            if (abs(($pcEnd - $hashCondon{$featureName}->{'start'})) <= $stopCodonSize) {
              $hashStats{'feature'}->{'startCodon'}++;
            }else{
              $hashStats{'feature'}->{$feature}++;
              $hashStats{'count'}->{$featureName}->{$feature}++;
            }
          }else{
            if (abs(($pcEnd - $hashCondon{$featureName}->{'stop'})) <= ($stopCodonSize - 1)) {
              $hashStats{'feature'}->{'startCodon'}++;
            }else{
              $hashStats{'feature'}->{$feature}++;
              $hashStats{'count'}->{$featureName}->{$feature}++;
            }
          }
        }else{
          $hashStats{'feature'}->{$feature}++;
          $hashStats{'count'}->{$featureName}->{$feature}++;
        }
      }else{
        $hashStats{'feature'}->{$feature}++;
        $hashStats{'count'}->{$featureName}->{$feature}++;
      }
    }else{
      $hashStats{'feature'}->{$feature}++;
      $hashStats{'count'}->{$featureName}->{$feature}++;
    }
  }else{
    next;
  }
}
close(INTERSECT);

if (defined($bed12File)) {
  `rm -f $tempCenter $intersectTemp $bed6File`;
}else{
  `rm -f $tempCenter $intersectTemp`;
}

my %hashFeatureLength;
foreach my $featureName (keys %{$hashStats{'count'}}) {
  if (exists($hashStats{'count'}->{$featureName})) {
    my @eachFeatures = keys %{$hashStats{'count'}->{$featureName}};
    foreach my $eachFeature(@eachFeatures) {
      $hashFeatureLength{$eachFeature}->{'count'}++;
      $hashFeatureLength{$eachFeature}->{'totalLength'} += $hashNameFeatureLength{$featureName}->{$eachFeature};
    }
  }
}

foreach my $feature (keys %hashFeatureLength) {
  my $number = $hashFeatureLength{$feature}->{'count'};
  my $totalLength = $hashFeatureLength{$feature}->{'totalLength'};
  $hashFeatureLength{$feature}->{'average'} = $totalLength / $number;
}
$hashFeatureLength{'startCodon'}->{'average'} = $startCodonSize;
$hashFeatureLength{'stopCodon'}->{'average'} = $stopCodonSize;

open (OUT, ">$output") or die "can't open $output, $!\n";
foreach my $eachFeature (@features) {
  next if (! exists($hashStats{'feature'}->{$eachFeature}));
  my $featureCount = $hashStats{'feature'}->{$eachFeature};
  my $featureAveLength = $hashFeatureLength{$eachFeature}->{'average'};
  my $enrichment = sprintf("%.2f", ($featureCount / $featureAveLength));
  print OUT join ("\t", ($eachFeature, $featureCount, $enrichment)), "\n";
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

This script is designed for computing the distribution of regions of transcripts.

=head1 CAUTIONS

Defaults: -f=('5utr,cds,stopCodon,3utr'), --sort=(false), -strand=(false).

=head1 SYNOPSIS

regionDistribution.pl [options] --input [file] -bed12 [file] -o [file]

 Options:
    -bed12                The annotation bed12 file
    -bed6                 The annotation bed6 file
    -cover                The region coverage of 5'UTR or 3'UTR as startCondon or stop codon
    -size                 Region size for star or stop codon region
    -sort                 Sort bed before using bedtools
    -strand               Only report overlap on the same strand
    -f | --feature        The region feature
    -i | --input          The input bed6 file
    -o | --output         The output file
    -h | --help           Brief help message
    -verbose              Print more warnings
    -man                  Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

regionDistribution.pl [options] --input [file] -bed12 [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will compute the distribution of regions of transcripts.

=cut
