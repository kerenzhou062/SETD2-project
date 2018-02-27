#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for extracting specific gene types from bed12 based on gtf.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);


my $help   = 0;
my $man    = 0;
my $bedFile;
my $genome = 'hg19';
my $output;
my $gff3;
my $gtf;
my @geneTypes;
my $verbose;

GetOptions (
  "b|bed=s{1,1}"                     =>\$bedFile,
  "g|genome=s{1,1}"                  =>\$genome,
  "o|output=s{1,1}"                  =>\$output,
  "t|type=s{1,}"                     =>\@geneTypes,
  "gff3=s{1,1}"                      =>\$gff3,
  "gtf=s{1,1}"                       =>\$gtf,
  "verbose"                          =>\$verbose,
  "h|help"                           =>\$help,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;


my %hashGenomeType;
$hashGenomeType{"hg19"} = {
  "mRNA"       => ["protein_coding"],
  "lncRNA"     => ["processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding",  "sense_intronic", "sense_overlapping", "Retained_intron", "TEC", "known_ncrna", "macro_lncRNA", "bidirectional_promoter_lncrna"],
  "sRNA"       => ["snRNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "misc_RNA", "miRNA", "ribozyme", "sRNA", "scaRNA", "vaultRNA"],
  "pseudogene" => ["Mt_tRNA_pseudogene", "tRNA_pseudogene", "snoRNA_pseudogene", "snRNA_pseudogene", "scRNA_pseudogene", "rRNA_pseudogene", "misc_RNA_pseudogene", "miRNA_pseudogene", "pseudogene", "processed_pseudogene", "polymorphic_pseudogene", "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene"],
  "all"        => 1
};
$hashGenomeType{"hg38"} = {
  "mRNA"       => ["protein_coding"],
  "lncRNA"     => ["processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding",  "sense_intronic", "sense_overlapping", "Retained_intron", "TEC", "known_ncrna", "macro_lncRNA", "bidirectional_promoter_lncrna"],
  "sRNA"       => ["snRNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "misc_RNA", "miRNA", "ribozyme", "sRNA", "scaRNA", "vaultRNA"],
  "pseudogene" => ["Mt_tRNA_pseudogene", "tRNA_pseudogene", "snoRNA_pseudogene", "snRNA_pseudogene", "scRNA_pseudogene", "rRNA_pseudogene", "misc_RNA_pseudogene", "miRNA_pseudogene", "pseudogene", "processed_pseudogene", "polymorphic_pseudogene", "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene"],
  "all"        => 1
};

my %hashExtGeneType;
my $allFlag;
foreach my $geneType(@geneTypes) {
  if ($geneType ne 'all') {
    if ( !exists($hashGenomeType{$genome}) ){
      if (defined($verbose)) {
        my $genomes = join ", ", keys %hashGenomeType;
        print STDERR "Invalid --genome ($genomes).\n";
      }
      exit 0;
    }
    if (!exists($hashGenomeType{$genome}->{$geneType})){
      if (defined($verbose)) {
        my $genomeTypes = join ", ", keys %{$hashGenomeType{$genome}};
        print STDERR "Invalid -geneType ($genomeTypes).\n";
      }
      exit 0;
    }
    my @ExtGeneTypes = @{$hashGenomeType{$genome}->{$geneType}};
    foreach my $ExtGeneType (@ExtGeneTypes) {
      $hashExtGeneType{$ExtGeneType}++;
    }
  }else{
    $allFlag = 1;
  }
}

my $annotation;
if (defined($gtf)) {
  $annotation = $gtf;
}elsif (defined($gff3)) {
  $annotation = $gff3;
}else{
  print STDERR "Unkown annotation file!\n";
  exit;
}

my %hashTxIdInfo;
open (IN, "<$annotation") or die "cannot open file: $annotation $!";
while(my $line=<IN>) {
  if ($line =~ m/^#/) {
    next;
  }
  chomp($line);
  my @lineContents       = split /\t/, $line;
  my $type       = $lineContents[2]; # gene, transcript or exon
  my $start = $lineContents[3] - 1;
  my $end = $lineContents[4];
  my $attribute  = $lineContents[8];
  next if ($type ne 'transcript');

  my %hashAttribute;
  if (defined($gff3)) {
    my @attributes = split (";", $attribute);
    foreach (@attributes) {
      my @item                     = split ("=", $_);
      my $item_name                = $item[0];
      my $item_value               = $item[1];
      $hashAttribute{$item_name} = $item_value;
    }
  }else{
    my @attributes = split ("; ", $attribute);
    foreach (@attributes) {
      my @item                     = split (" ", $_);
      my $item_name                = $item[0];
      my $item_value               = $item[1];
      $item_value =~ s/(^"|"$)//g;
      $hashAttribute{$item_name} = $item_value;
    }
  }
  my $txID     = $hashAttribute{'transcript_id'};
  my $txType   = $hashAttribute{'transcript_type'};
  my $txName   = $hashAttribute{'transcript_name'};
  my $geneID   = $hashAttribute{'gene_id'};
  my $geneType = $hashAttribute{'gene_type'};
  my $geneName = $hashAttribute{'gene_name'};
  $txID        =~ s/\.\d+$//;# remove the version of transcript_id
  if (!defined($allFlag)) {
    next if (!exists($hashExtGeneType{$txType}));
    if ($txType eq 'protein_coding') {
      next if ( ($hashAttribute{'tag'} eq 'cds_start_NF') or ($hashAttribute{'tag'} eq 'cds_end_NF'));
    }
  }
  $hashTxIdInfo{$txID} = [$geneID, $geneName, $geneType, $txName, $txType];
}
close(IN);

open (IN, "<$bedFile") or die "cannot open file: $bedFile $!";
open (OUT, ">$output") or die "cannot open file: $output $!";
while (my $line=<IN>) {
  chomp($line);
  my @lineContents = split ("\t", $line);
  my $txID  = $lineContents[3];
  my $txIDNoVersion;
  ($txIDNoVersion = $txID) =~ s/\.\d+$//;
  next if (! exists($hashTxIdInfo{$txIDNoVersion}));
  my ($geneID, $geneName, $geneType, $txName, $txType) = @{$hashTxIdInfo{$txIDNoVersion}};
  $lineContents[3] = join ":", ($geneID, $geneName, $geneType, $txID, $txName, $txType);
  print OUT join ("\t", @lineContents), "\n";
}
close(IN);
close(OUT);

`sort -t \$'\t' -k 1,1V -k 2,2n $output -o $output`;

################# Abbreviations for this script #################
#
# txID = transcript id
# txType = transcript type
# txName = transcript name
# txType = transcript type
# txIDNoVersion = transcript id without version (eg. ENST00000618323)
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for extracting specific gene types from bed12 based on gtf(gff3).

=head1 CAUTIONS

If "-g" is false, script will set as 'hg19';

=head1 SYNOPSIS

extTypeAnnotation.pl [options] -t mRNA -b [bed] -gtf [gtf] -o [file]

 Options:
    -b | --bed         The bed file
    -g | --genome      The genome build
    -o | --output      The output file name
    -t | --type        The the extracted geneType (mRNA, lncRNA, sRNA, pseudogene, all)
    -h | --help        Brief help message
    -gtf               The annotation file in gtf format
    -gff3              The annotation file in gff3 format
    -verbose           Print more error message
    -man               Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

extTypeAnnotation.pl [options] -t mRNA -b [bed] -gtf [gtf] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will extract specific gene types from bed12 based on gtf(gff3).

=cut