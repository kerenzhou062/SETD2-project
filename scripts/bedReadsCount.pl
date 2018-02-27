#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
###develop by K.R.Chow, designed for computing read counts of bed.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $input;
my @bamFiles;
my $rpmFlag;
my $output = './binResult.txt';

GetOptions (
  "i|input=s{1,1}"                   =>\$input,
  "bam=s{1,}"                        =>\@bamFiles,
  "rpm"                              =>\$rpmFlag,
  "o|output=s{1,1}"                  =>\$output,
  "h|help"                           =>\$help,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-Verbose => 2) if $man;

my @requeres = ('bedtools', 'featureCounts');
&Requere(@requeres);

my $outputPath  = &FileFolderPath($output);

###### generating temp bed from input bed file ######
my $bedType = &BedTypeJudge($input);
my $tempAnno = $outputPath .'/'. rand(1000000) .'.anno.tmp';
my %hashIdLine;
if ($bedType eq 'bed12') {
  %hashIdLine = &bedToGtf($input, $tempAnno);
}else{
  %hashIdLine = &bedToSaf($input, $tempAnno);
}
my $bams = join(" ", @bamFiles);

my $temp2 = $outputPath .'/'. rand(1000000) .'.featureCount.tmp';

if ($bedType eq 'bed12') {
  `featureCounts -T 10 -t exon -g transcript_id -a $tempAnno -F GTF -s 2 -o $temp2 $bams > /dev/null 2>&1`;
}

my $totalMappedReads;
open (SUMMARY, "<$temp2.summary") or die "cannot open file!";
<SUMMARY>;
while(my $line=<SUMMARY>){
  chomp($line);
  my @lineContents = split("\t", $line);
  for (my $i = 1; $i < scalar(@lineContents); $i++) {
    $totalMappedReads  += $lineContents[$i];
  }
}

my %hashIdExp;
open(COUNT, "<$temp2") or die "cannot open file:$temp2";
<COUNT>;
<COUNT>;
while (my $line=<COUNT>) {
  chomp($line);
  my @lineContents = split("\t", $line);
  my $id = $lineContents[0];
  my @readCounts = @lineContents[6..$#lineContents];
  my $sum = 0;
  my $average = 0;
  for (my $i = 0; $i < 3; $i++) {
    $sum += $readCounts[$i];
  }
  for (my $i = 3; $i < 6; $i++) {
    $sum -= $readCounts[$i];
  }
  my $expression;
  if (defined($rpmFlag)) {
    $expression = $sum / $totalMappedReads * 1e6;
  }else{
    $expression = $sum / 3;
  }
  $hashIdExp{$id} = $expression;
}

my $count = 1;
open(IN, "<$input") or die "cannot open file:$input";
open(OUT, ">$output") or die "cannot open file:$output";
while(my $line=<IN>) {
  chomp($line);
  my @lineContents = split("\t", $line);
  my $geneID = $lineContents[3];
  my $txID = $geneID . "-". $count;
  my $readCount = $hashIdExp{$txID};
  $lineContents[4] = $readCount;
  print OUT join("\t", @lineContents), "\n";
  $count++;
}
close(IN);
close(OUT);
`rm -f $tempAnno $temp2 $temp2.summary`;

##############  Subroutines   ################
##############                ################
##############  Subroutines   ################

sub FileFolderPath{
  my ($file,undef) = @_;
  my @temp = split (/\/|\\/, $file);
  if (scalar(@temp) == 1) {
    return "./";
  }else{
    splice (@temp, -1 , 1);
    return join ("/", @temp);
  }
}

sub Requere {
  my @requeres = @_;
  my $requereNum = 0;
  my $errorRequere;
  my %hashRequere;
  for my $path ( split /:/, $ENV{PATH} ) {
    foreach my $requere (@requeres) {
      if ( -f "$path/$requere" && -x _ ) {
        $requereNum++;
        $hashRequere{$requere}++;
      }
    }
  }
  if ($requereNum < scalar(@requeres)) {
    foreach my $requere (@requeres) {
      if (!exists($hashRequere{$requere})) {
        $errorRequere .= "$requere,";
      }
    }
    print STDERR "Requered softwares do not exists ($errorRequere)!\n";
    exit;
  }
}

sub BedTypeJudge{
  my ($bedFile, undef) = @_;
  open (BED, "<$bedFile") or die "Cannot open bed file: $bedFile, $!\n";
  my $aveCount;
  my $sum = 0;
  my $count = 1;
  my $lineCount = 0;
  while (my $line=<BED>) {
    next if ($line =~ /^#/);
    if ($count > 0 and $count < 30) {
      my @lineContents = split ("\t", $line);
      $sum += scalar(@lineContents);
      $lineCount++;
    }else{
      if ($count >= 30) {
        last;
      }
    }
    $count++;
  }
  close(BED);

  $aveCount = $sum / $lineCount;

  if ($aveCount == int($aveCount)) {
    if ($aveCount < 3) {
      return 'N';
    }elsif($aveCount == 3) {
      return 'bed3';
    }elsif ($aveCount < 6) {
      return 'bed3+';
    }elsif ($aveCount == 6) {
      return 'bed6';
    }elsif ($aveCount < 12) {
      return 'bed6+';
    }elsif ($aveCount == 12) {
      return 'bed12';
    }else{
      return 'bed12+';
    }
  }else{
    return 'N';
  }
}

sub bedToGtf{
  my ($bedFile, $tempGtf) = @_;
  open (BED, "<$bedFile") or die "Cannot open file: $bedFile, $!\n";
  open (GTF, ">$tempGtf") or die "Cannot open file: $tempGtf, $!\n";
  my %hashTemp;
  my $count = 1;
  while (my $line=<BED>) {
    chomp($line);
    my @lineContents = split("\t", $line);
    my $chr = $lineContents[0];
    my $start = $lineContents[1];
    my $end = $lineContents[2];
    my $geneID = $lineContents[3];
    my $strand = $lineContents[5];
    chop($lineContents[-2]);
    my @lengths = split(",", $lineContents[-2]);
    my @starts = split(",", $lineContents[-1]);
    my $txID = $geneID ."-". $count;
    $hashTemp{$txID} = $count;
    print GTF join("\t", ($chr, 'ExomePeak', 'transcript', $start + 1, $end, '.', $strand, '.', "transcript_id \"$txID\";")), "\n";
    for (my $i = 0; $i < scalar(@starts); $i++) {
      my $exonNum = $i + 1;
      my $exonStart = $start + $starts[$i];
      my $exonEnd = $exonStart + $lengths[$i];
      print GTF join("\t", ($chr, 'ExomePeak', 'exon', $exonStart + 1, $exonEnd, '.', $strand, '.', "transcript_id \"$txID\"; exon_number $exonNum")), "\n";
    }
    $count++;
  }
  close(BED);
  close(GTF);
  return %hashTemp;
}

################# Abbreviations for this script #################
#
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for computing read counts of bed.
Requere: featureCounts

=head1 CAUTIONS

Defaults: -rpm=(false)

=head1 SYNOPSIS

bedReadsCount.pl [options] --input [file] -bam [file file] -o [file]

 Options:
    -bam                  The input bams
    -rpm                  Calculate rpm for each position(false by default)
    -i | --input          The input bed file
    -o | --output         The The out put file
    -h | --help           Brief help message
    -man                  Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

bedReadsCount.pl [options] --input [file] -bed12 [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will compute the distribution patterns on transcripts.

=cut
