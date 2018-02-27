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
my $features = '5utr,cds,3utr';
my $input;
my $bed12File;
my $bed6File;
my $output = './binResult.txt';
my $binType = 'constant';
my $binSize = 100;
my $distance = 0;
my $type = 'percentage';
my $peakMethod = 'center';
my $totalPeak;
my $mappedFlag;
my $smooth = 'move';
my $sort;
my $strandFlag;
my $span = 5;
my $score = 1;
my $Verbose;

GetOptions (
  "d|distance=s{1,1}"                =>\$distance,
  "f|feature=s{1,1}"                 =>\$features,
  "i|input=s{1,1}"                   =>\$input,
  "o|output=s{1,1}"                  =>\$output,
  "s|size=s{1,1}"                    =>\$binSize,
  "t|type=s{1,1}"                    =>\$type,
  "bed12=s{1,1}"                     =>\$bed12File,
  "bed6=s{1,1}"                      =>\$bed6File,
  "binType=s{1,1}"                   =>\$binType,
  "method=s{1,1}"                    =>\$peakMethod,
  "total=s{1,1}"                     =>\$totalPeak,
  "mapped"                           =>\$mappedFlag,
  "smooth=s{1,1}"                    =>\$smooth,
  "span=s{1,1}"                      =>\$span,
  "sort"                             =>\$sort,
  "strand"                           =>\$strandFlag,
  "score=s{1,1}"                     =>\$score,
  "Verbose"                          =>\$Verbose,
  "h|help"                           =>\$help,
  "man"                              =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-Verbose => 2) if $man;

my @requeres = ('bedtools');
&Requere(@requeres);

if (defined($Verbose)) {
  if (defined($bed12File) and defined($bed6File)) {
    print STDERR "-bed6 will be ignore!\n";
  }
}

if ( ($type ne 'percentage') and ($type ne 'count') ) {
  $type = 'percentage';
  &Verbose("Unkown -t. 'percentage' is used!\n");
}


my $outputPath  = &FileFolderPath($output);

###### generating temp bed from input bed file ######
my $bedType = &BedTypeJudge($input);
my $tempBed = $outputPath .'/'. rand(1000000) .'.temp.bed';
my $totalPeakNum;
if ($peakMethod eq 'center') {
  if ($bedType eq 'N') {
    print STDERR "--input is not a bed like file!\n";
    exit;
  }else{
    my $temp2 = "$tempBed.tmp";
    if ($bedType eq 'bed12' or $bedType eq 'bed12+') {
      $totalPeakNum = &bedToTemp($bedType, $input, $tempBed, $score, 0);
    }else{
      $totalPeakNum = &bedToTemp($bedType, $input, $temp2, $score, 1);
      &bedToCenter($temp2, $tempBed, $score);
      `rm -f $temp2`;
    }
  }
}elsif($peakMethod eq 'coverage'){
  if ($bedType eq 'N') {
    print STDERR "--input is not a bed like file!\n";
    exit;
  }else{
    $totalPeakNum = &bedToTemp($bedType, $input, $tempBed, $score, 1);
  }
}else{
  print STDERR "Unkown -method!\n";
  exit;
}

###### generating annotation feature hash ######
if (defined($bed12File)) {
  my @minusRequeres = ('bed12ToBed6.pl');
  &Requere(@minusRequeres);
  $bed6File = $outputPath .'/'. rand(1000000) .'.bed6.tmp';
  `bed12ToBed6.pl -b $bed12File -o $bed6File`;
}

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

###calculating feature size ####
my @features = split (",", $features);
my %hashFeatureBinSize;
if ($binType ne 'average'){
  if ( $binType ne 'constant' ) {
    $binSize = 100;
    &Verbose("Unkown -binType. '100' is used for --size!\n");
  }
  for (my $i = 0; $i < scalar(@features); $i++) {
    $hashFeatureBinSize{$features[$i]} = $binSize;
  }
}else{
  my %hashFeatureSize;
  foreach my $featureName (keys %hashBed6) {
    foreach my $feature (keys %{$hashBed6{$featureName}}) {
      my $eachSize = $hashBed6{$featureName}->{$feature}->{'size'};
      $hashFeatureSize{$feature}->{'count'}++;
      $hashFeatureSize{$feature}->{'sumSize'} += $eachSize;
    }
  }
  my $fullFeatureLength = 0;
  foreach my $feature (@features) {
    my $count = $hashFeatureSize{$feature}->{'count'};
    my $sumSize = $hashFeatureSize{$feature}->{'sumSize'};
    $hashFeatureSize{$feature}->{'average'} = $sumSize / $count;
    $fullFeatureLength += $hashFeatureSize{$feature}->{'average'};
  }

  my $totalBinSize = $binSize * scalar(@features);
  foreach my $feature (@features) {
    my $featureLength = $hashFeatureSize{$feature}->{'average'};
    my $featurePercent = $featureLength / $fullFeatureLength;
    my $eachBinSize = sprintf ("%.0f", $featurePercent * $totalBinSize);
    $hashFeatureBinSize{$feature} = $eachBinSize;
  }
}

###interseting with annotation bed6 ####
my $intersectTemp = $outputPath .'/'. rand(1000000) .'_intersect.tmp';
if (defined($sort)) {
  `sort -t \$'\t' -k 1,1 -k 2,2n $tempBed -o $tempBed`;
  `sort -t \$'\t' -k 1,1 -k 2,2n $bed6File -o $bed6File`;
}
if (defined($strandFlag)) {
  `bedtools intersect -wb -s -a $tempBed -b $bed6File > $intersectTemp`;
}else{
  `bedtools intersect -wb -a $tempBed -b $bed6File > $intersectTemp`;
}


###calculating bin coverage ####
my %hashIntersectPeak;
open(INTERSECT,"<$intersectTemp") or die "can't open $intersectTemp, $!\n";
while (my $line=<INTERSECT>) {
  my (undef, undef, undef, $peakID, undef) = split ("\t", $line, 5);
  $hashIntersectPeak{$peakID}++;
}
close (INTERSECT);

my $totalBinCoverage = 0;
my %hashBin;
open(INTERSECT,"<$intersectTemp") or die "can't open $intersectTemp, $!\n";
while (my $line=<INTERSECT>) {
  chomp($line);
  my @lineContents = split ("\t", $line);
  my $start = $lineContents[1];
  my $end = $lineContents[2];
  my $peakID = $lineContents[3];
  my $peakCoverage = $lineContents[4];
  my $featureStart = $lineContents[7];
  my $featureEnd = $lineContents[8];
  my ($featureName, $feature, $order) = split (/\|/, $lineContents[9]);
  my $featureStrand = $lineContents[11];

  if ($distance != 0) {
    $start = int($start + $end) / 2;
    $end = $start + 1;
  }

  my $peakIntersetTime = $hashIntersectPeak{$peakID};
  my $averageCoverage = $peakCoverage/ $peakIntersetTime;
  my $featureBinSize = $hashFeatureBinSize{$feature};
  my $featureTotalSize = $hashBed6{$featureName}->{$feature}->{'size'}; #1 base

  for (my $j = ($start + 1); $j <= $end; $j++){
    my $shiftRe = 0; # 1 base
    my $bin; # 1 base
    for (my $i = 1; $i <= $order; $i++) {
      if ($i == $order) {
        $shiftRe += $j - $featureStart;
      }else{
        my $featureSize = $hashBed6{$featureName}->{$feature}->{'order'}->{$i};
        $shiftRe += $featureSize;
      }
    }

    if ($featureStrand eq '+') {
      $bin = ceil( $shiftRe / $featureTotalSize * $featureBinSize );
      $hashBin{'original'}->{$feature}->{$bin} += $averageCoverage;
    }else{
      $bin = ceil( ($featureTotalSize - $shiftRe + 1)/ $featureTotalSize * $featureBinSize );
      $hashBin{'original'}->{$feature}->{$bin} += $averageCoverage;
    }
    $totalBinCoverage += $averageCoverage;
  }
}
close(INTERSECT);

if (defined($bed12File)) {
  `rm -f $tempBed $intersectTemp $bed6File`;
}else{
  `rm -f $tempBed $intersectTemp`;
}

my $totalNum;
if (defined($totalPeak) and looks_like_number($totalPeak)) {
  $totalNum = $totalPeak;
}else{
  if (defined($mappedFlag)) {
    $totalNum = $totalBinCoverage;
  }else{
    if ($peakMethod eq 'coverage') {
      $totalNum = $totalBinCoverage;
    }else{
      $totalNum = $totalPeakNum;
    }
  }
}


###Normalizing bin coverage ####

my $spanWith = 0;
foreach my $feature(@features) {
  my $binShift = $spanWith + 1;
  my $eachBinSize = $hashFeatureBinSize{$feature};
  for (my $i = $binShift; $i < ($binShift + $eachBinSize); $i++){
    my $bin = $i - $binShift + 1;
    if (!exists($hashBin{'original'}->{$feature}->{$bin})) {
      if ($type eq 'percentage') {
        $hashBin{'percentage'}->{$i} = 0;
      }else{
        $hashBin{'count'}->{$i} = 0;
      }
    }else{
      if ($type eq 'percentage') {
        $hashBin{'percentage'}->{$i} = ($hashBin{'original'}->{$feature}->{$bin} / $totalNum) * 100;
      }else{
        $hashBin{'count'}->{$i} = $hashBin{'original'}->{$feature}->{$bin};
      }
    }
  }
  $spanWith += $eachBinSize;
}

###smooth bin data ####
if ($smooth ne 'no' and $smooth ne 'average' and $smooth ne 'move') {
  $smooth = 'move';
  if (defined($Verbose)) {
    print STDERR "Unkown -smooth, 'move' will be assigned!\n";
  }
}

my @ouputValueList;
$spanWith = 1;
foreach my $feature (@features) {
  my @binValueList;
  my $eachBinSize = $hashFeatureBinSize{$feature};
  for (my $i = $spanWith; $i < ($spanWith + $eachBinSize); $i++) {
    if ($type eq 'percentage') {
      push(@binValueList, $hashBin{'percentage'}->{$i});
    }else{
      push(@binValueList, $hashBin{'count'}->{$i});
    }
  }
  if ($smooth eq 'average') {
    push( @ouputValueList, &SmoothAverage( \@binValueList, $span ) );
  }elsif($smooth eq 'move'){
    push( @ouputValueList, &SmoothMove( \@binValueList, $span ) );
  }else{
    push(@ouputValueList, @binValueList);
  }
  $spanWith += $eachBinSize;
}

open (OUT, ">$output") or die "can't open $output, $!\n";
$spanWith = 0;
foreach my $feature (@features) {
  my $eachBinSize = $hashFeatureBinSize{$feature};
  for (my $i = $spanWith; $i < ($spanWith + $eachBinSize); $i++) {
    print OUT join ("\t", ($feature, $i + 1, $ouputValueList[$i])), "\n";
  }
  $spanWith += $eachBinSize;
}



##############  Subroutines   ################
##############                ################
##############  Subroutines   ################
sub SmoothMove{
  my ($valueListRef, $span) = @_;
  $span = int($span);
  my @movingSmoothValList;
  for (my $i = 0; $i < scalar(@$valueListRef); $i++) {
    my $sum         = 0;
    my $movingValue = 0;
    if ($i < ($span - 1)) {
      for (my $j = 0; $j <= $i; $j++) {
        $sum += $$valueListRef[$j];
      }
      $movingValue = $sum / ($i + 1);
    }else{
      for (my $j = ($i - $span + 1); $j < ($i + 1); $j++){
        $sum += $$valueListRef[$j];
      }
      $movingValue = $sum / $span;
    }
    $movingSmoothValList[$i] = $movingValue;
  }
  return @movingSmoothValList;
}

sub SmoothAverage{
  my ($valueListRef, $span) = @_;
  $span = ( int($span / 2) == ($span / 2) ) ? int($span) + 1 : int($span); # make the span be odd
  my $hSpan = int($span / 2);
  my $listLength = scalar(@$valueListRef);
  my @SmoothValList;

  for (my $i = 0; $i < $listLength; $i++) {
    my $sum         = 0;
    my $averageValue = 0;
    if ($i < $hSpan) {
      my $tSpan = $i * 2 + 1;
      for (my $j = 0; $j < $tSpan; $j++) {
        $sum += $$valueListRef[$j];
      }
      $averageValue = $sum / $tSpan;
    }elsif ($i >= $hSpan and $i < ($listLength - $hSpan)) {
      for (my $j = ($i - $hSpan); $j < ($i + $hSpan + 1); $j++) {
        $sum += $$valueListRef[$j];
      }
      $averageValue = $sum / $span;
    }else{
      my $tSpan = ($listLength - $i - 1) * 2 + 1;
      for (my $j = ($i - $tSpan); $j < ($listLength - 1); $j++) {
        $sum += $$valueListRef[$j];
      }
      $averageValue = $sum / $tSpan;
    }
    $SmoothValList[$i] = $averageValue;
  }
  return @SmoothValList;
}

sub Verbose{
  my $warningText = shift;
  if (defined($Verbose)) {
    print STDERR $warningText, "\n";
  }
}

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
    if ($count > 0 and $count < 20) {
      my @lineContents = split ("\t", $line);
      $sum += scalar(@lineContents);
      $lineCount++;
    }else{
      if ($count >= 20) {
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

sub bedToTemp{
  my ($bedType, $bedFile, $tempFile, $bedScore, $convertFlag) = @_;
  my $totalPeak = 0;
  open (BED, "<$bedFile") or die "Cannot open file: $bedFile, $!\n";
  open (TEMP, ">$tempFile") or die "Cannot open file: $tempFile, $!\n";
  if ($convertFlag == 0 && $peakMethod ne 'coverage') {
    while (my $line=<BED>) {
      chomp($line);
      my @lineContents = split("\t", $line);
      if (defined($distance)) {
        my $center = int( ($lineContents[1] + $lineContents[2]) / 2 );
        $lineContents[1] = ($center - $distance) < 0 ? 0 : ($center - $distance);
        $lineContents[2] = $center + $distance + 1;
      }
      $lineContents[3] = 'peak='. ($totalPeak + 1);
      if ($bedScore eq 'NA') {
        print TEMP join ("\t", @lineContents[0..5]), "\n";
      }else{
        if (looks_like_number($bedScore)) {
          $lineContents[4] = $bedScore;
        }else{
          $lineContents[4] = 1;
        }
      }
      print TEMP join ("\t", @lineContents[0..5]), "\n";
      $totalPeak++;
    }
  }else{
    if ($bedType eq 'bed12' or $bedType eq 'bed12+') {
      my $temp2 = "$tempFile.tmp";
      &Requere('bed12ToBed6.pl');
      `bed12ToBed6.pl -value original --repeat -average -b $input -o $temp2`;
      open (TEMP2, "<$temp2") or die "Cannot open file: $temp2, $!\n";
      while (my $line=<TEMP2>) {
        my @lineContents = split ("\t", $line);
        $lineContents[3] = "peak=". ($totalPeak + 1);
        if ($bedScore ne 'NA') {
          if (looks_like_number($bedScore)) {
            $lineContents[4] = $bedScore;
          }else{
            $lineContents[4] = 1;
          }
        }
        print TEMP join ("\t", @lineContents);
        $totalPeak++;
      }
      close (TEMP2);
      if ($bedScore eq 'NA') {
        ($totalPeak, undef) = split (/\s/, `wc -l $input`);
      }
      `rm -f $temp2`;
    }elsif($bedType eq 'bed3' or $bedType eq 'bed3+'){
      while (my $line=<BED>) {
        chomp($line);
        my @lineContents = split("\t", $line);
        if (defined($distance)) {
          my $center = int( ($lineContents[1] + $lineContents[2]) / 2 );
          $lineContents[1] = ($center - $distance) < 0 ? 0 : ($center - $distance);
          $lineContents[2] = $center + $distance + 1;
        }
        $lineContents[3] = 'peak='. ($totalPeak + 1);
        if (looks_like_number($bedScore)) {
          $lineContents[4] = $bedScore;
        }else{
          $lineContents[4] = 1;
        }
        $lineContents[5] = '+';
        print TEMP join ("\t", @lineContents), "\n";
        $totalPeak++;
      }
    }else{
      while (my $line=<BED>) {
        chomp($line);
        my @lineContents = split("\t", $line);
        if (defined($distance)) {
          my $center = int( ($lineContents[1] + $lineContents[2]) / 2 );
          $lineContents[1] = ($center - $distance) < 0 ? 0 : ($center - $distance);
          $lineContents[2] = $center + $distance + 1;
        }
        $lineContents[3] = 'peak='. ($totalPeak + 1);
        if ($bedScore ne 'NA') {
          if (looks_like_number($bedScore)) {
            $lineContents[4] = $bedScore;
          }else{
            $lineContents[4] = 1;
          }
        }
        print TEMP join ("\t", @lineContents[0..5]), "\n";
        $totalPeak++;
      }
    }
  }
  close(BED);
  close(TEMP);
  return $totalPeak;
}

sub bedToCenter{
  my ($bed6File, $centerFile, $bedScore) = @_;
  open (BED, "<$bed6File") or die "Cannot open file: $bed6File, $!\n";
  open (CENTER, ">$centerFile") or die "Cannot open file: $centerFile, $!\n";
  my $peakNumber = 1;
  while (my $line=<BED>) {
    chomp($line);
    my @lineContents = split("\t", $line);
    my $pcStart      = int( ($lineContents[1] + $lineContents[2]) / 2 );
    my $pcEnd        = $pcStart + 1;
    $lineContents[3] = "peak=". $peakNumber;
    if ($bedScore ne 'NA') {
      if (looks_like_number($bedScore)) {
          $lineContents[4] = $bedScore;
        }else{
          $lineContents[4] = 1;
        }
    }
    print CENTER join ("\t", ($lineContents[0], $pcStart, $pcEnd, @lineContents[3..5])), "\n";
    $peakNumber++;
  }
  close(BED);
  close(CENTER);
}

################# Abbreviations for this script #################
#
# bed6File (4thCol) = $geneID:$geneName:$geneType:$txID:$txName:$txType|cds|1
# pcEnd             = end position of a peak center
# shiftRe         = length shifts from the begining of a transcript region (5utr, cds, or 3utr)
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for computing the distribution patterns on transcripts
Requere: bedtools, bedCoverage.pl, bed12ToBed6.pl

=head1 CAUTIONS

Defaults: -smooth=(move), --size=(100), --distance=(undefined) -sort=(false), -span=(5), -strand=(false).

=head1 SYNOPSIS

bedBinDistribution.pl [options] --input [file] -bed12 [file] -o [file]

 Options:
    -bed12                The annotation bed12 file
    -bed6                 The annotation bed6 file
    -binType              bin type size ('constant'[default] or 'average')
    -mapped               Flag for using mapped peaks instead of total peaks
    -method               Method for calculating bin coverage ('center'[default] or 'coverage')
    -smooth               The smooth method ('no', 'average', or 'move'[default])
    -total                The total peak number for 'percentage' (inactivate by default)
    -sort                 Sort bed before using bedtools
    -span                 The span for smooth (int[default:5])
    -strand               Only report overlap on the same strand
    -score                The score value for bed peak ('int' [default=1] or 'NA')
    -d | --distance       Distance extended aroud peak center (works under -method 'coverage')
    -f | --feature        The bin features (5utr,cds,3utr[default])
    -i | --input          The input bed file
    -o | --output         The The out put file
    -s | --size           The bin size ('int'[default:100])
    -t | --type           The output type ('percentage'[default] or 'count')
    -h | --help           Brief help message
    -Verbose              Print more warnings
    -man                  Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

bedBinDistribution.pl [options] --input [file] -bed12 [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will compute the distribution patterns on transcripts.

=cut
