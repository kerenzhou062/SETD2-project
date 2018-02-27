#!/usr/bin/perl
use strict;
use warnings;
# z-score normalization, X* = X - mean / standard_deviation
die "perl $0 <HepG2_RNA-seq_RSEM_shCont_input.txt> rep1/rep2/rep3/mean <HepG2_shCont_gene_histone.txt> <HepG2_shCont_gene_m6a.txt> <join_exp-histone-m6A>" if (scalar(@ARGV)==0);

my $rsemFile   = $ARGV[0];
my $sample = $ARGV[1];
my $histoneFile = $ARGV[2];
my $m6aFile = $ARGV[3];
my $outputPrefix = $ARGV[4];

my $expIndex;
if ($sample eq 'rep1') {
  $expIndex = 2;
}elsif ($sample eq 'rep2') {
  $expIndex = 3;
}elsif ($sample eq 'rep3') {
  $expIndex = 4;
}elsif ($sample eq 'mean') {
  $expIndex = 5;
}else{
  print $sample, "\n";
  die "Should be 'rep1/rep2/rep3/mean'!\n";
}

my %hashGeneName;
my %hashGeneExp;
my @geneExp;
open (EXP, "<$rsemFile") or die "cannot open file:$!\n";
<EXP>;
while(my $line=<EXP>){
  chomp($line);
  my @lineCont = split("\t", $line);
  my $exp = $lineCont[$expIndex];
  next if ($exp == 0);
  push(@geneExp, $exp);
  my $geneID = $lineCont[0];
  my $geneName = $lineCont[1];
  $hashGeneExp{$geneID}->{'exp'} = $exp;
  $hashGeneName{$geneID} = $geneName;
}
close(EXP);

my $lowCutoff = &percentile(\@geneExp, 0.4);
my $highCutoff = &percentile(\@geneExp, 0.7);
foreach my $geneID (keys %hashGeneExp) {
  my $exp = $hashGeneExp{$geneID}->{'exp'};
  if ($exp >= $highCutoff) {
    $hashGeneExp{$geneID}->{'class'} = 'high';
  }elsif ($exp < $highCutoff and $exp > $lowCutoff) {
    $hashGeneExp{$geneID}->{'class'} = 'medium';
  }else{
    $hashGeneExp{$geneID}->{'class'} = 'low';
  }
}

my %hashGeneHistone;
my @histoneFoldEnri;
open (HISTONE, "<$histoneFile") or die "cannot open file:$!\n";
while (my $line=<HISTONE>) {
  chomp($line);
  my ($geneID, $foldEnri) = split("\t", $line);
  $hashGeneHistone{$geneID}->{'foldEnri'} = $foldEnri;
  push(@histoneFoldEnri, $foldEnri);
}

my $histoneMean = &Mean(\@histoneFoldEnri);
my $histoneSD = &standardDeviation(\@histoneFoldEnri);

foreach my $geneID (keys %hashGeneHistone) {
  my $foldEnri = $hashGeneHistone{$geneID}->{'foldEnri'};
  my $zscore = &Zscore($foldEnri, $histoneMean, $histoneSD);
  $hashGeneHistone{$geneID}->{'zscore'} = $zscore
}

my %hashGeneM6a;
my @m6aFoldEnri;
open (M6A, "<$m6aFile") or die "cannot open file:$!\n";
while (my $line=<M6A>) {
  chomp($line);
  my ($geneID, $foldEnri) = split("\t", $line);
  $hashGeneM6a{$geneID}->{'foldEnri'} = $foldEnri;
  push(@m6aFoldEnri, $foldEnri);
}

my $m6aMean = &Mean(\@m6aFoldEnri);
my $m6aSD = &standardDeviation(\@m6aFoldEnri);

foreach my $geneID (keys %hashGeneM6a) {
  my $foldEnri = $hashGeneM6a{$geneID}->{'foldEnri'};
  my $zscore = &Zscore($foldEnri, $m6aMean, $m6aSD);
  $hashGeneM6a{$geneID}->{'zscore'} = $zscore
}

my @shareGeneIDs;
foreach my $geneID (sort keys %hashGeneName) {
  if (exists($hashGeneHistone{$geneID}) and exists($hashGeneM6a{$geneID})){
    push (@shareGeneIDs, $geneID);
  }
}

my $foldEnriOutput = $outputPrefix .".txt";
my $foldEnriLogOutput = $outputPrefix ."_log2.txt";
my $foldEnriLogReshapgeOutput = $outputPrefix ."_log2_reshape.txt";
my $zscoreOutput = $outputPrefix ."_". "zscore.txt";
my $header = "geneID\tgeneName\texpClass\tm6A\tH3K36me3";

open(FOLDENRI, ">$foldEnriOutput") or die "cannot open file:$!\n";
open(FOLDENRILOG, ">$foldEnriLogOutput") or die "cannot open file:$!\n";
open(ZSCORE, ">$zscoreOutput") or die "cannot open file:$!\n";
open(FOLDENRILOGRESHAPE, ">$foldEnriLogReshapgeOutput") or die "cannot open file:$!\n";

print FOLDENRI "$header\n";
print FOLDENRILOG "$header\n";
print ZSCORE "$header\n";
print FOLDENRILOGRESHAPE "expClass\tType\tFoldLog2\n";

foreach my $geneID (@shareGeneIDs) {
  my $geneName = $hashGeneName{$geneID};
  my $expClass = $hashGeneExp{$geneID}->{'class'};
  my $histoneFoldEnri = $hashGeneHistone{$geneID}->{'foldEnri'};
  my $histoneFoldEnriLog = log($histoneFoldEnri) / log(2);
  my $m6aFoldEnri = $hashGeneM6a{$geneID}->{'foldEnri'};
  my $m6aFoldEnriLog = log($m6aFoldEnri) / log(2);
  my $histoneZscore = $hashGeneHistone{$geneID}->{'zscore'};
  my $m6aZscore = $hashGeneM6a{$geneID}->{'zscore'};
  print FOLDENRI join("\t", ($geneID, $geneName, $expClass, $m6aFoldEnri, $histoneFoldEnri)), "\n";
  print FOLDENRILOG join("\t", ($geneID, $geneName, $expClass, $m6aFoldEnriLog, $histoneFoldEnriLog)), "\n";
  print ZSCORE join("\t", ($geneID, $geneName, $expClass, $m6aZscore, $histoneZscore)), "\n";
  print FOLDENRILOGRESHAPE "$expClass\tm6A\t$m6aFoldEnriLog\n";
  print FOLDENRILOGRESHAPE "$expClass\tH3K36me3\t$histoneFoldEnriLog\n";
}

`sort -t \$'\t' -k1,1 -k2,2 $foldEnriLogReshapgeOutput -o $foldEnriLogReshapgeOutput`;






sub Mean{
  my ($data, undef) = @_; # $data
  my $sum = 0;
  for (my $i = 0; $i < scalar(@$data); $i++) {
    $sum += $data->[$i];
  }
  return ($sum / scalar(@$data));
}

sub standardDeviation{
  my ($data, undef) = @_; # $data
  my $sum = 0;
  my $mean = &Mean($data);
  for (my $i = 0; $i < scalar(@$data); $i++) {
    $sum += ($data->[$i] - $mean) ** 2;
  }
  return (sqrt($sum / scalar(@$data)));
}

sub Zscore{
  my ($value, $mean, $standardDeviation) = @_;
  my $zscore = ($value - $mean) / $standardDeviation;
  return $zscore;
}

sub percentile {
  my ($data, $percentile) = @_; # $data, 0.75
  my @sortData=sort {$a <=> $b} @$data;
  my $n = scalar(@sortData);
  my $percenNum;
  if ($n >= 4) {
    my $Ln = $n * $percentile;
    if ($Ln == int($Ln)) {
      $percenNum = ($sortData[$Ln -1] + $sortData[$Ln]) / 2;
    }else{
      $percenNum = $sortData[int($Ln)];
    }
  }elsif ($n == 3) {
    if ($percentile > 0.5) {
      $percenNum = ($sortData[2] - $sortData[1]) * ($percentile - 0.5) * 2 + $sortData[1];
    }elsif($percentile == 0.5){
      $percenNum = $sortData[1];
    }else{
      $percenNum = ($sortData[1] - $sortData[0]) * $percentile * 2 + $sortData[0];
    }
  }elsif ($n == 2) {
    $percenNum = ($sortData[1] - $sortData[0]) * $percentile + $sortData[0];
  }else{
    $percenNum = $sortData[0];
  }
  return $percenNum + 0;
}
