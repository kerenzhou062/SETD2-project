#!/usr/bin/perl
use strict;
use warnings;
### all the FC are log2(FC)
die "perl $0 <HepG2_RNA-seq_RSEM_shCont_input.txt> rep1/rep2/rep3/mean <HepG2_shCont_gene_histone.txt/HepG2_shCont_gene_m6a.txt> <join_exp-histone/join_exp-m6A>" if (scalar(@ARGV)==0);

my $rsemFile   = $ARGV[0];
my $sample = $ARGV[1];
my $modiFile = $ARGV[2];
my $outputPrefix = $ARGV[3];

my @samples = ("shCont", "shSETD2", "shMETTL14", "shMETTL3", "shWTAP");
my $sampleNum = 5;
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

my $lowCutoff = &percentile(\@geneExp, 0.25);
my $highCutoff = &percentile(\@geneExp, 0.75);
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

my @shareGeneIDs;
my %hashGeneModi;

open (MODI, "<$modiFile") or die "cannot open file:$!\n";
while (my $line=<MODI>) {
  chomp($line);
  my @lineCont = split("\t", $line);
  my $geneID = $lineCont[0];
  next if (! exists($hashGeneExp{$geneID}));
  push(@shareGeneIDs, $geneID);
  my @fondEnri = @lineCont[1..$#lineCont];
  $hashGeneModi{$geneID} = \@fondEnri;
}

my $foldEnriOutput = $outputPrefix .".txt";
my $foldEnriLogOutput = $outputPrefix ."_log2.txt";
my $header = "geneID\tgeneName\tfpkm\texpClass\t". join("\t", @samples);

open(FOLDENRILOG, ">$foldEnriLogOutput") or die "cannot open file:$!\n";

print FOLDENRILOG "$header\n";

foreach my $geneID (@shareGeneIDs) {
  my $geneName = $hashGeneName{$geneID};
  my $exp = $hashGeneExp{$geneID}->{'exp'};
  my $expClass = $hashGeneExp{$geneID}->{'class'};
  my @modiFoldEnriLog = @{$hashGeneModi{$geneID}};
  print FOLDENRILOG join("\t", ($geneID, $geneName, $exp, $expClass, @modiFoldEnriLog)), "\n";
}



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
