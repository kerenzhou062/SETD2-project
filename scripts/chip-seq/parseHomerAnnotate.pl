#!/usr/bin/perl
use strict;
use warnings;
die "perl $0 <transcript.bed6> <homer_annotatePeaks.txt> <gene_foldEnrichment.txt>" if (scalar(@ARGV)==0);

my $annotation = $ARGV[0];
my $homerAnno = $ARGV[1];
my $output = $ARGV[2];

my %hashTx;
my %hashGene;
open (IN, "<$annotation") or die "cannot open file:$!";
while (my $line=<IN>) {
  chomp($line);
  my @lineCont = split("\t", $line);
  my $chr = $lineCont[0];
  my $start = $lineCont[1];
  my $end = $lineCont[2];
  my $strand = $lineCont[5];
  my @txInfo = split(/\|/, $lineCont[3]);
  my $txID = $txInfo[0];
  my $geneID = $txInfo[-1];
  $hashTx{$txID}->{'start'} = $start;
  $hashTx{$txID}->{'end'} = $end;
  $hashTx{$txID}->{'strand'} = $strand;
  $hashTx{$txID}->{'gene'} = $geneID;
}
close(IN);

open (IN, "<$homerAnno") or die "cannot open file:$!";
<IN>;
while (my $line=<IN>) {
  chomp($line);
  my @lineCont = split("\t", $line);
  my @peakCoor = ($lineCont[2], $lineCont[3]);
  my $foldEnri = $lineCont[5];
  my $txID = $lineCont[10];
  my $geneID = $hashTx{$txID}->{'gene'};

  my @txCoor = ($hashTx{$txID}->{'start'}, $hashTx{$txID}->{'end'});
  my $overLenth = &CoorOverLen(\@peakCoor, \@txCoor);
  next if ($overLenth < -100);

  if (exists($hashGene{$geneID})){
    my $eFoldEnri = $hashGene{$geneID}->{'foldEnri'};
    if ($foldEnri > $eFoldEnri ) {
      $hashGene{$geneID}->{'foldEnri'} = $foldEnri;
      $hashGene{$geneID}->{'len'} = $overLenth;
    }else{
      if ($foldEnri == $eFoldEnri){
        my $eOverLenth = $hashGene{$geneID}->{'len'};
        if ($overLenth > $eOverLenth) {
          $hashGene{$geneID}->{'len'} = $overLenth;
        }
      }
    }
  }else{
    $hashGene{$geneID}->{'foldEnri'} = $foldEnri;
    $hashGene{$geneID}->{'len'} = $overLenth;
  }
}
close(IN);

open (OUT, ">$output") or die "cannot open file:$!";
foreach my $geneID (sort keys %hashGene){
  my $foldEnri = $hashGene{$geneID}->{'foldEnri'};
  my $logFoldEnri = log($foldEnri) / log(2);
  print OUT "$geneID\t$foldEnri\t$logFoldEnri\n";
}
close(OUT);

# calculated overlapped lengths of two peaks
sub CoorOverLen{
  my ($coor1, $coor2) = @_;
  my @coors = @$coor1;
  push(@coors, @$coor2);
  @coors = sort @coors;
  my $overlappedLen = $coors[2] - $coors[1];

  if (@$coor1 ~~ @coors[0..1]) {
    return(-$overlappedLen);
  }elsif (@$coor1 ~~ @coors[2..3]) {
    return(-$overlappedLen);
  }else{
    return $overlappedLen;
  }
}

