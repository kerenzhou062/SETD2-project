#!/usr/bin/perl
use strict;
use warnings;
die "perl $0 shCont_input_rep1.genes.results /data/zhoukr/ccle/CCLE_RNAseq_081117.rpkm.gct Hela_rep1 correlation.txt" if (scalar(@ARGV)==0);

my $RNASeq     = $ARGV[0];
my $ccleRNASeq = $ARGV[1];
my $RNASeqName = $ARGV[2];
my $output     = $ARGV[3];
my %hashRNASeqExp;
open (IN, "<$RNASeq") or die "cannot open file:$!";
<IN>;
while (my $line=<IN>) {
  chomp($line);
  my @lineContents = split("\t", $line);
  my $geneID       = $lineContents[0];
  $geneID          =~ s/\.\w+$//;
  my $exp          = $lineContents[6];
  next if ($exp < 1);
  $hashRNASeqExp{$geneID} = $exp;
}
close(IN);

my %hashCCLEExp;
my %hashCCLEGene;
open (IN, "<$ccleRNASeq") or die "cannot open file:$!";
<IN>;
<IN>;
my $header = <IN>;
chomp($header);
my @cellLineInfo = split("\t", $header);
splice(@cellLineInfo,0,2);
my $celllineNum = scalar(@cellLineInfo);
while (my $line=<IN>) {
  chomp($line);
  my @lineContents = split("\t", $line);
  my $geneID       = $lineContents[0];
  $geneID          =~ s/\.\w+$//;
  $hashCCLEGene{$geneID}++;
  splice(@lineContents,0,2);
  for (my $i = 0; $i < $celllineNum; $i++){
    $hashCCLEExp{$i}->{$geneID} = $lineContents[$i];
  }
}
close(IN);

my @geneIDs;
foreach my $geneID(sort %hashRNASeqExp) {
  if (exists($hashCCLEGene{$geneID})) {
    push (@geneIDs, $geneID);
  }
}

my @RNASeqExp = map{$hashRNASeqExp{$_}} @geneIDs;

my $tempFile = "./temp.tmp";
`rm -f $tempFile`;
open (OUT, ">$tempFile") or die "cannot open temp file :$!";
print OUT "$RNASeqName\t$RNASeqName\t";
print OUT join ("\t", @RNASeqExp);
print OUT "\n";

for (my $i = 0; $i < $celllineNum; $i++){
  my $celllineName = $cellLineInfo[$i];
  my @Exps = map{$hashCCLEExp{$i}->{$_}} @geneIDs;
  print OUT "$celllineName\t$celllineName\t";
  print OUT join ("\t", @Exps);
  print OUT "\n";
}

`coExpressionFDR -p 1 -q 1 -n 1 $tempFile > $output`;

`rm -f $tempFile`;
