#!/usr/bin/perl
use strict;
use warnings;
die "perl $0 /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/correlationStats" if (scalar(@ARGV)==0);

my $rsemFolder   = $ARGV[0];
my $outputFolder = $ARGV[1];
my $output = "$outputFolder/m6A-seq_sample_replicates_correlationStats.txt";
my @rsemFiles = split ("\n", `find $rsemFolder -type f -name "*.genes.results"`);
my %hashSampleGenExp;
foreach my $rsemFile(@rsemFiles) {
  my @contents = split (/\\|\//, $rsemFile);
  my ($target, $seqType, $rep, undef) = split (/_|\./, $contents[-1]);
  my $sample = "${target}_${seqType}_$rep";
  open (IN, "<$rsemFile") or die "cannot open file:$!";
  <IN>;
  while (my $line=<IN>) {
    chomp($line);
    my @lineContents = split("\t", $line);
    my $geneID       = $lineContents[0];
    my $exp          = $lineContents[5];
    next if ($exp < 1);
    $hashSampleGenExp{$sample}->{$geneID} = $exp;
  }
  close(IN);
}

my @samples = sort keys %hashSampleGenExp;
my $sampleNum = scalar(@samples);

my $tempFile = "$outputFolder/temp.tmp";
open (OUT, ">$output") or die "cannot open file:$!";
print OUT "target\ttype\trep\tpairedTarget\tpairedType\tpairRep\tgeneNum\tcoefficient-R\tp-val\n";

for (my $i = 0; $i < $sampleNum; $i++) {
  my $sample = $samples[$i];
  my ($target, $seqType, $rep) = split("_", $sample);
  my %hashRepsPair;
  my @pairs;
  my @allGeneIDs = sort keys %{$hashSampleGenExp{$sample}};
  for (my $j = ($i + 1); $j < $sampleNum; $j++) {
    my $pairSample = $samples[$j];
    my ($pairTarget, $pairSeqType, $pairRep) = split("_", $pairSample);
    my @overlapedGeneIDs;
    foreach my $geneID(@allGeneIDs) {
      next if (!exists($hashSampleGenExp{$pairSample}{$geneID}));
      push (@overlapedGeneIDs, $geneID);
    }
    my @repExp = map{$hashSampleGenExp{$sample}->{$_}} @overlapedGeneIDs;
    my @pairedExp = map{$hashSampleGenExp{$pairSample}->{$_}} @overlapedGeneIDs;
    open (TEMP, ">$tempFile") or die "cannot open file:$!";
    print TEMP join ("\t", ("${target}_${seqType}_${rep}", "${target}_${seqType}_${rep}", @repExp)), "\n";
    print TEMP join ("\t", ("${pairTarget}_${pairSeqType}_${pairRep}", "${pairTarget}_${pairSeqType}_${pairRep}", @pairedExp)), "\n";
    close(TEMP);
    my $result = (split ("\n", `coExpressionFDR -p 1 -q 1 -n 1 $tempFile`))[-1];
    my @contents = split ("\t", $result);
    @contents = ($target, $seqType, $rep, $pairTarget, $pairSeqType, $pairRep, $contents[5], $contents[6], $contents[7]);
    print OUT join("\t", @contents), "\n";
  }
}

close(OUT);
`rm -f $tempFile`;
