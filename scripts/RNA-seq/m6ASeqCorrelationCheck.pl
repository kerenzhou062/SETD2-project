#!/usr/bin/perl
use strict;
use warnings;
die "perl $0 input /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/RSEM /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/Hela/correlationStats" if (scalar(@ARGV)==0);

my $seqType      = $ARGV[0];
my $rsemFolder   = $ARGV[1];
my $outputFolder = $ARGV[2];
my $output = "$outputFolder/m6A-seq_${seqType}_correlationStats.txt";
my @rsemFiles = split ("\n", `find $rsemFolder -type f -name "*.genes.results" | grep "$seqType"`);
my %hashTargetRepFile;
foreach my $rsemFile(@rsemFiles) {
  my @contents = split (/\\|\//, $rsemFile);
  my ($target, undef, $rep, undef) = split (/_|\./, $contents[-1]);
  open (IN, "<$rsemFile") or die "cannot open file:$!";
  <IN>;
  while (my $line=<IN>) {
    chomp($line);
    my @lineContents = split("\t", $line);
    my $geneID       = $lineContents[0];
    my $exp          = $lineContents[5];
    next if ($exp < 1);
    $hashTargetRepFile{$target}->{$rep}->{$geneID} = $exp;
  }
  close(IN);
}
my @targets = sort keys %hashTargetRepFile;


my $tempFile = "$outputFolder/temp.tmp";
open (OUT, ">$output") or die "cannot open file:$!";
print OUT "sample\tpairedSample\tgeneNum\tcoefficient-R\tp-val\n";
foreach my $target(@targets) {
  my $repNum = keys %{$hashTargetRepFile{$target}};
  my %hashRepsPair;
  my @reps;
  for (my $i = 1; $i <= $repNum; $i++) {
    my $replicate = "rep".$i;
    push (@reps, $replicate);
    for (my $j = $i; $j <= $repNum; $j++){
      if (($repNum - $j) >= 1) {
        my $pairedRep = "rep".($j+1);
        if (exists($hashRepsPair{$replicate})) {
          push (@{$hashRepsPair{$replicate}}, $pairedRep);
        }else{
          $hashRepsPair{$replicate} = [];
          push (@{$hashRepsPair{$replicate}}, $pairedRep);
        }
      }
    }
  }
  @reps = sort keys %hashRepsPair;
  foreach my $rep(@reps) {
    my @repPairs = @{$hashRepsPair{$rep}};
    foreach my $repPair(@repPairs) {
      my @allGeneIDs = sort keys %{$hashTargetRepFile{$target}->{$rep}};
      my @overlapedGeneIDs;
      foreach my $geneID(@allGeneIDs) {
        next if (!exists($hashTargetRepFile{$target}->{$repPair}->{$geneID}));
        push (@overlapedGeneIDs, $geneID);
      }
      my @repExp = map{$hashTargetRepFile{$target}->{$rep}->{$_}} @overlapedGeneIDs;
      my @pairedExp = map{$hashTargetRepFile{$target}->{$repPair}->{$_}} @overlapedGeneIDs;
      open (TEMP, ">$tempFile") or die "cannot open file:$!";
      print TEMP join ("\t", ("${target}_${rep}", "${target}_${rep}", @repExp)), "\n";
      print TEMP join ("\t", ("${target}_${repPair}", "${target}_${repPair}", @pairedExp)), "\n";
      close(TEMP);
      my $result = (split ("\n", `coExpressionFDR -p 1 -q 1 -n 1 $tempFile`))[-1];
      my @contents = split ("\t", $result);
      @contents = ($contents[1], $contents[3], $contents[5], $contents[6], $contents[7]);
      print OUT join("\t", @contents), "\n";
    }
  }
}

close(OUT);
`rm -f $tempFile`;
