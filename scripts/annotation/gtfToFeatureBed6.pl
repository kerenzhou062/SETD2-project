#!/usr/bin/perl
use strict;
use warnings;

die "$0 <annotation.gtf> <transcript> <annotation.bed6>" if (scalar(@ARGV) == 0 or scalar(@ARGV) != 3);

my $gtf = $ARGV[0];
my $feature = $ARGV[1];
my $ouput = $ARGV[2];
open(IN, "<$gtf") or die "cannot open file: $!";
open(OUT, ">$ouput") or die "cannot open file: $!";
while (my $line=<IN>) {
  next if ($line =~ /^#/);
  chomp($line);
  my @linCont = split("\t", $line);
  my $type = $linCont[2];
  next if ($type ne $feature);
  my $chr = $linCont[0];
  next if ($chr =~ /^GL\d+/);
  my $start = $linCont[3] - 1;
  my $end = $linCont[4];
  my $strand = $linCont[6];
  my %hashAtt;
  my @attributes = split("; ", $linCont[-1]);
  foreach my $attribute (@attributes){
    my ($name, $value) = split(" ", $attribute);
    $value =~ s/"//g;
    $hashAtt{$name} = $value;
  }
  my $nameCol;
  if ($feature eq "gene"){
    $nameCol= join("|", ($hashAtt{"gene_id"}, $hashAtt{"gene_name"}, $hashAtt{"gene_type"}));
  }elsif ($feature eq "transcript"){
    $nameCol= join("|", ($hashAtt{"transcript_id"}, $hashAtt{"transcript_name"}, $hashAtt{"transcript_type"}, $hashAtt{"gene_id"}));
  }
  print OUT join("\t", ($chr, $start, $end, $nameCol, "0", $strand)), "\n";
}

`sort -t \$'\t' -k 1,1V -k 2,2n $ouput -o $ouput`;
