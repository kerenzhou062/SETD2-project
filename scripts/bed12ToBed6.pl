#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
###develop by K.R.Chow, designed for formatting bed12 to bed6.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $annotationFile;
my $bed12File;
my $output;
my $repeat;
my $average;
my $sort = "V";
my $value = 'span';
my $verbose;

GetOptions (
  "b|bed=s{1,1}"                          =>\$bed12File,
  "o|output=s{1,1}"                       =>\$output,
  "r|repeat"                              =>\$repeat,
  "average"                               =>\$average,
  "s|sort=s{1,1}"                         =>\$sort,
  "value=s{1,1}"                          =>\$value,
  "verbose"                               =>\$verbose,
  "h|help"                                =>\$help,
  "man"                                   =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

my $bedType = &BedTypeJudge($bed12File);
if ($bedType ne 'bed12' and $bedType ne 'bed12+') {
  print STDERR "The -b did not seem like bed12!\n";
  exit;
}

if ($value ne 'span' and $value ne 'original') {
  if (!looks_like_number($value)){
    print STDERR "-value should be 'span' or number!\n";
    exit;
  }
}


my %hashFeatureName;
open (BED12, "<$bed12File") or die "Cannot open file: $bed12File, $!\n";
open (OUT, ">$output") or die "Cannot open file: $output, $!\n";
while (my $line=<BED12>) {
  next if ($line =~ /^#/);
  chomp($line);
  my @lineContents    = split("\t", $line);
  my $chr             = $lineContents[0];
  my $start           = $lineContents[1];
  my $end             = $lineContents[2];
  my $name            = $lineContents[3];
  my $score          = $lineContents[4];
  my $strand          = $lineContents[5];
  my $cdsStart        = $lineContents[6];
  my $cdsEnd          = $lineContents[7];
  #my $blockCount     = $lineContents[9];
  my $blockLengthCol  = $lineContents[10]; $blockLengthCol =~ s/,$//;
  my $blockStartCol   = $lineContents[11]; $blockStartCol =~ s/,$//;
  my @blockLengthList = split(',', $blockLengthCol);
  my @blockStartList  = split(',', $blockStartCol);
  my $exonValue;
  if ($value ne 'span') {
    if ($value eq 'original') {
      $value = $score;
    }
    if (defined($average)) {
      $exonValue = $value / scalar(@blockLengthList);
    }else{
      $exonValue = $value;
    }
  }

  if (exists($hashFeatureName{$name})) {
    if (! defined($repeat)) {
      next;
    }
  }

  $hashFeatureName{$name}++;
  if (scalar(@blockLengthList) != scalar(@blockStartList)) {
    if (defined($verbose)){
      print STDERR "Block counts not match in bed12!\n";
    }
    next;
  }

  my $utr5Count = 0;
  my $cdsCount  = 0;
  my $utr3Count = 0;
  if ($cdsStart == $cdsEnd) {
    for (my $i = 0; $i < scalar(@blockLengthList); $i++) {
      my $blockStart = $start + $blockStartList[$i];
      my $blockEnd = $blockStart + $blockLengthList[$i];
      my $blockCount = $i + 1;
      if ($value eq 'span') {
        print OUT join ( "\t", ( $chr, $blockStart, $blockEnd, "$name|exon|$blockCount", $blockLengthList[$i], $strand ) ), "\n";
      }
    }
  }else{
    for (my $i = 0; $i < scalar(@blockLengthList); $i++) {
      my $blockStart = $start + $blockStartList[$i];
      my $blockEnd = $blockStart + $blockLengthList[$i];
      my $name4thCol;

      if ($cdsStart > $blockEnd) {
        if ($strand eq '+') {
          $utr5Count++;
          $name4thCol = "$name|5utr|$utr5Count";
        }else{
          $utr3Count++;
          $name4thCol = "$name|3utr|$utr3Count";
        }
        if ($value eq 'span') {
          print OUT join ( "\t", ( $chr, $blockStart, $blockEnd, $name4thCol, ($blockEnd - $blockStart), $strand ) ), "\n";
        }else{
          print OUT join ( "\t", ( $chr, $blockStart, $blockEnd, $name4thCol, $exonValue, $strand ) ), "\n";
        }
      }elsif ($cdsEnd < $blockStart) {
        if ($strand eq '+') {
          $utr3Count++;
          $name4thCol = "$name|3utr|$utr3Count";
        }else{
          $utr5Count++;
          $name4thCol = "$name|5utr|$utr5Count";
        }
        if ($value eq 'span') {
          print OUT join ( "\t", ( $chr, $blockStart, $blockEnd, $name4thCol, ($blockEnd - $blockStart), $strand ) ), "\n";
        }else{
          print OUT join ( "\t", ( $chr, $blockStart, $blockEnd, $name4thCol, $exonValue, $strand ) ), "\n";
        }
      }elsif ($cdsStart >= $blockStart and $cdsEnd <= $blockEnd) {
        if ($cdsStart > $blockStart) {
          if ($strand eq '+') {
            $utr5Count++;
            $name4thCol = "$name|5utr|$utr5Count";
          }else{
            $utr3Count++;
            $name4thCol = "$name|3utr|$utr3Count";
          }
          if ($value eq 'span') {
            print OUT join ( "\t", ( $chr, $blockStart, $cdsStart, $name4thCol, ($cdsStart - $blockStart), $strand ) ), "\n";
          }else{
            print OUT join ( "\t", ( $chr, $blockStart, $cdsStart, $name4thCol, $exonValue, $strand ) ), "\n";
          }
        }
        if ($cdsEnd < $blockEnd) {
          if ($strand eq '+') {
            $utr3Count++;
            $name4thCol = "$name|3utr|$utr3Count";
          }else{
            $utr5Count++;
            $name4thCol = "$name|5utr|$utr5Count";
          }
          if ($value eq 'span') {
            print OUT join ( "\t", ( $chr, $cdsEnd, $blockEnd, $name4thCol, ($blockEnd - $cdsEnd), $strand ) ), "\n";
          }else{
            print OUT join ( "\t", ( $chr, $cdsEnd, $blockEnd, $name4thCol, $exonValue, $strand ) ), "\n";
          }
        }
        $cdsCount++;
        if ($value eq 'span') {
          print OUT join ( "\t", ( $chr, $cdsStart, $cdsEnd, "$name|cds|$cdsCount", ($cdsEnd - $cdsStart), $strand ) ), "\n";
        }else{
          print OUT join ( "\t", ( $chr, $cdsStart, $cdsEnd, "$name|cds|$cdsCount", $exonValue, $strand ) ), "\n";
        }
      }elsif ($cdsStart >= $blockStart and $cdsEnd > $blockEnd) {
        if ($cdsStart > $blockStart) {
          if ($strand eq '+') {
            $utr5Count++;
            $name4thCol = "$name|5utr|$utr5Count";
          }else{
            $utr3Count++;
            $name4thCol = "$name|3utr|$utr3Count";
          }
          if ($value eq 'span') {
            print OUT join ( "\t", ( $chr, $blockStart, $cdsStart, $name4thCol, ($cdsStart - $blockStart), $strand ) ), "\n";
          }else{
            print OUT join ( "\t", ( $chr, $blockStart, $cdsStart, $name4thCol, $exonValue, $strand ) ), "\n";
          }
        }
        $cdsCount++;
        if ($value eq 'span') {
          print OUT join ( "\t", ( $chr, $cdsStart, $blockEnd, "$name|cds|$cdsCount", ($blockEnd - $cdsStart), $strand ) ), "\n";
        }else{
          print OUT join ( "\t", ( $chr, $cdsStart, $blockEnd, "$name|cds|$cdsCount", $exonValue, $strand ) ), "\n";
        }
      }elsif ($cdsStart < $blockStart and $cdsEnd <= $blockEnd) {
        if ($cdsEnd < $blockEnd) {
          if ($strand eq '+') {
            $utr3Count++;
            $name4thCol = "$name|3utr|$utr3Count";
          }else{
            $utr5Count++;
            $name4thCol = "$name|5utr|$utr5Count";
          }
          if ($value eq 'span') {
            print OUT join ( "\t", ( $chr, $cdsEnd, $blockEnd, $name4thCol, ($blockEnd - $cdsEnd), $strand ) ), "\n";
          }else{
            print OUT join ( "\t", ( $chr, $cdsEnd, $blockEnd, $name4thCol, $exonValue, $strand ) ), "\n";
          }
        }
        $cdsCount++;
        if ($value eq 'span') {
          print OUT join ( "\t", ( $chr, $blockStart, $cdsEnd, "$name|cds|$cdsCount", ($cdsEnd - $blockStart), $strand ) ), "\n";
        }else{
          print OUT join ( "\t", ( $chr, $blockStart, $cdsEnd, "$name|cds|$cdsCount", $exonValue, $strand ) ), "\n";
        }
      }else{
        $cdsCount++;
        if ($value eq 'span') {
          print OUT join ( "\t", ( $chr, $blockStart, $blockEnd, "$name|cds|$cdsCount", ($blockEnd - $blockStart), $strand ) ), "\n";
        }else{
          print OUT join ( "\t", ( $chr, $blockStart, $blockEnd, "$name|cds|$cdsCount", $exonValue, $strand ) ), "\n";
        }
      }
    }
  }
}
close(BED12);

if (defined($sort)) {
  if ($sort eq 'n') {
    `sort -t \$'\t' -k 1,1 -k 2,2n $output -o $output`;
  }elsif($sort eq 'V') {
    `sort -t \$'\t' -k 1,1V -k 2,2n $output -o $output`;
  }else{
    print STDERR "Error in -sort! Will sort the ouput with default.\n";
    `sort -t \$'\t' -k 1,1 -k 2,2n $output -o $output`;
  }
}

sub BedTypeJudge{
  my ($bedFile, undef) = @_;
  open (BED, "<$bedFile") or die "Cannot open bed file: $bedFile, $!\n";
  my $aveCount;
  my $sum = 0;
  my $lineCount = 1;
  while (my $line=<BED>) {
    next if ($line =~ /^#/);
    if ($lineCount >= 20 and $lineCount < 30) {
      my @lineContents = split ("\t", $line);
      $sum += scalar(@lineContents);
    }else{
      if ($lineCount >= 30) {
        last;
      }
    }
    $lineCount++;
  }
  close(BED);

  $aveCount = $sum / 10;

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

################# Abbreviations for this script #################
#
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for formatting bed12 to bed6.

=head1 CAUTIONS

Defaults: -s=(n), -repat=(false)
The 4th column: name|cds|cdsCount

=head1 SYNOPSIS

bed12Tobed6.pl [options] -b [file] -o [file]

 Options:
    -b | --bed                The input bed12 file
    -o | --output             The output file name
    -r | --repeat             Keep the repeat name record
    -s | --sort               Sort the ouput ('n'[default] or 'V')
    -average                  Average the input value (block count)
    -value                    Value for score column ('span' or int)
    -h | --help               Brief help message
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

bed12Tobed6.pl [options] -b [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will format bed12 to bed6.

=cut
