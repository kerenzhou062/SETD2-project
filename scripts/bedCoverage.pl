#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for calculating coverage of bed.

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $bedFile;
my $genomeSizeFile;
my $output;
my $strand;
my $split;
my $verbose;

GetOptions (
  "b|bed=s{1,1}"                          =>\$bedFile,
  "g|genome=s{1,1}"                       =>\$genomeSizeFile,
  "o|output=s{1,1}"                       =>\$output,
  "strand"                                =>\$strand,
  "split"                                 =>\$split,
  "verbose"                               =>\$verbose,
  "h|help"                                =>\$help,
  "man"                                   =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

my @requeres = ('bedtools');
&Requere(@requeres);

my $tempPath = &FileFolderPath($output);
my $tempFile = "$tempPath/bedCoverage". rand(1000) . rand(1000) .".tmp";

my $lineColNum = &BedToTemp($bedFile, $tempFile);

my $splitCommand = (defined($split) and ($lineColNum == 12)) ? ' -split' : '';

if (defined($strand)) {
  if ($lineColNum >= 6) {
    my $plusFile = $tempFile.".plus";
    my $minusFile = $tempFile.".minus";
    my $plusTemp = $tempFile.".plus.tmp";
    my $minusTemp = $tempFile.".minus.tmp";
    `awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\";}{if(\$6==\"+\"){print \$0;}}' $tempFile > $plusFile`;
    `awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\";}{if(\$6==\"-\"){print \$0;}}' $tempFile > $minusFile`;
    `bedtools genomecov -g $genomeSizeFile -i $plusFile -bg $splitCommand > $plusTemp`;
    `bedtools genomecov -g $genomeSizeFile -i $minusFile -bg $splitCommand > $minusTemp`;
    &TempToBed($plusTemp, $plusFile, '+');
    &TempToBed($minusTemp, $minusFile, '-');
    `cat $plusFile $minusFile | sort -t \$'\t' -k1,1 -k2,2n > $output`;
    `rm -f $tempFile $plusFile $plusTemp $minusFile $minusTemp`;
  }else{
    my $plusTemp = $tempFile.".plus";
    `bedtools genomecov -g $genomeSizeFile -i $tempFile -bg $splitCommand > $plusTemp`;
    &TempToBed($plusTemp, $output, '+');
    `rm -f $tempFile $plusTemp`;
  }
}else{
  my $plusTemp = $tempFile.".plus";
  `bedtools genomecov -g $genomeSizeFile -i $tempFile -bg $splitCommand > $plusTemp`;
  &TempToBed($plusTemp, $output, '+');
  `rm -f $tempFile $plusTemp`;
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

sub BedToTemp{
  my ($bed, $temp) = @_;
  open (BED, "<$bed") or die "Cannot open file: $bed, $!\n";
  open (TEMP, ">$temp") or die "Cannot open file: $temp, $!\n";
  my $lineColNum;
  while (my $line=<BED>) {
    next if ($line =~ /^#/);
    chomp($line);
    my @lineContents    = split("\t", $line);
    $lineColNum = scalar(@lineContents);
    if (($lineColNum == 3) or ($lineColNum == 6) or ($lineColNum == 12) ) {
      print TEMP "$line\n";
    }elsif ($lineColNum < 6) {
      $line = join ("\t", @lineContents[0..2]);
      print TEMP "$line\n";
    }elsif ($lineColNum < 12) {
      $line = join ("\t", @lineContents[0..5]);
      print TEMP "$line\n";
    }
  }
  close(BED);
  close(TEMP);
  `sort -t \$'\t' -k1,1 -k2,2n $temp -o $temp`;
  return $lineColNum;
}

sub TempToBed{
  my ($temp, $bed, $strand) = @_;
  open (TEMP, "<$temp") or die "Cannot open file: $temp, $!\n";
  open (BED, ">$bed") or die "Cannot open file: $bed, $!\n";
  my $count = 1;
  while (my $line=<TEMP>) {
    chomp ($line);
    my ($chr, $start, $end, $coverage) = split ("\t", $line);
    print BED join ("\t", ($chr, $start, $end, "posName=$count", $coverage, $strand)), "\n";
    $count++;
  }
  close (TEMP);
  close (BED);
  `sort -t \$'\t' -k1,1 -k2,2n $bed -o $bed`;
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

################# Abbreviations for this script #################
#
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for calculating coverage of bed.
Require: bedtools

=head1 CAUTIONS

Defaults: -s=(V), -repat=(false)
The 4th column: name|cds|cdsCount

=head1 SYNOPSIS

bedCoverage.pl [options] -b [file] -g [file] -o [file]

 Options:
    -b | --bed                The input bed file
    -g | --genome             The genome size file (sort by number)
    -o | --output             The output file name
    -strand                   Coverage on the same strand
    -split                    Split features on bed12
    -h | --help               Brief help message
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

bedCoverage.pl [options] -b [file] -g [file] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will format bed12 to bed6.

=cut
