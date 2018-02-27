#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
###develop by K.R.Chow, designed for formatting exomePeak (bed12) to summit (bed6).

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my @inputFiles;
my $output;
my $header = 'false';
my $interval = 0.01;
my $verbose;

GetOptions (
  "i|input=s{1,}"                         =>\@inputFiles,
  "o|output=s{1,1}"                       =>\$output,
  "interval=s{1,1}"                       =>\$interval,
  "header=s{1,1}"                         =>\$header,
  "verbose"                               =>\$verbose,
  "h|help"                                =>\$help,
  "man"                                   =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

if ($header =~ /^false$/i) {
  $header = 0;
}elsif ($header =~ /^true$/i) {
  $header = 1;
}else{
  &verbose("\"Invalid -header: $header! Will set -header to false\"");
  $header = 0;
}

my $intervalDecimal;
if ($interval =~ /^0\.0+?(\d)$/i) {
  $intervalDecimal = $1;
}else{
  &verbose("\"Invalid -interval: $interval! Will set -interval to $interval\"");
  $interval = $interval;
}

my $precision = 1;
while (($interval * (10 ** $precision)) < $intervalDecimal) {
  $precision++;
}

my @outputPath  = split (/\\|\//, $output);
splice (@outputPath, -1, 1);
if (scalar(@outputPath) == 0) {
  unshift (@outputPath , ".");
}

my $outputPath  = join "/", @outputPath;

my %hashNumCount;
my $min;
my $max;
for (my $i = 0; $i < scalar(@inputFiles); $i++) {
  open (IN, "<$inputFiles[$i]") or die "cannot open file $inputFiles[$i]: $!\n";
  if ($header) {
    <IN>;
  }
  while (my $line=<IN>) {
    chomp($line);
    my $number = sprintf ("%.${precision}f", $line);
    if (!defined($min)) {
      $min = int($number);
    }else{
      if ($min > int($number)) {
        $min = int($number);
      }
    }
    if (!defined($max)) {
      $max = $number;
    }else{
      if ($max < $number) {
        $max = $number;
      }
    }
    $hashNumCount{$i}->{'sum'}++;
    $hashNumCount{$i}->{'number'}->{$number}++;
  }
  close(IN);
}

my %hashCumulVal;
for (my $i = 0; $i < scalar(@inputFiles); $i++) {
  my $totalNum = $hashNumCount{$i}->{'sum'};
  my $cumulNum = 0;
  for (my $k = $min; $k < ($max + 5*$interval); $k = $k + $interval) {
    $k = sprintf ("%.${precision}f", $k);
    if (exists($hashNumCount{$i}->{'number'}->{$k})) {
      my $number = $hashNumCount{$i}->{'number'}->{$k};
      $cumulNum += $number;
    }
    $hashCumulVal{$k}->{$i} = sprintf ("%.4f", ($cumulNum / $totalNum));
  }
}

open (OUT, ">$output") or die "cannot open file $output: $!\n";
for (my $k = $min; $k < ($max + $interval); $k = $k + $interval) {
  $k = sprintf ("%.${precision}f", $k);
  print OUT $k;
  for (my $i = 0; $i < scalar(@inputFiles); $i++) {
    print OUT "\t", $hashCumulVal{$k}->{$i};
  }
  print OUT "\n";
}



sub verbose{
  my $warningText = shift;
  if (defined($verbose)) {
    print STDERR $warningText, "\n";
  }
}


################# Abbreviations for this script #################
#
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for calculating the fold change of control peak under different conditions.

=head1 CAUTIONS

Defaults: -header=(false), -interval=(0.01);

=head1 SYNOPSIS

cumulativePlot.pl [options] --input [file1 file2] -o [file]

 Options:

    -i | --input              input files (file1 file2...)
    -o | --output             The output file name
    -interval                 Interval between tick
    -header                   Whether the input files contain header
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

cumulativePlot.pl [options] --input [file1 file2] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will calculate the fold change of control peak under different conditions.

=cut
