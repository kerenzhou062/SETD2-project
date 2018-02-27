#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Pod::Usage;
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
###gffcompare.pl develop by K.R.Chow, designed for formatting exomePeak (bed12) to summit (bed6).

Getopt::Long::Configure(qw(posix_default no_ignore_case));  # allow the program accept the command-line without order.
pod2usage(1) if (scalar(@ARGV) == 0);

my $help   = 0;
my $man    = 0;
my $output;
my $fpkmCutoff = 1;
my $tpmCutoff = 1;
my @gtfFiles;
my $prefixName = 'gffcompare';
my $controlGtf;
my $uniqueFlag = 1;
my $verbose;

GetOptions (
  "c|control=s{1,1}"                      =>\$controlGtf,
  "o|output=s{1,1}"                       =>\$output,
  "f|fpkm=s{1,1}"                         =>\$fpkmCutoff,
  "t|tpm=s{1,1}"                          =>\$tpmCutoff,
  "g|gtf=s{1,}"                           =>\@gtfFiles,
  "p|prefix=s{1,1}"                       =>\$prefixName,
  "verbose"                               =>\$verbose,
  "h|help"                                =>\$help,
  "man"                                   =>\$man
) or pod2usage(2); pod2usage(1) if $help; pod2usage(-verbose => 2) if $man;

if ($fpkmCutoff <= 0) {
  &verbose("Invalid -pval: $fpkmCutoff! Will set --fpkm to 1");
  $fpkmCutoff = 1;
}

if (! -d "$output") {
  print STDERR "-o must be a folder!";
}


my %hashGene;
unshift(@gtfFiles, $controlGtf);
my @gtfNames;
foreach my $gtf(@gtfFiles) {
  (my $gtfName = $gtf) =~ s/.*(\/|\\)//g;
  $gtfName =~ s/\..*//g;
  push (@gtfNames, $gtfName);
  open (GTF, "<$gtf") or die "cannot open $gtf: $!";
  <GTF>;
  while (my $line=<GTF>) {
    chomp($line);
    my ($geneName, $txID, $classCode, $qryGeneID, $qryID, $numExons, $FPKM, $TPM, $cov, $length, $majorIsoID, $refMatchLen) = split ("\t", $line);
    next if ($geneName eq '-');
    $hashGene{$geneName}->{$gtfName}->{$txID}->{'classCode'}   = $classCode;
    $hashGene{$geneName}->{$gtfName}->{$txID}->{'FPKM'}        = $FPKM;
    $hashGene{$geneName}->{$gtfName}->{$txID}->{'TPM'}         = $TPM;
    $hashGene{$geneName}->{$gtfName}->{$txID}->{'length'}      = $length;
    $hashGene{$geneName}->{$gtfName}->{$txID}->{'refMatchLen'} = $refMatchLen;
  }
  close(GTF);
}

my $maxExpFile = "$output/$prefixName.maxExp.txt";
my @geneNames = sort keys %hashGene;
open (MAXEXP, ">$maxExpFile") or die "cannot open $maxExpFile, $!\n";
print MAXEXP "geneName";
foreach my $gtfName(@gtfNames) {
  print MAXEXP "\tTranscript_${gtfName}\tclassCode\tTPM\tsumTPM\tlength\trefMatchLen";
}
print MAXEXP "\n";

foreach my $geneName (@geneNames) {
  print MAXEXP "$geneName";
  foreach my $gtfName(@gtfNames) {
    if (!exists($hashGene{$geneName}->{$gtfName})) {
      print MAXEXP "\tNA\tNA\tNA\tNA\tNA";
    }else{
      my %hashTranscripExp;
      my @txIDs = sort keys %{$hashGene{$geneName}->{$gtfName}};
      foreach my $txID(@txIDs) {
        $hashTranscripExp{$txID}->{'FPKM'} = $hashGene{$geneName}->{$gtfName}->{$txID}->{'FPKM'};
        $hashTranscripExp{$txID}->{'TPM'}  = $hashGene{$geneName}->{$gtfName}->{$txID}->{'TPM'};
      }
      my ($maxTpmTx, $sumTpm) = &maxExpTranscript(\%hashTranscripExp, 'TPM');
      my $classCode           = $hashGene{$geneName}->{$gtfName}->{$maxTpmTx}->{'classCode'};
      my $TPM                 = $hashGene{$geneName}->{$gtfName}->{$maxTpmTx}->{'TPM'};
      my $length              = $hashGene{$geneName}->{$gtfName}->{$maxTpmTx}->{'length'};
      my $refMatchLen         = $hashGene{$geneName}->{$gtfName}->{$maxTpmTx}->{'refMatchLen'};
      print MAXEXP "\t$maxTpmTx\t$classCode\t$TPM\t$sumTpm\t$length\t$refMatchLen";
    }
  }
  print MAXEXP "\n";
}



################################################ subroutines #############################

sub maxExpTranscript { # return $min, $max
  my ($refData,$type) = @_;
  my @sort_data = sort { $refData->{$a}->{$type} <=> $refData->{$b}->{$type} or $a cmp $b } keys %{$refData};
  my $sum;
  map{$sum += $refData->{$_}->{$type};} keys %{$refData};
  return ($sort_data[-1], $sum);
}

sub verbose{
  my $warningText = shift;
  if (defined($verbose)) {
    print STDERR $warningText, "\n";
  }
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

################# Abbreviations for this script #################
#
#
################# Abbreviations for this script #################

__END__

=head1 NAME

This script is designed for calculating the fold change of control peak under different conditions.

=head1 CAUTIONS

Defaults: -p=(0.05), -f=(1), -fdr=(0.05), -log=(false), -FE=(false)

=head1 SYNOPSIS

gffcompare.pl [options] -x [file1 file2] -o [file]

 Options:

    -o | --output             The output file name
    -p | --pvalue             The cutoff p-value
    -s | --sort               Sort the ouput ('n' or 'V')
    -x | --xls                The input xls files(Cont.xls Treat1.xls Treat2.xls...)
    -h | --help               Brief help message
    -cutoff                   The cutoff value for zero fold_enrichment
    -f                        Minimum overlap required as a fraction of A
    -F                        Minimum overlap required as a fraction of B
    -FE                       Flag for report the original fold_enrichment
    -fdr                      The cutoff of FDR
    -fold                     The cutoff of fold_enrichment
    -log                      Log scale (postive number)
    -overlap                  Only report overlapped parts (0 0 1 1...)
    -keep                     Keep unique peaks
    -verbose                  Print more warnings
    -man                      Full documentation

=head1 OPTIONS

=over 4

=item B<-h|--help>

gffcompare.pl [options] -x [file1 file2] -o [file]

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will calculate the fold change of control peak under different conditions.

=cut
