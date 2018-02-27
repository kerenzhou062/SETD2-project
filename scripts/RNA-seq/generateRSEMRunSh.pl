#!/usr/bin/perl
use strict;
use warnings;
die "perl generateRSEMRunSh.pl /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq HepG2" if (scalar(@ARGV)==0);
my $fastqFilesPath = $ARGV[0];
my $cellLine = $ARGV[1];
my $absPath = `readlink -f $fastqFilesPath`;
chomp ($absPath);
my @fastqFiles;
if ($cellLine eq "HepG2" or $cellLine eq "Hela" ){
  @fastqFiles = split ("\n", `find $absPath -maxdepth 1 -type f -name "*rep*.fastq"`); #HepG2_m6A-seq_shWTAP_IP_rep1.fastq
 }else{
  @fastqFiles = split ("\n", `find $absPath -maxdepth 1 -type f -name "*.fastq"`); #HepG2_m6A-seq_shWTAP_IP_rep1.fastq
 }
@fastqFiles = sort @fastqFiles;
if (! -e "$cellLine/bash") {
    `mkdir $cellLine/bash`;
}

if (! -e "$cellLine/bash/log") {
    `mkdir $cellLine/bash/log`;
}

if (! -e "$cellLine/RSEM") {
    `mkdir $cellLine/RSEM`;
}

my $shFile = "$cellLine/bash/2_do_RSEM.sh";

open (OUT, ">$shFile") or die "cannot open file\n";
foreach my $fastqFile (@fastqFiles) {
    my $filePathAbs = (split ("\n",`readlink -f $fastqFile`))[0];
    my @temp = split (/\\|\//, $fastqFile);
    my $fileName = $temp[-1];
    @temp = split (/_|\./, $fileName);
    my $shTarget = $temp[2];
    my $typeRep = join ("_", ($temp[3], $temp[4]));
    my $bashFile = "$cellLine/bash/run_RSEM_${shTarget}_${typeRep}";
    open (BASH, ">$bashFile.sh") or die "cannot open file\n";
    print BASH "echo -e $fileName\n";
    print BASH "rsem-calculate-expression  --phred33-quals \\
                           --strandedness reverse \\
                           --bowtie2 \\
                           --output-genome-bam \\
                           -p 8 \\
                           $filePathAbs \\
                           /data/zhoukr/hhl_setd2_m6a/reference/rsem-index/bowtie2-hg19/hg19 \\
                           /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/$cellLine/RSEM/${shTarget}_${typeRep}\n";
    print OUT "nohup $bashFile.sh > $cellLine/bash/log/run_RSEM_${shTarget}_${typeRep}.log 2>&1 &\n";
    close (BASH);
    `chmod 744 $bashFile.sh`;
}
close OUT;

`chmod 744 $shFile`;
