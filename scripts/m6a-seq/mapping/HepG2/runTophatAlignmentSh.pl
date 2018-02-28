#!/usr/bin/perl
use strict;
use warnings;

my $absPath = `readlink -f ./`;
chomp ($absPath);
my @fastqFiles = split ("\n", `find $absPath -maxdepth 1 -type f -name "*rep*.fastq"`); #HepG2_m6A-seq_shWTAP_IP_rep1.fastq
@fastqFiles = sort @fastqFiles;
if (! -e "$absPath/bash") {
	`mkdir $absPath/bash`;
}

my $shFile = "$absPath/bash/2_do_tophat_for_m6a_peak_sort_and_rename_bam.sh";

open (OUT, ">$shFile") or die "cannot open file\n";
foreach my $fastqFile (@fastqFiles) {
	my $filePathAbs = (split ("\n",`readlink -f $fastqFile`))[0];
	my @temp = split (/\\|\//, $fastqFile);
	my $fileName = $temp[-1];
	@temp = split (/_|\./, $fileName);
	my $shTarget = $temp[2];
	my $typeRep = join ("_", ($temp[3], $temp[4]));
	my $outPath = "$absPath/$shTarget";
	`mkdir $outPath` unless (-e "$outPath");
	`mkdir $outPath/$typeRep` unless (-e "$outPath/$typeRep");
	`mkdir $outPath/${shTarget}_m6A` unless (-e "$outPath/${shTarget}_m6A");
	my $bashFile = "$absPath/bash/tophat_for_m6a_peak_${shTarget}_${typeRep}";
	open (BASH, ">$bashFile.sh") or die "cannot open file\n";
	print BASH "echo \"\\n\\n\\n\"\necho $fileName\n";
	print BASH "tophat -p 2 -g 1 --library-type=fr-firststrand -G /data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.v24lift37.annotation.gtf -o $outPath/$typeRep/ /public/genomes/human/hg19/bowtie2Index_UCSC/hg19 $filePathAbs\n";
	print BASH "rename $outPath/$typeRep/accepted_hits.bam $outPath/$typeRep/$fileName.bam $outPath/$typeRep/accepted_hits.bam\n";
	print BASH "mv $outPath/$typeRep/$fileName.bam $outPath/\n";
	print BASH "samtools sort --threads 16 -m 2G -O bam -o $outPath/$fileName.sorted.bam $outPath/$fileName.bam\n";
	print BASH "rm $outPath/$fileName.bam\n";
	print BASH "samtools index -b $outPath/$fileName.sorted.bam\n\n";
	print OUT "nohup $bashFile.sh > $absPath/bash/tophatLog/tophat_for_m6a_peak_${shTarget}_${typeRep}.alignment.log 2>&1 &\n";
	close (BASH);
	`chmod 744 $bashFile.sh`;
}
close OUT;

`chmod 744 $shFile`;
