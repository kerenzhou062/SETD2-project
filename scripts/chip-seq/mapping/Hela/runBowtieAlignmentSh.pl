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

my $shFile = "$absPath/bash/2_do_bowtie_and_sam_to_sort_bam.sh";

open (OUT, ">$shFile") or die "cannot open file\n";
foreach my $fastqFile (@fastqFiles) {
	my @temp = split (/\\|\//, $fastqFile);
	my $fileName = $temp[-1];
	@temp = split (/_|\./, $fileName);
	my $shTarget = $temp[2];
	my $type = $temp[3];
	my $typeRep = join ("_", ($temp[3], $temp[4]));
	my $outPath = "$absPath/$shTarget";
	`mkdir $outPath` unless (-e "$outPath");
	`mkdir $outPath/$type` unless (-e "$outPath/$type");
	`mkdir $outPath/$typeRep` unless (-e "$outPath/$typeRep");
	`mkdir $outPath/${shTarget}_macs2` unless (-e "$outPath/${shTarget}_macs2");
	`mkdir $outPath/${shTarget}_macs` unless (-e "$outPath/${shTarget}_macs");
	my $bashFile = "$absPath/bash/bowtie_for_chip_peak_${shTarget}_${typeRep}";
	open (BASH, ">$bashFile.sh") or die "cannot open file\n";
	print BASH "echo -e \"$fileName\n\"\n";

	print BASH "bowtie -t -p 16 -v 2 -m 1 --best --strata --sam /public/genomes/human/hg19/bowtieIndex_UCSC/hg19 $fastqFile $outPath/$typeRep/$fileName.sam\n";
	print BASH "samtools view -h -bS -F 4 --threads 32 $outPath/$typeRep/$fileName.sam -o $outPath/$typeRep/$fileName.bam\n";
	print BASH "rm $outPath/$typeRep/$fileName.sam\n";
	print BASH "samtools sort --threads 16 -m 2G -O bam -o $outPath/$typeRep/$fileName.sorted.bam $outPath/$typeRep/$fileName.bam\n";
	print BASH "rm $outPath/$typeRep/$fileName.bam\n";
	print BASH "samtools index -b $outPath/$typeRep/$fileName.sorted.bam\n\n";
	print OUT "nohup $bashFile.sh > $absPath/bash/bowtieLog/bowtie_for_chip_peak_${shTarget}_${typeRep}.alignment.log 2>&1 &\n";
	close (BASH);
	`chmod 744 $bashFile.sh`;
}
close OUT;

`chmod 744 $shFile`;
