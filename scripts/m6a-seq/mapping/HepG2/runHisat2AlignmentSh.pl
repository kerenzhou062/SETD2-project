#!/usr/bin/perl
use strict;
use warnings;

my $absPath = `readlink -f ./`;
chomp ($absPath);
my @fastqFiles = split ("\n", `find $absPath -maxdepth 1 -type f -name "*input_rep*.fastq"`); #HepG2_m6A-seq_shWTAP_IP_rep1.fastq
@fastqFiles = sort @fastqFiles;
if (! -e "$absPath/bash") {
	`mkdir $absPath/bash`;
}


my %hashTargetFile;
foreach my $fastqFile (@fastqFiles) {
	my @temp = split (/\\|\//, $fastqFile);
	my $fileName = $temp[-1];
	(my $fileBaseName = $fileName) =~ s/\.fastq//;
	@temp = split (/_|\./, $fileName);
	my $shTarget = $temp[2];
	$hashTargetFile{$shTarget}->{$fastqFile}->{'basename'} = $fileBaseName;
	$hashTargetFile{$shTarget}->{$fastqFile}->{'fileName'} = $fileName;
}

my $shFile = "$absPath/bash/2_do_hisat2_for_transcriptome_sort_and_rename_bam.sh";

open (OUT, ">$shFile") or die "cannot open file\n";
foreach my $shTarget (sort keys %hashTargetFile) {
	my $outPath = "$absPath/transcriptome_assembly/$shTarget";
	`mkdir -p $outPath` unless (-e "$outPath");
	my $bashFile = "$absPath/bash/hisat2_for_transcriptome_${shTarget}";
	my $transcriptGtf = '';
	open (BASH, ">$bashFile.sh") or die "cannot open file\n";
	my @files = sort keys %{$hashTargetFile{$shTarget}};
	foreach my $file (@files) {
		my $basename = $hashTargetFile{$shTarget}->{$file}->{'basename'};
		my $fileName = $hashTargetFile{$shTarget}->{$file}->{'fileName'};
		print BASH "echo \"$fileName\"\n";
		print BASH "hisat2 -p 3 --rna-strandness R --dta -x /data/zhoukr/reference/genome/Index/Hisat2_index/hg19 -q $file -S $outPath/$basename.sam\n";
		print BASH "echo -e \"transcripts alignment done. sam to index bam by using samtools...\\n\"\n";
		print BASH "samtools sort --threads 8 -O bam -o $outPath/$basename.sorted.bam $outPath/$basename.sam\n";
		print BASH "samtools index -b $outPath/$basename.sorted.bam\n";
		print BASH "echo -e \"stringtie: start assbemling transcriptome...\\n\"\n";
		print BASH "stringtie -p 8 --rf -G /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o $outPath/$basename.gtf $outPath/$basename.sorted.bam\n";
		print BASH "sed -i 's/STRG/$basename/g' $outPath/$basename.gtf\n";
		print BASH "rm $outPath/$basename.sam\n";
		$transcriptGtf .= "$outPath/$basename.gtf\n"
	}
	print BASH "echo \"$transcriptGtf\" > $outPath/mergeGtfList.txt\n";
	print BASH "echo -e \"stringtie: start merging the assembled transcriptome...\\n\"\n";
	print BASH "stringtie --merge -p 8 -o $absPath/transcriptome_assembly/${shTarget}_transcriptome_merge.gtf $outPath/mergeGtfList.txt\n";
	print BASH "echo -e \"gffcompare: start comparing to the reference gtf annotation...\\n\"\n";
	print BASH "gffcompare -r /data/zhoukr/reference/genome/gtf/gencode.v24lift37.annotation.gtf -o $outPath/${shTarget}_merged_compare $absPath/transcriptome_assembly/${shTarget}_transcriptome_merge.gtf\n";
	print BASH "echo -e \"jobs have been done.\\n\"\n";
	print OUT "nohup $bashFile.sh > $absPath/bash/hisat2Log/hisat2_for_transcriptome_${shTarget}.transcriptome.log 2>&1 &\n";
	close (BASH);
	`chmod 744 $bashFile.sh`;
}
close OUT;

`chmod 744 $shFile`;
