#!/usr/bin/perl
# primers.pl
# extract flanking region of break points and design primers by primer3, check number of blat hit for all primers
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, USA
# Email: lixing_yang@hms.harvard.edu
# i.e. perl primers.pl -i fusion_candidates.txt -o primers -F /db/hg19/hg19_fasta/ -P /opt/primer3/src/ -L /opt/blat/ -V 10.119.10.7 -T 17777 -p LY

use strict;
use Getopt::Std;
use Bio::DB::Fasta;

my %opts = (f=>500, s=>"20,23,25,27", m=>"50,60,65", g=>"40,60", n=>5, r=>0, q=>0);
getopts("i:o:p:c:f:s:m:g:n:r:q:F:P:L:V:T:h", \%opts);

my $inputfile = $opts{i};
my $outputfile = $opts{o};
my $prefix = $opts{p};
my $column_offset = $opts{c};
my $flank = $opts{f};
my @size = split (/,/, $opts{s});
my @tm = split (/,/, $opts{m});
my @gc = split (/,/, $opts{g});
my $primer_number = $opts{n};
my $rmask = $opts{r};
my $print_seq = $opts{q}; # print flanking seq
my $reference_path = $opts{F};
my $primer3_path = $opts{P};
my $primer3_command;
if (defined($opts{P}))
{
	$primer3_path .= '/' unless ($primer3_path =~ /\/$/);
	$primer3_command = $primer3_path.'primer3_core';
}
else
{
	$primer3_command = 'primer3_core';
}
my $blat_path = $opts{L};
my $gfClient_command;
if (defined($opts{L}))
{
	$blat_path .= '/' unless ($blat_path =~ /\/$/);
	$gfClient_command = $blat_path.'gfClient';
}
else
{
	$gfClient_command = 'gfClient';
}
my $blat_server = $opts{V};
my $blat_port = $opts{T};

&print_usage if (defined($opts{h}));
die "Please specify inputfile\n" unless (defined($inputfile));
die "Please specify outputfile\n" unless (defined($outputfile));
die "Please specify path of reference\n" unless (defined($reference_path));

my $db = Bio::DB::Fasta->new($reference_path);

my $newline;
my $i = 1;
open FILE, "<$inputfile";
open OUT, ">$outputfile";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split ('\t', $newline);
	print OUT "$newline\n";
	print STDERR "$newline\n";
	my ($seq, %primer) = &get_primer($data[3+$column_offset], $data[4+$column_offset], $data[5+$column_offset]);
	print OUT "$i.1\t$seq\n" if ($print_seq);
	foreach (@size)
	{
		my $size = $_;
		for (my $j = 0; $j < $primer_number; $j++)
		{
			print OUT "$prefix$i.1\t$primer{$size}{$j}{seq}\t$primer{$size}{$j}{pos}\t$primer{$size}{$j}{tm}\t$primer{$size}{$j}{gc}\t$primer{$size}{$j}{match}\t$primer{$size}{$j}{size}\n";
		}
	}
	my ($seq, %primer) = &get_primer($data[8+$column_offset], $data[9+$column_offset], $data[10+$column_offset]);
	print OUT "$i.2\t$seq\n" if ($print_seq);
	foreach (@size)
	{
		my $size = $_;
		for (my $j = 0; $j < $primer_number; $j++)
		{
			print OUT "$prefix$i.2\t$primer{$size}{$j}{seq}\t$primer{$size}{$j}{pos}\t$primer{$size}{$j}{tm}\t$primer{$size}{$j}{gc}\t$primer{$size}{$j}{match}\t$primer{$size}{$j}{size}\n";
		}
	}
	$i++;
}
close FILE;
close OUT;

system "rm primer.for" if (-e 'primer.for');
system "rm primer.out" if (-e 'primer.out');
system "rm primer.rev" if (-e 'primer.rev');
system "rm primer_input" if (-e 'primer_input');
system "rm out.psl" if (-e 'out.psl');
system "rm query.fa" if (-e 'query.fa');

sub get_primer
{
	my $chr = shift;
	my $pos = shift;
	my $primer_fr = shift;
	my %primer;
	my ($seq, $pick_f, $pick_r);
	if ($primer_fr == 1)
	{
		$seq = $db->seq($chr, $pos - $flank=>$pos);
		$pick_f = 1;
		$pick_r = 0;
	}
	if ($primer_fr == -1)
	{
		$seq = $db->seq($chr, $pos=>$pos + $flank);
		$pick_f = 0;
		$pick_r = 1;
	}
	#print "$chr\t$pos\t$primer_fr\n$seq\n";
	
	if ($rmask)
	{
		$seq =~ s/[atgc]/N/g;
	}
	foreach (@size)
	{
		my $size = $_;
		open PRIMER, ">primer_input";
		print PRIMER "SEQUENCE_ID=primer
SEQUENCE_TEMPLATE=$seq
PRIMER_TASK=pick_detection_primers
PRIMER_PICK_LEFT_PRIMER=$pick_f
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=$pick_r
PRIMER_NUM_RETURN=$primer_number
PRIMER_MIN_SIZE=$size[0]
PRIMER_OPT_SIZE=$size
PRIMER_MAX_SIZE=$size[-1]
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=200-300
PRIMER_MIN_TM=$tm[0]
PRIMER_OPT_TM=$tm[1]
PRIMER_MAX_TM=$tm[2]
PRIMER_MIN_GC=$gc[0]
PRIMER_MAX_GC=$gc[1]
P3_FILE_FLAG=1
PRIMER_EXPLAIN_FLAG=1
=
";
		close PRIMER;
		
		my $primer_outfile = 'primer.out';
		system "$primer3_command primer_input > $primer_outfile";
		my $j = 0;
		open TEMP, "<$primer_outfile";
		while ($newline = <TEMP>)
		{
			#print "$newline";
			chomp $newline;
			if ($newline =~ /SEQUENCE=/)
			{
				$primer{$size}{$j}{seq} = $';#'
			}
			if ($newline =~ /PRIMER_\w{4,5}_\d{1,3}=/)
			{
				$primer{$size}{$j}{pos} = $';#'
			}
			if ($newline =~ /TM=/)
			{
				$primer{$size}{$j}{tm} = $';#'
			}
			if ($newline =~ /GC_PERCENT=/)
			{
				$primer{$size}{$j}{gc} = $';#'
				if ($blat_server)
				{
					$primer{$size}{$j}{match} = &blat($primer{$size}{$j}{seq});
				}
				my @pos = split (/,/, $primer{$size}{$j}{pos});
				$primer{$size}{$j}{size} = abs($flank - $pos[0]);
				if ($primer_fr == -1)
				{
					$primer{$size}{$j}{size} = $pos[0];
				}
				#print "$size $j $primer{$size}{$j}{seq}\t$primer{$size}{$j}{pos}\t$primer{$size}{$j}{tm}\t$primer{$size}{$j}{gc}\t$primer{$size}{$j}{match}\t$primer{$size}{$j}{size}\n";
				$j++;
			}
		}
		close TEMP;
	}
	return ($seq, %primer);
}

sub blat
{
	my $query = shift;
	
	my $queryfile = 'query.fa';
	open QR, ">$queryfile";
	print QR ">test\n$query\n";
	close QR;
	
	system "$gfClient_command $blat_server $blat_port / query.fa out.psl -nohead >/dev/null";
	open BTEMP, "cat out.psl|wc|";
	my $newline = <BTEMP>;
	chomp $newline;
	$newline =~ / +(\d{1,10}) +/;
	my $match = $1;
	close BTEMP;
	return $match;
}


sub print_usage
{
	die "primers.pl [options]
	-i FILE	input file, required, fusions file generated by fusions.pl
	-o FILE	output file, required
	-p STR	prefix of primers
	-c INT	column offset, default 0, if you have sample name in column 1 for fusion events, set to 1
	-f INT	flanking region, default 500, the region to design primers in
	-s STR	primer sizes seperated by ",", default 20,23,25,27
	-m STR	primer min, opt, max Tm seperated by ",", default 50,60,65
	-m STR	primer min, max GC content seperated by ",", default 40,60
	-n INT	number of primers to design for each primer size, default 5
	-r INT	mask repeats, default 0
	-q INT	print out flanking sequences to design primers
	-F STR	/path/to/reference/fasta/files, path only, not the files, required, same as option F in meerkat.pl
	-P STR	/path/to/primer3_core, path only, not the command, no need to specify if primer3_core is in PATH
	-L STR	/path/to/blat, path only, not the command, no need to specify if blat is in PATH
	-V STR	blat server (e.g. run server as: gfServer start 10.11.240.76 17777 /reference/hg18/hg18.2bit -stepSize=5, the server name will be 10.11.240.76)
	-T STR	blat port (17777 in above example)
	-h help\n";
}
