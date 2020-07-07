#!/usr/bin/perl
# filter_fusion.pl
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, USA
# Email: lixing_yang@hms.harvard.edu

use strict;
use Getopt::Std;
my $version = '';

my %opts = ();
getopts("i:o:f:m:n:p:q:s:h", \%opts);

my $inputfile = $opts{i};
my $outputfile = $opts{o};
my $offset = $opts{f};
my $mp_cutoff = $opts{m};
my $mpf_cutoff = $opts{n};
my $min_cov = $opts{p};
my $max_cov = $opts{q};
my $intra_size_cutoff = $opts{s};

my $newline;
open FILE, "<$inputfile";
open TUMOR, ">$outputfile";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split ('\t', $newline);
	my @rp = split ('_', $data[-1]);
	my $remove;
	if ($mp_cutoff and $data[$offset+16] < $mp_cutoff)
	{
		$remove = 1;
	}
	if ($mpf_cutoff and $rp[0] < $mpf_cutoff)
	{
		$remove = 1;
	}
	if ($max_cov)
	{
		if (($rp[0]+$rp[2]) < $max_cov and ($rp[0]+$rp[3]) < $max_cov)
		{
			$remove = 1;
		}
	}
	if ($intra_size_cutoff)
	{
		if ($data[$offset+3] eq $data[$offset+8] and abs($data[$offset+9]-$data[$offset+4]) < $intra_size_cutoff)
		{
			$remove = 1;
		}
	}
	
	next if $remove;
	if ($data[$offset+3] eq $data[$offset+8])
	{
		if ($data[$offset+5] == 1 and $data[$offset+10] == -1)
		{
			if ($data[$offset+18] < 0)
			{
				$data[$offset+13] = 'del_ins';
			}
			else
			{
				$data[$offset+13] = 'del';
			}
		}
		if ($data[$offset+5] == -1 and $data[$offset+10] == 1)
		{
			$data[$offset+13] = 'tandem_dup';
		}
		if ($data[$offset+5] == 1 and $data[$offset+10] == 1)
		{
			$data[$offset+13] = 'invers_f';
		}
		if ($data[$offset+5] == -1 and $data[$offset+10] == -1)
		{
			$data[$offset+13] = 'invers_r';
		}
	}
	else
	{
		$data[$offset+13] = 'transl_inter';
	}
	$data[$offset+14] = 'NA';
	my $toprint = join ("\t", @data);
	print TUMOR "$toprint\n";
}
close FILE;
close TUMOR;

sub print_usage
{
	die "filter_fusions.pl [options]
	-i FILE	input file, required
	-o FILE	output file, required
	-f INT	column offset, if the first column is sample ID, set f as 1
	-m INT	number of supporting discordant read pair cutoff
	-n INT	number of full length supporting discordant read pair cutoff
	-p INT	minimun coverage cutoff
	-q INT	maximum coverage cutoff
	-s INT	size cutoff
	-h help\n";
}