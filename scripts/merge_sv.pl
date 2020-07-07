#!/usr/bin/perl
# merge_sv.pl
# merge SVs and remove redundant ones
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, USA
# Email: lixing_yang@hms.harvard.edu
# i.e. perl scripts/merge_sv.pl outfile infile1 infile2 ...

use strict;

my $outfile = $ARGV[0];
shift @ARGV;

my $tempfile = $outfile.'.temp';
my (%variants, $newline);

foreach (@ARGV)
{
	open FILE, "<$_";
	while ($newline = <FILE>)
	{
		$variants{$newline}++;
	}
	close FILE;
}

foreach my $key1 (keys %variants)
{
	foreach my $key2 (keys %variants)
	{
		next if ($key1 eq $key2);
		my @data1 = split (/\t/, $key1);
		my @data2 = split (/\t/, $key2);
		if ($data1[2] =~ /\// and $data1[2] =~ /$data2[2]/)
		{
			my @ids = split (/\//, $data1[2]);
			if ($ids[0] eq $data2[2] or $ids[1] eq $data2[2])
			{
				delete $variants{$key2};
			}
		}
	}
}


open OUT, ">$tempfile";
foreach my $key (keys %variants)
{
	print OUT "$key";
}
close OUT;

system "sort $tempfile -k 1,1 -k 6,6 -k 7,7n >$outfile";
system "rm $tempfile";
