#!/usr/bin/perl
# discon.pl
# count number of discordant and concordant read pairs for SVs called by Meerkat
# take Meerkat variant file as input
# add RP tag at the end
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, USA

use strict;
use Getopt::Std;

my %opts = (D=>0, d=>3);
getopts("i:o:D:B:C:I:K:Q:S:d:h", \%opts);

my $inputfile = $opts{i};
my $outputfile = $opts{o};
my $count_dis = $opts{D};
my $bamfile = $opts{B};
my $clusterfile = $opts{C};
my $isinfofile = $opts{I};
my $blackrgfile = $opts{K};
my $samtools_path = $opts{S};
my $sd_cutoff = $opts{d};
my $min_mapq = $opts{Q};
$clusterfile = undef if ($count_dis);

&print_usage unless (defined($inputfile) and defined($outputfile));
&print_usage if (defined($opts{h}));
die "bam file doesn't exist\n" unless (-e $bamfile);
die "isinfo file doesn't exist\n" unless (-e $isinfofile);

my $samtools_command;
if (defined($opts{S}))
{
	$samtools_path .= '/' unless ($samtools_path =~ /\/$/);
	$samtools_command = $samtools_path.'samtools';
}
else
{
	$samtools_command = 'samtools';
}

my ($newline, %is);
die "$isinfofile file not exist\n" unless (-e $isinfofile);
open FILE, "<$isinfofile";
while ($newline = <FILE>)
{
	chomp $newline;
	if ($newline =~ /Read length/)
	{
		my ($trash, $rg) = split ("\t", $newline);
		$newline = <FILE>;
		chomp $newline;
		$is{$rg}{'rl'} = $newline;
	}
	if ($newline =~ /Median/)
	{
		my ($trash, $rg) = split ("\t", $newline);
		$newline = <FILE>;
		chomp $newline;
		$is{$rg}{'median'} = $newline;
	}
	if ($newline =~ /Standard deviation/)
	{
		my ($trash, $rg) = split ("\t", $newline);
		$newline = <FILE>;
		chomp $newline;
		$is{$rg}{'sd'} = $newline;
		$is{'sdu'} = $newline if ($newline > $is{'sdu'});
		$is{$rg}{'isu'} = $is{$rg}{'median'} + $is{$rg}{'sd'}*$sd_cutoff;
		$is{$rg}{'isu'} = int($is{$rg}{'isu'}) + 1;
		$is{$rg}{'isd'} = $is{$rg}{'median'} - $is{$rg}{'sd'}*$sd_cutoff;
		$is{$rg}{'isd'} = int($is{$rg}{'isd'}) - 1;
		$is{'isu'} = $is{$rg}{'isu'} if ($is{'isu'} < $is{$rg}{'isu'});
	}
}
my $window_size = $is{'isu'}; # upstream and downstream of break points

my %blackrg;
open FILE, "<$blackrgfile";
while ($newline = <FILE>)
{
	chomp $newline;
	$blackrg{$newline} = 1;
}
close FILE;

my (@variants, %id);
open FILE, "<$inputfile";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split ('\t', $newline);
	push @variants, \@data;
	my @ids = split ("\/", $data[2]);
	foreach (@ids)
	{
		if ($_ =~ /_/)
		{
			$_ = $`;
		}
		$id{$_} = 1;
	}
	#print "@ids\n";
}
close FILE;

my %remap; # number of unmapped and softclipped pairs in each cluster
open FILE, "<$clusterfile";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split ('\t', $newline);
	if ($id{$data[0]} and $data[1] == 0)
	{
		if ($data[5] =~ /mu1|mu2|sc/)
		{
			$remap{$data[0]}{$`} = 1
		}
		else
		{
			$remap{$data[0]}{disc}++;
		}
	}
}
close FILE;

for my $id (keys %remap)
{
	my $i=0;
	for my $readname (keys %{$remap{$id}})
	{
		next if ($readname eq 'disc' or $readname eq 'remap');
		$i++;
	}
	$remap{$id}{remap} = $i;
	#print "$id\t$remap{$id}{disc}\t$remap{$id}{remap}\n";
}

open OUT, ">$outputfile";
foreach (@variants)
{
	my @variant = @$_;
	while (1)
	{
		if ($variant[-1] =~ /BP:/)
		{
			last;
		}
		else
		{
			delete $variant[-1];
		}
	}
	my $toprint = join ("\t", @variant);
	my @id = split ("\/", $variant[2]);
	foreach (@id)
	{
		if ($_ =~ /_/)
		{
			$_ = $`;
		}
	}
	my @disc = split ("\/", $variant[3]);
	if (defined($disc[2]))
	{
		$disc[0] = $disc[1];
		$disc[1] = $disc[2];
		undef($disc[2]);
	}
	
	if ($variant[0] eq 'del')
	{
		my $hm = $variant[9];
		my $concord1 = &ccd($variant[5], $variant[6], $hm);
		my $concord2 = &ccd($variant[5], $variant[7], $hm);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], 1, $variant[5], $variant[7], -1, $hm);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concord1".'_'."$concord2\n";
		#print "$d_total\n";
	}
	elsif ($variant[0] eq 'del_ins')
	{
		my $hm = -$variant[12];
		my $concord1 = &ccd($variant[5], $variant[6], $hm);
		my $concord2 = &ccd($variant[5], $variant[7], $hm);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], 1, $variant[5], $variant[7], -1, $hm);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concord1".'_'."$concord2\n";
		#print "$d_total\n";
	}
	elsif ($variant[0] eq 'tandem_dup')
	{
		my $hm = $variant[9];
		my $concord1 = &ccd($variant[5], $variant[6], $hm);
		my $concord2 = &ccd($variant[5], $variant[7], $hm);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], -1, $variant[5], $variant[7], 1, $hm);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concord1".'_'."$concord2\n";
	}
	elsif ($variant[0] eq 'invers_f')
	{
		my $hm = $variant[9];
		my $concord1 = &ccd($variant[5], $variant[6], $hm);
		my $concord2 = &ccd($variant[5], $variant[7], $hm);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], 1, $variant[5], $variant[7], 1, $hm);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concord1".'_'."$concord2\n";
	}
	elsif ($variant[0] eq 'invers_r')
	{
		my $hm = $variant[9];
		my $concord1 = &ccd($variant[5], $variant[6], $hm);
		my $concord2 = &ccd($variant[5], $variant[7], $hm);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], -1, $variant[5], $variant[7], -1, $hm);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concord1".'_'."$concord2\n";
	}
	elsif ($variant[0] eq 'invers')
	{
		my @hm = split ("\/", $variant[9]);
		my $concorda1 = &ccd($variant[5], $variant[6], $hm[0]);
		my $concorda2 = &ccd($variant[5], $variant[7], $hm[0]);
		my $concordb1 = &ccd($variant[5], $variant[6], $hm[1]);
		my $concordb2 = &ccd($variant[5], $variant[7], $hm[1]);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], 1, $variant[5], $variant[7], 1, $hm[0]);
			#print "$remap{$id[0]}{disc}\n";
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		unless (-e $clusterfile)
		{
			$remap{$id[1]}{disc} = &dis($variant[5], $variant[6], -1, $variant[5], $variant[7], -1, $hm[1]);
			#print "$remap{$id[1]}{disc}\n";
			unless ($count_dis)
			{
				$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
				($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
	}
	elsif ($variant[0] =~ /inssd/)
	{
		my @hm = split ("\/", $variant[14]);
		my $concorda1 = &ccd($variant[5], $variant[7], $hm[0]);
		my $concorda2 = &ccd($variant[9], $variant[11], $hm[0]);
		my $concordb1 = &ccd($variant[5], $variant[6], $hm[1]);
		my $concordb2 = &ccd($variant[9], $variant[10], $hm[1]);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[7], -1, $variant[9], $variant[11], 1, $hm[0]);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		unless (-e $clusterfile)
		{
			$remap{$id[1]}{disc} = &dis($variant[5], $variant[6], 1, $variant[9], $variant[10], -1, $hm[1]);
			unless ($count_dis)
			{
				$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
				($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
	}
	elsif ($variant[0] =~ /inssu/)
	{
		my @hm = split ("\/", $variant[14]);
		my $concorda1 = &ccd($variant[9], $variant[10], $hm[0]);
		my $concorda2 = &ccd($variant[5], $variant[6], $hm[0]);
		my $concordb1 = &ccd($variant[9], $variant[11], $hm[1]);
		my $concordb2 = &ccd($variant[5], $variant[7], $hm[1]);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[9], $variant[10], -1, $variant[5], $variant[6], 1, $hm[0]);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		unless (-e $clusterfile)
		{
			$remap{$id[1]}{disc} = &dis($variant[9], $variant[11], 1, $variant[5], $variant[7], -1, $hm[1]);
			unless ($count_dis)
			{
				$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
				($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
	}
	elsif ($variant[0] =~ /insod/)
	{
		my @hm = split ("\/", $variant[14]);
		my $concorda1 = &ccd($variant[5], $variant[6], $hm[0]);
		my $concorda2 = &ccd($variant[9], $variant[11], $hm[0]);
		my $concordb1 = &ccd($variant[5], $variant[7], $hm[1]);
		my $concordb2 = &ccd($variant[9], $variant[10], $hm[1]);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], 1, $variant[9], $variant[11], 1, $hm[0]);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		unless (-e $clusterfile)
		{
			$remap{$id[1]}{disc} = &dis($variant[5], $variant[7], -1, $variant[9], $variant[10], -1, $hm[1]);
			unless ($count_dis)
			{
				$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
				($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
	}
	elsif ($variant[0] =~ /insou/)
	{
		my @hm = split ("\/", $variant[14]);
		my $concorda1 = &ccd($variant[9], $variant[11], $hm[0]);
		my $concorda2 = &ccd($variant[5], $variant[6], $hm[0]);
		my $concordb1 = &ccd($variant[9], $variant[10], $hm[1]);
		my $concordb2 = &ccd($variant[5], $variant[7], $hm[1]);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[9], $variant[11], 1, $variant[5], $variant[6], 1, $hm[0]);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		unless (-e $clusterfile)
		{
			$remap{$id[1]}{disc} = &dis($variant[9], $variant[10], -1, $variant[5], $variant[7], -1, $hm[1]);
			unless ($count_dis)
			{
				$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
				($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
	}
	elsif ($variant[0] eq 'del_invers')
	{
		my @hm = split ("\/", $variant[15]);
		my $concorda1 = &ccd($variant[5], $variant[6], $hm[0]);
		my $concorda2 = &ccd($variant[9], $variant[11], $hm[0]);
		my $concordb1 = &ccd($variant[5], $variant[7], $hm[1]);
		my $concordb2 = &ccd($variant[9], $variant[10], $hm[1]);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], 1, $variant[9], $variant[11], 1, $hm[0]);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		unless (-e $clusterfile)
		{
			$remap{$id[1]}{disc} = &dis($variant[9], $variant[10], -1, $variant[5], $variant[7], -1, $hm[1]);
			unless ($count_dis)
			{
				$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
				($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
	}
	elsif ($variant[0] =~ /inss$/)
	{
		my @hm = split ("\/", $variant[13]);
		if ($variant[5] lt $variant[9])
		{
			my $concorda1 = &ccd($variant[5], $variant[6], $hm[0]);
			my $concorda2 = &ccd($variant[9], $variant[10], $hm[0]);
			my $concordb1 = &ccd($variant[5], $variant[7], $hm[1]);
			my $concordb2 = &ccd($variant[9], $variant[11], $hm[1]);
			unless (-e $clusterfile)
			{
				$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], 1, $variant[9], $variant[10], -1, $hm[0]);
				unless ($count_dis)
				{
					$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
					($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
				}
			}
			unless (-e $clusterfile)
			{
				$remap{$id[1]}{disc} = &dis($variant[5], $variant[7], -1, $variant[9], $variant[11], 1, $hm[1]);
				unless ($count_dis)
				{
					$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
					($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
				}
			}
			$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
		}
		else
		{
			my $concorda1 = &ccd($variant[5], $variant[7], $hm[0]);
			my $concorda2 = &ccd($variant[9], $variant[11], $hm[0]);
			my $concordb1 = &ccd($variant[5], $variant[6], $hm[1]);
			my $concordb2 = &ccd($variant[9], $variant[10], $hm[1]);
			unless (-e $clusterfile)
			{
				$remap{$id[0]}{disc} = &dis($variant[5], $variant[7], 1, $variant[9], $variant[11], -1, $hm[0]);
				unless ($count_dis)
				{
					$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
					($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
				}
			}
			unless (-e $clusterfile)
			{
				$remap{$id[1]}{disc} = &dis($variant[5], $variant[6], -1, $variant[9], $variant[10], 1, $hm[1]);
				unless ($count_dis)
				{
					$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
					($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
				}
			}
			$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
		}
	}
	elsif ($variant[0] =~ /inso$/)
	{
		my @hm = split ("\/", $variant[13]);
		if ($variant[5] lt $variant[9])
		{
			my $concorda1 = &ccd($variant[5], $variant[6], $hm[0]);
			my $concorda2 = &ccd($variant[9], $variant[11], $hm[0]);
			my $concordb1 = &ccd($variant[5], $variant[7], $hm[1]);
			my $concordb2 = &ccd($variant[9], $variant[10], $hm[1]);
			unless (-e $clusterfile)
			{
				$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], 1, $variant[9], $variant[11], 1, $hm[0]);
				unless ($count_dis)
				{
					$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
					($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
				}
			}
			unless (-e $clusterfile)
			{
				$remap{$id[1]}{disc} = &dis($variant[5], $variant[7], -1, $variant[9], $variant[10], -1, $hm[1]);
				unless ($count_dis)
				{
					$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
					($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
				}
			}
			$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
		}
		else
		{
			my $concorda1 = &ccd($variant[5], $variant[6], $hm[0]);
			my $concorda2 = &ccd($variant[9], $variant[11], $hm[0]);
			my $concordb1 = &ccd($variant[5], $variant[7], $hm[1]);
			my $concordb2 = &ccd($variant[9], $variant[10], $hm[1]);
			unless (-e $clusterfile)
			{
				$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], 1, $variant[9], $variant[11], 1, $hm[0]);
				unless ($count_dis)
				{
					$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
					($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
				}
			}
			unless (-e $clusterfile)
			{
				$remap{$id[1]}{disc} = &dis($variant[5], $variant[7], -1, $variant[9], $variant[10], -1, $hm[1]);
				unless ($count_dis)
				{
					$remap{$id[1]}{remap} = $disc[1] - $remap{$id[1]}{disc};
					($remap{$id[1]}{disc}, $remap{$id[1]}{remap}) = ('?', '?') if ($remap{$id[1]}{disc} == 0 or $remap{$id[1]}{remap} < 0);
				}
			}
			$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concorda1".'_'."$concorda2/$remap{$id[1]}{disc}".'_'."$remap{$id[1]}{remap}".'_'."$concordb1".'_'."$concordb2\n";
		}
	}
	elsif ($variant[0] eq 'transl_inter')
	{
		my $concord1 = &ccd($variant[5], $variant[6], $variant[11]);
		my $concord2 = &ccd($variant[8], $variant[9], $variant[11]);
		unless (-e $clusterfile)
		{
			$remap{$id[0]}{disc} = &dis($variant[5], $variant[6], $variant[7], $variant[8], $variant[9], $variant[10], $variant[11]);
			unless ($count_dis)
			{
				$remap{$id[0]}{remap} = $disc[0] - $remap{$id[0]}{disc};
				($remap{$id[0]}{disc}, $remap{$id[0]}{remap}) = ('?', '?') if ($remap{$id[0]}{disc} == 0 or $remap{$id[0]}{remap} < 0);
			}
		}
		$toprint .= "\tRP:$remap{$id[0]}{disc}".'_'."$remap{$id[0]}{remap}".'_'."$concord1".'_'."$concord2\n";
	}
	else
	{
		my $hm = $variant[9];
		my $concord1 = &ccd($variant[5], $variant[6], $hm);
		my $concord2 = &ccd($variant[5], $variant[7], $hm);
		$toprint .= "\tRP:$variant[3]".'_'."".'_'."$concord1".'_'."$concord2\n";
		#print "$d_total\n";
	}
	print OUT "$toprint";
}
close OUT;

# count number of concordant read pairs
sub ccd
{
	my $chr = shift;
	my $pos = shift;
	my $homology = shift;
	my $error_allowed = 5;
	
	my ($p1, $p2);
	$p1 = $pos - $window_size;
	$p2 = $pos + $window_size;
	my $region = $chr.':'.$p1.'-'.$p2;
	#print STDERR "$region\n";
	
	my $concord = 0;
	my $newline1;
	open TEMP, "$samtools_command view -X $bamfile $region|";
	while ($newline1 = <TEMP>)
	{
		chomp $newline1;
		my @tempdata = split ('\t', $newline1);
		#print "@tempdata\n";
		if (defined($min_mapq))
		{
			if ($tempdata[4] < $min_mapq)
			{
				next;	
			}
		}
		else
		{
			if ($newline1 =~ /XT/)
			{
				next unless ($newline1 =~ /XT:A:U/);
			}
		}
		my $rg = 'none';
		if ($newline1 =~ /RG:Z:(\S{1,50})/)
		{
			$rg = $1;
			next if ($blackrg{$rg});
		}
		next unless ($tempdata[6] eq '=');
		next if ($tempdata[8] > $is{$rg}{'isu'});
		#print "\t\t@tempdata\n";
		if ($tempdata[3] <= ($pos-abs($homology)/2-$error_allowed) and $tempdata[7] >= ($pos+abs($homology)/2+$error_allowed))
		{
			#print "$is{$rg}{'isu'}\t$homology\t@tempdata\n";
			$concord++;
		}
	}
	close TEMP;
	
	return $concord;
}

# count number of full length discordant read pairs
sub dis
{
	my $chr1 = shift;
	my $pos1 = shift;
	my $ori1 = shift;
	my $chr2 = shift;
	my $pos2 = shift;
	my $ori2 = shift;
	my $homology = shift;
	my $error_allowed = 5;
	
	my ($p1, $p2);
	$p1 = $pos1 - $window_size;
	$p2 = $pos1 + $window_size;
	my $region1 = $chr1.':'.$p1.'-'.$p2;
	my ($p3, $p4);
	$p3 = $pos2 - $window_size;
	$p4 = $pos2 + $window_size;
	my $region2 = $chr2.':'.$p3.'-'.$p4;
	#print STDERR "$region1\t$region2\n";
	
	my $discord = 0;
	my $newline1;
	open TEMP, "$samtools_command view -X $bamfile $region1|";
	while ($newline1 = <TEMP>)
	{
		chomp $newline1;
		my @tempdata = split ('\t', $newline1);
		#print "@tempdata\n";
		if (defined($min_mapq))
		{
			if ($tempdata[4] < $min_mapq)
			{
				next;	
			}
		}
		else
		{
			if ($newline1 =~ /XT/)
			{
				next unless ($newline1 =~ /XT:A:U/);
			}
		}
		my $rg = 'none';
		if ($newline1 =~ /RG:Z:(\S{1,50})/)
		{
			$rg = $1;
			next if ($blackrg{$rg});
		}
		next if ($tempdata[1] =~ /u|U/);
		
		my $strand = 1;
		$strand = -1 if ($tempdata[1] =~ /r/);
		my $mstrand = 1;
		$mstrand = -1 if ($tempdata[1] =~ /R/);
		my ($disc1, $disc2);
		if ($ori1 eq $strand and $ori2 eq $mstrand)
		{
			if ($chr1 eq $chr2)
			{
				if ($tempdata[6] eq '=')
				{
					if ($tempdata[8] > $is{$rg}{'isu'})
					{
						$disc1 = 1 if ($ori1 == 1 and $tempdata[3] <= ($pos1+abs($homology)/2+$error_allowed) and $tempdata[3] >= $pos1-$is{$rg}{'isu'}-$error_allowed);
						$disc1 = 1 if ($ori1 == -1 and $tempdata[3] >= ($pos1-abs($homology)/2-$error_allowed) and $tempdata[3] <= $pos1+$is{$rg}{'isu'}+$error_allowed);
						$disc2 = 1 if ($ori2 == 1 and $tempdata[7] <= ($pos2+abs($homology)/2+$error_allowed) and $tempdata[7] >= $pos2-$is{$rg}{'isu'}-$error_allowed);
						$disc2 = 1 if ($ori2 == -1 and $tempdata[7] >= ($pos2-abs($homology)/2-$error_allowed) and $tempdata[7] <= $pos2+$is{$rg}{'isu'}+$error_allowed);
						if ($disc1 and $disc2)
						{
							$discord++;
						}
					}
					unless ($ori1 == 1 and $ori2 == -1)
					{
						$disc1 = 1 if ($ori1 == 1 and $tempdata[3] <= ($pos1+abs($homology)/2+$error_allowed) and $tempdata[3] >= $pos1-$is{$rg}{'isu'}-$error_allowed);
						$disc1 = 1 if ($ori1 == -1 and $tempdata[3] >= ($pos1-abs($homology)/2-$error_allowed) and $tempdata[3] <= $pos1+$is{$rg}{'isu'}+$error_allowed);
						$disc2 = 1 if ($ori2 == 1 and $tempdata[7] <= ($pos2+abs($homology)/2+$error_allowed) and $tempdata[7] >= $pos2-$is{$rg}{'isu'}-$error_allowed);
						$disc2 = 1 if ($ori2 == -1 and $tempdata[7] >= ($pos2-abs($homology)/2-$error_allowed) and $tempdata[7] <= $pos2+$is{$rg}{'isu'}+$error_allowed);
						#print "@tempdata\n$pos1\t$ori1\t$pos2\t$ori2\t$disc1\t$disc2\n";
						if ($disc1 and $disc2)
						{
							$discord++;
						}
					}
				}
			}
			else
			{
				if ($tempdata[6] eq $chr2)
				{
					$disc1 = 1 if ($ori1 == 1 and $tempdata[3] <= ($pos1+abs($homology)/2+$error_allowed) and $tempdata[3] >= $pos1-$is{$rg}{'isu'}-$error_allowed);
					$disc1 = 1 if ($ori1 == -1 and $tempdata[3] >= ($pos1-abs($homology)/2-$error_allowed) and $tempdata[3] <= $pos1+$is{$rg}{'isu'}+$error_allowed);
					$disc2 = 1 if ($ori2 == 1 and $tempdata[7] <= ($pos2+abs($homology)/2+$error_allowed) and $tempdata[7] >= $pos2-$is{$rg}{'isu'}-$error_allowed);
					$disc2 = 1 if ($ori2 == -1 and $tempdata[7] >= ($pos2-abs($homology)/2-$error_allowed) and $tempdata[7] <= $pos2+$is{$rg}{'isu'}+$error_allowed);
					if ($disc1 and $disc2)
					{
						$discord++;
					}
				}
			}
		}
	}
	close TEMP;
	
	return $discord;
}

sub print_usage
{
	die "discon.pl [options]
	-i FILE	input file, required
	-o FILE	output file, required
	-D INT	count number of discordant pairs from bam, [0/1], turn this function on for genotyping, default 0
	-B FILE	bam file, required
	-C FILE	cluster file generated by Meerkat, required
	-I FILE	isinfo file from Meerkat run, required
	-K FILE	file name of read group to be ignored, one read group ID per line, same as R option in meerkat.pl
	-S STR	/path/to/samtools, path only, not the command, no need to specify if samtools is in PATH
	-d FLT	standard deviation cutoff to call discordant read pairs, default 3
	-Q INT	minimum mapping quality for reads to be used, default 0
	-h help\n";
}