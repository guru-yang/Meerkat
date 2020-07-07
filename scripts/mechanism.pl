#!/usr/bin/perl
# mechanism.pl
# assign mechanism to each variant
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
# Email: lixing_yang@hms.harvard.edu
# Input: sorted bam file, repeat file
# Output:
# TEI: transposable element insertion, single or multiple TE insertion.
# TEA: Alternative TE, usually a deletion with insertion event is called, deletion is a TE in reference genome and insertion is the same type of TE in reference genome. It should be arisen from sequence divergence of TE.
# VNTR: variable number of tandem repeat, deletion or insertion of satellite repeat, simple repeat or low complexity repeat.
# NAHR: non-allellic homologous recombination, >100bp homology.
# alt-EJ: alternative end joining, 3-100bp homology.
# NHEJ: non-homologous end joining, 0-2bp homology or 1-10bp insertion at deletion break point.
# FoSTeS: fork stalling and template switching, template switch, >10bp insertion at deletion break points.
# NA: unclassified.
# i.e. perl scripts/mechanism.pl -R /db/hg18/rmsk-hg18.txt -b

use strict;
use Getopt::Std;
use FindBin '$Bin';
my $version = 'v.0.174';

my %opts = (o=>1, t=>100000, z=>1000000000);
getopts("b:o:z:R:h", \%opts);

my $bamfile = $opts{b};
my $include_other = $opts{o};
my $rmskfile = $opts{R};
my $te_size_max = $opts{t};
my $sv_size_cutoff = $opts{z};
my $del_ins_size_cutoff_d = 0.8; # size ratio of del and ins in del_ins events
my $del_ins_size_cutoff_u = 1.2;
my $ovl = 0.8; # overlap of a predicted events and a annotated TE

&print_usage unless (defined($opts{b}) and defined($opts{R}));
&print_usage if (defined($opts{h}));
my $time0 = time;
my $local0 = localtime($time0);
print STDERR "$local0 Mechanism $version started\n";

my $prefix = $bamfile;
if ($prefix =~ /.bam$/)
{
	$prefix = $`;
}
if ($prefix =~ /.sorted$/)
{
	$prefix = $`;
}
my $intra_refine_type = $prefix.'.intra.refined.typ.sorted';
my $inter_refine_type = $prefix.'.inter.refined.typ.sorted';
my $outfile = $prefix.'.variants';
#$intra_refine_type = 'sim/mechanism/chr10.intra';# need to remove
#$inter_refine_type = 'sim/mechanism/chr10.inter';# need to remove

my (%te, %tei, %sr, %sri);
my $newline;
open RMSK, "<$rmskfile";
while ($newline = <RMSK>)
{
	chomp $newline;
	my @data = split (/\t/, $newline);
	if ($data[11] eq 'LINE' or $data[11] eq 'SINE' or $data[11] eq 'LTR' or $data[11] eq 'DNA')
	{
		$tei{$data[5]} = 0 unless ($tei{$data[5]});
		$te{$data[5]}[$tei{$data[5]}][0] = $data[6]; # start
		$te{$data[5]}[$tei{$data[5]}][1] = $data[7]; # end
		$te{$data[5]}[$tei{$data[5]}][2] = $data[11];# class
		$te{$data[5]}[$tei{$data[5]}][3] = $data[10];# name
		#$te{$data[5]}[$tei{$data[5]}][2] = $data[9]; # strand
		#print "$te{$data[5]}[$tei{$data[5]}][0]\t$te{$data[5]}[$tei{$data[5]}][1]\t$te{$data[5]}[$tei{$data[5]}][2]\t$te{$data[5]}[$tei{$data[5]}][3]\t$te{$data[5]}[$tei{$data[5]}][4]\t$te{$data[5]}[$tei{$data[5]}][5]\n";
		$tei{$data[5]}++;
	}
	if ($include_other and $data[11] eq 'Other')
	{
		$tei{$data[5]} = 0 unless ($tei{$data[5]});
		$te{$data[5]}[$tei{$data[5]}][0] = $data[6]; # start
		$te{$data[5]}[$tei{$data[5]}][1] = $data[7]; # end
		$te{$data[5]}[$tei{$data[5]}][2] = $data[11];# class
		$te{$data[5]}[$tei{$data[5]}][3] = $data[10];# name
		#$te{$data[5]}[$tei{$data[5]}][2] = $data[9]; # strand
		#print "$te{$data[5]}[$tei{$data[5]}][0]\t$te{$data[5]}[$tei{$data[5]}][1]\t$te{$data[5]}[$tei{$data[5]}][2]\t$te{$data[5]}[$tei{$data[5]}][3]\t$te{$data[5]}[$tei{$data[5]}][4]\t$te{$data[5]}[$tei{$data[5]}][5]\n";
		$tei{$data[5]}++;
	}
	if ($data[11] eq 'Satellite' or $data[11] eq 'Simple_repeat' or $data[11] eq 'Low_complexity')
	{
		$sri{$data[5]} = 0 unless ($sri{$data[5]});
		$sr{$data[5]}[$sri{$data[5]}][0] = $data[6]; # start
		$sr{$data[5]}[$sri{$data[5]}][1] = $data[7]; # end
		$sr{$data[5]}[$sri{$data[5]}][2] = $data[11];# class
		#$sr{$data[5]}[$sri{$data[5]}][3] = $data[12];# family
		#$sr{$data[5]}[$sri{$data[5]}][4] = $data[9]; # strand
		#$sr{$data[5]}[$sri{$data[5]}][5] = $data[10];# name
		#print "$sr{$data[5]}[$sri{$data[5]}][0]\t$sr{$data[5]}[$sri{$data[5]}][1]\t$sr{$data[5]}[$sri{$data[5]}][2]\t$sr{$data[5]}[$sri{$data[5]}][3]\t$sr{$data[5]}[$sri{$data[5]}][4]\t$sr{$data[5]}[$sri{$data[5]}][5]\n";
		$sri{$data[5]}++;
	}
}
close RMSK;

my @variant;
my $i = 0;
open FILE, "<$intra_refine_type";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split (/\t/, $newline);
	$variant[$i] = \@data;
	#print "@{$variant[$i]}\n";
	$i++;
}
close FILE;
open FILE, "<$inter_refine_type";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split (/\t/, $newline);
	$variant[$i] = \@data;
	$i++;
}
close FILE;

open OUT, ">$outfile";
LOOP: foreach (@variant)
{
	my @data = @$_;
	my $type = $data[0];
	my ($mechanism, @bp_annotation);
	#print "@data\n";
	# call mechanism
	if ($type eq 'del')
	{
		# large events
		if ($data[7] > $sv_size_cutoff)
		{
			;
		}
		else
		{
			my ($te_class, $te_name) = &tei($data[4], $data[5], $data[6]);
			$mechanism = 'TEI_'.$te_class.'_'.$te_name if ($te_class);
			$mechanism = 'TEI_complex' if ($te_class eq 'complex');
			goto LEND if ($mechanism);
			my $vntr = &vntr($data[4], $data[5], $data[6]);
			$mechanism = 'VNTR' if ($vntr);
			goto LEND if ($mechanism);
			$mechanism = 'NAHR' if ($data[8] > 100);
			goto LEND if ($mechanism);
			$mechanism = 'alt-EJ' if ($data[8] >= 2 and $data[8] <= 100);
			goto LEND if ($mechanism);
			$mechanism = 'NHEJ';
			goto LEND if ($mechanism);
		}
	}
	elsif ($type =~ /ins/)
	{
		# large events
		if ($data[7] > $sv_size_cutoff or $data[11] > $sv_size_cutoff)
		{
			;
		}
		else
		{
			# deletion with insertion inside the break points
			if ($type =~ /del/)
			{
				# known source of insertion
				if ($data[8] ne '-')
				{
					my ($te_class_d, $te_name_d) = &tei($data[4], $data[5], $data[6]);
					my ($te_class_i, $te_name_i) = &tei($data[8], $data[9], $data[10]);
					$mechanism = 'TEA_'.$te_class_d.'_'.$te_class_i if ($te_class_d and $te_class_i);
					goto LEND if ($mechanism);
					$mechanism = 'TEI_'.$te_class_i.'_'.$te_name_i if ($te_class_i);
					$mechanism = 'TEI_complex' if ($te_class_i eq 'complex');
					goto LEND if ($mechanism);
					my $vntr_d = &vntr($data[4], $data[5], $data[6]);
					my $vntr_i = &vntr($data[8], $data[9], $data[10]);
					$mechanism = 'VNTR' if ($vntr_d and $vntr_i);
					goto LEND if ($mechanism);
					$mechanism = 'NHEJ' if ($data[11] <= 10);
					goto LEND if ($mechanism);
				}
				# unknown source of insertion
				else
				{
					my ($te_class, $te_name) = &tei($data[4], $data[5], $data[6]);
					$mechanism = 'TEI_'.$te_class.'_'.$te_name if ($te_class);
					$mechanism = 'TEI_complex' if ($te_class eq 'complex');
					goto LEND if ($mechanism);
					my $vntr = &vntr($data[4], $data[5], $data[6]);
					$mechanism = 'VNTR' if ($vntr);
					goto LEND if ($mechanism);
					$mechanism = 'NHEJ' if ($data[11] <= 10);
					goto LEND if ($mechanism);
				}
				$mechanism = 'FoSTeS' if ($data[11] > 10);
				goto LEND if ($mechanism);
			}
			# insertion
			else
			{
				my ($te_class_i, $te_name_i) = &tei($data[8], $data[9], $data[10]);
				$mechanism = 'TEI_'.$te_class_i.'_'.$te_name_i if ($te_class_i);
				$mechanism = 'TEI_complex' if ($te_class_i eq 'complex');
				goto LEND if ($mechanism);
				my $vntr_i = &vntr($data[8], $data[9], $data[10]);
				$mechanism = 'VNTR' if ($vntr_i);
				goto LEND if ($mechanism);
			}
		}
	}
	elsif ($type =~ /invers/)
	{
		# del_invers
		if ($type =~ /del/)
		{
			$mechanism = 'FoSTeS';
		}
		# invers_f, invers_r and invers
		else
		{
			;
		}
	}
	elsif ($type eq 'tandem_dup')
	{
		my $vntr = &vntr($data[4], $data[5], $data[6]);
		$mechanism = 'VNTR' if ($vntr);
		goto LEND if ($mechanism);
	}
	elsif ($type eq 'transl_inter')
	{
		;
	}
LEND:	$mechanism = 'NA' unless ($mechanism);
	shift(@data);
	unshift(@data, $mechanism);
	unshift(@data, $type);
	
	# annotate break points
	if ($type eq 'del' or $type eq 'tandem_dup' or $type eq 'invers' or $type eq 'invers_f' or $type eq 'invers_r')
	{
		if (&sr_overlap($data[5], $data[6]))
		{
			$bp_annotation[0] = 'SR';
		}
		elsif (&te_overlap($data[5], $data[6]))
		{
			$bp_annotation[0] = 'TE';
		}
		else
		{
			$bp_annotation[0] = '';
		}
		if (&sr_overlap($data[5], $data[7]))
		{
			$bp_annotation[1] = 'SR';
		}
		elsif (&te_overlap($data[5], $data[7]))
		{
			$bp_annotation[1] = 'TE';
		}
		else
		{
			$bp_annotation[1] = '';
		}
	}
	elsif ($type =~ /ins/ or $type eq 'del_invers')
	{
		# known source of insertion
		if ($data[9] ne '-')
		{
			if (&sr_overlap($data[5], $data[6]))
			{
				$bp_annotation[0] = 'SR';
			}
			elsif (&te_overlap($data[5], $data[6]))
			{
				$bp_annotation[0] = 'TE';
			}
			else
			{
				$bp_annotation[0] = '';
			}
			if (&sr_overlap($data[5], $data[7]))
			{
				$bp_annotation[1] = 'SR';
			}
			elsif (&te_overlap($data[5], $data[7]))
			{
				$bp_annotation[1] = 'TE';
			}
			else
			{
				$bp_annotation[1] = '';
			}
			if (&sr_overlap($data[9], $data[10]))
			{
				$bp_annotation[2] = 'SR';
			}
			elsif (&te_overlap($data[9], $data[10]))
			{
				$bp_annotation[2] = 'TE';
			}
			else
			{
				$bp_annotation[2] = '';
			}
			if (&sr_overlap($data[9], $data[11]))
			{
				$bp_annotation[3] = 'SR';
			}
			elsif (&te_overlap($data[9], $data[11]))
			{
				$bp_annotation[3] = 'TE';
			}
			else
			{
				$bp_annotation[3] = '';
			}
		}
		# unknown source of insertion
		else
		{
			if (&sr_overlap($data[5], $data[6]))
			{
				$bp_annotation[0] = 'SR';
			}
			elsif (&te_overlap($data[5], $data[6]))
			{
				$bp_annotation[0] = 'TE';
			}
			else
			{
				$bp_annotation[0] = '';
			}
			if (&sr_overlap($data[5], $data[7]))
			{
				$bp_annotation[1] = 'SR';
			}
			elsif (&te_overlap($data[5], $data[7]))
			{
				$bp_annotation[1] = 'TE';
			}
			else
			{
				$bp_annotation[1] = '';
			}
		}
	}
	elsif ($type eq 'transl_inter')
	{
		if (&sr_overlap($data[5], $data[6]))
		{
			$bp_annotation[0] = 'SR';
		}
		elsif (&te_overlap($data[5], $data[6]))
		{
			$bp_annotation[0] = 'TE';
		}
		else
		{
			$bp_annotation[0] = '';
		}
		if (&sr_overlap($data[8], $data[9]))
		{
			$bp_annotation[1] = 'SR';
		}
		elsif (&te_overlap($data[8], $data[9]))
		{
			$bp_annotation[1] = 'TE';
		}
		else
		{
			$bp_annotation[1] = '';
		}
	}
	my $bp_annotation = join ("_", @bp_annotation);
	$bp_annotation = 'BP:'.$bp_annotation;
	push @data, $bp_annotation;
	my $toprint = join ("\t", @data);
	print OUT "$toprint\n";
}
close OUT;

my $time1 = time;
my $local1 = localtime($time1);
my $difference = $time1 - $time0;
my $seconds    =  $difference % 60;
$difference = ($difference - $seconds) / 60;
my $minutes    =  $difference % 60;
$difference = ($difference - $minutes) / 60;
print STDERR "$local1 Finished\n";
print STDERR "Time used: $difference:$minutes:$seconds\n";

# overlap satellite repeat, simple repeat, low complexity repeat
sub sr_overlap
{
	my $chr = shift;
	my $bp = shift;
	foreach (@{$sr{$chr}})
	{
		next if ($$_[0] < $bp and $$_[1] < $bp);
		last if ($$_[0] > $bp and $$_[1] > $bp);
		if ($bp >= $$_[0] and $bp <= $$_[1])
		{
			return 1;
		}
	}
	return 0;
}

# overlap TE
sub te_overlap
{
	my $chr = shift;
	my $bp = shift;
	foreach (@{$te{$chr}})
	{
		next if ($$_[0] < $bp and $$_[1] < $bp);
		last if ($$_[0] > $bp and $$_[1] > $bp);
		if ($bp >= $$_[0] and $bp <= $$_[1])
		{
			return 1;
		}
	}
	return 0;	
}

sub vntr
{
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	return 0 if (abs($end-$start) > $te_size_max);
	my @overlap;
	foreach (@{$sr{$chr}})
	{
		my $sr = $_;
		next if ($$sr[0] < $end and $$sr[1] < $end and $$sr[0] < $start and $$sr[1] < $start);
		last if ($$sr[0] > $end and $$sr[1] > $end and $$sr[0] > $start and $$sr[1] > $start);
		if (&covered($start, $end, $$sr[0], $$sr[1]))
		{
			my ($overlap1, $overlap2) = &overlap($start, $end, $$sr[0], $$sr[1]);
			if (abs($overlap2-$overlap1)/(abs($end-$start)+0.5)>=$ovl)
			{
				return 1;
			}
		}
	}
	return 0;	
}

# TE insertion only
sub tei
{
	my $chr = shift;
	my $start = shift;
	my $end = shift;
	return 0 if (abs($end-$start) > $te_size_max);
	my @overlap;
	foreach (@{$te{$chr}})
	{
		my $te = $_;
		next if ($$te[0] < $end and $$te[1] < $end and $$te[0] < $start and $$te[1] < $start);
		last if ($$te[0] > $end and $$te[1] > $end and $$te[0] > $start and $$te[1] > $start);
		if (&covered($start, $end, $$te[0], $$te[1]))
		{
			my ($overlap1, $overlap2) = &overlap($start, $end, $$te[0], $$te[1]);
			if (abs($overlap2-$overlap1)/(abs($$te[1]-$$te[0])+0.5)>=$ovl)
			{
				if (abs($overlap2-$overlap1)/(abs($end-$start)+0.5)>=$ovl)
				{
					return ($$te[2], $$te[3]);
				}
				else
				{
					push @overlap, $overlap2-$overlap1;
				}
			}
		}
	}
	if ($overlap[0])
	{
		my $overlap;
		foreach (@overlap)
		{
			$overlap += $_;
		}
		if ($overlap/(abs($end-$start)+0.5)>=$ovl)
		{
			return ('complex', 1);
		}
	}
	return 0;	
}

sub overlap
{
	my $a1 = shift;
	my $a2 = shift;
	my $b1 = shift;
	my $b2 = shift;
	return ($a1, $a2) unless ($b1 and $b2);
	($a1, $a2) = sort { $a <=> $b } ($a1, $a2);
	($b1, $b2) = sort { $a <=> $b } ($b1, $b2);
	if (($a1 >= $b1 and $a1 <= $b2) or ($a2 >= $b1 and $a2 <= $b2) or ($b1 >= $a1 and $b1 <= $a2) or ($b2 >= $a1 and $b2 <= $a2))
	{
		my $a = ($a1>$b1)?$a1:$b1;
		my $b = ($a2>$b2)?$b2:$a2;
		return ($a, $b);
	}
	else
	{
		return (0, 0);
	}
}

sub covered
{
	my $a1 = shift;
	my $a2 = shift;
	my $b1 = shift;
	my $b2 = shift;
	return 0 if ($a1 =~ /\D/ or $a2 =~ /\D/ or $b1 =~ /\D/ or $b2 =~ /\D/);
	($a1, $a2) = sort{ $a <=> $b } ($a1, $a2);
	($b1, $b2) = sort { $a <=> $b } ($b1, $b2);
	if (($a1 >= $b1 and $a1 <= $b2) or ($a2 >= $b1 and $a2 <= $b2) or ($b1 >= $a1 and $b1 <= $a2) or ($b2 >= $a1 and $b2 <= $a2))
	{
		return 1;
	}
}

sub print_usage
{
	die "mechanism.pl [options]
	-b FILE	sorted and indexed bam file, required
	-o INT	[0/1], include rmsk type \"Other\" in TE, default 1
	-t INT	max size of TE, default 100,000
	-z INT	size limit of SVs to be processed, default 1,000,000,000
	-R STR	/path/to/repeat_mask_file, required, can be downloaded from UCSC
	-h help\n";
}
