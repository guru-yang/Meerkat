#!/usr/bin/perl
# meerkat.pl
# Identify SVs by read pairs and split reads, give precise break points by local alignment
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
# Email: lixing_yang@hms.harvard.edu
# System requirements:
# samtools 0.1.5 or above
# BWA 0.5.7 or above
# NCBI blast 2.2.10 or above
# i.e. perl scripts/meerkat.pl -F /db/hg18/hg18_fasta/ -W /opt/bwa/ -B /opt/blast/bin/ -b

use strict;
use Getopt::Std;
use Bio::DB::Fasta;
use FindBin '$Bin';
my $version = 'v.0.189';

my %opts = (k=>1, d=>3, p=>2, o=>0, q=>1, z=>1000000000, s=>20, m=>1, a=>1, u=>0, Q=>0, l=>1, t=>1, R=>'null', P=>'all');
getopts("b:k:d:c:p:o:q:z:s:m:a:u:Q:g:f:l:t:R:F:S:W:B:P:n:h", \%opts);

my $bamfile = $opts{b};
my $blacklist = $opts{k};
my $sd_cutoff_disc = $opts{d};
my $sd_cutoff_cl = $opts{c};
my $support_mps = $opts{p};
my $support_mpf = $opts{o};
my $support_reads = $opts{q};
my $sv_size_cutoff = $opts{z};
my $cut_sr = $opts{s};
my $remove_dup = $opts{m};
my $ad_align = $opts{a};
my $use_all_align = $opts{u};
my $min_mapq = $opts{Q};
my $alt_map_max = $opts{g};
my $alt_map_max_clip = $opts{f};
my $clip = $opts{l};
my $threads_bwa = $opts{t};
my $black_rg = $opts{R};
my $reference_path = $opts{F};
my $samtools_path = $opts{S};
my $bwa_path = $opts{W};
my $blastall_path = $opts{B};
my $step = $opts{P};
my $line = $opts{n};
$ad_align = 0 if $use_all_align;
$sd_cutoff_cl = $sd_cutoff_disc unless (defined($opts{c}));
if (defined($opts{h}))
{
	&print_usage;
	die;
}
die "bam file not specified\n" unless (defined($bamfile));
if (defined($step))
{
	die "wrong step specified\nmust be one of the following: all|dc|cl|mpd|alg|srd|rf, default all\n" unless ($step =~ /all|dc|cl|mpd|alg|srd|ft|rf/);
}
#	control progress
my $discordant =1 if ($step eq 'all' or $step eq 'dc'); # extract discordant read pairs
my $cluster = 	1 if ($step eq 'all' or $step =~ /cl/); # generate discordant read pair clusters
my $mpd = 	1 if ($step eq 'all' or $step eq 'mpd'); # call candidates from discordant clusters
my $alg = 	1 if ($step eq 'all' or $step =~ /alg/); # construct search space based on mpd candidates and align split reads to search space
my $srd = 	1 if ($step eq 'all' or $step eq 'srd'); # call SVs from pseudo split reads
my $filter =	1 if ($step eq 'all' or $step eq 'srd' or $step eq 'ft'); # filter outputs
my $blast_bp = 	1 if ($step eq 'all' or $step eq 'rf'); # blast break points reads to break points regions
my $cl =        0;
$cl =           1 if ($step eq 'cl1');
$cl =           2 if ($step eq 'cl2');
$cl =           3 if ($step eq 'cl3');
my $algs =        0;
$algs =           1 if ($step eq 'alg1');
$algs =           2 if ($step eq 'alg2');

die "path to reference fasta files not specified\n" if ($alg and !(defined($opts{F})));

my $prefix = substr ($bamfile, 0, -4);
if ($prefix =~ /.sorted$/)
{
	$prefix = $`;
}
my $blacklistfile = $prefix.'.blacklist.gz' if ($opts{k});

my ($newline, %is);
# parse median and standard deviation of insert size from file
my $isinfofile = $prefix.'.isinfo';
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
		if ($is{'rlu'})
		{
			$is{'rlu'} = $is{$rg}{'rl'} if ($is{$rg}{'rl'} > $is{'rlu'});
		}
		else
		{
			$is{'rlu'} = $is{$rg}{'rl'};
		}
		if ($is{'rld'})
		{
			$is{'rld'} = $is{$rg}{'rl'} if ($is{$rg}{'rl'} < $is{'rld'});
		}
		else
		{
			$is{'rld'} = $is{$rg}{'rl'};
		}
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
		$is{$rg}{'isu'} = $is{$rg}{'median'} + $is{$rg}{'sd'}*$sd_cutoff_cl;
		$is{$rg}{'isu'} = int($is{$rg}{'isu'}) + 1;
		$is{$rg}{'isd'} = $is{$rg}{'median'} - $is{$rg}{'sd'}*$sd_cutoff_cl;
		$is{$rg}{'isd'} = int($is{$rg}{'isd'}) - 1;
		if ($is{'isu'})
		{
			$is{'isu'} = $is{$rg}{'isu'} if ($is{$rg}{'isu'} > $is{'isu'});
		}
		else
		{
			$is{'isu'} = $is{$rg}{'isu'};
		}
		#print "$rg\t$is{$rg}{'rl'}\t$is{$rg}{'median'}\t$is{$rg}{'sd'}\t$is{$rg}{'isu'}\t$is{$rg}{'isd'}\n";
	}
}

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
my $bwa_command;
if (defined($opts{W}))
{
	$bwa_path .= '/' unless ($bwa_path =~ /\/$/);
	$bwa_command = $bwa_path.'bwa';
}
else
{
	$bwa_command = 'bwa';
}
my $formatdb_command;
if (defined($opts{B}))
{
	$blastall_path .= '/' unless ($blastall_path =~ /\/$/);
	$formatdb_command = $blastall_path.'formatdb';
}
else
{
	$formatdb_command = 'formatdb';
}
my $blastall_command;
if (defined($opts{B}))
{
	$blastall_path .= '/' unless ($blastall_path =~ /\/$/);
	$blastall_command = $blastall_path.'blastall';
}
else
{
	$blastall_command = 'blastall';
}
$Bin =~ /scripts/;
my $bin_path = $`.'bin/';
my $dre_command = $bin_path.'dre';
my $scluster_command = $bin_path.'sclus';

my $reference_db = Bio::DB::Fasta->new($reference_path) if ($alg or $blast_bp);

my $time0 = time;
my $local0 = localtime($time0);
print STDERR "$local0 Meerkat $version started\n";

&discord ($bamfile, $prefix, $blacklistfile, $clip, $dre_command, $sd_cutoff_disc, $isinfofile, $samtools_command, $remove_dup, $black_rg) if ($discordant);
&cluster ($prefix, $blacklistfile, \%is, $cut_sr, $ad_align, $alt_map_max, $alt_map_max_clip, $clip, $samtools_command, $cl, $scluster_command) if ($cluster);
&mpd (\%is, $prefix, $support_mps, $support_mpf) if ($mpd);
my $time = time;
my $local = localtime($time);
print STDERR "$local called events by read pairs\n" if ($mpd);
&alg($line, $prefix, $bwa_command, $reference_db, \%is, $cut_sr, $threads_bwa) if ($alg);
my $time = time;
my $local = localtime($time);
print STDERR "$local aligned split reads\n" if ($alg);
&srd($line, $prefix, \%is, $cut_sr, $sv_size_cutoff, $support_reads) if ($srd);
my $time = time;
my $local = localtime($time);
print STDERR "$local confirmed events by split reads\n" if ($srd);
&filter($prefix) if ($filter);
my $time = time;
my $local = localtime($time);
print STDERR "$local filtered events\n" if ($srd);
&blast_bp($line, $prefix, $reference_db, \%is, $cut_sr, $formatdb_command, $blastall_command) if ($blast_bp);
my $time = time;
my $local = localtime($time);
print STDERR "$local refined break points by local alignments\n" if ($blast_bp);

my $time1 = time;
my $difference = $time1 - $time0;
my $seconds    =  $difference % 60;
$difference = ($difference - $seconds) / 60;
my $minutes    =  $difference % 60;
$difference = ($difference - $minutes) / 60;
print STDERR "Time used: $difference:$minutes:$seconds\n";


#	blast break points reads to break points regions
sub blast_bp
{
#	file format of prefix.sr.intra.refined
#	del, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, homology sizes
#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion and insertion
#	del_invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of inversion (2 col), inversion size, distance of inverstion and deletion (2 col), homology at break points
#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, rchr (donor), ange of insertion (2 col), insert size, distance of deletion and insertion, homology at break points
#	invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size, homology at break points
#	invers_*, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size, homology at break points
#	tandem_dup, cluster id, number of supporting read pairs, number of supporting split reads, chr, tandem duplication boundary 1, tandem duplication boundary 2, tandem duplication size, homology at break points

#	file format of prefix.sr.inter.refine
#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr of insertion donor, range of insertion (2 col), insert size, homology at break points
#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of insert site, insert site size, chr of insertion donor, range of insertion (2 col), insert size, homology at break points
#	transl_inter, cluster id, number of supporting read pairs, number of supporting split reads, chr of 1st cluster, boundary of 1st cluster, orientation of 1st cluster, chr of 2nd cluster, boundary of 2nd cluster, orientation of 2nd cluster, homology at break points
	my $line = shift;
	my $prefix = shift;
	my $db = shift;
	my $ref_is = shift;
	my $cut_sr = shift;
	my $formatdb_command = shift;
	my $blastall_command = shift;
	my $unmappedfile = $prefix.'.unmapped.fq.gz';
	my $softclipfile = $prefix.'.softclips.fq.gz';
	my $intra_filter = $prefix.'.sr.intra.filtered';
	my $inter_filter = $prefix.'.sr.inter.filtered';
	my $intra_refine = $prefix.'.intra.refined';
	my $inter_refine = $prefix.'.inter.refined';
	my $blastdir = $prefix.'.blast/';
	my $srdir = $prefix.'.sr/';
	system "rm $blastdir*.*" if (-e $blastdir);
	system "rm $srdir*.*" if (-e $srdir);
	system "mkdir $blastdir" unless (-e $blastdir);
	system "mkdir $srdir" unless (-e $srdir);
	my $bp_readsfile = $prefix.'.bp_reads';
	
#	generate query sequences
	my (%readslist, %readsdetail);
#	%readslist: reads list for each break point pairs. $readslist{readname} = breakpoint pairs name
#	%readsdetail: seq of each read. $readsdetail{breakpoint name}{readname} = seq
	my $ref_readslist = \%readslist;
	my $ref_readsdetail = \%readsdetail;
	open BPREADSRF, "<$bp_readsfile";
	while ($newline = <BPREADSRF>)
	{
		chomp $newline;
		my $name = $newline;
		$newline = <BPREADSRF>;
		chomp $newline;
		my @reads_name = split (/\t/, $newline);
		foreach (@reads_name)
		{
			$readslist{$_} = $name;
		}
	}
	close BPREADSRF;
	
	open(UMFQRF, "gunzip -c $unmappedfile |");
	while(<UMFQRF>)
	{
		chomp;
		my $readname = substr($_, 1);
		$newline = <UMFQRF>;
		chomp $newline;
		my $length = length($newline);
		$readname .= '_'.$length;
		if (defined($readslist{$readname}))
		{
			$readsdetail{$readslist{$readname}}{$readname} = $newline;
			#print "$readslist{$readname}\t$readname\t$readsdetail{$readslist{$readname}}{$readname}\n";
		}
		$newline = <UMFQRF>;
		$newline = <UMFQRF>;
	}
	close UMFQRF;
	
	open(SCFQRF, "gunzip -c $softclipfile |");
	while(<SCFQRF>)
	{
		chomp;
		my $readname = substr($_, 1);
		$newline = <SCFQRF>;
		chomp $newline;
		my $length = length($newline);
		$readname .= '_'.$length;
		if (defined($readslist{$readname}))
		{
			$readsdetail{$readslist{$readname}}{$readname} = $newline;
			#print "$readslist{$readname}\t$readname\t$readsdetail{$readslist{$readname}}{$readname}\n";
		}
		$newline = <SCFQRF>;
		$newline = <SCFQRF>;
	}
	close SCFQRF;
	
	my $n;
	for (my $i=0;$i<100;$i++)
	{
		$n .= 'N';
	}
	
#	generate sbjct sequences for intra chr events
	my $i = 0;
	my @inss;
	open SRINTRAFTRF, "<$intra_filter";
	open INTRAREFINE, ">$intra_refine";
	while ($newline = <SRINTRAFTRF>)
	{
		chomp $newline;
		$i++;
		next if (defined($line) and $i < $line);
		my @data = split (/\t/, $newline);
		my @mpd_id = split (/\//, $data[1]);
		my $ref_data = \@data;
		if ($data[0] eq 'del')
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $position3 = $data[6] - $$ref_is{'rlu'};
			my $position4 = $data[6] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position3 => $position4);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[6];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast, $right_bound_blast);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
			}
			$left_bound_blast = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast = $data[6] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $eventsizedel = $right_bound_blast - $left_bound_blast - 1;
			if ($eventsizedel > 0)
			{
				if ($$ref_bps[2]>0)
				{
					print INTRAREFINE "del_ins\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast\t$right_bound_blast\t$eventsizedel\t-\t-\t-\t$$ref_bps[2]\n";
				}
				else
				{
					my $homology = -$$ref_bps[2];
					$homology = 0 unless ($$ref_bps[2]);
					print INTRAREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast\t$right_bound_blast\t$eventsizedel\t$homology\n";
				}
			}
		}
		if ($data[0] =~ 'inssu')
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $position3 = $data[9] - $$ref_is{'rlu'};
			my $position4 = $data[9] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position3 => $position4);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[9];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 0, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast1, $right_bound_blast1);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
			}
			$left_bound_blast1 = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast1 = $data[9] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology1 = -$$ref_bps[2];
			$homology1 = 0 unless ($$ref_bps[2]);
			
			my $position1 = $data[6] - $$ref_is{'rlu'};
			my $position2 = $data[6] + $$ref_is{'rlu'};
			my $position3 = $data[10] - $$ref_is{'rlu'};
			my $position4 = $data[10] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position3 => $position4);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[6].'__'.$data[4].'__'.$data[10];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[1], $name, $seq, $ref_readsdetail, 0, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast2, $right_bound_blast2);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast2 = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast2 = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
			}
			$left_bound_blast2 = $data[6] unless ($$ref_bps[0]);
			$right_bound_blast2 = $data[10] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology2 = -$$ref_bps[2];
			$homology2 = 0 unless ($$ref_bps[2]);
			my $eventsizedel = $left_bound_blast2 - $left_bound_blast1 - 1;
			my $eventsizeins = $right_bound_blast2 - $right_bound_blast1 + 1;
			my $distance = $left_bound_blast1 - $right_bound_blast2;
			if ($eventsizedel <= 10)
			{
				print INTRAREFINE "inssu\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast1\t$right_bound_blast2\t$eventsizeins\t$distance\t$homology1/$homology2\n";
			}
			else
			{
				print INTRAREFINE "del_inssu\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast1\t$right_bound_blast2\t$eventsizeins\t$distance\t$homology1/$homology2\n";
			}
		}
		if ($data[0] =~ 'inssd')
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $position3 = $data[9] - $$ref_is{'rlu'};
			my $position4 = $data[9] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position3 => $position4);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[9];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[1], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast1, $right_bound_blast1);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
			}
			#print "$ref_bps\n";
			$left_bound_blast1 = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast1 = $data[9] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology2 = -$$ref_bps[2];
			$homology2 = 0 unless ($$ref_bps[2]);
			
			my $position1 = $data[6] - $$ref_is{'rlu'};
			my $position2 = $data[6] + $$ref_is{'rlu'};
			my $position3 = $data[10] - $$ref_is{'rlu'};
			my $position4 = $data[10] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position3 => $position4);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[6].'__'.$data[4].'__'.$data[10];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast2, $right_bound_blast2);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast2 = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast2 = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
			}
			#print "$ref_bps\n";
			$left_bound_blast2 = $data[6] unless ($$ref_bps[0]);
			$right_bound_blast2 = $data[10] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology1 = -$$ref_bps[2];
			$homology1 = 0 unless ($$ref_bps[2]);
			my $eventsizedel = $left_bound_blast2 - $left_bound_blast1 - 1;
			my $eventsizeins = $right_bound_blast2 - $right_bound_blast1 + 1;
			my $distance = $right_bound_blast1 - $left_bound_blast2;
			if ($eventsizedel <= 10)
			{
				print INTRAREFINE "inssd\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast1\t$right_bound_blast2\t$eventsizeins\t$distance\t$homology1/$homology2\n";
			}
			else
			{
				print INTRAREFINE "del_inssd\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast1\t$right_bound_blast2\t$eventsizeins\t$distance\t$homology1/$homology2\n";
			}
		}
		if ($data[0] =~ 'insou')
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $position3 = $data[10] - $$ref_is{'rlu'};
			my $position4 = $data[10] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position4 => $position3);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[10];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 0, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast1, $right_bound_blast1);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			}
			$left_bound_blast1 = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast1 = $data[10] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology1 = -$$ref_bps[2];
			$homology1 = 0 unless ($$ref_bps[2]);
			
			my $position1 = $data[6] - $$ref_is{'rlu'};
			my $position2 = $data[6] + $$ref_is{'rlu'};
			my $position3 = $data[9] - $$ref_is{'rlu'};
			my $position4 = $data[9] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position4 => $position3);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[6].'__'.$data[4].'__'.$data[9];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[1], $name, $seq, $ref_readsdetail, 0, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast2, $right_bound_blast2);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast2 = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast2 = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			}
			$left_bound_blast2 = $data[6] unless ($$ref_bps[0]);
			$right_bound_blast2 = $data[9] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology2 = -$$ref_bps[2];
			$homology2 = 0 unless ($$ref_bps[2]);
			my $eventsizedel = $left_bound_blast2 - $left_bound_blast1 - 1;
			my $eventsizeins = $right_bound_blast1 - $right_bound_blast2 + 1;
			my $distance = $left_bound_blast1 - $right_bound_blast1;
			if ($eventsizedel <= 10)
			{
				if ($eventsizedel <= 2 and abs($eventsizeins) <= 2)
				{
					if ($right_bound_blast1 > $left_bound_blast2)
					{
						my $eventsizeinv = $right_bound_blast1 - $left_bound_blast2;
						print INTRAREFINE "invers\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast2\t$right_bound_blast1\t$eventsizeinv\t$homology1/$homology2\n";
					}
					else
					{
						my $eventsizeinv = $left_bound_blast1 - $right_bound_blast2;
						print INTRAREFINE "invers\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$right_bound_blast2\t$left_bound_blast1\t$eventsizeinv\t$homology1/$homology2\n";
					}
				}
				else
				{
					print INTRAREFINE "insou\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast2\t$right_bound_blast1\t$eventsizeins\t$distance\t$homology1/$homology2\n";
				}
			}
			else
			{
				print INTRAREFINE "del_insou\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast2\t$right_bound_blast1\t$eventsizeins\t$distance\t$homology1/$homology2\n";
			}
		}
		if ($data[0] =~ 'insod')
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $position3 = $data[10] - $$ref_is{'rlu'};
			my $position4 = $data[10] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position4 => $position3);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[10];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast1, $right_bound_blast1);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			}
			$left_bound_blast1 = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast1 = $data[10] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology1 = -$$ref_bps[2];
			$homology1 = 0 unless ($$ref_bps[2]);
			
			my $position1 = $data[6] - $$ref_is{'rlu'};
			my $position2 = $data[6] + $$ref_is{'rlu'};
			my $position3 = $data[9] - $$ref_is{'rlu'};
			my $position4 = $data[9] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position4 => $position3);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[6].'__'.$data[4].'__'.$data[9];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[1], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast2, $right_bound_blast2);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast2 = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast2 = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			}
			$left_bound_blast2 = $data[6] unless ($$ref_bps[0]);
			$right_bound_blast2 = $data[9] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology2 = -$$ref_bps[2];
			$homology2 = 0 unless ($$ref_bps[2]);
			my $eventsizedel = $left_bound_blast2 - $left_bound_blast1 - 1;
			my $eventsizeins = $right_bound_blast1 - $right_bound_blast2 + 1;
			my $distance = $right_bound_blast2 - $left_bound_blast2;
			if ($eventsizedel <= 10)
			{
				if ($eventsizedel <= 2 and abs($eventsizeins) <= 2)
				{
					if ($right_bound_blast1 > $left_bound_blast2)
					{
						my $eventsizeinv = $right_bound_blast1 - $left_bound_blast2;
						print INTRAREFINE "invers\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast2\t$right_bound_blast1\t$eventsizeinv\t$homology1/$homology2\n";
					}
					else
					{
						my $eventsizeinv = $left_bound_blast1 - $right_bound_blast2;
						print INTRAREFINE "invers\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$right_bound_blast2\t$left_bound_blast1\t$eventsizeinv\t$homology1/$homology2\n";
					}
				}
				else
				{
					print INTRAREFINE "insod\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast2\t$right_bound_blast1\t$eventsizeins\t$distance\t$homology1/$homology2\n";
				}
			}
			else
			{
				print INTRAREFINE "del_insod\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast2\t$right_bound_blast1\t$eventsizeins\t$distance\t$homology1/$homology2\n";
			}
		}
		if ($data[0] eq 'invers')
		{
			my $position1 = $data[5] - ($$ref_is{'rlu'} + $cut_sr);
			my $position2 = $data[5] + ($$ref_is{'rlu'} + $cut_sr);
			my $position3 = $data[6] - ($$ref_is{'rlu'} + $cut_sr);
			my $position4 = $data[6] + ($$ref_is{'rlu'} + $cut_sr);
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position4 => $position3);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[6];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast1, $right_bound_blast1, $left_bound_blast2, $right_bound_blast2);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position1 + $$ref_bps[1] - 1;
				$left_bound_blast2 = $position1 + $$ref_bps2[0] - 1;
				$right_bound_blast2 = $position1 + $$ref_bps2[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
				$left_bound_blast2 = $position1 + $$ref_bps2[0] - 1;
				$right_bound_blast2 = $position4 - ($$ref_bps2[1] - ($position2-$position1+100)) + 2;
			}
			#print "@{$ref_bps}\n@{$ref_bps2}\n";
			$left_bound_blast1 = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast1 = $data[6] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology = -$$ref_bps[2];
			$homology = 0 unless ($$ref_bps[2]);
			if ($$ref_bps2[1] and abs($$ref_bps2[0]-$$ref_bps[0])>10 and abs($$ref_bps2[1]-$$ref_bps[1])>10)
			{
				if ($left_bound_blast1 < $left_bound_blast2)
				{
					my $eventsizedel = $right_bound_blast2 - $left_bound_blast1 - 1;
					my $eventsizeinv = $right_bound_blast1 - $left_bound_blast2;
					my $distance1 = $left_bound_blast2 - $left_bound_blast1;
					my $distance2 = $right_bound_blast2 - $right_bound_blast1;
					print INTRAREFINE "del_invers\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$right_bound_blast2\t$eventsizedel\t$data[4]\t$left_bound_blast2\t$right_bound_blast1\t$eventsizeinv\t$distance1\t$distance2\t$homology\n";
				}
				else
				{
					my $eventsizedel = $right_bound_blast1 - $left_bound_blast2 - 1;
					my $eventsizeinv = $right_bound_blast2 - $left_bound_blast1;
					my $distance1 = $left_bound_blast1 - $left_bound_blast2;
					my $distance2 = $right_bound_blast1 - $right_bound_blast2;
					print INTRAREFINE "del_invers\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast2\t$right_bound_blast1\t$eventsizedel\t$data[4]\t$left_bound_blast1\t$right_bound_blast2\t$eventsizeinv\t$distance1\t$distance2\t$homology\n";
				}
			}
			else
			{
				my $eventsizeinv = $right_bound_blast1 - $left_bound_blast1;
				print INTRAREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$right_bound_blast1\t$eventsizeinv\t$homology\n";
			}
		}
		if ($data[0] eq 'del_invers')
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $position3 = $data[10] - $$ref_is{'rlu'};
			my $position4 = $data[10] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position4 => $position3);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[10];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast1, $right_bound_blast1);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast1 = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			}
			$left_bound_blast1 = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast1 = $data[10] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology1 = -$$ref_bps[2];
			$homology1 = 0 unless ($$ref_bps[2]);
			
			my $position1 = $data[6] - $$ref_is{'rlu'};
			my $position2 = $data[6] + $$ref_is{'rlu'};
			my $position3 = $data[9] - $$ref_is{'rlu'};
			my $position4 = $data[9] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position4 => $position3);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position3 => $position2);
			}
			my $name = $data[4].'__'.$data[6].'__'.$data[4].'__'.$data[9];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[1], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast2, $right_bound_blast2);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position3 + $$ref_bps[1] - 1;
				$right_bound_blast2 = $position3 + $$ref_bps[0] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast2 = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			}
			$left_bound_blast2 = $data[6] unless ($$ref_bps[0]);
			$right_bound_blast2 = $data[9] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology2 = -$$ref_bps[2];
			$homology2 = 0 unless ($$ref_bps[2]);
			my $eventsizedel = $left_bound_blast2 - $left_bound_blast1 - 1;
			my $eventsizeinv = $right_bound_blast1 - $right_bound_blast2;
			my $distance1 = $right_bound_blast2 - $left_bound_blast1;
			my $distance2 = $left_bound_blast2 - $right_bound_blast1;
			#print "$left_bound_blast1\t$left_bound_blast2\t$right_bound_blast1\t$right_bound_blast2\n$eventsizedel\t$eventsizeinv\t$distance1\t$distance2\n";
			if ($distance1<=10 and $distance2<=10)
			{
				print INTRAREFINE "invers\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$homology1/$homology2\n";
			}
			elsif ($eventsizedel <= 10 and $eventsizeinv <= 10)
			{
				my $eventsizeinv = $right_bound_blast1 - $left_bound_blast1;
				print INTRAREFINE "invers\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$right_bound_blast1\t$eventsizeinv\t$homology1/$homology2\n";
			}
			elsif ($distance1 < -10)
			{
				my $eventsizedel = $left_bound_blast2 - $right_bound_blast1 - 1;
				my $eventsizeins = $left_bound_blast1 - $right_bound_blast2 + 1;
				my $distance = $right_bound_blast1 - $left_bound_blast1;
				if ($eventsizedel <= 10)
				{
					print INTRAREFINE "insou\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$right_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast2\t$left_bound_blast1\t$eventsizeins\t$distance\t$homology1/$homology2\n";
				}
				else
					{
					print INTRAREFINE "del_insou\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$right_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast2\t$left_bound_blast1\t$eventsizeins\t$distance\t$homology1/$homology2\n";
				}
			}
			elsif ($distance2 < -10)
			{
				my $eventsizedel = $right_bound_blast2 - $left_bound_blast1 - 1;
				my $eventsizeins = $right_bound_blast1 - $left_bound_blast2 + 1;
				my $distance = $left_bound_blast2 - $right_bound_blast2;
				if ($eventsizedel <= 10)
				{
					print INTRAREFINE "insod\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$right_bound_blast2\t$eventsizedel\t$data[4]\t$left_bound_blast2\t$right_bound_blast1\t$eventsizeins\t$distance\t$homology1/$homology2\n";
				}
				else
				{
					print INTRAREFINE "del_insod\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$right_bound_blast2\t$eventsizedel\t$data[4]\t$left_bound_blast2\t$right_bound_blast1\t$eventsizeins\t$distance\t$homology1/$homology2\n";
				}
			}
			else
			{
				print INTRAREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[4]\t$right_bound_blast2\t$right_bound_blast1\t$eventsizeinv\t$distance1\t$distance2\t$homology1/$homology2\n";
			}
		}
		if ($data[0] eq 'tandem_dup')
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $position3 = $data[6] - $$ref_is{'rlu'};
			my $position4 = $data[6] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position3 => $position4);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[6];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast, $right_bound_blast);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
			}
			$left_bound_blast = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast = $data[6] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology = -$$ref_bps[2];
			$homology = 0 unless ($$ref_bps[2]);
			my $eventsizetransl = $right_bound_blast - $left_bound_blast;
			print INTRAREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast\t$right_bound_blast\t$eventsizetransl\t$homology\n";
		}
		if ($data[0] eq 'invers_f')
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $position3 = $data[6] - $$ref_is{'rlu'};
			my $position4 = $data[6] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position4 => $position3);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[6];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast, $right_bound_blast);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			}
			$left_bound_blast = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast = $data[6] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology = -$$ref_bps[2];
			$homology = 0 unless ($$ref_bps[2]);
			my $eventsizeinv = $right_bound_blast - $left_bound_blast;
			print INTRAREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast\t$right_bound_blast\t$eventsizeinv\t$homology\n";
		}
		if ($data[0] eq 'invers_r')
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $position3 = $data[6] - $$ref_is{'rlu'};
			my $position4 = $data[6] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			$seq .= $db->seq($data[4], $position4 => $position3);
			if (&covered($position1, $position2, $position3, $position4))
			{
				$seq = $db->seq($data[4], $position1 => $position4);
			}
			my $name = $data[4].'__'.$data[5].'__'.$data[4].'__'.$data[6];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2, $left_bound_blast, $right_bound_blast);
			if (&covered($position1, $position2, $position3, $position4))
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast = $position1 + $$ref_bps[1] - 1;
			}
			else
			{
				($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				$left_bound_blast = $position1 + $$ref_bps[0] - 1;
				$right_bound_blast = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			}
			$left_bound_blast = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast = $data[6] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology = -$$ref_bps[2];
			$homology = 0 unless ($$ref_bps[2]);
			my $eventsizeinv = $right_bound_blast - $left_bound_blast;
			print INTRAREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast\t$right_bound_blast\t$eventsizeinv\t$homology\n";
		}
		last if (defined($line) and $i == $line);
	}
	close SRINTRAFTRF;	
	close INTRAREFINE;
	
#	generate sbjct sequences for inter chr events
	my $i = 0;
	open SRINTERFTRF, "<$inter_filter";
	open INTERREFINE, ">$inter_refine";
	while ($newline = <SRINTERFTRF>)
	{
		chomp $newline;
		$i++;
		next if (defined($line) and $i < $line);
		my @data = split (/\t/, $newline);
		my @mpd_id = split (/\//, $data[1]);
		my $ref_data = \@data;
		if ($data[0] =~ /inss/)
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			my $position3 = $data[9] - $$ref_is{'rlu'};
			my $position4 = $data[9] + $$ref_is{'rlu'};
			$seq .= $db->seq($data[8], $position3 => $position4);
			my $name = $data[4].'__'.$data[5].'__'.$data[8].'__'.$data[9];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
			my $left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
			my $right_bound_blast1 = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
			$left_bound_blast1 = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast1 = $data[9] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology1 = -$$ref_bps[2];
			$homology1 = 0 unless ($$ref_bps[2]);
			
			my $position1 = $data[6] - $$ref_is{'rlu'};
			my $position2 = $data[6] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			my $position3 = $data[10] - $$ref_is{'rlu'};
			my $position4 = $data[10] + $$ref_is{'rlu'};
			$seq .= $db->seq($data[8], $position3 => $position4);
			my $name = $data[4].'__'.$data[6].'__'.$data[8].'__'.$data[10];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[1], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
			my $left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
			my $right_bound_blast2 = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
			$left_bound_blast2 = $data[6] unless ($$ref_bps[0]);
			$right_bound_blast2 = $data[10] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology2 = -$$ref_bps[2];
			$homology2 = 0 unless ($$ref_bps[2]);
			my $eventsizedel = $left_bound_blast2 - $left_bound_blast1 - 1;
			my $eventsizeins = $right_bound_blast2 - $right_bound_blast1 + 1;
			if ($eventsizedel <= 10)
			{
				print INTERREFINE "inss\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[8]\t$right_bound_blast1\t$right_bound_blast2\t$eventsizeins\t$homology1/$homology2\n";
			}
			else
			{
				print INTERREFINE "del_inss\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[8]\t$right_bound_blast1\t$right_bound_blast2\t$eventsizeins\t$homology1/$homology2\n";
			}
		}
		if ($data[0] =~ /inso/)
		{
			my $position1 = $data[5] - $$ref_is{'rlu'};
			my $position2 = $data[5] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			my $position3 = $data[10] - $$ref_is{'rlu'};
			my $position4 = $data[10] + $$ref_is{'rlu'};
			$seq .= $db->seq($data[8], $position4 => $position3);
			my $name = $data[4].'__'.$data[5].'__'.$data[8].'__'.$data[10];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
			my $left_bound_blast1 = $position1 + $$ref_bps[0] - 1;
			my $right_bound_blast1 = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			$left_bound_blast1 = $data[5] unless ($$ref_bps[0]);
			$right_bound_blast1 = $data[10] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology1 = -$$ref_bps[2];
			$homology1 = 0 unless ($$ref_bps[2]);
			
			my $position1 = $data[6] - $$ref_is{'rlu'};
			my $position2 = $data[6] + $$ref_is{'rlu'};
			my $seq = $db->seq($data[4], $position1 => $position2);
			$seq .= $n;
			my $position3 = $data[9] - $$ref_is{'rlu'};
			my $position4 = $data[9] + $$ref_is{'rlu'};
			$seq .= $db->seq($data[8], $position4 => $position3);
			my $name = $data[4].'__'.$data[6].'__'.$data[8].'__'.$data[9];
			my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[1], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
			my ($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
			my $left_bound_blast2 = $position1 + $$ref_bps[0] - 1;
			my $right_bound_blast2 = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
			$left_bound_blast2 = $data[6] unless ($$ref_bps[0]);
			$right_bound_blast2 = $data[9] unless ($$ref_bps[1]);
			$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
			my $homology2 = -$$ref_bps[2];
			$homology2 = 0 unless ($$ref_bps[2]);
			my $eventsizedel = $left_bound_blast2 - $left_bound_blast1 - 1;
			my $eventsizeins = $right_bound_blast1 - $right_bound_blast2 + 1;
			if ($eventsizedel <= 10)
			{
				print INTERREFINE "inso\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[8]\t$right_bound_blast2\t$right_bound_blast1\t$eventsizeins\t$homology1/$homology2\n";
			}
			else
			{
				print INTERREFINE "del_inso\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast1\t$left_bound_blast2\t$eventsizedel\t$data[8]\t$right_bound_blast2\t$right_bound_blast1\t$eventsizeins\t$homology1/$homology2\n";
			}
		}
		if ($data[0] eq 'transl_inter')
		{
			if ($data[6] == 1 and $data[9] == -1)
			{
				my $position1 = $data[5] - $$ref_is{'rlu'};
				my $position2 = $data[5] + $$ref_is{'rlu'};
				my $seq = $db->seq($data[4], $position1 => $position2);
				$seq .= $n;
				my $position3 = $data[8] - $$ref_is{'rlu'};
				my $position4 = $data[8] + $$ref_is{'rlu'};
				$seq .= $db->seq($data[7], $position3 => $position4);
				my $name = $data[4].'__'.$data[5].'__'.$data[7].'__'.$data[8];
				my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
				my ($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				my $left_bound_blast = $position1 + $$ref_bps[0] - 1;
				my $right_bound_blast = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
				$left_bound_blast = $data[5] unless ($$ref_bps[0]);
				$right_bound_blast = $data[8] unless ($$ref_bps[1]);
				$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
				my $homology = -$$ref_bps[2];
				$homology = 0 unless ($$ref_bps[2]);
				print INTERREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast\t$data[6]\t$data[7]\t$right_bound_blast\t$data[9]\t$homology\n";
			}
			if ($data[6] == -1 and $data[9] == 1)
			{
				my $position1 = $data[5] - $$ref_is{'rlu'};
				my $position2 = $data[5] + $$ref_is{'rlu'};
				my $seq = $db->seq($data[4], $position1 => $position2);
				$seq .= $n;
				my $position3 = $data[8] - $$ref_is{'rlu'};
				my $position4 = $data[8] + $$ref_is{'rlu'};
				$seq .= $db->seq($data[7], $position3 => $position4);
				my $name = $data[4].'__'.$data[5].'__'.$data[7].'__'.$data[8];
				my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
				my ($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				my $left_bound_blast = $position1 + $$ref_bps[0] - 1;
				my $right_bound_blast = $position3 + ($$ref_bps[1] - ($position2-$position1+100)) - 2;
				$left_bound_blast = $data[5] unless ($$ref_bps[0]);
				$right_bound_blast = $data[8] unless ($$ref_bps[1]);
				$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
				my $homology = -$$ref_bps[2];
				$homology = 0 unless ($$ref_bps[2]);
				print INTERREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast\t$data[6]\t$data[7]\t$right_bound_blast\t$data[9]\t$homology\n";
			}
			if ($data[6] == 1 and $data[9] == 1)
			{
				my $position1 = $data[5] - $$ref_is{'rlu'};
				my $position2 = $data[5] + $$ref_is{'rlu'};
				my $seq = $db->seq($data[4], $position1 => $position2);
				$seq .= $n;
				my $position3 = $data[8] - $$ref_is{'rlu'};
				my $position4 = $data[8] + $$ref_is{'rlu'};
				$seq .= $db->seq($data[7], $position4 => $position3);
				my $name = $data[4].'__'.$data[5].'__'.$data[7].'__'.$data[8];
				my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
				my ($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				my $left_bound_blast = $position1 + $$ref_bps[0] - 1;
				my $right_bound_blast = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
				$left_bound_blast = $data[5] unless ($$ref_bps[0]);
				$right_bound_blast = $data[8] unless ($$ref_bps[1]);
				$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
				my $homology = -$$ref_bps[2];
				$homology = 0 unless ($$ref_bps[2]);
				print INTERREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast\t$data[6]\t$data[7]\t$right_bound_blast\t$data[9]\t$homology\n";
			}
			if ($data[6] == -1 and $data[9] == -1)
			{
				my $position1 = $data[5] - $$ref_is{'rlu'};
				my $position2 = $data[5] + $$ref_is{'rlu'};
				my $seq = $db->seq($data[4], $position1 => $position2);
				$seq .= $n;
				my $position3 = $data[8] - $$ref_is{'rlu'};
				my $position4 = $data[8] + $$ref_is{'rlu'};
				$seq .= $db->seq($data[7], $position4 => $position3);
				my $name = $data[4].'__'.$data[5].'__'.$data[7].'__'.$data[8];
				my $bltoutfile = &run_blast($blastdir, $ref_data, $mpd_id[0], $name, $seq, $ref_readsdetail, 1, $formatdb_command, $blastall_command);
				my ($ref_bps, $ref_bps2) = &parse_blast($bltoutfile, $ref_is, $cut_sr, $ref_readslist, $ref_readsdetail, $support_reads);
				my $left_bound_blast = $position1 + $$ref_bps[0] - 1;
				my $right_bound_blast = $position4 - ($$ref_bps[1] - ($position2-$position1+100)) + 2;
				$left_bound_blast = $data[5] unless ($$ref_bps[0]);
				$right_bound_blast = $data[8] unless ($$ref_bps[1]);
				$$ref_bps[2] = -1000 unless ($$ref_bps[0]);
				my $homology = -$$ref_bps[2];
				$homology = 0 unless ($$ref_bps[2]);
				print INTERREFINE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$left_bound_blast\t$data[6]\t$data[7]\t$right_bound_blast\t$data[9]\t$homology\n";
			}
		}
		last if (defined($line) and $i == $line);
	}
	close SRINTERFTRF;
	close INTERREFINE;
	
	my $intra_refine_coord = $prefix.'.intra.refined.crd.sorted';
	my $inter_refine_coord = $prefix.'.inter.refined.crd.sorted';
	my $intra_refine_type = $prefix.'.intra.refined.typ.sorted';
	my $inter_refine_type = $prefix.'.inter.refined.typ.sorted';
	#system "sort $intra_refine -k 5,5 -k 6,6n > $intra_refine_coord";
	#system "sort $inter_refine -k 5,5 -k 6,6n > $inter_refine_coord";
	system "sort $intra_refine -k 1,1 -k 5,5 -k 6,6n > $intra_refine_type";
	system "sort $inter_refine -k 1,1 -k 5,5 -k 6,6n > $inter_refine_type";
	system "rm $intra_refine";
	system "rm $inter_refine";
	system "rm -rf $blastdir";
}

#	extract break points from blast output
sub parse_blast
{
	my $bltoutfile = shift;
	my $ref_is = shift;
	my $cut_sr = shift;
	my $ref_readslist = shift;
	my $ref_readsdetail = shift;
	my $support_reads = shift;
	my $bp1 = $$ref_is{'rlu'};
	my $bp2 = $$ref_is{'rlu'}*3+100;
	my $srfile = $bltoutfile;
	$srfile =~ s/blast/sr/;
	$srfile =~ s/bltout/srout/;
	my $sbjfafile = $bltoutfile;
	$sbjfafile =~ s/bltout/sbj.fa/;
	
	open SBJFA, "<$sbjfafile";
	$newline = <SBJFA>;
	$newline = <SBJFA>;
	chomp $newline;
	my @ref_seq = split (/N+/, $newline);
	#print "$ref_seq[0]\n$ref_seq[1]\n";
	open SROUT, ">$srfile";
	
	my (%freq_s1, %freq_s2, %freq_dist, %alg, $lastquery, $i);
#	%freq_s1: frequency of break points on one region
#	%freq_s2: frequency of break points on the other region
#	%freq_dist: frequency of query break points distance, positive: unmatched bases, negative: homology
#	%alg: alignment details
#	$alg{readname}{i}: reference to ith alignment detail
	open BLTOUT, "<$bltoutfile";
	while ($newline = <BLTOUT>)
	{
		chomp $newline;
		my ($query, $subject, $identity, $hitlength, $mismatch, $gap, $qrstart, $qrend, $sbjstart, $sbjend, $evalue, $score);
		my @data = split (/\t/, $newline);
		if ($lastquery eq $data[0])
		{
			$i++;
		}
		else
		{
			$i = 0;
			$lastquery = $data[0];
		}
		@{$alg{$data[0]}{$i}} = @data;
#		print "$i\t@{$alg{$data[0]}{$i}}\n";
	}
	close BLTOUT;
	
	my $ref_seq_ori_det1 = 0; # determine the orientation of 1st ref seq, 0 not yet determined, 1 determined already
	my $ref_seq_ori_det2 = 0;
	my $ref_seq_ori1 = 1; # orientation of ref seq, 1 forward, 0 reverse
	my $ref_seq_ori2 = 1;
	my @aligned_reads;
	foreach my $key (keys %alg)
	{
		my $value = $alg{$key};
		next unless (defined($alg{$key}{1}[0]));
		# only 2 entries
		unless (defined($alg{$key}{2}[0]))
		{
			# not process reads not entirely used (such read should belong to other break point, caused by overlapping regions produced by small events)
			next if (($alg{$key}{1}[6] >= $alg{$key}{0}[6] and $alg{$key}{1}[6] <= $alg{$key}{0}[7] and $alg{$key}{1}[7] >= $alg{$key}{0}[6] and $alg{$key}{1}[7] <= $alg{$key}{0}[7]) or ($alg{$key}{0}[6] >= $alg{$key}{1}[6] and $alg{$key}{0}[6] <= $alg{$key}{1}[7] and $alg{$key}{0}[7] >= $alg{$key}{1}[6] and $alg{$key}{0}[7] <= $alg{$key}{1}[7]));
			# skip if not the entire read is used
			#if ($alg{$key}{1}[6] > $alg{$key}{0}[6])
			#{
			#	next if (($alg{$key}{1}[6]-$alg{$key}{0}[7])>5 and !($is_del));
			#}
			#else
			#{
			#	next if (($alg{$key}{0}[6]-$alg{$key}{1}[7])>5 and !($is_del));
			#}
			# 1st half read on top
			if ($alg{$key}{0}[7] > $alg{$key}{0}[6] and $alg{$key}{0}[7] < $alg{$key}{1}[7])
			{
				if ($alg{$key}{0}[9] > $alg{$key}{1}[8])
				{
					$freq_s2{$alg{$key}{0}[9]}++;
					$freq_s1{$alg{$key}{1}[8]}++;
				}
				else
				{
					$freq_s1{$alg{$key}{0}[9]}++;
					$freq_s2{$alg{$key}{1}[8]}++;
				}
				my $dist = $alg{$key}{1}[6] - $alg{$key}{0}[7] - 1;
				$freq_dist{$dist}++;
				
				unless ($ref_seq_ori_det1)
				{
					if (($alg{$key}{0}[8] < $alg{$key}{1}[8] and $alg{$key}{0}[8] > $alg{$key}{0}[9]) or ($alg{$key}{0}[8] > $alg{$key}{1}[8] and $alg{$key}{0}[8] < $alg{$key}{0}[9]))
					{
						($ref_seq[0] = reverse $ref_seq[0]) =~ tr/gatcGATC/ctagCTAG/;
						$ref_seq_ori1 = 0;
					}
					$ref_seq_ori_det1 = 1;
				}
				unless ($ref_seq_ori_det2)
				{
					if (($alg{$key}{0}[8] < $alg{$key}{1}[8] and $alg{$key}{1}[8] > $alg{$key}{1}[9]) or ($alg{$key}{0}[8] > $alg{$key}{1}[8] and $alg{$key}{1}[8] < $alg{$key}{1}[9]))
					{
						($ref_seq[1] = reverse $ref_seq[1]) =~ tr/gatcGATC/ctagCTAG/;
						$ref_seq_ori2 = 0;
					}
					$ref_seq_ori_det2 = 1;
				}
				my ($space, $num_space, $trim_size);
				my $aligned_read = $$ref_readsdetail{$$ref_readslist{$alg{$key}{0}[0]}}{$alg{$key}{0}[0]};
				if ($alg{$key}{0}[8] < $alg{$key}{1}[8])
				{
					if ($ref_seq_ori1)
					{
						$num_space = ($alg{$key}{0}[8]>$alg{$key}{0}[9]?$alg{$key}{0}[9]:$alg{$key}{0}[8]) - 1;
					}
					else
					{
						$num_space = length($ref_seq[0]) - ($alg{$key}{0}[8]>$alg{$key}{0}[9]?$alg{$key}{0}[8]:$alg{$key}{0}[9]);
					}
					$trim_size = 1 - ($alg{$key}{0}[6]>$alg{$key}{0}[7]?$alg{$key}{0}[7]:$alg{$key}{0}[6]);
				}
				else
				{
					if ($ref_seq_ori1)
					{
						$num_space = ($alg{$key}{1}[8]>$alg{$key}{1}[9]?$alg{$key}{1}[9]:$alg{$key}{1}[8]) - 1;
						($aligned_read = reverse $aligned_read) =~ tr/gatcGATC/ctagCTAG/;
					}
					else
					{
						$num_space = length($ref_seq[0]) - ($alg{$key}{1}[8]>$alg{$key}{1}[9]?$alg{$key}{1}[8]:$alg{$key}{1}[9]);
						($aligned_read = reverse $aligned_read) =~ tr/gatcGATC/ctagCTAG/;
					}
					$trim_size = 1 - ($alg{$key}{0}[6]>$alg{$key}{0}[7]?$alg{$key}{0}[7]:$alg{$key}{0}[6]);
				}
				if ($trim_size*$num_space == 0)
				{
					$aligned_read = substr ($aligned_read, -$trim_size);
				}
				else
				{
					$aligned_read = substr ($aligned_read, -$trim_size);
					for (my $i=0; $i<-$trim_size; $i++)
					{
						$aligned_read = ' '.$aligned_read;
					}
				}
				for (my $i=0; $i<$num_space; $i++)
				{
					$aligned_read = ' '.$aligned_read;
				}
				push @aligned_reads, $aligned_read if ($aligned_read);
				#print "$alg{$key}{0}[0]\t$num_space\t$trim_size\t$ref_seq_ori1\t$ref_seq_ori2\n";
				#print "$alg{$key}{1}[6]\t$alg{$key}{0}[7]\t$dist\n";
			}
			# 2nd half read on top
			if ($alg{$key}{1}[7] > $alg{$key}{1}[6] and $alg{$key}{1}[7] < $alg{$key}{0}[7])
			{
				if ($alg{$key}{1}[9] > $alg{$key}{0}[8])
				{
					$freq_s2{$alg{$key}{1}[9]}++;
					$freq_s1{$alg{$key}{0}[8]}++;
				}
				else
				{
					$freq_s1{$alg{$key}{1}[9]}++;
					$freq_s2{$alg{$key}{0}[8]}++;
				}
				my $dist = $alg{$key}{0}[6] - $alg{$key}{1}[7] - 1;
				$freq_dist{$dist}++;
				
				unless ($ref_seq_ori_det1)
				{
					if (($alg{$key}{1}[8] < $alg{$key}{0}[8] and $alg{$key}{1}[8] > $alg{$key}{1}[9]) or ($alg{$key}{1}[8] > $alg{$key}{0}[8] and $alg{$key}{1}[8] < $alg{$key}{1}[9]))
					{
						($ref_seq[0] = reverse $ref_seq[0]) =~ tr/gatcGATC/ctagCTAG/;
						$ref_seq_ori1 = 0;
					}
					$ref_seq_ori_det1 = 1;
				}
				unless ($ref_seq_ori_det2)
				{
					if (($alg{$key}{1}[8] < $alg{$key}{0}[8] and $alg{$key}{0}[8] > $alg{$key}{0}[9]) or ($alg{$key}{1}[8] > $alg{$key}{0}[8] and $alg{$key}{0}[8] < $alg{$key}{0}[9]))
					{
						($ref_seq[1] = reverse $ref_seq[1]) =~ tr/gatcGATC/ctagCTAG/;
						$ref_seq_ori2 = 0;
					}
					$ref_seq_ori_det2 = 1;
				}
				my ($space, $num_space, $trim_size);
				my $aligned_read = $$ref_readsdetail{$$ref_readslist{$alg{$key}{1}[0]}}{$alg{$key}{1}[0]};
				if ($alg{$key}{1}[8] < $alg{$key}{0}[8])
				{
					if ($ref_seq_ori1)
					{
						$num_space = ($alg{$key}{1}[8]>$alg{$key}{1}[9]?$alg{$key}{1}[9]:$alg{$key}{1}[8]) - 1;
					}
					else
					{
						$num_space = length($ref_seq[0]) - ($alg{$key}{1}[8]>$alg{$key}{1}[9]?$alg{$key}{1}[8]:$alg{$key}{1}[9]);
					}
					$trim_size = 1 - ($alg{$key}{1}[6]>$alg{$key}{1}[7]?$alg{$key}{1}[7]:$alg{$key}{1}[6]);
				}
				else
				{
					if ($ref_seq_ori1)
					{
						$num_space = ($alg{$key}{0}[8]>$alg{$key}{0}[9]?$alg{$key}{0}[9]:$alg{$key}{0}[8]) - 1;
						($aligned_read = reverse $aligned_read) =~ tr/gatcGATC/ctagCTAG/;
					}
					else
					{
						$num_space = length($ref_seq[0]) - ($alg{$key}{0}[8]>$alg{$key}{0}[9]?$alg{$key}{0}[8]:$alg{$key}{0}[9]);
						($aligned_read = reverse $aligned_read) =~ tr/gatcGATC/ctagCTAG/;
					}
					$trim_size = 1 - ($alg{$key}{1}[6]>$alg{$key}{1}[7]?$alg{$key}{1}[7]:$alg{$key}{1}[6]);
				}
				if ($trim_size*$num_space == 0)
				{
					$aligned_read = substr ($aligned_read, -$trim_size);
				}
				else
				{
					$aligned_read = substr ($aligned_read, -$trim_size);
					for (my $i=0; $i<-$trim_size; $i++)
					{
						$aligned_read = ' '.$aligned_read;
					}
				}
				for (my $i=0; $i<$num_space; $i++)
				{
					$aligned_read = ' '.$aligned_read;
				}
				push @aligned_reads, $aligned_read if ($aligned_read);
				#print "$alg{$key}{0}[0]\t$num_space\t$trim_size\t$ref_seq_ori1\t$ref_seq_ori2\n";
				#print "$alg{$key}{0}[6]\t$alg{$key}{1}[7]\t$dist\n";
			}
		}
		# more than 2 entries
		else
		{
			# get top 2 hits
			my %alg_top2;
			my $maxk = keys %{$alg{$key}};
			# initialize first entry
			@{$alg_top2{$key}{0}} = @{$alg{$key}{0}};
			delete $alg{$key}{0};
			# get best hit for first entry
			for (my $k=1;$k<$maxk;$k++)
			{
				if ($alg_top2{$key}{0}[6] == $alg{$key}{$k}[6] and $alg_top2{$key}{0}[7] == $alg{$key}{$k}[7])
				{
					# forward strand
					if ($alg{$key}{$k}[9] > $alg{$key}{$k}[8])
					{
						# first region
						if (abs($alg{$key}{$k}[9] - $bp1) < abs($alg{$key}{$k}[8] - $bp2))
						{
							if (abs($alg{$key}{$k}[8]-$bp1)<abs($alg_top2{$key}{0}[8]-$bp1))
							{
								@{$alg_top2{$key}{0}} = @{$alg{$key}{$k}};
							}
						}
						# second region
						else
						{
							if (abs($alg{$key}{$k}[9]-$bp2)<abs($alg_top2{$key}{0}[9]-$bp2))
							{
								@{$alg_top2{$key}{0}} = @{$alg{$key}{$k}};
							}
						}
					}
					# reverse strand
					else
					{
						# first region
						if (abs($alg{$key}{$k}[9] - $bp1) < abs($alg{$key}{$k}[8] - $bp2))
						{
							if (abs($alg{$key}{$k}[9]-$bp1)<abs($alg_top2{$key}{0}[9]-$bp1))
							{
								@{$alg_top2{$key}{0}} = @{$alg{$key}{$k}};
							}
						}
						# second region
						else
						{
							if (abs($alg{$key}{$k}[8]-$bp2)<abs($alg_top2{$key}{0}[8]-$bp2))
							{
								@{$alg_top2{$key}{0}} = @{$alg{$key}{$k}};
							}
						}
					}
					delete $alg{$key}{$k};
				}
				else
				{
					last;
				}
			}
			# initialize second entry
			my $startk;
			for (my $k=1;$k<$maxk;$k++)
			{
				if ($alg{$key}{$k}[0])
				{
					next if (($alg{$key}{$k}[6] >= $alg_top2{$key}{0}[6] and $alg{$key}{$k}[6] <= $alg_top2{$key}{0}[7] and $alg{$key}{$k}[7] >= $alg_top2{$key}{0}[6] and $alg{$key}{$k}[7] <= $alg_top2{$key}{0}[7]) or ($alg_top2{$key}{0}[6] >= $alg{$key}{$k}[6] and $alg_top2{$key}{0}[6] <= $alg{$key}{$k}[7] and $alg_top2{$key}{0}[7] >= $alg{$key}{$k}[6] and $alg_top2{$key}{0}[7] <= $alg{$key}{$k}[7]));
					@{$alg_top2{$key}{1}} = @{$alg{$key}{$k}};
					delete $alg{$key}{$k};
					$startk = $k + 1;
					last;
				}
			}
			#print "@{$alg_top2{$key}{0}}\n@{$alg_top2{$key}{1}}\n";
			# get best hit for second entry
			for (my $k=$startk;$k<$maxk;$k++)
			{
				if ($alg_top2{$key}{1}[6] == $alg{$key}{$k}[6] and $alg_top2{$key}{1}[7] == $alg{$key}{$k}[7])
				{
					# forward strand
					if ($alg{$key}{$k}[9] > $alg{$key}{$k}[8])
					{
						# first region
						if (abs($alg{$key}{$k}[9] - $bp1) < abs($alg{$key}{$k}[8] - $bp2))
						{
							if (abs($alg{$key}{$k}[8]-$bp1)<abs($alg_top2{$key}{1}[8]-$bp1))
							{
								@{$alg_top2{$key}{1}} = @{$alg{$key}{$k}};
							}
						}
						# second region
						else
						{
							if (abs($alg{$key}{$k}[9]-$bp2)<abs($alg_top2{$key}{1}[9]-$bp2))
							{
								@{$alg_top2{$key}{1}} = @{$alg{$key}{$k}};
							}
						}
					}
					# reverse strand
					else
					{
						# first region
						if (abs($alg{$key}{$k}[9] - $bp1) < abs($alg{$key}{$k}[8] - $bp2))
						{
							if (abs($alg{$key}{$k}[9]-$bp1)<abs($alg_top2{$key}{1}[9]-$bp1))
							{
								@{$alg_top2{$key}{1}} = @{$alg{$key}{$k}};
							}
						}
						# second region
						else
						{
							if (abs($alg{$key}{$k}[8]-$bp2)<abs($alg_top2{$key}{1}[8]-$bp2))
							{
								@{$alg_top2{$key}{1}} = @{$alg{$key}{$k}};
							}
						}
					}
					delete $alg{$key}{$k};
				}
				else
				{
					last;
				}
			}
			#print "@{$alg_top2{$key}{0}}\n@{$alg_top2{$key}{1}}\n";
			
			# skip if not the entire read is used
			#if ($alg_top2{$key}{1}[6] > $alg_top2{$key}{0}[6])
			#{
			#	next if (($alg_top2{$key}{1}[6]-$alg_top2{$key}{0}[7])>5 and !($is_del));
			#}
			#else
			#{
			#	next if (($alg_top2{$key}{0}[6]-$alg_top2{$key}{1}[7])>5 and !($is_del));
			#}
			
			# 1st half read on top
			if ($alg_top2{$key}{0}[7] > $alg_top2{$key}{0}[6] and $alg_top2{$key}{0}[7] < $alg_top2{$key}{1}[7])
			{
				if ($alg_top2{$key}{0}[9] > $alg_top2{$key}{1}[8])
				{
					$freq_s2{$alg_top2{$key}{0}[9]}++;
					$freq_s1{$alg_top2{$key}{1}[8]}++;
				}
				else
				{
					$freq_s1{$alg_top2{$key}{0}[9]}++;
					$freq_s2{$alg_top2{$key}{1}[8]}++;
				}
				my $dist = $alg_top2{$key}{1}[6] - $alg_top2{$key}{0}[7] - 1;
				$freq_dist{$dist}++;
				
				unless ($ref_seq_ori_det1)
				{
					if (($alg_top2{$key}{0}[8] < $alg_top2{$key}{1}[8] and $alg_top2{$key}{0}[8] > $alg_top2{$key}{0}[9]) or ($alg_top2{$key}{0}[8] > $alg_top2{$key}{1}[8] and $alg_top2{$key}{0}[8] < $alg_top2{$key}{0}[9]))
					{
						($ref_seq[0] = reverse $ref_seq[0]) =~ tr/gatcGATC/ctagCTAG/;
						$ref_seq_ori1 = 0;
					}
					$ref_seq_ori_det1 = 1;
				}
				unless ($ref_seq_ori_det2)
				{
					if (($alg_top2{$key}{0}[8] < $alg_top2{$key}{1}[8] and $alg_top2{$key}{1}[8] > $alg_top2{$key}{1}[9]) or ($alg_top2{$key}{0}[8] > $alg_top2{$key}{1}[8] and $alg_top2{$key}{1}[8] < $alg_top2{$key}{1}[9]))
					{
						($ref_seq[1] = reverse $ref_seq[1]) =~ tr/gatcGATC/ctagCTAG/;
						$ref_seq_ori2 = 0;
					}
					$ref_seq_ori_det2 = 1;
				}
				my ($space, $num_space, $trim_size);
				my $aligned_read = $$ref_readsdetail{$$ref_readslist{$alg_top2{$key}{0}[0]}}{$alg_top2{$key}{0}[0]};
				if ($alg_top2{$key}{0}[8] < $alg_top2{$key}{1}[8])
				{
					if ($ref_seq_ori1)
					{
						$num_space = ($alg_top2{$key}{0}[8]>$alg_top2{$key}{0}[9]?$alg_top2{$key}{0}[9]:$alg_top2{$key}{0}[8]) - 1;
					}
					else
					{
						$num_space = length($ref_seq[0]) - ($alg_top2{$key}{0}[8]>$alg_top2{$key}{0}[9]?$alg_top2{$key}{0}[8]:$alg_top2{$key}{0}[9]);
					}
					$trim_size = 1 - ($alg_top2{$key}{0}[6]>$alg_top2{$key}{0}[7]?$alg_top2{$key}{0}[7]:$alg_top2{$key}{0}[6]);
				}
				else
				{
					if ($ref_seq_ori1)
					{
						$num_space = ($alg_top2{$key}{1}[8]>$alg_top2{$key}{1}[9]?$alg_top2{$key}{1}[9]:$alg_top2{$key}{1}[8]) - 1;
						($aligned_read = reverse $aligned_read) =~ tr/gatcGATC/ctagCTAG/;
					}
					else
					{
						$num_space = length($ref_seq[0]) - ($alg_top2{$key}{1}[8]>$alg_top2{$key}{1}[9]?$alg_top2{$key}{1}[8]:$alg_top2{$key}{1}[9]);
						($aligned_read = reverse $aligned_read) =~ tr/gatcGATC/ctagCTAG/;
					}
					$trim_size = 1 - ($alg_top2{$key}{0}[6]>$alg_top2{$key}{0}[7]?$alg_top2{$key}{0}[7]:$alg_top2{$key}{0}[6]);
				}
				if ($trim_size*$num_space == 0)
				{
					$aligned_read = substr ($aligned_read, -$trim_size);
				}
				else
				{
					$aligned_read = substr ($aligned_read, -$trim_size);
					for (my $i=0; $i<-$trim_size; $i++)
					{
						$aligned_read = ' '.$aligned_read;
					}
				}
				for (my $i=0; $i<$num_space; $i++)
				{
					$aligned_read = ' '.$aligned_read;
				}
				push @aligned_reads, $aligned_read if ($aligned_read);
				#print "$alg_top2{$key}{0}[0]\t$num_space\t$trim_size\t$ref_seq_ori1\t$ref_seq_ori2\n";
				#print "$alg_top2{$key}{1}[6]\t$alg_top2{$key}{0}[7]\t$dist\n";
			}
			# 2nd half read on top
			if ($alg_top2{$key}{1}[7] > $alg_top2{$key}{1}[6] and $alg_top2{$key}{1}[7] < $alg_top2{$key}{0}[7])
			{
				if ($alg_top2{$key}{1}[9] > $alg_top2{$key}{0}[8])
				{
					$freq_s2{$alg_top2{$key}{1}[9]}++;
					$freq_s1{$alg_top2{$key}{0}[8]}++;
				}
				else
				{
					$freq_s1{$alg_top2{$key}{1}[9]}++;
					$freq_s2{$alg_top2{$key}{0}[8]}++;
				}
				my $dist = $alg_top2{$key}{0}[6] - $alg_top2{$key}{1}[7] - 1;
				$freq_dist{$dist}++;
				
				unless ($ref_seq_ori_det1)
				{
					if (($alg_top2{$key}{1}[8] < $alg_top2{$key}{0}[8] and $alg_top2{$key}{1}[8] > $alg_top2{$key}{1}[9]) or ($alg_top2{$key}{1}[8] > $alg_top2{$key}{0}[8] and $alg_top2{$key}{1}[8] < $alg_top2{$key}{1}[9]))
					{
						($ref_seq[0] = reverse $ref_seq[0]) =~ tr/gatcGATC/ctagCTAG/;
						$ref_seq_ori1 = 0;
					}
					$ref_seq_ori_det1 = 1;
				}
				unless ($ref_seq_ori_det2)
				{
					if (($alg_top2{$key}{1}[8] < $alg_top2{$key}{0}[8] and $alg_top2{$key}{0}[8] > $alg_top2{$key}{0}[9]) or ($alg_top2{$key}{1}[8] > $alg_top2{$key}{0}[8] and $alg_top2{$key}{0}[8] < $alg_top2{$key}{0}[9]))
					{
						($ref_seq[1] = reverse $ref_seq[1]) =~ tr/gatcGATC/ctagCTAG/;
						$ref_seq_ori2 = 0;
					}
					$ref_seq_ori_det2 = 1;
				}
				my ($space, $num_space, $trim_size);
				my $aligned_read = $$ref_readsdetail{$$ref_readslist{$alg_top2{$key}{1}[0]}}{$alg_top2{$key}{1}[0]};
				if ($alg_top2{$key}{1}[8] < $alg_top2{$key}{0}[8])
				{
					if ($ref_seq_ori1)
					{
						$num_space = ($alg_top2{$key}{1}[8]>$alg_top2{$key}{1}[9]?$alg_top2{$key}{1}[9]:$alg_top2{$key}{1}[8]) - 1;
					}
					else
					{
						$num_space = length($ref_seq[0]) - ($alg_top2{$key}{1}[8]>$alg_top2{$key}{1}[9]?$alg_top2{$key}{1}[8]:$alg_top2{$key}{1}[9]);
					}
					$trim_size = 1 - ($alg_top2{$key}{1}[6]>$alg_top2{$key}{1}[7]?$alg_top2{$key}{1}[7]:$alg_top2{$key}{1}[6]);
				}
				else
				{
					if ($ref_seq_ori1)
					{
						$num_space = ($alg_top2{$key}{0}[8]>$alg_top2{$key}{0}[9]?$alg_top2{$key}{0}[9]:$alg_top2{$key}{0}[8]) - 1;
						($aligned_read = reverse $aligned_read) =~ tr/gatcGATC/ctagCTAG/;
					}
					else
					{
						$num_space = length($ref_seq[0]) - ($alg_top2{$key}{0}[8]>$alg_top2{$key}{0}[9]?$alg_top2{$key}{0}[8]:$alg_top2{$key}{0}[9]);
						($aligned_read = reverse $aligned_read) =~ tr/gatcGATC/ctagCTAG/;
					}
					$trim_size = 1 - ($alg_top2{$key}{1}[6]>$alg_top2{$key}{1}[7]?$alg_top2{$key}{1}[7]:$alg_top2{$key}{1}[6]);
				}
				if ($trim_size*$num_space == 0)
				{
					$aligned_read = substr ($aligned_read, -$trim_size);
				}
				else
				{
					$aligned_read = substr ($aligned_read, -$trim_size);
					for (my $i=0; $i<-$trim_size; $i++)
					{
						$aligned_read = ' '.$aligned_read;
					}
				}
				for (my $i=0; $i<$num_space; $i++)
				{
					$aligned_read = ' '.$aligned_read;
				}
				push @aligned_reads, $aligned_read if ($aligned_read);
				#print "$alg_top2{$key}{0}[0]\t$num_space\t$trim_size\t$ref_seq_ori1\t$ref_seq_ori2\n";
				#print "$alg_top2{$key}{0}[6]\t$alg_top2{$key}{1}[7]\t$dist\n";
			}
		}
	}
		
	my (@results, @results2);
	my $fqi = 0;
	foreach my $i (sort{$freq_s1{$b} <=> $freq_s1{$a}} keys %freq_s1)
	{
		$fqi++;
		if ($fqi == 1)
		{
			push @results, $i;
			next;
		}
		next if (abs($results[0]-$i)<=2);
		if ($freq_s1{$i} >= $support_reads)
		{
			push @results2, $i;
		}
		#print "s1 $i\t$freq_s1{$i}\n";
		last;
	}
	my $fqi = 0;
	foreach my $i (sort{$freq_s2{$b} <=> $freq_s2{$a}} keys %freq_s2)
	{
		$fqi++;
		next if (abs($results[1]-$i)<=2);
		if ($fqi == 1)
		{
			push @results, $i;
			next;
		}
		if ($freq_s2{$i} >= $support_reads)
		{
			push @results2, $i;
		}
		#print "s2 $i\t$freq_s2{$i}\n";
		last;
	}
	my @raw_results = @results;
	#print "@results\n";
	
	# adjust break points if query break points having overlap or gap
	# adjust break point
	foreach my $i (sort{$freq_dist{$b} <=> $freq_dist{$a}} keys %freq_dist)
	{
		#print "$i\t@results\n";
		if ($i<0)
		{
			if ($i%2)
			{
				$results[0] -= -$i/2+0.5;
				$results[1] += -$i/2-0.5;
			}
			else
			{
				$results[0] -= -$i/2;
				$results[1] += -$i/2;
			}
		}
		last;
	}
	# get largest homolgy
	my $max_freq;
	foreach my $i (sort{$freq_dist{$b} <=> $freq_dist{$a}} keys %freq_dist)
	{
		$max_freq = $freq_dist{$i};
		last;
	}
	my $ins_size; # positive number for insertion, negative number for homology
	foreach my $i (sort{$freq_dist{$b} <=> $freq_dist{$a}} keys %freq_dist)
	{
		#print "$i\t@results\n";
		if (!(defined($ins_size)) or $i < $ins_size)
		{
			$ins_size = $i;
		}
		last if ($freq_dist{$i} < $max_freq);
	}
	@results = sort {$a <=> $b} @results;
	@results2 = sort {$a <=> $b} @results2;
	push @results, int($ins_size);
	
	print SROUT "$ref_seq[0]\n\n";
	foreach (@aligned_reads)
	{
		print SROUT "$_\n";
	}
	my $pos1 = $raw_results[0] + 1;
	$pos1 = length($ref_seq[0]) - $pos1 + 3 unless ($ref_seq_ori1);
	my $pos2 = $raw_results[1] - 100 - length($ref_seq[0]);
	$pos2 = length($ref_seq[1]) - $pos2 + 1 unless ($ref_seq_ori2);
	$pos2 = $pos2 - $results[2];
	my $trim = $pos2 - $pos1;
	if ($trim > 0)
	{
		$ref_seq[1] = substr ($ref_seq[1], $trim);
	}
	else
	{
		for (my $i=0;$i<-$trim;$i++)
		{
			$ref_seq[1] = ' '.$ref_seq[1];
		}
	}
	print SROUT "\n$ref_seq[1]\n";
	close SROUT;
	#print "$bltoutfile\n@results\n@results2\n";
	#print "$ref_seq_ori1\t$ref_seq_ori2\t$pos1\t$pos2\t$trim\n@aligned_reads\n";
	return (\@results, \@results2);
}

#	run blast on intra chr deletions
sub run_blast
{
	my $blastdir = shift;
	my $ref_data = shift;
	my $mpd_id = shift; 
	my $name = shift;
	my $seq = shift;
	my $ref_readsdetail = shift;
	my $forward = shift;
	my $formatdb_command = shift;
	my $blastall_command = shift;
	
	my $sbjctfile = $blastdir.$mpd_id.'.sbj.fa';
	open SBJCTSEQ, ">$sbjctfile";
	print SBJCTSEQ ">$name\n$seq\n";
	close SBJCTSEQ;
			
	my $queryfile = $blastdir.$mpd_id.'.query.fa';
	my $bltoutfile = $blastdir.$mpd_id.'.bltout';
	open QRSEQ, ">$queryfile";
	if ($$ref_readsdetail{$name})
	{
		foreach my $key (keys %{$$ref_readsdetail{$name}})
		{
			my $value = $$ref_readsdetail{$name}{$key};
			print QRSEQ">$key\n$value\n";
		}
	}
	# flip the order of break points
	else
	{
		my @temp_bp = split (/__/, $name);
		($temp_bp[0], $temp_bp[1], $temp_bp[2], $temp_bp[3]) = ($temp_bp[2], $temp_bp[3], $temp_bp[0], $temp_bp[1]);
		my $alt_name = join ("__", @temp_bp);
		foreach my $key (keys %{$$ref_readsdetail{$alt_name}})
		{
			my $value = $$ref_readsdetail{$alt_name}{$key};
			print QRSEQ">$key\n$value\n";
		}
	}
	close QRSEQ;
	system "$formatdb_command -i $sbjctfile -p F";
	system "$blastall_command -p blastn -d $sbjctfile -i $queryfile -o $bltoutfile -m 8 -F F";
	return $bltoutfile;
}

#	filter outputs
sub filter
{
	my $prefix = shift;
	my $srintra_out = $prefix.'.sr.intra.out';
	my $srintra_filter = $prefix.'.sr.intra.filtered';
	my $srinter_out = $prefix.'.sr.inter.out';
	my $srinter_filter = $prefix.'.sr.inter.filtered';
	
	my (%data, %cluster_event_map, %eventsize);
	# $data{event id} = complete event string
	# $cluster_event_map{$cluster_id} = array of event ids
	# $eventsize{event id} = event size
	my $ftai=0;
	open SRINTRAOUTFT, "<$srintra_out";
	while ($newline = <SRINTRAOUTFT>)
	{
		chomp $newline;
		my @data = split (/\t/, $newline);
		if ($data[1] =~ /\//)
		{
			my ($a, $b) = ($`, $');#'
			if ($a =~ /_/)
			{
				$a = $`;
			}
			if ($b =~ /_/)
			{
				$b = $`;
			}
			push @{$cluster_event_map{$a}}, $ftai;
			push @{$cluster_event_map{$b}}, $ftai;
			$eventsize{$ftai} += abs($data[7]);
			$eventsize{$ftai} += abs($data[11]);
		}
		else
		{
			my $a = $data[1];
			if ($a =~ /_/)
			{
				$a = $`;
			}
			push @{$cluster_event_map{$a}}, $ftai;
			$eventsize{$ftai} += abs($data[7]);
		}
		$data{$ftai} = $newline;
		$ftai++;
	}
	close SRINTRAOUTFT;
	
	foreach my $cluster_id (%cluster_event_map)
	{
		if ($cluster_event_map{$cluster_id})
		{
			if (@{$cluster_event_map{$cluster_id}} > 1)
			{
				#print "$cluster_id\t@{$cluster_event_map{$cluster_id}}\n";
				my ($select_id, @delete_id);
				$select_id = $cluster_event_map{$cluster_id}[0];
				foreach (@{$cluster_event_map{$cluster_id}})
				{
					if ($eventsize{$_}<$eventsize{$select_id})
					{
						$select_id = $_;
					}
				}
				foreach (@{$cluster_event_map{$cluster_id}})
				{
					 if ($_ ne $select_id)
					 {
					 	push @delete_id, $_;
					 }
				}
				foreach (@delete_id)
				{
					my @current_event = split (/\t/, $data{$_});
					if ($current_event[1] =~ /\//)
					{
						if ($current_event[0] =~ /inssd/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my ($srd1, $srd2) = split (/\//, $current_event[3]);
							my @temp_event;
							$temp_event[0] = 'tandem_dup';
							$temp_event[1] = $cl_id1;
							$temp_event[2] = $mpd1;
							$temp_event[3] = $srd1;
							$temp_event[4] = $current_event[4];
							$temp_event[5] = $current_event[6];
							$temp_event[6] = $current_event[10];
							$temp_event[7] = $temp_event[6] - $temp_event[5];
							my $modified_event = join("\t", @temp_event);
							$data{$_} = $modified_event;
							#print "$data{$_}\n";
						}
						if ($current_event[0] =~ /inssu/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my ($srd1, $srd2) = split (/\//, $current_event[3]);
							my @temp_event;
							$temp_event[0] = 'tandem_dup';
							$temp_event[1] = $cl_id1;
							$temp_event[2] = $mpd1;
							$temp_event[3] = $srd1;
							$temp_event[4] = $current_event[4];
							$temp_event[5] = $current_event[9];
							$temp_event[6] = $current_event[5];
							$temp_event[7] = $temp_event[6] - $temp_event[5];
							my $modified_event = join("\t", @temp_event);
							$data{$_} = $modified_event;
							#print "$modified_event\n";
						}
						if ($current_event[0] =~ /insod/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my ($srd1, $srd2) = split (/\//, $current_event[3]);
							my @temp_event;
							$temp_event[0] = 'invers_f';
							$temp_event[1] = $cl_id1;
							$temp_event[2] = $mpd1;
							$temp_event[3] = $srd1;
							$temp_event[4] = $current_event[4];
							$temp_event[5] = $current_event[5];
							$temp_event[6] = $current_event[10];
							$temp_event[7] = $temp_event[6] - $temp_event[5];
							my $modified_event = join("\t", @temp_event);
							$data{$_} = $modified_event;
							#print "$modified_event\n";
						}
						if ($current_event[0] =~ /insou/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my ($srd1, $srd2) = split (/\//, $current_event[3]);
							my @temp_event;
							$temp_event[0] = 'invers_f';
							$temp_event[1] = $cl_id1;
							$temp_event[2] = $mpd1;
							$temp_event[3] = $srd1;
							$temp_event[4] = $current_event[4];
							$temp_event[5] = $current_event[10];
							$temp_event[6] = $current_event[5];
							$temp_event[7] = $temp_event[6] - $temp_event[5];
							my $modified_event = join("\t", @temp_event);
							$data{$_} = $modified_event;
							#print "$modified_event\n";
						}
						if ($current_event[0] eq 'del_invers')
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my ($srd1, $srd2) = split (/\//, $current_event[3]);
							my @temp_event;
							$temp_event[0] = 'invers_f';
							$temp_event[1] = $cl_id1;
							$temp_event[2] = $mpd1;
							$temp_event[3] = $srd1;
							$temp_event[4] = $current_event[4];
							$temp_event[5] = $current_event[5];
							$temp_event[6] = $current_event[10];
							$temp_event[7] = $temp_event[6] - $temp_event[5];
							my $modified_event = join("\t", @temp_event);
							$data{$_} = $modified_event;
							#print "$modified_event\n";
						}
					}
					else
					{
						delete $data{$_};
					}
				}
			}
		}
	}
	for (my $i=0;$i<=$ftai;$i++)
	{
		my @dataa = split (/\t/, $data{$i});
		delete $dataa[1];
		delete $dataa[2];
		my $stringa = join("__", @dataa);
		for (my $j=0;$j<=$ftai;$j++)
		{
			next if ($i == $j);
			my @datab = split (/\t/, $data{$j});
			delete $datab[1];
			delete $datab[2];
			my $stringb = join("__", @datab);
			if ($stringa eq $stringb)
			{
				delete $data{$j};
			}
		}
	}
	
	open SRINTRAFILTER, ">$srintra_filter";
	for (my $i=0;$i<=$ftai;$i++)
	{
		print SRINTRAFILTER "$data{$i}\n" if ($data{$i});
	}
	close SRINTRAFILTER;
	
	my (%data, %cluster_event_map, %eventsize);
	# $cluster_event_map{$cluster_id} = array of event ids
	# $eventsize{event id} = event size
	my $ftei=0;
	open SRINTEROUTFT, "<$srinter_out";
	while ($newline = <SRINTEROUTFT>)
	{
		chomp $newline;
		my @data = split (/\t/, $newline);
		if ($data[1] =~ /\//)
		{
			my ($a, $b) = ($`, $');#'
			if ($a =~ /_/)
			{
				$a = $`;
			}
			if ($b =~ /_/)
			{
				$b = $`;
			}
			push @{$cluster_event_map{$a}}, $ftei;
			push @{$cluster_event_map{$b}}, $ftei;
			$eventsize{$ftei} += abs($data[7]);
			$eventsize{$ftei} += abs($data[11]);
		}
		else
		{
			my $a = $data[1];
			if ($a =~ /_/)
			{
				$a = $`;
			}
			push @{$cluster_event_map{$a}}, $ftei;
			$eventsize{$ftei} = 1000000000;
		}
		$data{$ftei} = $newline;
		$ftei++;
	}
	close SRINTEROUTFT;
	
	foreach my $cluster_id (%cluster_event_map)
	{
		if ($cluster_event_map{$cluster_id})
		{
			if (@{$cluster_event_map{$cluster_id}} > 1)
			{
				#print "$cluster_id\t@{$cluster_event_map{$cluster_id}}\n";
				my ($select_id, @delete_id);
				$select_id = $cluster_event_map{$cluster_id}[0];
				foreach (@{$cluster_event_map{$cluster_id}})
				{
					if ($eventsize{$_}<$eventsize{$select_id})
					{
						$select_id = $_;
					}
				}
				foreach (@{$cluster_event_map{$cluster_id}})
				{
					 if ($_ ne $select_id)
					 {
					 	push @delete_id, $_;
					 }
				}
				foreach (@delete_id)
				{
					my @current_event = split (/\t/, $data{$_});
					if ($current_event[1] =~ /\//)
					{
						if ($current_event[0] =~ /inss/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my ($srd1, $srd2) = split (/\//, $current_event[3]);
							my @temp_event;
							$temp_event[0] = 'transl_inter';
							$temp_event[1] = $cl_id1;
							$temp_event[2] = $mpd1;
							$temp_event[3] = $srd1;
							if ($current_event[4] lt $current_event[8])
							{
								$temp_event[4] = $current_event[4];
								$temp_event[5] = $current_event[5];
								$temp_event[6] = 1;
								$temp_event[7] = $current_event[8];
								$temp_event[8] = $current_event[9];
								$temp_event[9] = -1;
							}
							else
							{
								$temp_event[4] = $current_event[8];
								$temp_event[5] = $current_event[10];
								$temp_event[6] = 1;
								$temp_event[7] = $current_event[4];
								$temp_event[8] = $current_event[6];
								$temp_event[9] = -1;
							}
							my $modified_event = join("\t", @temp_event);
							$data{$_} = $modified_event;
							#print "inss $modified_event\n";
						}
						if ($current_event[0] =~ /inso/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my ($srd1, $srd2) = split (/\//, $current_event[3]);
							my @temp_event;
							$temp_event[0] = 'transl_inter';
							$temp_event[1] = $cl_id1;
							$temp_event[2] = $mpd1;
							$temp_event[3] = $srd1;
							if ($current_event[4] lt $current_event[8])
							{
								$temp_event[4] = $current_event[4];
								$temp_event[5] = $current_event[5];
								$temp_event[6] = 1;
								$temp_event[7] = $current_event[8];
								$temp_event[8] = $current_event[10];
								$temp_event[9] = 1;
							}
							else
							{
								$temp_event[4] = $current_event[8];
								$temp_event[5] = $current_event[10];
								$temp_event[6] = 1;
								$temp_event[7] = $current_event[4];
								$temp_event[8] = $current_event[5];
								$temp_event[9] = 1;
							}
							my $modified_event = join("\t", @temp_event);
							$data{$_} = $modified_event;
							#print "inso $modified_event\n";
						}
					}
					else
					{
						delete $data{$_};
					}
				}
			}
		}
	}
	for (my $i=0;$i<=$ftei;$i++)
	{
		my @dataa = split (/\t/, $data{$i});
		delete $dataa[1];
		delete $dataa[2];
		my $stringa = join("__", @dataa);
		for (my $j=0;$j<=$ftei;$j++)
		{
			next if ($i == $j);
			my @datab = split (/\t/, $data{$j});
			delete $datab[1];
			delete $datab[2];
			my $stringb = join("__", @datab);
			if ($stringa eq $stringb)
			{
				delete $data{$j};
			}
		}
	}
	
	open SRINTERFILTER, ">$srinter_filter";
	for (my $i=0;$i<=$ftei;$i++)
	{
		print SRINTERFILTER "$data{$i}\n" if ($data{$i});
	}
	close SRINTERFILTER;
	
	my $srintra_coord = $prefix.'.sr.intra.filtered.crd.sorted';
	my $srinter_coord = $prefix.'.sr.inter.filtered.crd.sorted';
	#system "sort $srintra_filter -k 5,5 -k 6,6n > $srintra_coord";
	#system "sort $srinter_filter -k 5,5 -k 6,6n > $srinter_coord";
}

#	cluster discordant split reads in one region (for inversion only)
sub sr_cluster_1
{
	my $ref_discord_sr = shift;
	my $ref_is = shift;
	my $cut_sr = shift;
	my $cl2 = shift;
	my $sr_dis_cutoff = $$ref_is{'rlu'} - $cut_sr;
	my %cluster;
	my $count = 0;
	
	for my $a (@{$ref_discord_sr})
	{
		my @data1 = split (/\t/, $a);
		my $readname = $data1[0];
		my $seqid = $data1[2];
		my $start = $data1[3];
		my $strand = 1;
		$strand = -1 if ($data1[1] =~ /r/);
		my $mseqid = $data1[6];
		my $mstart = $data1[7];
		my $mstrand = 1;
		$mstrand = -1 if ($data1[1] =~ /R/);
		my $isize = $data1[8];
		next unless ($mseqid eq '=');
		#print "$start\t$mstart\t$strand\t$mstrand\t$isize\n";
		
		$readname =~ /_(\d{2,3})$/;
		my $readlength = $1;
		my ($incluster, $topush);
SR1:		for (my $i=$count-5;$i<=$count;$i++)
		{
			if ($cluster{$i}[0])
			{
				my @data2 = split (/\t/, $cluster{$i}[-1]);
				my $readnamet = $data2[0];
				my $seqidt = $data2[2];
				my $startt = $data2[3];
				my $strandt = 1;
				$strandt = -1 if ($data2[1] =~ /r/);
				my $mseqidt = $data2[6];
				my $mstartt = $data2[7];
				my $mstrandt = 1;
				$mstrandt = -1 if ($data2[1] =~ /R/);
				my $isizet = $data2[8];
				$readnamet =~ /_(\d{2,3})$/;
				my $readlengtht = $1;
				$readnamet =~ /_(\d{2,3})$/;
				my $readlengtht = $1;
				#my $temp = abs(abs($start-$startt+$mstart-$mstartt)-abs($readlength-$readlengtht)); print "$temp\tabs(abs($start-$startt+$mstart-$mstartt)-abs($readlength-$readlengtht))\n";
				if (abs(abs($start-$startt+$mstart-$mstartt)-abs($readlength-$readlengtht)) <= 2 and abs($startt-$start)<$$ref_is{'rlu'} and abs($mstartt-$mstart)<$$ref_is{'rlu'})
				{
					$topush = $i;
					$incluster = 1;
					last SR1;
				}
			}
		}
		unless ($incluster)
		{
			$count++;
			push @{$cluster{$count}}, $a;
		}
		if (defined($topush))
		{
			push @{$cluster{$topush}}, $a;
		}
	}
	
	my %clustersize;
	for (my $i=1;$i<=$count;$i++)
	{
		$clustersize{$i} = @{$cluster{$i}} if ($cluster{$i});
	}
	
	my (@boundary1, @boundary2, @support_sr, %bpread);
	my $k = 0;
	foreach my $i (sort { $clustersize {$b} <=> $clustersize {$a}} keys %clustersize )
	{
		my (@start, @mstart, @isize);
		foreach (@{$cluster{$i}})
		{
			my @data1 = split (/\t/, $_);
			my $readname = $data1[0];
			my $start = $data1[3];
			my $mstart = $data1[7];
			#print "$i $start $mstart\n";
			push @start, $start;
			push @mstart, $mstart;
			push @{$bpread{$k}}, $readname;
		}
		@start = sort { $a <=> $b } @start;
		@mstart = sort { $a <=> $b } @mstart;
		#print "$start[0]\t$start[-1]\t$mstart[0]\t$mstart[-1]\n";
		if (($start[-1]-$start[0])<$$ref_is{'rlu'} and ($mstart[-1]-$mstart[0])<$$ref_is{'rlu'})
		{
			my $support_sr = @{$cluster{$i}};
			if ($cl2)
			{
				push @boundary1, $start[0];
				push @boundary1, $start[-1];
				push @boundary2, $mstart[0];
				push @boundary2, $mstart[-1];
				push @support_sr, $support_sr;
			}
			else
			{
				return ($start[0], $start[-1], $mstart[0], $mstart[-1], $support_sr, \%bpread);
			}
		}
		else
		{
			$bpread{$k} = undef;
			next;
		}
		last if ($k);
		$k++;
	}
	if ($cl2)
	{
		return (\@boundary1, \@boundary2, \@support_sr, \%bpread);
	}
	return 0;
}

#	cluster discordant split reads
sub sr_cluster
{
	my $ref_discord_sr = shift;
	my $ref_is = shift;
	my $cut_sr = shift;
	my $cl2 = shift; # there are 2 clusters in this region
	my $internal_break = shift; # 0 no break, 1 break on left side, 2 break on right side
	my %cluster;
	my $count = 0;
	
	for my $a (@{$ref_discord_sr})
	{
		my @data1 = split (/\t/, $a);
		my $readname = $data1[0];
		my $seqid = $data1[2];
		my $start = $data1[3];
		my $strand = 1;
		$strand = -1 if ($data1[1] =~ /r/);
		my $mseqid = $data1[6];
		my $mstart = $data1[7];
		my $mstrand = 1;
		$mstrand = -1 if ($data1[1] =~ /R/);
		my $isize = $data1[8];
		next unless ($mseqid eq '=');
		
		$readname =~ /_(\d{2,3})$/;
		my $readlength = $1;
		my ($incluster, $topush);
SR:		for (my $i=$count-5;$i<=$count;$i++)
		{
			foreach (@{$cluster{$i}})
			{
				my @data2 = split (/\t/, $_);
				my $readnamet = $data2[0];
				my $seqidt = $data2[2];
				my $startt = $data2[3];
				my $strandt = 1;
				$strandt = -1 if ($data2[1] =~ /r/);
				my $mseqidt = $data2[6];
				my $mstartt = $data2[7];
				my $mstrandt = 1;
				$mstrandt = -1 if ($data2[1] =~ /R/);
				my $isizet = $data2[8];
				$readnamet =~ /_(\d{2,3})$/;
				my $readlengtht = $1;
				if ($readlength>=$readlengtht and abs(($isize-$isizet)-($readlength-$readlengtht)) <= 2 and abs($startt-$start)<$$ref_is{'rlu'} and abs($mstartt-$mstart)<$$ref_is{'rlu'})
				{
					$topush = $i;
					#print "($isize-$isizet-($readlength-$readlengtht))\n";
					$incluster = 1;
					last SR;
				}
				if ($readlength<$readlengtht and abs(($isizet-$isize)-($readlengtht-$readlength)) <= 2 and abs($startt-$start)<$$ref_is{'rlu'} and abs($mstartt-$mstart)<$$ref_is{'rlu'})
				{
					$topush = $i;
					#print "($isizet-$isize-($readlengtht-$readlength))\n";
					$incluster = 1;
					last SR;
				}
			}
		}
		unless ($incluster)
		{
			$count++;
			push @{$cluster{$count}}, $a;
		}
		if (defined($topush))
		{
			push @{$cluster{$topush}}, $a;
		}
	}
	
	my %clustersize;
	for (my $i=1;$i<=$count;$i++)
	{
		$clustersize{$i} = @{$cluster{$i}};
	}

	# split 1 cluster into 2, split at largest distance read
	if ($internal_break)
	{
		my $tosplit;
		$tosplit = 1 unless ($clustersize{2});
		if ($tosplit)
		{		
			my $clusteri;
			foreach my $i (sort { $clustersize {$b} <=> $clustersize {$a}} keys %clustersize )
			{
				$clusteri = $i;
				last;
			}
			my (%cluster_temp, $laststart, $lastmstart, $max_dis, $max_name, $max_pos);
			if ($internal_break == 1)
			{
				# find largest distance
				foreach (@{$cluster{$clusteri}})
				{
					my @data1 = split (/\t/, $_);
					my $readname = $data1[0];
					my $start = $data1[3];
					#print "$readname\t$start\n";
					if ($laststart)
					{
						my $distance = abs($start-$laststart);
						if ($distance > $max_dis)
						{
							$max_dis = $distance;
							$max_name = $readname;
							$max_pos = $start;
						}
					}
					$laststart = $start;
				}
				#print "$max_name\t$max_dis\n";
				# split one cluster into two
				my $cluster2;
				foreach (@{$cluster{$clusteri}})
				{
					my @data1 = split (/\t/, $_);
					my $readname = $data1[0];
					my $start = $data1[3];
					#print "$readname ";
					if ($readname eq $max_name and $start == $max_pos)
					{
						$cluster2 = 1;
					}
					if ($cluster2)
					{
						push @{$cluster_temp{2}}, $_;
						#print "2 $_\n";
					}
					else
					{
						push @{$cluster_temp{1}}, $_;
						#print "1 $_\n";
					}
				}
				%cluster = %cluster_temp;
			}
			if ($internal_break == 2)
			{
				# sort by mstart
				my (%sort_mstart, %unsort_mstart);
				foreach (@{$cluster{$clusteri}})
				{
					my @data1 = split (/\t/, $_);
					my $readname = $data1[0];
					my $mstart = $data1[7];
					$unsort_mstart{$readname}{'ob'} = $_;
					$unsort_mstart{$readname}{'T'} = $mstart;
				}
				my $t = 0;
				foreach my $readname (sort { $unsort_mstart{$a}{'T'} <=> $unsort_mstart{$b}{'T'}} keys %unsort_mstart )
				{
					$sort_mstart{$t} = $unsort_mstart{$readname}{'ob'};
					$t++;
				}
				# find largest distance
				for (my $t=0;$sort_mstart{$t};$t++)
				{
					my @data1 = split (/\t/, $_);
					my $readname = $data1[0];
					my $mstart = $data1[7];
					if ($lastmstart)
					{
						my $distance = abs($mstart-$lastmstart);
						if ($distance > $max_dis)
						{
							$max_dis = $distance;
							$max_name = $readname;
							$max_pos = $mstart;
						}
					}
					$lastmstart = $mstart;
				}
				# split one cluster into two
				my $cluster2;
				for (my $t=0;$sort_mstart{$t};$t++)
				{
					my @data1 = split (/\t/, $_);
					my $readname = $data1[0];
					my $mstart = $data1[7];
					if ($readname eq $max_name and $mstart == $max_pos)
					{
						$cluster2 = 1;
					}
					if ($cluster2)
					{
						push @{$cluster_temp{2}}, $sort_mstart{$t};
					}
					else
					{
						push @{$cluster_temp{1}}, $sort_mstart{$t};
					}
				}
				%cluster = %cluster_temp;
			}
			
			# update clustersize
			%clustersize = undef;
			for (my $i=1;$i<=2;$i++)
			{
				$clustersize{$i} = @{$cluster{$i}} if ($cluster{$i});
			}
		}
	}
	
	my (@boundary1, @boundary2, @support_sr, %bpread);
	my $k = 0;
	foreach my $i (sort { $clustersize {$b} <=> $clustersize {$a}} keys %clustersize )
	{
		my (@start, @mstart);
		foreach (@{$cluster{$i}})
		{
			my @data1 = split (/\t/, $_);
			my $readname = $data1[0];
			my $start = $data1[3];
			my $mstart = $data1[7];
			#print "$i $start $mstart\n";
			push @start, $start;
			push @mstart, $mstart;
			push @{$bpread{$k}}, $readname;
		}
		@start = sort { $a <=> $b } @start;
		@mstart = sort { $a <=> $b } @mstart;
		#print "$start[0]\t$start[-1]\t$mstart[0]\t$mstart[-1]\n";
		if (($start[-1]-$start[0])<$$ref_is{'rlu'} and ($mstart[-1]-$mstart[0])<$$ref_is{'rlu'})
		{
			my $support_sr = @{$cluster{$i}};
			if ($cl2)
			{
				push @boundary1, $start[0];
				push @boundary1, $start[-1];
				push @boundary2, $mstart[0];
				push @boundary2, $mstart[-1];
				push @support_sr, $support_sr;
			}
			else
			{
				return ($start[0], $start[-1], $mstart[0], $mstart[-1], $support_sr, \%bpread);
			}
		}
		else
		{
			$bpread{$k} = undef;
			next;
		}
		last if ($k);
		$k++;
	}
	if ($cl2)
	{
		return (\@boundary1, \@boundary2, \@support_sr, \%bpread);
	}
	return 0;
}

#	call all SVs based on split reads
sub srd
{
#	file format of prefix.sr.intra.out
#	del, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size
#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion and insertion
#	del_invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of inversion (2 col), inversion size, distance of inverstion and deletion (2 col)
#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion and insertion
#	invers, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size
#	invers_*, cluster id, number of supporting read pairs, number of supporting split reads, chr, inversion left boundary, inversion right boundary, inversion size
#	tandem_dup, cluster id, number of supporting read pairs, number of supporting split reads, chr, tandem duplication boundary 1, tandem duplication boundary 2, tandem duplication size

#	file format of prefix.sr.inter.out
#	del_ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of deletion (2 col), deletion size, chr of insertion donor, range of insertion (2 col), insert size
#	ins*, cluster id, number of supporting read pairs, number of supporting split reads, chr, range of insert site, insert site size, chr of insertion donor, range of insertion (2 col), insert size
#	transl_inter, cluster id, number of supporting read pairs, number of supporting split reads, chr of 1st cluster, boundary of 1st cluster, orientation of 1st cluster, chr of 2nd cluster, boundary of 2nd cluster, orientation of 2nd cluster

#	$left_bound_sr1: left bound in primary cluster
#	$right_bound_sr1: right bound in primary cluster
#	$left_bound_sr2: left bound in secondary cluster
#	$right_bound_sr2: right bound in secondary cluster
#	primary cluster: cluster to initiate call in mpd
#	secondary cluster: cluster to accompany primary cluster in mpd

	my $line = shift;
	my $prefix = shift;
	my $ref_is = shift;
	my $cut_sr = shift;
	my $sv_size_cutoff = shift;
	my $support_reads = shift;
	my $bpinfofile = $prefix.'.bp.info';
	my $mpintra_outfile = $prefix.'.mp.intra.out';
	my $mpinter_outfile = $prefix.'.mp.inter.out';
	my $srintra_outfile = $prefix.'.sr.intra.out';
	my $srinter_outfile = $prefix.'.sr.inter.out';
	my $bp_readsfile = $prefix.'.bp_reads';
	my $sr_sortbam = $prefix.'.sr.sorted.bam';
	system "rm $bp_readsfile" if (-e $bp_readsfile);
	
	my ($newline, $newline1, %cluster_region);
	# %cluster_region: index of correponding region for each cluster
	# $cluster_region{$cluster_id} = region name
	open BPDTSRD, "<$bpinfofile";
	while ($newline = <BPDTSRD>)
	{
		chomp $newline;
		my @data = split (/\t/, $newline);
		my @cl = split (/\//, $data[0]);
		foreach (@cl)
		{
			$cluster_region{$_} = $data[1];
			#print "$_\t$cluster_region{$_}\n";
		}
	}
	close BPDTSRD;
	
	my ($i, @result_sr, @unsupport_del);
	my $sri = 0; # index in @result_sr
	my $usdi = 0; # index in @unsupport_del
	open MPINTRASRD, "<$mpintra_outfile";
	open BPREAD, ">$bp_readsfile";
	while ($newline = <MPINTRASRD>)
	{
		chomp $newline;
		$i++;
		next if (defined($line) and $i < $line);
		last if (defined($line) and $i > $line);
		my @data = split (/\t/, $newline);
		my @cl = split (/\//, $data[1]);
		my (@discord_sr1, @discord_sr2);
		my ($ref_discord_sr1, $ref_discord_sr2) = (\@discord_sr1, \@discord_sr2);
		my %bp_window;
		# the window in candidate region where a break point is located
		# $bp_window{$data[4]} = 1/2, 1 left window, 2 right window
		if ($data[0] eq 'del')
		{
			next unless ($cluster_region{$cl[0]});
			my ($trash, $position1, $position2, $trash, $position3, $position4) = split (/__/, $cluster_region{$cl[0]});
			if (&covered($data[4], $data[4], $position1, $position2))
			{
				$bp_window{$data[4]} = 1;
			}
			if (&covered($data[4], $data[4], $position3, $position4))
			{
				$bp_window{$data[4]} = 2;
			}
			if (&covered($data[5], $data[5], $position1, $position2))
			{
				$bp_window{$data[5]} = 1;
			}
			if (&covered($data[5], $data[5], $position3, $position4))
			{
				$bp_window{$data[5]} = 2;
			}
			#print "$cluster_region{$cl[0]}\n$data[4]\t$bp_window{$data[4]}\n$data[5]\t$bp_window{$data[5]}\n";
			my @alignments;
			open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
			while ($newline1 = <SAM>)
			{
				chomp $newline1;
				push @alignments, $newline1;
			}
			close SAM;
			if ($bp_window{$data[4]} ne $bp_window{$data[5]})
			{
				@discord_sr1 = undef;
				for my $a (@alignments)
				{
					my @data1 = split (/\t/, $a);
					my $start = $data1[3];
					my $strand = 1;
					$strand = -1 if ($data1[1] =~ /r/);
					my $mseqid = $data1[6];
					my $mstart = $data1[7];
					my $mstrand = 1;
					$mstrand = -1 if ($data1[1] =~ /R/);
					my $isize = $data1[8];
					if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand)
					{
						#print "$start $mstart\t$isize\t$strand $mstrand\n";
						push @discord_sr1, $a;
					}
				}
				if (@discord_sr1 >= $support_reads)
				{
					my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 0);
					#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
					# large deletion with small close insertion
					if ($$ref_support_sr[1] >= $support_reads)
					{
						# inssu
						if (&covered($$ref_boundary1[0], $$ref_boundary1[1], $$ref_boundary1[2], $$ref_boundary1[3]) and (($$ref_boundary2[0]>=($position2-$position1+101) and $$ref_boundary2[2]<=($position2-$position1+101) and $$ref_boundary1[0]<=($position2-$position1+101) and $$ref_boundary1[2]<=($position2-$position1+101)) or ($$ref_boundary2[2]>=($position2-$position1+101) and $$ref_boundary2[0]<=($position2-$position1+101) and $$ref_boundary1[0]<=($position2-$position1+101) and $$ref_boundary1[2]<=($position2-$position1+101))))
						{
							if ($$ref_boundary2[0]<$$ref_boundary2[2])
							{
								my $left_bound_sr1 = $position1 + $$ref_boundary1[0];
								my $right_bound_sr1 = $position1 + $$ref_boundary2[1] + $cut_sr;
								my $left_bound_sr2 = $position1 + $$ref_boundary1[3] + $cut_sr;
								my $right_bound_sr2 = $$ref_boundary2[2] - ($position2-$position1+101) + $position3;
								my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
								my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
									if ($del_size > 10)
									{
										$result_sr[$sri][0] = 'del_inssu';
									}
									else
									{
										$result_sr[$sri][0] = 'inssu';
									}
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $right_bound_sr1;
									$result_sr[$sri][6] = $right_bound_sr2;
									$result_sr[$sri][7] = $del_size;
									$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
									$sri++;
								}
							}
							else
							{
								my $left_bound_sr1 = $position1 + $$ref_boundary1[2];
								my $right_bound_sr1 = $position1 + $$ref_boundary2[3] + $cut_sr;
								my $left_bound_sr2 = $position1 + $$ref_boundary1[1] + $cut_sr;
								my $right_bound_sr2 = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
								my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
								my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{1}});
									my $reads1 = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
									if ($del_size > 10)
									{
										$result_sr[$sri][0] = 'del_inssu';
									}
									else
									{
										$result_sr[$sri][0] = 'inssu';
									}
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $right_bound_sr1;
									$result_sr[$sri][6] = $right_bound_sr2;
									$result_sr[$sri][7] = $del_size;
									$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
									$sri++;
								}
							}
						}
						# inssd
						elsif (&covered($$ref_boundary2[0], $$ref_boundary2[1], $$ref_boundary2[2], $$ref_boundary2[3]) and (($$ref_boundary1[0]<=($position2-$position1+101) and $$ref_boundary1[2]>=($position2-$position1+101) and $$ref_boundary2[0]>=($position2-$position1+101) and $$ref_boundary2[2]>=($position2-$position1+101)) or ($$ref_boundary1[2]<=($position2-$position1+101) and $$ref_boundary1[0]>=($position2-$position1+101) and $$ref_boundary2[0]>=($position2-$position1+101) and $$ref_boundary2[2]>=($position2-$position1+101))))
						{
							if ($$ref_boundary1[0]<$$ref_boundary1[2])
							{
								my $left_bound_sr1 = $$ref_boundary1[2] - ($position2-$position1+101) + $position3;
								my $right_bound_sr1 = $$ref_boundary2[3] - ($position2-$position1+101) + $position3 + $cut_sr;
								my $left_bound_sr2 = $position1 + $$ref_boundary1[1] + $cut_sr;
								my $right_bound_sr2 = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
								my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
								my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
									if ($del_size > 10)
									{
										$result_sr[$sri][0] = 'del_inssd';
									}
									else
									{
										$result_sr[$sri][0] = 'inssd';
									}
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr2;
									$result_sr[$sri][6] = $left_bound_sr1;
									$result_sr[$sri][7] = $del_size;
									$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
									$sri++;
								}
							}
							else
							{
								my $left_bound_sr1 = $$ref_boundary1[0] - ($position2-$position1+101) + $position3;
								my $right_bound_sr1 = $$ref_boundary2[1] - ($position2-$position1+101) + $position3 + $cut_sr;
								my $left_bound_sr2 = $position1 + $$ref_boundary1[3] + $cut_sr;
								my $right_bound_sr2 = $$ref_boundary2[2] - ($position2-$position1+101) + $position3;
								my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
								my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{1}});
									my $reads1 = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
									if ($del_size > 10)
									{
										$result_sr[$sri][0] = 'del_inssd';
									}
									else
									{
										$result_sr[$sri][0] = 'inssd';
									}
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr2;
									$result_sr[$sri][6] = $left_bound_sr1;
									$result_sr[$sri][7] = $del_size;
									$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
									$sri++;
								}
							}
						}
						# del
						else
						{
							# re-process deletion to require direction
							@discord_sr1 = undef;
							for my $a (@alignments)
							{
								my @data1 = split (/\t/, $a);
								my $start = $data1[3];
								my $strand = 1;
								$strand = -1 if ($data1[1] =~ /r/);
								my $mseqid = $data1[6];
								my $mstart = $data1[7];
								my $mstrand = 1;
								$mstrand = -1 if ($data1[1] =~ /R/);
								my $isize = $data1[8];
								my $flags;
								$flags = 'FIRST_MATE' if ($data1[1] =~ /1/);
								$flags = 'SECOND_MATE' if ($data1[1] =~ /2/);
								if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start < ($position2-$position1) and $mstart >= ($position2-$position1+100))
								{
									if (($flags =~ /FIRST_MATE/ and $strand == 1) or ($flags =~ /SECOND_MATE/ and $strand == -1))
									{
										#print "$start $mstart\t$isize\t$strand $mstrand\n";
										push @discord_sr1, $a;
									}
								}
							}
							my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 0);
							#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
							
							my ($left_bound_sr, $right_bound_sr);
							if ($$ref_boundary1[1] < ($position2-$position1+101) and $$ref_boundary2[0] >= ($position2-$position1+101))
							{
								$left_bound_sr = $position1 + $$ref_boundary1[1] + $cut_sr;
								$right_bound_sr = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
							}
							if ($$ref_boundary1[3] < ($position2-$position1+101) and $$ref_boundary2[2] >= ($position2-$position1+101))
							{
								$left_bound_sr = $position1 + $$ref_boundary1[3] + $cut_sr;
								$right_bound_sr = $$ref_boundary2[2] - ($position2-$position1+101) + $position3;
							}
							my $del_size = $right_bound_sr - $left_bound_sr - 1;
							if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr and $right_bound_sr)
							{
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = $$ref_support_sr[0];
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr;
								$result_sr[$sri][6] = $right_bound_sr;
								$result_sr[$sri][7] = "$del_size\n";
								$sri++;
							}
							else
							{
								$unsupport_del[$usdi][0] = $data[0];
								$unsupport_del[$usdi][1] = $data[1];
								$unsupport_del[$usdi][2] = $data[2];
								$unsupport_del[$usdi][3] = $data[3];
								$unsupport_del[$usdi][4] = $data[4];
								$unsupport_del[$usdi][5] = $data[5];
								$unsupport_del[$usdi][6] = $data[6];
								$unsupport_del[$usdi][7] = "$data[7]\n";
								$usdi++;
							}
						}
					}
					# del
					else
					{
						# re-process deletion to require direction
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							my $flags;
							$flags = 'FIRST_MATE' if ($data1[1] =~ /1/);
							$flags = 'SECOND_MATE' if ($data1[1] =~ /2/);
							if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start < ($position2-$position1) and $mstart >= ($position2-$position1+100))
							{
								if (($flags =~ /FIRST_MATE/ and $strand == 1) or ($flags =~ /SECOND_MATE/ and $strand == -1))
								{
									#print "$start $mstart\t$isize\t$strand $mstrand\n";
									push @discord_sr1, $a;
								}
							}
						}
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 0);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
							
						my $left_bound_sr = $position1 + $$ref_boundary1[1] + $cut_sr if ($$ref_boundary1[1] < ($position2-$position1+101));
						my $right_bound_sr = $$ref_boundary2[0] - ($position2-$position1+101) + $position3 if ($$ref_boundary2[0] >= ($position2-$position1+101));
						$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
						$right_bound_sr = $data[5] unless ($$ref_boundary2[0]);
						my $del_size = $right_bound_sr - $left_bound_sr - 1;
						if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr and $right_bound_sr)
						{
							my $reads = join("\t", @{$$ref_bpread{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = $data[0];
							$result_sr[$sri][1] = $data[1];
							$result_sr[$sri][2] = $data[2];
							$result_sr[$sri][3] = $$ref_support_sr[0];
							$result_sr[$sri][4] = $data[3];
							$result_sr[$sri][5] = $left_bound_sr;
							$result_sr[$sri][6] = $right_bound_sr;
							$result_sr[$sri][7] = "$del_size\n";
							$sri++;
						}
						else
						{
							$unsupport_del[$usdi][0] = $data[0];
							$unsupport_del[$usdi][1] = $data[1];
							$unsupport_del[$usdi][2] = $data[2];
							$unsupport_del[$usdi][3] = $data[3];
							$unsupport_del[$usdi][4] = $data[4];
							$unsupport_del[$usdi][5] = $data[5];
							$unsupport_del[$usdi][6] = $data[6];
							$unsupport_del[$usdi][7] = "$data[7]\n";
							$usdi++;
						}
					}
				}
				else
				{
					$unsupport_del[$usdi][0] = $data[0];
					$unsupport_del[$usdi][1] = $data[1];
					$unsupport_del[$usdi][2] = $data[2];
					$unsupport_del[$usdi][3] = $data[3];
					$unsupport_del[$usdi][4] = $data[4];
					$unsupport_del[$usdi][5] = $data[5];
					$unsupport_del[$usdi][6] = $data[6];
					$unsupport_del[$usdi][7] = "$data[7]\n";
					$usdi++;
				}
			}
			else
			{
				if ($bp_window{$data[4]} == 1 and $bp_window{$data[5]} == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						my $flags;
						$flags = 'FIRST_MATE' if ($data1[1] =~ /1/);
						$flags = 'SECOND_MATE' if ($data1[1] =~ /2/);
						if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start < ($position2-$position1) and $mstart < ($position2-$position1))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand $flags\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 0);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						if ($$ref_support_sr[1] >= $support_reads)
						{
							# inssu
							if (&covered($$ref_boundary1[0], $$ref_boundary1[1], $$ref_boundary1[2], $$ref_boundary1[3]))
							{
								if ($$ref_boundary2[0]<$$ref_boundary2[2])
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[0];
									my $right_bound_sr1 = $position1 + $$ref_boundary2[1] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr2 = $position1 + $$ref_boundary2[2];
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssu';
										}
										else
										{
											$result_sr[$sri][0] = 'inssu';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[2];
									my $right_bound_sr1 = $position1 + $$ref_boundary2[3] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr2 = $position1 + $$ref_boundary2[0];
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssu';
										}
										else
										{
											$result_sr[$sri][0] = 'inssu';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
										$sri++;
									}
								}
							}
							# inssd
							elsif (&covered($$ref_boundary2[0], $$ref_boundary2[1], $$ref_boundary2[2], $$ref_boundary2[3]))
							{
								if ($$ref_boundary1[0]<$$ref_boundary1[2])
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[2];
									my $right_bound_sr1 = $position1 + $$ref_boundary2[3] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr2 = $position1 + $$ref_boundary2[0];
									my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssd';
										}
										else
										{
											$result_sr[$sri][0] = 'inssd';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr2;
										$result_sr[$sri][6] = $left_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[0];
									my $right_bound_sr1 = $position1 + $$ref_boundary2[1] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr2 = $position1 + $$ref_boundary2[2];
									my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssd';
										}
										else
										{
											$result_sr[$sri][0] = 'inssd';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr2;
										$result_sr[$sri][6] = $left_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
								}
							}
							else
							{
								my $left_bound_sr = $position1 + $$ref_boundary1[1] + $cut_sr;
								my $right_bound_sr = $position1 + $$ref_boundary2[0];
								$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
								$right_bound_sr = $data[5] unless ($$ref_boundary2[0]);
								my $del_size = $right_bound_sr - $left_bound_sr - 1;
								if ($$ref_support_sr[0] >= $support_reads)
								{
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = $data[0];
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = $$ref_support_sr[0];
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$del_size\n";
									$sri++;
								}
								else
								{
									$unsupport_del[$usdi][0] = $data[0];
									$unsupport_del[$usdi][1] = $data[1];
									$unsupport_del[$usdi][2] = $data[2];
									$unsupport_del[$usdi][3] = $data[3];
									$unsupport_del[$usdi][4] = $data[4];
									$unsupport_del[$usdi][5] = $data[5];
									$unsupport_del[$usdi][6] = $data[6];
									$unsupport_del[$usdi][7] = "$data[7]\n";
									$usdi++;
								}
							}
						}
						else
						{
							my $left_bound_sr = $position1 + $$ref_boundary1[1] + $cut_sr;
							my $right_bound_sr = $position1 + $$ref_boundary2[0];
							$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
							$right_bound_sr = $data[5] unless ($$ref_boundary2[0]);
							my $del_size = $right_bound_sr - $left_bound_sr - 1;
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = $$ref_support_sr[0];
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr;
								$result_sr[$sri][6] = $right_bound_sr;
								$result_sr[$sri][7] = "$del_size\n";
								$sri++;
							}
							else
							{
								$unsupport_del[$usdi][0] = $data[0];
								$unsupport_del[$usdi][1] = $data[1];
								$unsupport_del[$usdi][2] = $data[2];
								$unsupport_del[$usdi][3] = $data[3];
								$unsupport_del[$usdi][4] = $data[4];
								$unsupport_del[$usdi][5] = $data[5];
								$unsupport_del[$usdi][6] = $data[6];
								$unsupport_del[$usdi][7] = "$data[7]\n";
								$usdi++;
							}
						}
					}
					else
					{
						$unsupport_del[$usdi][0] = $data[0];
						$unsupport_del[$usdi][1] = $data[1];
						$unsupport_del[$usdi][2] = $data[2];
						$unsupport_del[$usdi][3] = $data[3];
						$unsupport_del[$usdi][4] = $data[4];
						$unsupport_del[$usdi][5] = $data[5];
						$unsupport_del[$usdi][6] = $data[6];
						$unsupport_del[$usdi][7] = "$data[7]\n";
						$usdi++;
					}
				}
				# both break points in window 2
				else
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start >= ($position2-$position1+100) and $mstart >= ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 0);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						if ($$ref_support_sr[1] >= $support_reads)
						{
							# inssu
							if (&covered($$ref_boundary1[0], $$ref_boundary1[1], $$ref_boundary1[2], $$ref_boundary1[3]))
							{
								if ($$ref_boundary2[0]<$$ref_boundary2[2])
								{
									my $left_bound_sr1 = $$ref_boundary1[0] - ($position2-$position1+101) + $position3;
									my $right_bound_sr1 = $$ref_boundary2[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr2 = $$ref_boundary2[2] - ($position2-$position1+101) + $position3;
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssu';
										}
										else
										{
											$result_sr[$sri][0] = 'inssu';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $$ref_boundary1[2] - ($position2-$position1+101) + $position3;
									my $right_bound_sr1 = $$ref_boundary2[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr2 = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssu';
										}
										else
										{
											$result_sr[$sri][0] = 'inssu';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
										$sri++;
									}
								}
							}
							# inssd
							elsif (&covered($$ref_boundary2[0], $$ref_boundary2[1], $$ref_boundary2[2], $$ref_boundary2[3]))
							{
								if ($$ref_boundary1[0]<$$ref_boundary1[2])
								{
									my $left_bound_sr1 = $$ref_boundary1[2] - ($position2-$position1+101) + $position3;
									my $right_bound_sr1 = $$ref_boundary2[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr2 = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
									my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssd';
										}
										else
										{
											$result_sr[$sri][0] = 'inssd';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr2;
										$result_sr[$sri][6] = $left_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $$ref_boundary1[0] - ($position2-$position1+101) + $position3;
									my $right_bound_sr1 = $$ref_boundary2[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr2 = $$ref_boundary2[2] - ($position2-$position1+101) + $position3;
									my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssd';
										}
										else
										{
											$result_sr[$sri][0] = 'inssd';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr2;
										$result_sr[$sri][6] = $left_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
								}
							}
							else
							{
								my $left_bound_sr = $$ref_boundary1[1] - ($position2-$position1+101) + $position3 + $cut_sr;
								my $right_bound_sr = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
								$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
								$right_bound_sr = $data[5] unless ($$ref_boundary2[0]);
								my $del_size = $right_bound_sr - $left_bound_sr - 1;
								if ($$ref_support_sr[0] >= $support_reads)
								{
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = $data[0];
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = $$ref_support_sr[0];
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$del_size\n";
									$sri++;
								}
								else
								{
									$unsupport_del[$usdi][0] = $data[0];
									$unsupport_del[$usdi][1] = $data[1];
									$unsupport_del[$usdi][2] = $data[2];
									$unsupport_del[$usdi][3] = $data[3];
									$unsupport_del[$usdi][4] = $data[4];
									$unsupport_del[$usdi][5] = $data[5];
									$unsupport_del[$usdi][6] = $data[6];
									$unsupport_del[$usdi][7] = "$data[7]\n";
									$usdi++;
								}
							}
						}
						else
						{
							my $left_bound_sr = $$ref_boundary1[1] - ($position2-$position1+101) + $position3 + $cut_sr;
							my $right_bound_sr = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
							$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
							$right_bound_sr = $data[5] unless ($$ref_boundary2[0]);
							my $del_size = $right_bound_sr - $left_bound_sr - 1;
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = $$ref_support_sr[0];
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr;
								$result_sr[$sri][6] = $right_bound_sr;
								$result_sr[$sri][7] = "$del_size\n";
								$sri++;
							}
							else
							{
								$unsupport_del[$usdi][0] = $data[0];
								$unsupport_del[$usdi][1] = $data[1];
								$unsupport_del[$usdi][2] = $data[2];
								$unsupport_del[$usdi][3] = $data[3];
								$unsupport_del[$usdi][4] = $data[4];
								$unsupport_del[$usdi][5] = $data[5];
								$unsupport_del[$usdi][6] = $data[6];
								$unsupport_del[$usdi][7] = "$data[7]\n";
								$usdi++;
							}
						}
					}
					else
					{
						$unsupport_del[$usdi][0] = $data[0];
						$unsupport_del[$usdi][1] = $data[1];
						$unsupport_del[$usdi][2] = $data[2];
						$unsupport_del[$usdi][3] = $data[3];
						$unsupport_del[$usdi][4] = $data[4];
						$unsupport_del[$usdi][5] = $data[5];
						$unsupport_del[$usdi][6] = $data[6];
						$unsupport_del[$usdi][7] = "$data[7]\n";
						$usdi++;
					}
				}
			}
		}
		if ($data[0] =~ /inssd/)
		{
			next unless ($cluster_region{$cl[0]} and $cluster_region{$cl[1]});
			my ($trash, $positiona1, $positiona2, $trash, $positiona3, $positiona4, $orientation) = split (/__/, $cluster_region{$cl[0]});
			my ($trash, $positionb1, $positionb2, $trash, $positionb3, $positionb4, $orientation) = split (/__/, $cluster_region{$cl[1]});
			if (&covered($data[5], $data[5], $positiona1, $positiona2))
			{
				$bp_window{$data[5]} = 1;
			}
			if (&covered($data[5], $data[5], $positiona3, $positiona4))
			{
				$bp_window{$data[5]} = 2;
			}
			if (&covered($data[9], $data[9], $positiona1, $positiona2))
			{
				$bp_window{$data[9]} = 1;
			}
			if (&covered($data[9], $data[9], $positiona3, $positiona4))
			{
				$bp_window{$data[9]} = 2;
			}
			if (&covered($data[4], $data[4], $positionb1, $positionb2))
			{
				$bp_window{$data[4]} = 1;
			}
			if (&covered($data[4], $data[4], $positionb3, $positionb4))
			{
				$bp_window{$data[4]} = 2;
			}
			if (&covered($data[8], $data[8], $positionb1, $positionb2))
			{
				$bp_window{$data[8]} = 1;
			}
			if (&covered($data[8], $data[8], $positionb3, $positionb4))
			{
				$bp_window{$data[8]} = 2;
			}
			#print "$cluster_region{$cl[0]}\n$cluster_region{$cl[1]}\n$data[4]\t$bp_window{$data[4]}\n$data[5]\t$bp_window{$data[5]}\n$data[8]\t$bp_window{$data[8]}\n$data[9]\t$bp_window{$data[9]}\n";
			if ($cluster_region{$cl[0]} eq $cluster_region{$cl[1]})
			{
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				@discord_sr1 = undef;
				for my $a (@alignments)
				{
					my @data1 = split (/\t/, $a);
					my $start = $data1[3];
					my $strand = 1;
					$strand = -1 if ($data1[1] =~ /r/);
					my $mseqid = $data1[6];
					my $mstart = $data1[7];
					my $mstrand = 1;
					$mstrand = -1 if ($data1[1] =~ /R/);
					my $isize = $data1[8];
					if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand)
					{
						#print "$start $mstart\t$isize\t$strand $mstrand\n";
						push @discord_sr1, $a;
					}
				}
				if (@discord_sr1 >= $support_reads)
				{
					my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 1);
					#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
					if ($$ref_boundary1[2] > $$ref_boundary1[0])
					{
						my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
						$left_bound_sr1 = $positiona1 + $$ref_boundary1[2] if ($bp_window{$data[5]} == 1 and $$ref_boundary1[2] < ($positiona2-$positiona1+101));
						$left_bound_sr1 = $$ref_boundary1[2] - ($positiona2-$positiona1+101) + $positiona3 if ($bp_window{$data[5]} == 2 and $$ref_boundary1[2] >= ($positiona2-$positiona1+101));
						$right_bound_sr1 = $positiona1 + $$ref_boundary2[3] + $cut_sr if ($bp_window{$data[9]} == 1 and $$ref_boundary2[3] < ($positiona2-$positiona1+101));
						$right_bound_sr1 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($bp_window{$data[9]} == 2 and $$ref_boundary2[3] >= ($positiona2-$positiona1+101));
						$left_bound_sr2 = $positiona1 + $$ref_boundary1[1] + $cut_sr if ($bp_window{$data[4]} == 1 and $$ref_boundary1[1] < ($positiona2-$positiona1+101));
						$left_bound_sr2 = $$ref_boundary1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($bp_window{$data[4]} == 2 and $$ref_boundary1[1] >= ($positiona2-$positiona1+101));
						$right_bound_sr2 = $positiona1 + $$ref_boundary2[0] if ($bp_window{$data[8]} == 1 and $$ref_boundary2[0] < ($positiona2-$positiona1+101));
						$right_bound_sr2 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3 if ($bp_window{$data[8]} == 2 and $$ref_boundary2[0] >= ($positiona2-$positiona1+101));
						$left_bound_sr1 = $data[5] unless ($$ref_boundary1[2]);
						$right_bound_sr1 = $data[9] unless ($$ref_boundary2[3]);
						$left_bound_sr2 = $data[4] unless ($$ref_boundary1[1]);
						$right_bound_sr2 = $data[8] unless ($$ref_boundary2[0]);
						my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
						my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
						if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $left_bound_sr2 and $right_bound_sr1 and $right_bound_sr2)
						{
							my $reads0 = join("\t", @{$$ref_bpread{0}});
							my $reads1 = join("\t", @{$$ref_bpread{1}});
							print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
							$result_sr[$sri][0] = $data[0];
							$result_sr[$sri][1] = $data[1];
							$result_sr[$sri][2] = $data[2];
							$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
							$result_sr[$sri][4] = $data[3];
							$result_sr[$sri][5] = $left_bound_sr2;
							$result_sr[$sri][6] = $left_bound_sr1;
							$result_sr[$sri][7] = $del_size;
							$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
							$sri++;
						}
						else
						{
							my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
							my ($mpd1, $mpd2) = split (/\//, $data[2]);
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $left_bound_sr = $left_bound_sr1;
								my $right_bound_sr = $right_bound_sr1;
								my $reads = join("\t", @{$$ref_bpread{0}});
								if ($$ref_bpread{1})
								{
									$reads = join("\t", @{$$ref_bpread{1}});
								}
								else
								{
									$left_bound_sr = $positiona1 + $$ref_boundary1[0];
									$right_bound_sr = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
								}
								if ($left_bound_sr and $right_bound_sr)
								{
									my $size = $right_bound_sr - $left_bound_sr;
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "tandem_dup\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
							if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
							{
								my $left_bound_sr = $left_bound_sr2;
								my $right_bound_sr = $right_bound_sr2;
								my $size = $right_bound_sr - $left_bound_sr - 1;
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "del\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
						
					}
					else
					{
						my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
						$left_bound_sr1 = $positiona1 + $$ref_boundary1[0] if ($bp_window{$data[5]} == 1 and $$ref_boundary1[0] < ($positiona2-$positiona1+101));
						$left_bound_sr1 = $$ref_boundary1[0] - ($positiona2-$positiona1+101) + $positiona3 if ($bp_window{$data[5]} == 2 and $$ref_boundary1[0] >= ($positiona2-$positiona1+101));
						$right_bound_sr1 = $positiona1 + $$ref_boundary2[1] + $cut_sr if ($bp_window{$data[9]} == 1 and $$ref_boundary2[1] < ($positiona2-$positiona1+101));
						$right_bound_sr1 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($bp_window{$data[9]} == 2 and $$ref_boundary2[1] >= ($positiona2-$positiona1+101));
						$left_bound_sr2 = $positiona1 + $$ref_boundary1[3] + $cut_sr if ($bp_window{$data[4]} == 1 and $$ref_boundary1[3] < ($positiona2-$positiona1+101));
						$left_bound_sr2 = $$ref_boundary1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($bp_window{$data[4]} == 2 and $$ref_boundary1[3] >= ($positiona2-$positiona1+101));
						$right_bound_sr2 = $positiona1 + $$ref_boundary2[2] if ($bp_window{$data[8]} == 1 and $$ref_boundary2[2] < ($positiona2-$positiona1+101));
						$right_bound_sr2 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3 if ($bp_window{$data[8]} == 2 and $$ref_boundary2[2] >= ($positiona2-$positiona1+101));
						$left_bound_sr1 = $data[5] unless ($$ref_boundary1[0]);
						$right_bound_sr1 = $data[9] unless ($$ref_boundary2[1]);
						$left_bound_sr2 = $data[4] unless ($$ref_boundary1[3]);
						$right_bound_sr2 = $data[8] unless ($$ref_boundary2[2]);
						my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
						my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
						if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $left_bound_sr2 and $right_bound_sr1 and $right_bound_sr2)
						{
							my $reads0 = join("\t", @{$$ref_bpread{1}});
							my $reads1 = join("\t", @{$$ref_bpread{0}});
							print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
							$result_sr[$sri][0] = $data[0];
							$result_sr[$sri][1] = $data[1];
							$result_sr[$sri][2] = $data[2];
							$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
							$result_sr[$sri][4] = $data[3];
							$result_sr[$sri][5] = $left_bound_sr2;
							$result_sr[$sri][6] = $left_bound_sr1;
							$result_sr[$sri][7] = $del_size;
							$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
							$sri++;
						}
						else
						{
							my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
							my ($mpd1, $mpd2) = split (/\//, $data[2]);
							if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
							{
								my $left_bound_sr = $left_bound_sr1;
								my $right_bound_sr = $right_bound_sr1;
								my $size = $right_bound_sr - $left_bound_sr;
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "tandem_dup\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $left_bound_sr = $left_bound_sr2;
								my $right_bound_sr = $right_bound_sr2;
								my $reads = join("\t", @{$$ref_bpread{0}});
								if ($$ref_bpread{1})
								{
									$reads = join("\t", @{$$ref_bpread{1}});
								}
								else
								{
									$left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									$right_bound_sr = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
								}
								if ($left_bound_sr and $right_bound_sr)
								{
									my $size = $right_bound_sr - $left_bound_sr - 1;
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "del\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
						}
					}
				}
			}
			else
			{
				my (@src_return1, $left_bound_sr1, $right_bound_sr1, @src_return2, $left_bound_sr2, $right_bound_sr2);
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[5]} == 1 and $bp_window{$data[9]} == 2)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
						$left_bound_sr1 = $positiona1 + $src_return1[0];
						$right_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
						$left_bound_sr1 = $data[5] unless ($src_return1[0]);
						$right_bound_sr1 = $data[9] unless ($src_return1[3]);
					}
				}
				else
				{
					if ($bp_window{$data[5]} == 1 and $bp_window{$data[9]} == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart < ($positiona2-$positiona1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr1 = $positiona1 + $src_return1[0];
							$right_bound_sr1 = $positiona1 + $src_return1[3] + $cut_sr;
							$left_bound_sr1 = $data[5] unless ($src_return1[0]);
							$right_bound_sr1 = $data[9] unless ($src_return1[3]);
						}
					}
					if ($bp_window{$data[5]} == 2 and $bp_window{$data[9]} == 2)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start > ($positiona2-$positiona1+100) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr1 = $src_return1[0] - ($positiona2-$positiona1+101) + $positiona3;
							$right_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$left_bound_sr1 = $data[5] unless ($src_return1[0]);
							$right_bound_sr1 = $data[9] unless ($src_return1[3]);
						}
					}
				}
				
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[1]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[4]} == 1 and $bp_window{$data[8]} == 2)
				{
					@discord_sr2 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr2, $a;
						}
					}
					if (@discord_sr2 >= $support_reads)
					{
						@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
						$left_bound_sr2 = $positionb1 + $src_return2[1] + $cut_sr;
						$right_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
						$left_bound_sr2 = $data[4] unless ($src_return2[1]);
						$right_bound_sr2 = $data[8] unless ($src_return2[2]);
					}
				}
				else
				{
					if ($bp_window{$data[4]} == 1 and $bp_window{$data[8]} == 1)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start < ($positionb2-$positionb1) and $mstart < ($positionb2-$positionb1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr2 = $positionb1 + $src_return2[1] + $cut_sr;
							$right_bound_sr2 = $positionb1 + $src_return2[2];
							$left_bound_sr2 = $data[4] unless ($src_return2[1]);
							$right_bound_sr2 = $data[8] unless ($src_return2[2]);
						}
					}
					if ($bp_window{$data[4]} == 2 and $bp_window{$data[8]} == 2)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start > ($positionb2-$positionb1+100) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr2 = $src_return2[1] - ($positionb2-$positionb1+101) + $positionb3 + $cut_sr;
							$right_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
							$left_bound_sr2 = $data[4] unless ($src_return2[1]);
							$right_bound_sr2 = $data[8] unless ($src_return2[2]);
						}
					}
					
				}
				
				my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
				my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
				if ($src_return1[4] >= $support_reads and $src_return2[4] >= $support_reads and $left_bound_sr1 and $left_bound_sr2 and $right_bound_sr1 and $right_bound_sr2)
				{
					my $reads0 = join("\t", @{$src_return2[5]{0}});
					my $reads1 = join("\t", @{$src_return1[5]{0}});
					print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
					$result_sr[$sri][0] = $data[0];
					$result_sr[$sri][1] = $data[1];
					$result_sr[$sri][2] = $data[2];
					$result_sr[$sri][3] = "$src_return1[4]/$src_return2[4]";
					$result_sr[$sri][4] = $data[3];
					$result_sr[$sri][5] = $left_bound_sr2;
					$result_sr[$sri][6] = $left_bound_sr1;
					$result_sr[$sri][7] = $del_size;
					$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
					$sri++;
				}
				else
				{
					my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
					my ($mpd1, $mpd2) = split (/\//, $data[2]);
					if ($src_return1[4] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
					{
						my $left_bound_sr = $left_bound_sr1;
						my $right_bound_sr = $right_bound_sr1;
						my $size = $right_bound_sr - $left_bound_sr;
						my $reads = join("\t", @{$src_return1[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "tandem_dup\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
					if ($src_return2[4] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
					{
						my $left_bound_sr = $left_bound_sr2;
						my $right_bound_sr = $right_bound_sr2;
						my $size = $right_bound_sr - $left_bound_sr - 1;
						my $reads = join("\t", @{$src_return2[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "del\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
				}
			}
		}
		if ($data[0] =~ /inssu/)
		{
			next unless ($cluster_region{$cl[0]} and $cluster_region{$cl[1]});
			my ($trash, $positiona1, $positiona2, $trash, $positiona3, $positiona4, $orientation) = split (/__/, $cluster_region{$cl[0]});
			my ($trash, $positionb1, $positionb2, $trash, $positionb3, $positionb4, $orientation) = split (/__/, $cluster_region{$cl[1]});
			if (&covered($data[8], $data[8], $positiona1, $positiona2))
			{
				$bp_window{$data[8]} = 1;
			}
			if (&covered($data[8], $data[8], $positiona3, $positiona4))
			{
				$bp_window{$data[8]} = 2;
			}
			if (&covered($data[4], $data[4], $positiona1, $positiona2))
			{
				$bp_window{$data[4]} = 1;
			}
			if (&covered($data[4], $data[4], $positiona3, $positiona4))
			{
				$bp_window{$data[4]} = 2;
			}
			if (&covered($data[9], $data[9], $positionb1, $positionb2))
			{
				$bp_window{$data[9]} = 1;
			}
			if (&covered($data[9], $data[9], $positionb3, $positionb4))
			{
				$bp_window{$data[9]} = 2;
			}
			if (&covered($data[5], $data[5], $positionb1, $positionb2))
			{
				$bp_window{$data[5]} = 1;
			}
			if (&covered($data[5], $data[5], $positionb3, $positionb4))
			{
				$bp_window{$data[5]} = 2;
			}
			#print "$cluster_region{$cl[0]}\n$cluster_region{$cl[1]}\n$data[4]\t$bp_window{$data[4]}\n$data[5]\t$bp_window{$data[5]}\n$data[8]\t$bp_window{$data[8]}\n$data[9]\t$bp_window{$data[9]}\n";
			if ($cluster_region{$cl[0]} eq $cluster_region{$cl[1]})
			{
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				@discord_sr1 = undef;
				for my $a (@alignments)
				{
					my @data1 = split (/\t/, $a);
					my $start = $data1[3];
					my $strand = 1;
					$strand = -1 if ($data1[1] =~ /r/);
					my $mseqid = $data1[6];
					my $mstart = $data1[7];
					my $mstrand = 1;
					$mstrand = -1 if ($data1[1] =~ /R/);
					my $isize = $data1[8];
					if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand)
					{
						#print "$readname\t$start $mstart\t$isize\t$strand $mstrand\n";
						push @discord_sr1, $a;
					}
				}
				if (@discord_sr1 >= $support_reads)
				{
					my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 2);
					#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n@{$$ref_bpread{0}}\n@{$$ref_bpread{1}}\n";
					if ($$ref_boundary2[2] > $$ref_boundary2[0])
					{
						my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
						$left_bound_sr1 = $positiona1 + $$ref_boundary1[0] if ($bp_window{$data[8]} == 1 and $$ref_boundary1[0] < ($positiona2-$positiona1+101));
						$left_bound_sr1 = $$ref_boundary1[0] - ($positiona2-$positiona1+101) + $positiona3 if ($bp_window{$data[8]} == 2 and $$ref_boundary1[0] >= ($positiona2-$positiona1+101));
						$right_bound_sr1 = $positiona1 + $$ref_boundary2[1] + $cut_sr if ($bp_window{$data[4]} == 1 and $$ref_boundary2[1] < ($positiona2-$positiona1+101));
						$right_bound_sr1 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($bp_window{$data[4]} == 2 and $$ref_boundary2[1] >= ($positiona2-$positiona1+101));
						$left_bound_sr2 = $positiona1 + $$ref_boundary1[3] + $cut_sr if ($bp_window{$data[9]} == 1 and $$ref_boundary1[3] < ($positiona2-$positiona1+101));
						$left_bound_sr2 = $$ref_boundary1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($bp_window{$data[9]} == 2 and $$ref_boundary1[3] >= ($positiona2-$positiona1+101));
						$right_bound_sr2 = $positiona1 + $$ref_boundary2[2] if ($bp_window{$data[5]} == 1 and $$ref_boundary2[2] < ($positiona2-$positiona1+101));
						$right_bound_sr2 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3 if ($bp_window{$data[5]} == 2 and $$ref_boundary2[2] >= ($positiona2-$positiona1+101));
						$left_bound_sr1 = $data[8] unless ($$ref_boundary1[0]);
						$right_bound_sr1 = $data[4] unless ($$ref_boundary2[1]);
						$left_bound_sr2 = $data[9] unless ($$ref_boundary1[3]);
						$right_bound_sr2 = $data[5] unless ($$ref_boundary2[2]);
						my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
						my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
						if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $left_bound_sr2 and $right_bound_sr1 and $right_bound_sr2)
						{
							my $reads0 = join("\t", @{$$ref_bpread{0}});
							my $reads1 = join("\t", @{$$ref_bpread{1}});
							print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
							$result_sr[$sri][0] = $data[0];
							$result_sr[$sri][1] = $data[1];
							$result_sr[$sri][2] = $data[2];
							$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
							$result_sr[$sri][4] = $data[3];
							$result_sr[$sri][5] = $right_bound_sr1;
							$result_sr[$sri][6] = $right_bound_sr2;
							$result_sr[$sri][7] = $del_size;
							$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
							$sri++;
						}
						else
						{
							my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
							my ($mpd1, $mpd2) = split (/\//, $data[2]);
							if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
							{
								my $left_bound_sr = $left_bound_sr1;
								my $right_bound_sr = $right_bound_sr1;
								my $size = $right_bound_sr - $left_bound_sr;
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "tandem_dup\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
							if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
							{
								my $left_bound_sr = $left_bound_sr2;
								my $right_bound_sr = $right_bound_sr2;
								my $size = $right_bound_sr - $left_bound_sr - 1;
								my $reads = join("\t", @{$$ref_bpread{1}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "del\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
					}
					else
					{
						my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
						$left_bound_sr1 = $positiona1 + $$ref_boundary1[2] if ($bp_window{$data[8]} == 1 and $$ref_boundary1[2] < ($positiona2-$positiona1+101));
						$left_bound_sr1 = $$ref_boundary1[2] - ($positiona2-$positiona1+101) + $positiona3 if ($bp_window{$data[8]} == 2 and $$ref_boundary1[2] >= ($positiona2-$positiona1+101));
						$right_bound_sr1 = $positiona1 + $$ref_boundary2[3] + $cut_sr if ($bp_window{$data[4]} == 1 and $$ref_boundary2[3] < ($positiona2-$positiona1+101));
						$right_bound_sr1 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($bp_window{$data[4]} == 2 and $$ref_boundary2[3] >= ($positiona2-$positiona1+101));
						$left_bound_sr2 = $positiona1 + $$ref_boundary1[1] + $cut_sr if ($bp_window{$data[9]} == 1 and $$ref_boundary1[1] < ($positiona2-$positiona1+101));
						$left_bound_sr2 = $$ref_boundary1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($bp_window{$data[9]} == 2 and $$ref_boundary1[1] >= ($positiona2-$positiona1+101));
						$right_bound_sr2 = $positiona1 + $$ref_boundary2[0] if ($bp_window{$data[5]} == 1 and $$ref_boundary2[0] < ($positiona2-$positiona1+101));
						$right_bound_sr2 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3 if ($bp_window{$data[5]} == 2 and $$ref_boundary2[0] >= ($positiona2-$positiona1+101));
						$left_bound_sr1 = $data[8] unless ($$ref_boundary1[2]);
						$right_bound_sr1 = $data[4] unless ($$ref_boundary2[3]);
						$left_bound_sr2 = $data[9] unless ($$ref_boundary1[1]);
						$right_bound_sr2 = $data[5] unless ($$ref_boundary2[0]);
						my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
						my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
						if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $left_bound_sr2 and $right_bound_sr1 and $right_bound_sr2)
						{
							my $reads0 = join("\t", @{$$ref_bpread{1}});
							my $reads1 = join("\t", @{$$ref_bpread{0}});
							print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
							$result_sr[$sri][0] = $data[0];
							$result_sr[$sri][1] = $data[1];
							$result_sr[$sri][2] = $data[2];
							$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
							$result_sr[$sri][4] = $data[3];
							$result_sr[$sri][5] = $right_bound_sr1;
							$result_sr[$sri][6] = $right_bound_sr2;
							$result_sr[$sri][7] = $del_size;
							$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
							$sri++;
						}
						else
						{
							my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
							my ($mpd1, $mpd2) = split (/\//, $data[2]);
							if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
							{
								my $left_bound_sr = $left_bound_sr1;
								my $right_bound_sr = $right_bound_sr1;
								my $size = $right_bound_sr - $left_bound_sr;
								my $reads = join("\t", @{$$ref_bpread{1}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "tandem_dup\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
							if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
							{
								my $left_bound_sr = $left_bound_sr2;
								my $right_bound_sr = $right_bound_sr2;
								my $size = $right_bound_sr - $left_bound_sr - 1;
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "del\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
					}
				}
			}
			else
			{
				my (@src_return1, $left_bound_sr1, $right_bound_sr1, @src_return2, $left_bound_sr2, $right_bound_sr2);
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[8]} == 1 and $bp_window{$data[4]} == 2)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
						$left_bound_sr1 = $positiona1 + $src_return1[0];
						$right_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
						$left_bound_sr1 = $data[8] unless ($src_return1[0]);
						$right_bound_sr1 = $data[4] unless ($src_return1[3]);
					}
				}
				else
				{
					if ($bp_window{$data[8]} == 1 and $bp_window{$data[4]} == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart < ($positiona2-$positiona1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr1 = $positiona1 + $src_return1[0];
							$right_bound_sr1 = $positiona1 + $src_return1[3] + $cut_sr;
							$left_bound_sr1 = $data[8] unless ($src_return1[0]);
							$right_bound_sr1 = $data[4] unless ($src_return1[3]);
						}
					}
					if ($bp_window{$data[8]} == 2 and $bp_window{$data[4]} == 2)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start > ($positiona2-$positiona1+100) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr1 = $src_return1[0] - ($positiona2-$positiona1+101) + $positiona3;
							$right_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$left_bound_sr1 = $data[8] unless ($src_return1[0]);
							$right_bound_sr1 = $data[4] unless ($src_return1[3]);
						}
					}
				}
							
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[1]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[9]} == 1 and $bp_window{$data[5]} == 2)
				{
					@discord_sr2 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr2, $a;
						}
					}
					if (@discord_sr2 >= $support_reads)
					{
						@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
						$left_bound_sr2 = $positionb1 + $src_return2[1] + $cut_sr;
						$right_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
						$left_bound_sr2 = $data[9] unless ($src_return2[1]);
						$right_bound_sr2 = $data[5] unless ($src_return2[2]);
					}
				}
				else
				{
					if ($bp_window{$data[9]} == 1 and $bp_window{$data[5]} == 1)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start < ($positionb2-$positionb1) and $mstart < ($positionb2-$positionb1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr2 = $positionb1 + $src_return2[1] + $cut_sr;
							$right_bound_sr2 = $positionb1 + $src_return2[2];
							$left_bound_sr2 = $data[9] unless ($src_return2[1]);
							$right_bound_sr2 = $data[5] unless ($src_return2[2]);
						}
					}
					if ($bp_window{$data[9]} == 2 and $bp_window{$data[5]} == 2)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start > ($positionb2-$positionb1+100) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr2 = $src_return2[1] - ($positionb2-$positionb1+101) + $positionb3 + $cut_sr;
							$right_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
							$left_bound_sr2 = $data[9] unless ($src_return2[1]);
							$right_bound_sr2 = $data[5] unless ($src_return2[2]);
						}
					}
				}
				
				my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
				my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
				if ($src_return1[4] >= $support_reads and $src_return2[4] >= $support_reads and $left_bound_sr1 and $left_bound_sr2 and $right_bound_sr1 and $right_bound_sr2)
				{
					my $reads0 = join("\t", @{$src_return1[5]{0}});
					my $reads1 = join("\t", @{$src_return2[5]{0}});
					print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
					$result_sr[$sri][0] = $data[0];
					$result_sr[$sri][1] = $data[1];
					$result_sr[$sri][2] = $data[2];
					$result_sr[$sri][3] = "$src_return1[4]/$src_return2[4]";
					$result_sr[$sri][4] = $data[3];
					$result_sr[$sri][5] = $right_bound_sr1;
					$result_sr[$sri][6] = $right_bound_sr2;
					$result_sr[$sri][7] = $del_size;
					$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
					$sri++;
				}
				else
				{
					my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
					my ($mpd1, $mpd2) = split (/\//, $data[2]);
					if ($src_return1[4] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
					{
						my $left_bound_sr = $left_bound_sr1;
						my $right_bound_sr = $right_bound_sr1;
						my $size = $right_bound_sr - $left_bound_sr;
						my $reads = join("\t", @{$src_return1[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "tandem_dup\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
					if ($src_return2[4] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
					{
						my $left_bound_sr = $left_bound_sr2;
						my $right_bound_sr = $right_bound_sr2;
						my $size = $right_bound_sr - $left_bound_sr - 1;
						my $reads = join("\t", @{$src_return2[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "del\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
				}
			}
		}
		if ($data[0] =~ /insod/)
		{
			next unless ($cluster_region{$cl[0]} and $cluster_region{$cl[1]});
			my ($trash, $positiona1, $positiona2, $trash, $positiona3, $positiona4, $orientationa) = split (/__/, $cluster_region{$cl[0]});
			my ($trash, $positionb1, $positionb2, $trash, $positionb3, $positionb4, $orientationb) = split (/__/, $cluster_region{$cl[1]});
			if (&covered($data[4], $data[4], $positiona1, $positiona2))
			{
				$bp_window{$data[4]} = 1;
			}
			if (&covered($data[4], $data[4], $positiona3, $positiona4))
			{
				$bp_window{$data[4]} = 2;
			}
			if (&covered($data[9], $data[9], $positiona1, $positiona2))
			{
				$bp_window{$data[9]} = 1;
			}
			if (&covered($data[9], $data[9], $positiona3, $positiona4))
			{
				$bp_window{$data[9]} = 2;
			}
			if (&covered($data[5], $data[5], $positionb1, $positionb2))
			{
				$bp_window{$data[5]} = 1;
			}
			if (&covered($data[5], $data[5], $positionb3, $positionb4))
			{
				$bp_window{$data[5]} = 2;
			}
			if (&covered($data[8], $data[8], $positionb1, $positionb2))
			{
				$bp_window{$data[8]} = 1;
			}
			if (&covered($data[8], $data[8], $positionb3, $positionb4))
			{
				$bp_window{$data[8]} = 2;
			}
			#print "$cluster_region{$cl[0]}\n$cluster_region{$cl[1]}\n$data[4]\t$bp_window{$data[4]}\n$data[5]\t$bp_window{$data[5]}\n$data[8]\t$bp_window{$data[8]}\n$data[9]\t$bp_window{$data[9]}\n";
			if ($cluster_region{$cl[0]} eq $cluster_region{$cl[1]})
			{
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[4]} ne $bp_window{$data[9]} or $bp_window{$data[5]} ne $bp_window{$data[8]})
				{
					if ($orientationa == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "1 $start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
					}
					@discord_sr2 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand)
						{
							#print "2 $start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr2, $a;
						}
					}
					# small deletion, large insertion, close, 3 break points in one window
					if (@discord_sr1 >= $support_reads and @discord_sr2 >= $support_reads and $bp_window{$data[5]} == 1 and $bp_window{$data[8]} == 1)
					{
						my @src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
						my @src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
						#print "@src_return1\n@src_return2\n";
						my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
						$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr if ($src_return1[1] < ($positiona2-$positiona1+101));
						$right_bound_sr1 = $positiona4 - ($src_return1[2] - ($positiona2-$positiona1+101)) if ($src_return1[2] >= ($positiona2-$positiona1+101));
						$left_bound_sr2 = $positiona1 + $src_return2[0] if ($src_return2[0] < ($positiona2-$positiona1+101));
						$right_bound_sr2 = $positiona1 + $src_return2[2] if ($src_return2[2] < ($positiona2-$positiona1+101));
						$left_bound_sr1 = $data[4] unless ($src_return1[1]);
						$right_bound_sr1 = $data[9] unless ($src_return1[2]);
						$left_bound_sr2 = $data[5] unless ($src_return2[0]);
						$right_bound_sr2 = $data[8] unless ($src_return2[2]);
						my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
						my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
						if ($src_return1[4] >= $support_reads and $src_return2[4] >= $support_reads and $left_bound_sr1 and $left_bound_sr2 and $right_bound_sr1 and $right_bound_sr2)
						{
							my $reads0 = join("\t", @{$src_return1[5]{0}});
							my $reads1 = join("\t", @{$src_return2[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
							$result_sr[$sri][0] = $data[0];
							$result_sr[$sri][1] = $data[1];
							$result_sr[$sri][2] = $data[2];
							$result_sr[$sri][3] = "$src_return1[4]/$src_return2[4]";
							$result_sr[$sri][4] = $data[3];
							$result_sr[$sri][5] = $left_bound_sr1;
							$result_sr[$sri][6] = $left_bound_sr2;
							$result_sr[$sri][7] = $del_size;
							$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
							$sri++;
						}
						else
						{
							my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
							my ($mpd1, $mpd2) = split (/\//, $data[2]);
							if ($src_return1[4] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
							{
								my $left_bound_sr = $left_bound_sr1;
								my $right_bound_sr = $right_bound_sr1;
								my $size = $right_bound_sr - $left_bound_sr;
								my $reads = join("\t", @{$src_return1[5]{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
							if ($src_return2[4] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
							{
								#my $left_bound_sr = $left_bound_sr2;
								#my $right_bound_sr = $right_bound_sr2;
								#my $size = $right_bound_sr - $left_bound_sr;
								#my $reads = join("\t", @{$src_return2[5]{0}});
								#print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								#$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								#$sri++;
							}
						}
					}
					# small deletion, small insertion, far, 2 break points in one window, the other 2 in the other window
					elsif (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 1);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						if ($$ref_boundary1[2] > $$ref_boundary1[0])
						{
							my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
							$left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr if ($$ref_boundary1[1] < ($positiona2-$positiona1+101));
							$right_bound_sr1 = $positiona4 - ($$ref_boundary2[0] - ($positiona2-$positiona1+101)) if ($$ref_boundary2[0] >= ($positiona2-$positiona1+101));
							$left_bound_sr2 = $positiona1 + $$ref_boundary1[2] if ($$ref_boundary1[2] < ($positiona2-$positiona1+101));
							$right_bound_sr2 = $positiona4 - ($$ref_boundary2[3] - ($positiona2-$positiona1+101) + $cut_sr) if ($$ref_boundary2[3] >= ($positiona2-$positiona1+101));
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[1]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[0]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[2]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[3]);
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $left_bound_sr2 and $right_bound_sr1 and $right_bound_sr2)
							{
								my $reads0 = join("\t", @{$$ref_bpread{0}});
								my $reads1 = join("\t", @{$$ref_bpread{1}});
								print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr1;
								$result_sr[$sri][6] = $left_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
								{
									my $left_bound_sr = $left_bound_sr1;
									my $right_bound_sr = $right_bound_sr1;
									my $size = $right_bound_sr - $left_bound_sr;
									my $reads = join("\t",@{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
								if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
								{
									my $left_bound_sr = $left_bound_sr2;
									my $right_bound_sr = $right_bound_sr2;
									my $size = $right_bound_sr - $left_bound_sr;
									my $reads = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
						}
						else
						{
							my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
							$left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr if ($$ref_boundary1[3] < ($positiona2-$positiona1+101));
							$right_bound_sr1 = $positiona4 - ($$ref_boundary2[2] - ($positiona2-$positiona1+101)) if ($$ref_boundary2[2] >= ($positiona2-$positiona1+101));
							$left_bound_sr2 = $positiona1 + $$ref_boundary1[0] if ($$ref_boundary1[0] < ($positiona2-$positiona1+101));
							$right_bound_sr2 = $positiona4 - ($$ref_boundary2[1] - ($positiona2-$positiona1+101) + $cut_sr) if ($$ref_boundary2[1] >= ($positiona2-$positiona1+101));
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[3]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[2]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[0]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[1]);
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
							{
								my $reads0 = join("\t", @{$$ref_bpread{1}});
								my $reads1 = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr1;
								$result_sr[$sri][6] = $left_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
								{
									my $left_bound_sr = $left_bound_sr1;
									my $right_bound_sr = $right_bound_sr1;
									my $size = $right_bound_sr - $left_bound_sr;
									my $reads = join("\t",@{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
								if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
								{
									my $left_bound_sr = $left_bound_sr2;
									my $right_bound_sr = $right_bound_sr2;
									my $size = $right_bound_sr - $left_bound_sr;
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
						}
					}
					# large deletion, small insertion, close, bp_info 0, or small deletion, small insertion, far, bp_info 0
					elsif (@discord_sr2 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 1);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						if ($$ref_boundary1[2] > $$ref_boundary1[0])
						{
							my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
							$left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr if ($$ref_boundary1[1] < ($positiona2-$positiona1+101));
							$right_bound_sr1 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($$ref_boundary2[1] >= ($positiona2-$positiona1+101));
							$left_bound_sr2 = $$ref_boundary1[2] - ($positiona2-$positiona1+101) + $positiona3 if ($$ref_boundary1[2] >= ($positiona2-$positiona1+101));
							$right_bound_sr2 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3 if ($$ref_boundary2[2] >= ($positiona2-$positiona1+101));
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[1]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[1]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[2]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[2]);
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
							{
								my $reads0 = join("\t", @{$$ref_bpread{0}});
								my $reads1 = join("\t", @{$$ref_bpread{1}});
								print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr1;
								$result_sr[$sri][6] = $left_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
								{
									my $left_bound_sr = $left_bound_sr1;
									my $right_bound_sr = $right_bound_sr1;
									my $size = $right_bound_sr - $left_bound_sr;
									my $reads = join("\t",@{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
								if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
								{
									#my $left_bound_sr = $left_bound_sr2;
									#my $right_bound_sr = $right_bound_sr2;
									#my $size = $right_bound_sr - $left_bound_sr;
									#my $reads = join("\t", @{$$ref_bpread{1}});
									#print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									#$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									#$sri++;
								}
							}
						}
						else
						{
							my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
							$left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr if ($$ref_boundary1[3] < ($positiona2-$positiona1+101));
							$right_bound_sr1 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr if ($$ref_boundary2[3] >= ($positiona2-$positiona1+101));
							$left_bound_sr2 = $$ref_boundary1[0] - ($positiona2-$positiona1+101) + $positiona3 if ($$ref_boundary1[0] >= ($positiona2-$positiona1+101));
							$right_bound_sr2 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3 if ($$ref_boundary2[0] >= ($positiona2-$positiona1+101));
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[3]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[3]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[0]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[0]);
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
							{
								my $reads0 = join("\t", @{$$ref_bpread{1}});
								my $reads1 = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr1;
								$result_sr[$sri][6] = $left_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
								{
									my $left_bound_sr = $left_bound_sr1;
									my $right_bound_sr = $right_bound_sr1;
									my $size = $right_bound_sr - $left_bound_sr;
									my $reads = join("\t",@{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
								if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
								{
									#my $left_bound_sr = $left_bound_sr2;
									#my $right_bound_sr = $right_bound_sr2;
									#my $size = $right_bound_sr - $left_bound_sr;
									#my $reads = join("\t", @{$$ref_bpread{0}});
									#print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									#$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									#$sri++;
								}
							}
						}
					}
				}
				# small deletion, small insertion, close, 4 break points in one window
				else
				{
					if ($bp_window{$data[4]} == 1 and $bp_window{$data[5]} == 1 and $bp_window{$data[8]} == 1 and $bp_window{$data[9]} == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1+101) and $mstart < ($positiona2-$positiona1+101))
							{
								#print "$i $start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1,);
							#print "@{$ref_boundary1}\n@{$ref_boundary2}\n";
							# 2 cluster event
							if ($$ref_support_sr[1])
							{
								if ($$ref_boundary1[1]<$$ref_boundary1[3])
								{
									my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
									$left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									$right_bound_sr1 = $positiona1 + $$ref_boundary2[1] + $cut_sr;
									$left_bound_sr2 = $positiona1 + $$ref_boundary1[2];
									$right_bound_sr2 = $positiona1 + $$ref_boundary2[2];
									$left_bound_sr1 = $data[4] unless ($$ref_boundary1[1]);
									$right_bound_sr1 = $data[9] unless ($$ref_boundary2[1]);
									$left_bound_sr2 = $data[5] unless ($$ref_boundary1[2]);
									$right_bound_sr2 = $data[8] unless ($$ref_boundary2[2]);
									my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = $data[0];
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $left_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t",@{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{1}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
								else
								{
									my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
									$left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr if ($$ref_boundary1[3] < ($positiona2-$positiona1+101));
									$right_bound_sr1 = $positiona1 + $$ref_boundary2[3] + $cut_sr if ($$ref_boundary2[3] < ($positiona2-$positiona1+101));
									$left_bound_sr2 = $positiona1 + $$ref_boundary1[0] if ($$ref_boundary1[0] < ($positiona2-$positiona1+101));
									$right_bound_sr2 = $positiona1 + $$ref_boundary2[0] if ($$ref_boundary2[0] < ($positiona2-$positiona1+101));
									$left_bound_sr1 = $data[4] unless ($$ref_boundary1[3]);
									$right_bound_sr1 = $data[9] unless ($$ref_boundary2[3]);
									$left_bound_sr2 = $data[5] unless ($$ref_boundary1[0]);
									$right_bound_sr2 = $data[8] unless ($$ref_boundary2[0]);
									my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = $data[0];
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $left_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $reads = join("\t", @{$$ref_bpread{1}});
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
							}
							elsif ($$ref_support_sr[0] >= $support_reads)
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								my $left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr if ($$ref_boundary1[1] < ($positiona2-$positiona1+101));
								my $right_bound_sr = $positiona1 + $$ref_boundary2[1] + $cut_sr if ($$ref_boundary2[1] < ($positiona2-$positiona1+101));
								my $reads = join("\t", @{$$ref_bpread{0}});
								if ($left_bound_sr and $right_bound_sr)
								{
									my $size = $right_bound_sr - $left_bound_sr;
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
						}
					}
					# 4 break points in window 2
					else
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start >= ($positiona2-$positiona1+101) and $mstart >= ($positiona2-$positiona1+101))
							{
								#print "$i $start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1,);
							#print "@{$ref_boundary1}\n@{$ref_boundary2}\n";
							# 2 cluster event
							if ($$ref_support_sr[1])
							{
								if ($$ref_boundary1[1]<$$ref_boundary1[3])
								{
									my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
									$left_bound_sr1 = $$ref_boundary1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									$right_bound_sr1 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									$left_bound_sr2 = $$ref_boundary1[2] - ($positiona2-$positiona1+101) + $positiona3;
									$right_bound_sr2 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
									$left_bound_sr1 = $data[4] unless ($$ref_boundary1[1]);
									$right_bound_sr1 = $data[9] unless ($$ref_boundary2[1]);
									$left_bound_sr2 = $data[5] unless ($$ref_boundary1[2]);
									$right_bound_sr2 = $data[8] unless ($$ref_boundary2[2]);
									my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = $data[0];
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $left_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t",@{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{1}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
								else
								{
									my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
									$left_bound_sr1 = $$ref_boundary1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									$right_bound_sr1 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3;
									$left_bound_sr2 = $$ref_boundary1[0] - ($positiona2-$positiona1+101) + $positiona3;
									$right_bound_sr2 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
									$left_bound_sr1 = $data[4] unless ($$ref_boundary1[3]);
									$right_bound_sr1 = $data[9] unless ($$ref_boundary2[3]);
									$left_bound_sr2 = $data[5] unless ($$ref_boundary1[0]);
									$right_bound_sr2 = $data[8] unless ($$ref_boundary2[0]);
									my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = $data[0];
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $left_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $reads = join("\t", @{$$ref_bpread{1}});
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
							}
							elsif ($$ref_support_sr[0] >= $support_reads)
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								my $left_bound_sr = $$ref_boundary1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
								my $right_bound_sr = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
								my $reads = join("\t", @{$$ref_bpread{0}});
								if ($left_bound_sr and $right_bound_sr)
								{
									my $size = $right_bound_sr - $left_bound_sr;
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
						}
					}
				}
			}
			else
			{
				my (@src_return1, $left_bound_sr1, $right_bound_sr1, @src_return2, $left_bound_sr2, $right_bound_sr2);
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[4]} ne $bp_window{$data[9]})
				{
					if ($orientationa == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
							$right_bound_sr1 = $positiona4 - ($src_return1[2] - ($positiona2-$positiona1+101));
							$left_bound_sr1 = $data[4] unless ($src_return1[1]);
							$right_bound_sr1 = $data[9] unless ($src_return1[2]);
						}
					}
					else
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
							#print "@src_return1\n";
							$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
							$right_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$left_bound_sr1 = $data[4] unless ($src_return1[1]);
							$right_bound_sr1 = $data[9] unless ($src_return1[3]);
						}
					}
				}
				else
				{
					if ($bp_window{$data[4]} == 1 and $bp_window{$data[9]} == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart < ($positiona2-$positiona1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
							#print "@src_return1\n";
							$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
							$right_bound_sr1 = $positiona1 + $src_return1[3] + $cut_sr;
							$left_bound_sr1 = $data[4] unless ($src_return1[1]);
							$right_bound_sr1 = $data[9] unless ($src_return1[3]);
						}
					}
					if ($bp_window{$data[4]} == 2 and $bp_window{$data[9]} == 2)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start > ($positiona2-$positiona1+100) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
							#print "@src_return1\n";
							$left_bound_sr1 = $src_return1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$right_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$left_bound_sr1 = $data[4] unless ($src_return1[1]);
							$right_bound_sr1 = $data[9] unless ($src_return1[3]);
						}
					}
				}
				
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[1]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[5]} ne $bp_window{$data[8]})
				{
					if ($orientationb == 1)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>100 and $strand == $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
							#print "@src_return2\n";
							$left_bound_sr2 = $positionb1 + $src_return2[0];
							$right_bound_sr2 = $positionb4 - ($src_return2[3] - ($positionb2-$positionb1+101) + $cut_sr);
							$left_bound_sr2 = $data[5] unless ($src_return1[0]);
							$right_bound_sr2 = $data[8] unless ($src_return1[3]);
						}
					}
					else
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
							#print "@src_return2\n";
							$left_bound_sr2 = $positionb1 + $src_return2[0];
							$right_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
							$left_bound_sr2 = $data[5] unless ($src_return2[0]);
							$right_bound_sr2 = $data[8] unless ($src_return2[2]);
						}
					}
				}
				else
				{
					if ($bp_window{$data[5]} == 1 and $bp_window{$data[8]} == 1)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positionb2-$positionb1) and $mstart < ($positionb2-$positionb1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
							#print "@src_return2\n";
							$left_bound_sr2 = $positionb1 + $src_return2[0];
							$right_bound_sr2 = $positionb1 + $src_return2[2];
							$left_bound_sr2 = $data[5] unless ($src_return2[0]);
							$right_bound_sr2 = $data[8] unless ($src_return2[2]);
						}
					}
					if ($bp_window{$data[5]} == 2 and $bp_window{$data[8]} == 2)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start > ($positionb2-$positionb1+100) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
							#print "@src_return2\n";
							$left_bound_sr2 = $src_return2[0] - ($positionb2-$positionb1+101) + $positionb3;
							$right_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
							$left_bound_sr2 = $data[5] unless ($src_return2[0]);
							$right_bound_sr2 = $data[8] unless ($src_return2[2]);
						}
					}
				}
				
				my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
				my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
				if ($src_return1[4] >= $support_reads and $src_return2[4] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
				{
					my $reads0 = join("\t", @{$src_return1[5]{0}});
					my $reads1 = join("\t", @{$src_return2[5]{0}});
					print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads1\n";
					$result_sr[$sri][0] = $data[0];
					$result_sr[$sri][1] = $data[1];
					$result_sr[$sri][2] = $data[2];
					$result_sr[$sri][3] = "$src_return1[4]/$src_return2[4]";
					$result_sr[$sri][4] = $data[3];
					$result_sr[$sri][5] = $left_bound_sr1;
					$result_sr[$sri][6] = $left_bound_sr2;
					$result_sr[$sri][7] = $del_size;
					$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
					$sri++;
				}
				else
				{
					my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
					my ($mpd1, $mpd2) = split (/\//, $data[2]);
					if ($src_return1[4] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
					{
						my $left_bound_sr = $left_bound_sr1;
						my $right_bound_sr = $right_bound_sr1;
						my $size = $right_bound_sr - $left_bound_sr;
						my $reads = join("\t", @{$src_return1[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
					if ($src_return2[4] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
					{
						my $left_bound_sr = $left_bound_sr2;
						my $right_bound_sr = $right_bound_sr2;
						my $size = $right_bound_sr - $left_bound_sr;
						my $reads = join("\t", @{$src_return2[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
				}
			}
		}
		if ($data[0] =~ /insou/)
		{
			next unless ($cluster_region{$cl[0]} and $cluster_region{$cl[1]});
			my ($trash, $positiona1, $positiona2, $trash, $positiona3, $positiona4, $orientationa) = split (/__/, $cluster_region{$cl[0]});
			my ($trash, $positionb1, $positionb2, $trash, $positionb3, $positionb4, $orientationb) = split (/__/, $cluster_region{$cl[1]});
			if (&covered($data[9], $data[9], $positiona1, $positiona2))
			{
				$bp_window{$data[9]} = 1;
			}
			if (&covered($data[9], $data[9], $positiona3, $positiona4))
			{
				$bp_window{$data[9]} = 2;
			}
			if (&covered($data[4], $data[4], $positiona1, $positiona2))
			{
				$bp_window{$data[4]} = 1;
			}
			if (&covered($data[4], $data[4], $positiona3, $positiona4))
			{
				$bp_window{$data[4]} = 2;
			}
			if (&covered($data[8], $data[8], $positionb1, $positionb2))
			{
				$bp_window{$data[8]} = 1;
			}
			if (&covered($data[8], $data[8], $positionb3, $positionb4))
			{
				$bp_window{$data[8]} = 2;
			}
			if (&covered($data[5], $data[5], $positionb1, $positionb2))
			{
				$bp_window{$data[5]} = 1;
			}
			if (&covered($data[5], $data[5], $positionb3, $positionb4))
			{
				$bp_window{$data[5]} = 2;
			}
			#print "$cluster_region{$cl[0]}\n$cluster_region{$cl[1]}\n$data[4]\t$bp_window{$data[4]}\n$data[5]\t$bp_window{$data[5]}\n$data[8]\t$bp_window{$data[8]}\n$data[9]\t$bp_window{$data[9]}\n";
			if ($cluster_region{$cl[0]} eq $cluster_region{$cl[1]})
			{
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[4]} ne $bp_window{$data[9]} or $bp_window{$data[5]} ne $bp_window{$data[8]})
				{
					if ($orientationa == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "1 $start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
					}
					@discord_sr2 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand)
						{
							#print "2 $start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr2, $a;
						}
					}
					# small deletion, large insertion, close, 3 break points in one window
					if (@discord_sr1 >= $support_reads and @discord_sr2 >= $support_reads and $bp_window{$data[4]} == 2 and $bp_window{$data[9]} == 2)
					{
						my @src_return2 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
						my @src_return1 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
						#print "@src_return1\n@src_return2\n";
						my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
						$left_bound_sr1 = $positiona4 - ($src_return1[2] - ($positiona2-$positiona1+101)) if ($src_return1[2] >= ($positiona2-$positiona1+101));
						$right_bound_sr1 = $positiona4 - ($src_return1[0] - ($positiona2-$positiona1+101)) if ($src_return1[0] >= ($positiona2-$positiona1+101));
						$left_bound_sr2 = $positiona1 + $src_return2[0] if ($src_return2[0] < ($positiona2-$positiona1+101));
						$right_bound_sr2 = $positiona4 - ($src_return2[3] - ($positiona2-$positiona1+101) + $cut_sr) if ($src_return2[3] >= ($positiona2-$positiona1+101));
						$left_bound_sr1 = $data[9] unless ($src_return1[2]);
						$right_bound_sr1 = $data[4] unless ($src_return1[0]);
						$left_bound_sr2 = $data[8] unless ($src_return2[0]);
						$right_bound_sr2 = $data[5] unless ($src_return2[3]);
						my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
						my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
						if ($src_return1[4] >= $support_reads and $src_return2[4] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
						{
							my $reads0 = join("\t", @{$src_return1[5]{0}});
							my $reads1 = join("\t", @{$src_return2[5]{0}});
							print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
							$result_sr[$sri][0] = $data[0];
							$result_sr[$sri][1] = $data[1];
							$result_sr[$sri][2] = $data[2];
							$result_sr[$sri][3] = "$src_return1[4]/$src_return2[4]";
							$result_sr[$sri][4] = $data[3];
							$result_sr[$sri][5] = $right_bound_sr1;
							$result_sr[$sri][6] = $right_bound_sr2;
							$result_sr[$sri][7] = $del_size;
							$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
							$sri++;
						}
						else
						{
							my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
							my ($mpd1, $mpd2) = split (/\//, $data[2]);
							if ($src_return1[4] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
							{
								#my $left_bound_sr = $left_bound_sr1;
								#my $right_bound_sr = $right_bound_sr1;
								#my $size = $right_bound_sr - $left_bound_sr;
								#my $reads = join("\t", @{$src_return1[5]{0}});
								#print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								#$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								#$sri++;
							}
							if ($src_return2[4] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
							{
								my $left_bound_sr = $left_bound_sr2;
								my $right_bound_sr = $right_bound_sr2;
								my $size = $right_bound_sr - $left_bound_sr;
								my $reads = join("\t", @{$src_return2[5]{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
					}
					# small deletion, small insertion, far, 2 break points in one window, the other 2 in the other window
					elsif (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 2);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						if ($$ref_boundary2[2] > $$ref_boundary2[0])
						{
							my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
							$left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr if ($$ref_boundary1[3] < ($positiona2-$positiona1+101));
							$right_bound_sr1 = $positiona4 - ($$ref_boundary2[2] - ($positiona2-$positiona1+101)) if ($$ref_boundary2[2] >= ($positiona2-$positiona1+101));
							$left_bound_sr2 = $positiona1 + $$ref_boundary1[0] if ($$ref_boundary1[0] < ($positiona2-$positiona1+101));
							$right_bound_sr2 = $positiona4 - ($$ref_boundary2[1] - ($positiona2-$positiona1+101) + $cut_sr) if ($$ref_boundary2[1] >= ($positiona2-$positiona1+101));
							$left_bound_sr1 = $data[9] unless ($$ref_boundary1[3]);
							$right_bound_sr1 = $data[4] unless ($$ref_boundary2[2]);
							$left_bound_sr2 = $data[8] unless ($$ref_boundary1[0]);
							$right_bound_sr2 = $data[5] unless ($$ref_boundary2[1]);
							my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
							my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
							{
								my $reads0 = join("\t", @{$$ref_bpread{1}});
								my $reads1 = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $right_bound_sr1;
								$result_sr[$sri][6] = $right_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[0] >= $support_reads)
								{
									my $left_bound_sr = $left_bound_sr1;
									my $right_bound_sr = $right_bound_sr1;
									my $reads = join("\t", @{$$ref_bpread{0}});
									if ($$ref_bpread{1})
									{
										$reads = join("\t", @{$$ref_bpread{1}});
									}
									else
									{
										$left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
										$right_bound_sr = $positiona4 - ($$ref_boundary2[0] - ($positiona2-$positiona1+101));
									}
									if ($left_bound_sr and $right_bound_sr)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
								}
								if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
								{
									my $left_bound_sr = $left_bound_sr2;
									my $right_bound_sr = $right_bound_sr2;
									my $size = $right_bound_sr - $left_bound_sr;
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
						}
						else
						{
							my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
							$left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr if ($$ref_boundary1[1] < ($positiona2-$positiona1+101));
							$right_bound_sr1 = $positiona4 - ($$ref_boundary2[0] - ($positiona2-$positiona1+101)) if ($$ref_boundary2[0] >= ($positiona2-$positiona1+101));
							$left_bound_sr2 = $positiona1 + $$ref_boundary1[2] if ($$ref_boundary1[2] < ($positiona2-$positiona1+101));
							$right_bound_sr2 = $positiona4 - ($$ref_boundary2[3] - ($positiona2-$positiona1+101) + $cut_sr) if ($$ref_boundary2[3] >= ($positiona2-$positiona1+101));
							$left_bound_sr1 = $data[9] unless ($$ref_boundary1[1]);
							$right_bound_sr1 = $data[4] unless ($$ref_boundary2[0]);
							$left_bound_sr2 = $data[8] unless ($$ref_boundary1[2]);
							$right_bound_sr2 = $data[5] unless ($$ref_boundary2[3]);
							my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
							my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
							{
								my $reads0 = join("\t", @{$$ref_bpread{0}});
								my $reads1 = join("\t", @{$$ref_bpread{1}});
								print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $right_bound_sr1;
								$result_sr[$sri][6] = $right_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
								{
									my $left_bound_sr = $left_bound_sr1;
									my $right_bound_sr = $right_bound_sr1;
									my $size = $right_bound_sr - $left_bound_sr;
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
								if ($$ref_support_sr[0] >= $support_reads)
								{
									my $left_bound_sr = $left_bound_sr2;
									my $right_bound_sr = $right_bound_sr2;
									my $reads = join("\t", @{$$ref_bpread{0}});
									if ($$ref_bpread{1})
									{
										$reads = join("\t", @{$$ref_bpread{1}});
									}
									else
									{
										$left_bound_sr = $positiona1 + $$ref_boundary1[0];
										$right_bound_sr = $positiona4 - ($$ref_boundary2[1] - ($positiona2-$positiona1+101) + $cut_sr);
									}
									if ($left_bound_sr and $right_bound_sr)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
								}
							}
						}
					}
					# large deletion, small insertion, close, 3 break points in one window
					elsif (@discord_sr2 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 1);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						if ($$ref_boundary2[2] > $$ref_boundary2[0])
						{
							my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
							$left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr if ($$ref_boundary1[1] < ($positiona2-$positiona1+101));
							$right_bound_sr1 = $positiona1 + $$ref_boundary2[1] + $cut_sr if ($$ref_boundary2[1] < ($positiona2-$positiona1+101));
							$left_bound_sr2 = $positiona1 + $$ref_boundary1[2] if ($$ref_boundary1[2] < ($positiona2-$positiona1+101));
							$right_bound_sr2 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3 if ($$ref_boundary2[2] >= ($positiona2-$positiona1+101));
							$left_bound_sr1 = $data[9] unless ($$ref_boundary1[1]);
							$right_bound_sr1 = $data[4] unless ($$ref_boundary2[1]);
							$left_bound_sr2 = $data[8] unless ($$ref_boundary1[2]);
							$right_bound_sr2 = $data[5] unless ($$ref_boundary2[2]);
							my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
							my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
							{
								my $reads0 = join("\t", @{$$ref_bpread{1}});
								my $reads1 = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $right_bound_sr1;
								$result_sr[$sri][6] = $right_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[0] >= $support_reads)
								{
									#my $left_bound_sr = $left_bound_sr1;
									#my $right_bound_sr = $right_bound_sr1;
									#my $reads = join("\t", @{$$ref_bpread{0}});
									#if ($$ref_bpread{1})
									#{
									#	$reads = join("\t", @{$$ref_bpread{1}});
									#}
									#else
									#{
									#	$left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									#	$right_bound_sr = $positiona4 - ($$ref_boundary2[0] - ($positiona2-$positiona1+101));
									#}
									#if ($left_bound_sr and $right_bound_sr)
									#{
									#	my $size = $right_bound_sr - $left_bound_sr;
									#	print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									#	$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									#	$sri++;
									#}
								}
								if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
								{
									my $left_bound_sr = $left_bound_sr2;
									my $right_bound_sr = $right_bound_sr2;
									my $size = $right_bound_sr - $left_bound_sr;
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
						}
						else
						{
							my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
							$left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr if ($$ref_boundary1[3] < ($positiona2-$positiona1+101));
							$right_bound_sr1 = $positiona1 + $$ref_boundary2[3] + $cut_sr if ($$ref_boundary2[3] < ($positiona2-$positiona1+101));
							$left_bound_sr2 = $positiona1 + $$ref_boundary1[0] if ($$ref_boundary1[0] < ($positiona2-$positiona1+101));
							$right_bound_sr2 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3 if ($$ref_boundary2[0] >= ($positiona2-$positiona1+101));
							$left_bound_sr1 = $data[9] unless ($$ref_boundary1[3]);
							$right_bound_sr1 = $data[4] unless ($$ref_boundary2[3]);
							$left_bound_sr2 = $data[8] unless ($$ref_boundary1[0]);
							$right_bound_sr2 = $data[5] unless ($$ref_boundary2[0]);
							my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
							my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
							{
								my $reads0 = join("\t", @{$$ref_bpread{0}});
								my $reads1 = join("\t", @{$$ref_bpread{1}});
								print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $right_bound_sr1;
								$result_sr[$sri][6] = $right_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
								{
									#my $left_bound_sr = $left_bound_sr1;
									#my $right_bound_sr = $right_bound_sr1;
									#my $size = $right_bound_sr - $left_bound_sr;
									#my $reads = join("\t", @{$$ref_bpread{0}});
									#print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									#$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									#$sri++;
								}
								if ($$ref_support_sr[0] >= $support_reads)
								{
									my $left_bound_sr = $left_bound_sr2;
									my $right_bound_sr = $right_bound_sr2;
									my $reads = join("\t", @{$$ref_bpread{0}});
									if ($$ref_bpread{1})
									{
										$reads = join("\t", @{$$ref_bpread{1}});
									}
									else
									{
										$left_bound_sr = $positiona1 + $$ref_boundary1[0];
										$right_bound_sr = $positiona4 - ($$ref_boundary2[1] - ($positiona2-$positiona1+101) + $cut_sr);
									}
									if ($left_bound_sr and $right_bound_sr)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
								}
							}
						}
					}
				}
				# small deletion, small insertion, close, 4 break points in one window
				else
				{
					if ($bp_window{$data[4]} == 1 and $bp_window{$data[5]} == 1 and $bp_window{$data[8]} == 1 and $bp_window{$data[9]} == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart < ($positiona2-$positiona1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
							#print "@{$ref_boundary1}\n@{$ref_boundary2}\n@{$ref_support_sr}\n";
							# 2 cluster event
							if ($$ref_support_sr[1])
							{
								if ($$ref_boundary2[1]<$$ref_boundary2[3])
								{
									my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
									$left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									$right_bound_sr1 = $positiona1 + $$ref_boundary2[1] + $cut_sr;
									$left_bound_sr2 = $positiona1 + $$ref_boundary1[2];
									$right_bound_sr2 = $positiona1 + $$ref_boundary2[2];
									$left_bound_sr1 = $data[9] unless ($$ref_boundary1[1]);
									$right_bound_sr1 = $data[4] unless ($$ref_boundary2[1]);
									$left_bound_sr2 = $data[8] unless ($$ref_boundary1[2]);
									$right_bound_sr2 = $data[5] unless ($$ref_boundary2[2]);
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = $data[0];
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $reads = join("\t", @{$$ref_bpread{1}});
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
								else
								{
									my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
									$left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									$right_bound_sr1 = $positiona1 + $$ref_boundary2[3] + $cut_sr;
									$left_bound_sr2 = $positiona1 + $$ref_boundary1[0];
									$right_bound_sr2 = $positiona1 + $$ref_boundary2[0];
									$left_bound_sr1 = $data[9] unless ($$ref_boundary1[3]);
									$right_bound_sr1 = $data[4] unless ($$ref_boundary2[3]);
									$left_bound_sr2 = $data[8] unless ($$ref_boundary1[0]);
									$right_bound_sr2 = $data[5] unless ($$ref_boundary2[0]);
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = $data[0];
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $reads = join("\t", @{$$ref_bpread{1}});
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
							}
							elsif ($$ref_support_sr[0] >= $support_reads)
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								my $left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
								my $right_bound_sr = $positiona1 + $$ref_boundary2[1] + $cut_sr;
								my $reads = join("\t", @{$$ref_bpread{0}});
								if ($left_bound_sr and $right_bound_sr)
								{
									my $size = $right_bound_sr - $left_bound_sr;
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
						}
					}
					# 4 break points in window 2
					else
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start >= ($positiona2-$positiona1+100) and $mstart >= ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
							#print "@{$ref_boundary1}\n@{$ref_boundary2}\n@{$ref_support_sr}\n";
							# 2 cluster event
							if ($$ref_support_sr[1])
							{
								if ($$ref_boundary2[1]<$$ref_boundary2[3])
								{
									my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
									$left_bound_sr1 = $$ref_boundary1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									$right_bound_sr1 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									$left_bound_sr2 = $$ref_boundary1[2] - ($positiona2-$positiona1+101) + $positiona3;
									$right_bound_sr2 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
									$left_bound_sr1 = $data[9] unless ($$ref_boundary1[1]);
									$right_bound_sr1 = $data[4] unless ($$ref_boundary2[1]);
									$left_bound_sr2 = $data[8] unless ($$ref_boundary1[2]);
									$right_bound_sr2 = $data[5] unless ($$ref_boundary2[2]);
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = $data[0];
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $reads = join("\t", @{$$ref_bpread{1}});
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
								else
								{
									my ($left_bound_sr1, $right_bound_sr1, $left_bound_sr2, $right_bound_sr2);
									$left_bound_sr1 = $$ref_boundary1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									$right_bound_sr1 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									$left_bound_sr2 = $$ref_boundary1[0] - ($positiona2-$positiona1+101) + $positiona3;
									$right_bound_sr2 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
									$left_bound_sr1 = $data[9] unless ($$ref_boundary1[3]);
									$right_bound_sr1 = $data[4] unless ($$ref_boundary2[3]);
									$left_bound_sr2 = $data[8] unless ($$ref_boundary1[0]);
									$right_bound_sr2 = $data[5] unless ($$ref_boundary2[0]);
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = $data[0];
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $reads = join("\t", @{$$ref_bpread{1}});
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
							}
							elsif ($$ref_support_sr[0] >= $support_reads)
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								my $left_bound_sr = $$ref_boundary1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
								my $right_bound_sr = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
								my $reads = join("\t", @{$$ref_bpread{0}});
								if ($left_bound_sr and $right_bound_sr)
								{
									my $size = $right_bound_sr - $left_bound_sr;
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
									$sri++;
								}
							}
						}
					}
				}
			}
			else
			{
				my (@src_return1, $left_bound_sr1, $right_bound_sr1, @src_return2, $left_bound_sr2, $right_bound_sr2);
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[4]} ne $bp_window{$data[9]})
				{
					if ($orientationa == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr1 = $positiona1 + $src_return1[1];
							$right_bound_sr1 = $positiona4 - ($src_return1[2] - ($positiona2-$positiona1+101) + $cut_sr);
							$left_bound_sr1 = $data[9] unless ($src_return1[1]);
							$right_bound_sr1 = $data[4] unless ($src_return1[2]);
						}
					}
					else
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
							#print "@src_return1\n";
							$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
							$right_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$left_bound_sr1 = $data[9] unless ($src_return1[1]);
							$right_bound_sr1 = $data[4] unless ($src_return1[3]);
						}
					}
				}
				else
				{
					if ($bp_window{$data[4]} == 1 and $bp_window{$data[9]} == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart < ($positiona2-$positiona1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
							#print "@src_return1\n";
							$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
							$right_bound_sr1 = $positiona1 + $src_return1[3] + $cut_sr;
							$left_bound_sr1 = $data[9] unless ($src_return1[1]);
							$right_bound_sr1 = $data[4] unless ($src_return1[3]);
						}
					}
					if ($bp_window{$data[4]} == 2 and $bp_window{$data[9]} == 2)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start > ($positiona2-$positiona1+100) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
							#print "@src_return1\n";
							$left_bound_sr1 = $src_return1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$right_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$left_bound_sr1 = $data[9] unless ($src_return1[1]);
							$right_bound_sr1 = $data[4] unless ($src_return1[3]);
						}
					}
				}
				
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[1]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[5]} ne $bp_window{$data[8]})
				{
					if ($orientationb == 1)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>100 and $strand == $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr2 = $positionb1 + $src_return2[0] + $cut_sr;
							$right_bound_sr2 = $positionb4 - ($src_return2[3] - ($positionb2-$positionb1+101));
							$left_bound_sr2 = $data[8] unless ($src_return2[0]);
							$right_bound_sr2 = $data[5] unless ($src_return2[3]);
						}
					}
					else
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
							#print "@src_return2\n";
							$left_bound_sr2 = $positionb1 + $src_return2[0];
							$right_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
							$left_bound_sr2 = $data[8] unless ($src_return2[0]);
							$right_bound_sr2 = $data[5] unless ($src_return2[2]);
						}
					}
				}
				else
				{
					if ($bp_window{$data[5]} == 1 and $bp_window{$data[8]} == 1)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positionb2-$positionb1) and $mstart < ($positionb2-$positionb1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
							#print "@src_return2\n";
							$left_bound_sr2 = $positionb1 + $src_return2[0];
							$right_bound_sr2 = $positionb1 + $src_return2[2];
							$left_bound_sr2 = $data[8] unless ($src_return2[0]);
							$right_bound_sr2 = $data[5] unless ($src_return2[2]);
						}
					}
					if ($bp_window{$data[5]} == 2 and $bp_window{$data[8]} == 2)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start > ($positionb2-$positionb1+100) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
							#print "@src_return2\n";
							$left_bound_sr2 = $src_return2[0] - ($positionb2-$positionb1+101) + $positionb3;
							$right_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
							$left_bound_sr2 = $data[8] unless ($src_return2[0]);
							$right_bound_sr2 = $data[5] unless ($src_return2[2]);
						}
					}
				}
					
				my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
				my $ins_size = $left_bound_sr1 - $left_bound_sr2 + 1;
				if ($src_return1[4] >= $support_reads and $src_return2[4] >= $support_reads and $left_bound_sr1 and $right_bound_sr1 and $left_bound_sr2 and $right_bound_sr2)
				{
					my $reads0 = join("\t", @{$src_return1[5]{0}});
					my $reads1 = join("\t", @{$src_return2[5]{0}});
					print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
					$result_sr[$sri][0] = $data[0];
					$result_sr[$sri][1] = $data[1];
					$result_sr[$sri][2] = $data[2];
					$result_sr[$sri][3] = "$src_return1[4]/$src_return2[4]";
					$result_sr[$sri][4] = $data[3];
					$result_sr[$sri][5] = $right_bound_sr1;
					$result_sr[$sri][6] = $right_bound_sr2;
					$result_sr[$sri][7] = $del_size;
					$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$left_bound_sr1\t$ins_size\n";
					$sri++;
				}
				else
				{
					my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
					my ($mpd1, $mpd2) = split (/\//, $data[2]);
					if ($src_return1[4] >= $support_reads and $left_bound_sr1 and $right_bound_sr1)
					{
						my $left_bound_sr = $left_bound_sr1;
						my $right_bound_sr = $right_bound_sr1;
						my $size = $right_bound_sr - $left_bound_sr;
						my $reads = join("\t", @{$src_return1[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
					if ($src_return2[4] >= $support_reads and $left_bound_sr2 and $right_bound_sr2)
					{
						my $left_bound_sr = $left_bound_sr2;
						my $right_bound_sr = $right_bound_sr2;
						my $size = $right_bound_sr - $left_bound_sr;
						my $reads = join("\t", @{$src_return2[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
				}
			}
		}
		if ($data[0] eq 'invers' or $data[0] eq 'del_invers')
		{
			next unless ($cluster_region{$cl[0]} and $cluster_region{$cl[1]});
			my ($trash, $positiona1, $positiona2, $trash, $positiona3, $positiona4, $orientationa) = split (/__/, $cluster_region{$cl[0]});
			my ($trash, $positionb1, $positionb2, $trash, $positionb3, $positionb4, $orientationb) = split (/__/, $cluster_region{$cl[1]});
			if (&covered($data[4], $data[4], $positiona1, $positiona2))
			{
				$bp_window{$data[4]} = 1;
			}
			if (&covered($data[4], $data[4], $positiona3, $positiona4))
			{
				$bp_window{$data[4]} = 2;
			}
			if (&covered($data[8], $data[8], $positiona1, $positiona2))
			{
				$bp_window{$data[8]} = 1;
			}
			if (&covered($data[8], $data[8], $positiona3, $positiona4))
			{
				$bp_window{$data[8]} = 2;
			}
			if (&covered($data[5], $data[5], $positionb1, $positionb2))
			{
				$bp_window{$data[5]} = 1;
			}
			if (&covered($data[5], $data[5], $positionb3, $positionb4))
			{
				$bp_window{$data[5]} = 2;
			}
			if (&covered($data[9], $data[9], $positionb1, $positionb2))
			{
				$bp_window{$data[9]} = 1;
			}
			if (&covered($data[9], $data[9], $positionb3, $positionb4))
			{
				$bp_window{$data[9]} = 2;
			}
			#print "$cluster_region{$cl[0]}\n$cluster_region{$cl[1]}\n$data[4]\t$bp_window{$data[4]}\n$data[5]\t$bp_window{$data[5]}\n$data[8]\t$bp_window{$data[8]}\n$data[9]\t$bp_window{$data[9]}\n";
			if ($cluster_region{$cl[0]} eq $cluster_region{$cl[1]})
			{
				if ($bp_window{$data[4]} ne $bp_window{$data[8]} or $bp_window{$data[5]} ne $bp_window{$data[9]})
				{
					my @alignments;
					open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
					while ($newline1 = <SAM>)
					{
						chomp $newline1;
						push @alignments, $newline1;
					}
					close SAM;
					if ($orientationa == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 0);
							#print "@{$ref_boundary1}\t\n@{$ref_boundary2}\n@{$ref_support_sr}\n";
							if (abs($$ref_boundary1[2]-$$ref_boundary1[1])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[0]-$$ref_boundary2[3])<($$ref_is{'rlu'}-$cut_sr) or abs($$ref_boundary1[0]-$$ref_boundary1[3])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[2]-$$ref_boundary2[1])<($$ref_is{'rlu'}-$cut_sr))
							{
								my $left_bound_sr = $positiona1 + int(($$ref_boundary1[1] + $$ref_boundary1[2] + $cut_sr)/2);
								my $right_bound_sr = $positiona4 - (int(($$ref_boundary2[1] + $$ref_boundary2[2] + $cut_sr)/2) - ($positiona2-$positiona1+101));
								my $inv_size = $right_bound_sr - $left_bound_sr;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads0\t$reads1\n";
									$result_sr[$sri][0] = 'invers';
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$inv_size\n";
									$sri++;
								}
								else
								{
									my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
									my ($mpd1, $mpd2) = split (/\//, $data[2]);
									if ($$ref_support_sr[0] >= $support_reads)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										my $reads = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
									if ($$ref_support_sr[1] >= $support_reads)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										my $reads = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
								}
							}
							else
							{
								if ($$ref_boundary1[2]>$$ref_boundary1[0])
								{
									my $left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr1 = $positiona4 - ($$ref_boundary2[3] - ($positiona2-$positiona1+101));
									my $left_bound_sr2 = $positiona1 + $$ref_boundary1[2];
									my $right_bound_sr2 = $positiona4 - (($$ref_boundary2[0] + $cut_sr) - ($positiona2-$positiona1+101));
									my $del_size = $right_bound_sr1 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr2 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr1 - $right_bound_sr2;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{0}});
										$reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr2\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr1;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{1}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
								else
								{
									my $left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr1 = $positiona4 - ($$ref_boundary2[1] - ($positiona2-$positiona1+101));
									my $left_bound_sr2 = $positiona1 + $$ref_boundary1[0];
									my $right_bound_sr2 = $positiona4 - (($$ref_boundary2[2] + $cut_sr) - ($positiona2-$positiona1+101));
									my $del_size = $right_bound_sr1 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr2 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr1 - $right_bound_sr2;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{1}});
										$reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr2\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr2;
											my $reads = join("\t", @{$$ref_bpread{0}});
											if ($$ref_bpread{1})
											{
												$reads = join("\t", @{$$ref_bpread{1}});
											}
											else
											{
												$left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
											}
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr1;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
							}
						}
					}
					else
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
							#print "@{$ref_boundary1}\t\n@{$ref_boundary2}\n@{$ref_support_sr}\n";
							if (abs($$ref_boundary1[2]-$$ref_boundary1[1])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[0]-$$ref_boundary2[3])<($$ref_is{'rlu'}-$cut_sr) or abs($$ref_boundary1[0]-$$ref_boundary1[3])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[2]-$$ref_boundary2[1])<($$ref_is{'rlu'}-$cut_sr))
							{
								my $left_bound_sr = $positiona1 + int(($$ref_boundary1[1] + $$ref_boundary1[2] + $cut_sr)/2);
								my $right_bound_sr = int(($$ref_boundary2[1] + $$ref_boundary2[2] + $cut_sr)/2 - ($positiona2-$positiona1+101)) + $positiona3;
								my $inv_size = $right_bound_sr - $left_bound_sr;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads0\t$reads1\n";
									$result_sr[$sri][0] = 'invers';
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$inv_size\n";
									$sri++;
								}
								else
								{
									my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
									my ($mpd1, $mpd2) = split (/\//, $data[2]);
									if ($$ref_support_sr[0] >= $support_reads)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										my $reads = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
									if ($$ref_support_sr[1] >= $support_reads)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										my $reads = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
								}
							}
							else
							{
								if ($$ref_boundary1[2]>$$ref_boundary1[0])
								{
									my $left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr1 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3;
									my $left_bound_sr2 = $positiona1 + $$ref_boundary1[2];
									my $right_bound_sr2 = $$ref_boundary2[0] + $cut_sr - ($positiona2-$positiona1+101) + $positiona3;
									my $del_size = $right_bound_sr1 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr2 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr1 - $right_bound_sr2;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{0}});
										$reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr2\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{1}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
								else
								{
									my $left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr1 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3;
									my $left_bound_sr2 = $positiona1 + $$ref_boundary1[0];
									my $right_bound_sr2 = $$ref_boundary2[2] + $cut_sr - ($positiona2-$positiona1+101) + $positiona3;
									my $del_size = $right_bound_sr1 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr2 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr1 - $right_bound_sr2;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{1}});
										$reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr2\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $reads = join("\t", @{$$ref_bpread{0}});
											if ($$ref_bpread{1})
											{
												$reads = join("\t", @{$$ref_bpread{1}});
											}
											else
											{
												$left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
											}
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
							}
						}
					}
				}
				# 4 break points in one window
				else
				{
					my @alignments;
					open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
					while ($newline1 = <SAM>)
					{
						chomp $newline1;
						push @alignments, $newline1;
					}
					close SAM;
					if ($bp_window{$data[4]} == 1 and $bp_window{$data[5]} == 1 and $bp_window{$data[8]} == 1 and $bp_window{$data[9]} == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart < ($positiona2-$positiona1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
							#print "@{$ref_boundary1}\n@{$ref_boundary2}\n@{$ref_support_sr}\n";
							if (abs($$ref_boundary1[2]-$$ref_boundary1[1])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[2]-$$ref_boundary2[1])<($$ref_is{'rlu'}-$cut_sr) or abs($$ref_boundary1[0]-$$ref_boundary1[3])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[0]-$$ref_boundary2[3])<($$ref_is{'rlu'}-$cut_sr))
							{
								my $left_bound_sr = $positiona1 + $$ref_boundary1[2];
								my $right_bound_sr = $positiona1 + $$ref_boundary2[2];
								my $inv_size = $right_bound_sr - $left_bound_sr;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads0\t$reads1\n";
									$result_sr[$sri][0] = 'invers';
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$inv_size\n";
									$sri++;
								}
								else
								{
									my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
									my ($mpd1, $mpd2) = split (/\//, $data[2]);
									if ($$ref_support_sr[0] >= $support_reads)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										my $reads = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
									if ($$ref_support_sr[1] >= $support_reads)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										my $reads = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
								}
							}
							else
							{
								if ($$ref_boundary1[2]>$$ref_boundary1[0])
								{
									my $left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr1 = $positiona1 + $$ref_boundary2[0] + $cut_sr;
									my $left_bound_sr2 = $positiona1 + $$ref_boundary1[2];
									my $right_bound_sr2 = $positiona1 + $$ref_boundary2[3];
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{0}});
										$reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{1}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
								else
								{
									my $left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr1 = $positiona1 + $$ref_boundary2[2] + $cut_sr;
									my $left_bound_sr2 = $positiona1 + $$ref_boundary1[0];
									my $right_bound_sr2 = $positiona1 + $$ref_boundary2[1];
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{1}});
										$reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $reads = join("\t", @{$$ref_bpread{0}});
											if ($$ref_bpread{1})
											{
												$reads = join("\t", @{$$ref_bpread{1}});
											}
											else
											{
												$left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
												$right_bound_sr = $positiona1 + $$ref_boundary2[0] + $cut_sr;
											}
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
							}
						}
					}
					if ($bp_window{$data[4]} == 2 and $bp_window{$data[5]} == 2 and $bp_window{$data[8]} == 2 and $bp_window{$data[9]} == 2)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start > ($positiona2-$positiona1+100) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
							#print "@{$ref_boundary1}\n@{$ref_boundary2}\n@{$ref_support_sr}\n";
							if (abs($$ref_boundary1[2]-$$ref_boundary1[1])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[2]-$$ref_boundary2[1])<($$ref_is{'rlu'}-$cut_sr) or abs($$ref_boundary1[0]-$$ref_boundary1[3])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[0]-$$ref_boundary2[3])<($$ref_is{'rlu'}-$cut_sr))
							{
								my $left_bound_sr = $$ref_boundary1[2] - ($positiona2-$positiona1+101) + $positiona3;
								my $right_bound_sr = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
								my $inv_size = $right_bound_sr - $left_bound_sr;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads0\t$reads1\n";
									$result_sr[$sri][0] = 'invers';
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$inv_size\n";
									$sri++;
								}
								else
								{
									my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
									my ($mpd1, $mpd2) = split (/\//, $data[2]);
									if ($$ref_support_sr[0] >= $support_reads)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										my $reads = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
									if ($$ref_support_sr[1] >= $support_reads)
									{
										my $size = $right_bound_sr - $left_bound_sr;
										my $reads = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
										$sri++;
									}
								}
							}
							else
							{
								if ($$ref_boundary1[2]>$$ref_boundary1[0])
								{
									my $left_bound_sr1 = $$ref_boundary1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									my $right_bound_sr1 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[2] - ($positiona2-$positiona1+101) + $positiona3;
									my $right_bound_sr2 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3;
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{0}});
										$reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{1}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
								else
								{
									my $left_bound_sr1 = $$ref_boundary1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									my $right_bound_sr1 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[0] - ($positiona2-$positiona1+101) + $positiona3;
									my $right_bound_sr2 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3;
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{1}});
										$reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
									else
									{
										my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
										my ($mpd1, $mpd2) = split (/\//, $data[2]);
										if ($$ref_support_sr[0] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr1;
											my $right_bound_sr = $right_bound_sr1;
											my $reads = join("\t", @{$$ref_bpread{0}});
											if ($$ref_bpread{1})
											{
												$reads = join("\t", @{$$ref_bpread{1}});
											}
											else
											{
												$left_bound_sr = $$ref_boundary1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
												$right_bound_sr = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
											}
											my $size = $right_bound_sr - $left_bound_sr;
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
										if ($$ref_support_sr[1] >= $support_reads)
										{
											my $left_bound_sr = $left_bound_sr2;
											my $right_bound_sr = $right_bound_sr2;
											my $size = $right_bound_sr - $left_bound_sr;
											my $reads = join("\t", @{$$ref_bpread{0}});
											print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
											$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
											$sri++;
										}
									}
								}
							}
						}
					}
				}
			}
			else
			{
				my (@src_return1, $left_bound_sr1, $right_bound_sr1, @src_return2, $left_bound_sr2, $right_bound_sr2);
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[4]} ne $bp_window{$data[8]})
				{
					if ($orientationa == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
							$right_bound_sr1 = $positiona4 - ($src_return1[2] - ($positiona2-$positiona1+101));
						}
					}
					else
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
							$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
							$right_bound_sr1 = $src_return1[2] - ($positiona2-$positiona1+101) + $positiona3;
						}
					}
				}
				else
				{
					if ($bp_window{$data[4]} == 1 and $bp_window{$data[8]} == 1)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart < ($positiona2-$positiona1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
							#print "@src_return1\n";
							$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
							$right_bound_sr1 = $positiona1 + $src_return1[3] + $cut_sr;
						}
					}
					if ($bp_window{$data[4]} == 2 and $bp_window{$data[8]} == 2)
					{
						@discord_sr1 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start > ($positiona2-$positiona1+100) and $mstart > ($positiona2-$positiona1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr1, $a;
							}
						}
						if (@discord_sr1 >= $support_reads)
						{
							@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
							#print "@src_return1\n";
							$left_bound_sr1 = $src_return1[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$right_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
						}
					}
				}
				
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[1]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($bp_window{$data[5]} ne $bp_window{$data[9]})
				{
					if ($orientationa == 1)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>100 and $strand == $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
							$left_bound_sr2 = $positionb1 + $src_return2[0];
							$right_bound_sr2 = $positionb4 - ($src_return2[3] - ($positionb2-$positionb1+101) + $cut_sr);
						}
					}
					else
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
							$left_bound_sr2 = $positionb1 + $src_return2[0];
							$right_bound_sr2 = $src_return2[3] - ($positionb2-$positionb1+101) + $positionb3 + $cut_sr;
						}
					}
				}
				else
				{
					if ($bp_window{$data[5]} == 1 and $bp_window{$data[9]} == 1)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start < ($positionb2-$positionb1) and $mstart < ($positionb2-$positionb1))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
							#print "@src_return2\n";
							$left_bound_sr2 = $positionb1 + $src_return2[0];
							$right_bound_sr2 = $positionb1 + $src_return2[2];
						}
					}
					if ($bp_window{$data[5]} == 2 and $bp_window{$data[9]} == 2)
					{
						@discord_sr2 = undef;
						for my $a (@alignments)
						{
							my @data1 = split (/\t/, $a);
							my $start = $data1[3];
							my $strand = 1;
							$strand = -1 if ($data1[1] =~ /r/);
							my $mseqid = $data1[6];
							my $mstart = $data1[7];
							my $mstrand = 1;
							$mstrand = -1 if ($data1[1] =~ /R/);
							my $isize = $data1[8];
							if ($isize>0 and $strand != $mstrand and $start > ($positionb2-$positionb1+100) and $mstart > ($positionb2-$positionb1+100))
							{
								#print "$start $mstart\t$isize\t$strand $mstrand\n";
								push @discord_sr2, $a;
							}
						}
						if (@discord_sr2 >= $support_reads)
						{
							@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
							#print "@src_return2\n";
							$left_bound_sr2 = $src_return2[0] - ($positionb2-$positionb1+101) + $positionb3;
							$right_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
						}
					}
				}
				
				my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
				my $inv_size = $right_bound_sr1 - $left_bound_sr2;
				my $distance1 = $left_bound_sr2 - $left_bound_sr1;
				my $distance2 = $right_bound_sr2 - $right_bound_sr1;
				if ($src_return1[4] >= $support_reads and $src_return2[4] >= $support_reads)
				{
					my $reads0 = join("\t", @{$src_return1[5]{0}});
					my $reads1 = join("\t", @{$src_return2[5]{0}});
					print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
					$result_sr[$sri][0] = $data[0];
					$result_sr[$sri][1] = $data[1];
					$result_sr[$sri][2] = $data[2];
					$result_sr[$sri][3] = "$src_return1[4]/$src_return2[4]";
					$result_sr[$sri][4] = $data[3];
					$result_sr[$sri][5] = $left_bound_sr1;
					$result_sr[$sri][6] = $right_bound_sr2;
					$result_sr[$sri][7] = $del_size;
					$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
					$sri++;
				}
				else
				{
					my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
					my ($mpd1, $mpd2) = split (/\//, $data[2]);
					if ($src_return1[4] >= $support_reads)
					{
						my $left_bound_sr = $left_bound_sr1;
						my $right_bound_sr = $right_bound_sr1;
						my $size = $right_bound_sr - $left_bound_sr;
						my $reads = join("\t", @{$src_return1[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "invers_f\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
					if ($src_return2[4] >= $support_reads)
					{
						my $left_bound_sr = $left_bound_sr2;
						my $right_bound_sr = $right_bound_sr2;
						my $size = $right_bound_sr - $left_bound_sr;
						my $reads = join("\t", @{$src_return2[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "invers_r\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
				}
			}
		}
		if ($data[0] eq 'tandem_dup')
		{
			next unless ($cluster_region{$cl[0]});
			my ($trash, $position1, $position2, $trash, $position3, $position4, $orientation) = split (/__/, $cluster_region{$cl[0]});
			if (&covered($data[4], $data[4], $position1, $position2))
			{
				$bp_window{$data[4]} = 1;
			}
			if (&covered($data[4], $data[4], $position3, $position4))
			{
				$bp_window{$data[4]} = 2;
			}
			if (&covered($data[5], $data[5], $position1, $position2))
			{
				$bp_window{$data[5]} = 1;
			}
			if (&covered($data[5], $data[5], $position3, $position4))
			{
				$bp_window{$data[5]} = 2;
			}
			#print "$cluster_region{$cl[0]}\n$data[4]\t$bp_window{$data[4]}\n$data[5]\t$bp_window{$data[5]}\n";
			my @alignments;
			open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
			while ($newline1 = <SAM>)
			{
				chomp $newline1;
				push @alignments, $newline1;
			}
			close SAM;
			if ($bp_window{$data[4]} ne $bp_window{$data[5]})
			{
				@discord_sr1 = undef;
				for my $a (@alignments)
				{
					my @data1 = split (/\t/, $a);
					my $start = $data1[3];
					my $strand = 1;
					$strand = -1 if ($data1[1] =~ /r/);
					my $mseqid = $data1[6];
					my $mstart = $data1[7];
					my $mstrand = 1;
					$mstrand = -1 if ($data1[1] =~ /R/);
					my $isize = $data1[8];
					if ($isize>100 and $strand == $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
					{
						#print "$start $mstart\t$isize\t$strand $mstrand\n";
						push @discord_sr1, $a;
					}
				}
				if (@discord_sr1 >= $support_reads)
				{
					my @src_return = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
					my $left_bound_sr = $position1 + $src_return[0];
					my $right_bound_sr = $src_return[3] - ($position2-$position1+101) + $position3 + $cut_sr;
					$left_bound_sr = $data[4] unless ($src_return[0]);
					$right_bound_sr = $data[5] unless ($src_return[3]);
					my $size = $right_bound_sr - $left_bound_sr;
					if ($src_return[4] >= $support_reads)
					{
						my $reads = join("\t", @{$src_return[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
						$sri++;
					}
				}
			}
			else
			{
				if ($bp_window{$data[4]} == 1 and $bp_window{$data[5]} == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start < ($position2-$position1) and $mstart < ($position2-$position1))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 0);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						# 2 cluster event
						if ($$ref_support_sr[1] >= $support_reads)
						{
							# inssu
							if (&covered($$ref_boundary1[0], $$ref_boundary1[1], $$ref_boundary1[2], $$ref_boundary1[3]))
							{
								if ($$ref_boundary2[0]<$$ref_boundary2[2])
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[0];
									my $right_bound_sr1 = $position1 + $$ref_boundary2[1] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr2 = $position1 + $$ref_boundary2[2];
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssu';
										}
										else
										{
											$result_sr[$sri][0] = 'inssu';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[2];
									my $right_bound_sr1 = $position1 + $$ref_boundary2[3] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr2 = $position1 + $$ref_boundary2[0];
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssu';
										}
										else
										{
											$result_sr[$sri][0] = 'inssu';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
										$sri++;
									}
								}
							}
							# inssd
							if (&covered($$ref_boundary2[0], $$ref_boundary2[1], $$ref_boundary2[2], $$ref_boundary2[3]))
							{
								if ($$ref_boundary1[0]<$$ref_boundary1[2])
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[2];
									my $right_bound_sr1 = $position1 + $$ref_boundary2[3] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr2 = $position1 + $$ref_boundary2[0];
									my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssd';
										}
										else
										{
											$result_sr[$sri][0] = 'inssd';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr2;
										$result_sr[$sri][6] = $left_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[0];
									my $right_bound_sr1 = $position1 + $$ref_boundary2[1] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr2 = $position1 + $$ref_boundary2[2];
									my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssd';
										}
										else
										{
											$result_sr[$sri][0] = 'inssd';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr2;
										$result_sr[$sri][6] = $left_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
								}
							}
						}
						# single cluster
						else
						{
							my $left_bound_sr = $position1 + $$ref_boundary1[0];
							my $right_bound_sr = $position1 + $$ref_boundary2[1] + $cut_sr;
							$left_bound_sr = $data[4] unless ($$ref_boundary1[0]);
							$right_bound_sr = $data[5] unless ($$ref_boundary2[1]);
							my $size = $right_bound_sr - $left_bound_sr;
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
					}
				}
				if ($bp_window{$data[4]} == 2 and $bp_window{$data[5]} == 2)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>($$ref_is{rlu}-$cut_sr+5) and $strand == $mstrand and $start > ($position2-$position1+100) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 0);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						# 2 cluster event
						if ($$ref_support_sr[1] >= $support_reads)
						{
							# inssu
							if (&covered($$ref_boundary1[0], $$ref_boundary1[1], $$ref_boundary1[2], $$ref_boundary1[3]))
							{
								if ($$ref_boundary2[0]<$$ref_boundary2[2])
								{
									my $left_bound_sr1 = $$ref_boundary1[0] - ($position2-$position1+101) + $position3;
									my $right_bound_sr1 = $$ref_boundary2[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr2 = $$ref_boundary2[2] - ($position2-$position1+101) + $position3;
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssu';
										}
										else
										{
											$result_sr[$sri][0] = 'inssu';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $$ref_boundary1[2] - ($position2-$position1+101) + $position3;
									my $right_bound_sr1 = $$ref_boundary2[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr2 = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
									my $del_size = $right_bound_sr2 - $right_bound_sr1 - 1;
									my $ins_size = $left_bound_sr2 - $left_bound_sr1 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$right_bound_sr1\__$data[3]\__$left_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssu';
										}
										else
										{
											$result_sr[$sri][0] = 'inssu';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $right_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr1\t$left_bound_sr2\t$ins_size\n";
										$sri++;
									}
								}
							}
							# inssd
							if (&covered($$ref_boundary2[0], $$ref_boundary2[1], $$ref_boundary2[2], $$ref_boundary2[3]))
							{
								if ($$ref_boundary1[0]<$$ref_boundary1[2])
								{
									my $left_bound_sr1 = $$ref_boundary1[2] - ($position2-$position1+101) + $position3;
									my $right_bound_sr1 = $$ref_boundary2[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr2 = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
									my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{0}});
										my $reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssd';
										}
										else
										{
											$result_sr[$sri][0] = 'inssd';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr2;
										$result_sr[$sri][6] = $left_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $$ref_boundary1[0] - ($position2-$position1+101) + $position3;
									my $right_bound_sr1 = $$ref_boundary2[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr2 = $$ref_boundary2[2] - ($position2-$position1+101) + $position3;
									my $del_size = $left_bound_sr1 - $left_bound_sr2 - 1;
									my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my $reads0 = join("\t", @{$$ref_bpread{1}});
										my $reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr2\__$data[3]\__$right_bound_sr2\n$reads0\n$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads1\n";
										if ($del_size > 10)
										{
											$result_sr[$sri][0] = 'del_inssd';
										}
										else
										{
											$result_sr[$sri][0] = 'inssd';
										}
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr2;
										$result_sr[$sri][6] = $left_bound_sr1;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
										$sri++;
									}
								}
							}
						}
						# single cluster
						else
						{
							my $left_bound_sr = $$ref_boundary1[0] - ($position2-$position1+101) + $position3;
							my $right_bound_sr = $$ref_boundary2[1] - ($position2-$position1+101) + $position3 + $cut_sr;
							$left_bound_sr = $data[4] unless ($$ref_boundary1[0]);
							$right_bound_sr = $data[5] unless ($$ref_boundary2[1]);
							my $size = $right_bound_sr - $left_bound_sr;
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
					}
				}
			}
		}
		if ($data[0] eq 'invers_f')
		{
			next unless ($cluster_region{$cl[0]});
			my ($trash, $position1, $position2, $trash, $position3, $position4, $orientation) = split (/__/, $cluster_region{$cl[0]});
			if (&covered($data[4], $data[4], $position1, $position2))
			{
				$bp_window{$data[4]} = 1;
			}
			if (&covered($data[4], $data[4], $position3, $position4))
			{
				$bp_window{$data[4]} = 2;
			}
			if (&covered($data[5], $data[5], $position1, $position2))
			{
				$bp_window{$data[5]} = 1;
			}
			if (&covered($data[5], $data[5], $position3, $position4))
			{
				$bp_window{$data[5]} = 2;
			}
			#print "$cluster_region{$cl[0]}\n$data[4]\t$bp_window{$data[4]}\n$data[5]\t$bp_window{$data[5]}\n";
			my @alignments;
			open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
			while ($newline1 = <SAM>)
			{
				chomp $newline1;
				push @alignments, $newline1;
			}
			close SAM;
			if ($bp_window{$data[4]} ne $bp_window{$data[5]})
			{
				if ($orientation == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my @src_return = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
						my $left_bound_sr = $position1 + $src_return[1] + $cut_sr;
						my $right_bound_sr = $position4 - ($src_return[2] - ($position2-$position1+101));
						$left_bound_sr = $data[4] unless ($src_return[1]);
						$right_bound_sr = $data[5] unless ($src_return[2]);
						my $size = $right_bound_sr - $left_bound_sr;
						if ($src_return[4] >= $support_reads)
						{
							my $reads = join("\t", @{$src_return[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
							$sri++;
						}
					}
				}
				else
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my @src_return = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
						my $left_bound_sr = $position1 + $src_return[1] + $cut_sr;
						my $right_bound_sr = $src_return[2] - ($position2-$position1+101) + $position3;
						$left_bound_sr = $data[4] unless ($src_return[1]);
						$right_bound_sr = $data[5] unless ($src_return[2]);
						my $size = $right_bound_sr - $left_bound_sr;
						if ($src_return[4] >= $support_reads)
						{
							my $reads = join("\t", @{$src_return[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
							$sri++;
						}
					}
				}
			}
			else
			{
				if ($bp_window{$data[4]} == 1 and $bp_window{$data[5]} == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start < ($position2-$position1) and $mstart < ($position2-$position1))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						# 2 cluster event
						if ($$ref_support_sr[1] >= $support_reads)
						{
							if (abs($$ref_boundary1[2]-$$ref_boundary1[1])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[2]-$$ref_boundary2[1])<($$ref_is{'rlu'}-$cut_sr) or abs($$ref_boundary1[0]-$$ref_boundary1[3])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[0]-$$ref_boundary2[3])<($$ref_is{'rlu'}-$cut_sr))
							{
								my $left_bound_sr = $position1 + $$ref_boundary1[2];
								my $right_bound_sr = $position1 + $$ref_boundary2[2];
								my $inv_size = $right_bound_sr - $left_bound_sr;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads0\t$reads1\n";
									$result_sr[$sri][0] = 'invers';
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$inv_size\n";
									$sri++;
								}
							}
							else
							{
								if ($$ref_boundary1[2]>$$ref_boundary1[0])
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr1 = $position1 + $$ref_boundary2[0] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[2];
									my $right_bound_sr2 = $position1 + $$ref_boundary2[3];
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{0}});
										$reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr1 = $position1 + $$ref_boundary2[2] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[0];
									my $right_bound_sr2 = $position1 + $$ref_boundary2[1];
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{1}});
										$reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
								}
							}
						}
						# single cluster
						else
						{
							my $left_bound_sr = $position1 + $$ref_boundary1[1] + $cut_sr;
							my $right_bound_sr = $position1 + $$ref_boundary2[1] + $cut_sr;
							$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
							$right_bound_sr = $data[5] unless ($$ref_boundary2[1]);
							my $size = $right_bound_sr - $left_bound_sr;
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
					}
				}
				if ($bp_window{$data[4]} == 2 and $bp_window{$data[5]} == 2)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start > ($position2-$position1+100) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						# 2 cluster event
						if ($$ref_support_sr[1] >= $support_reads)
						{
							if (abs($$ref_boundary1[2]-$$ref_boundary1[1])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[2]-$$ref_boundary2[1])<($$ref_is{'rlu'}-$cut_sr) or abs($$ref_boundary1[0]-$$ref_boundary1[3])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[0]-$$ref_boundary2[3])<($$ref_is{'rlu'}-$cut_sr))
							{
								my $left_bound_sr = $$ref_boundary1[2] - ($position2-$position1+101) + $position3;
								my $right_bound_sr = $$ref_boundary2[2] - ($position2-$position1+101) + $position3;
								my $inv_size = $right_bound_sr - $left_bound_sr;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads0\t$reads1\n";
									$result_sr[$sri][0] = 'invers';
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$inv_size\n";
									$sri++;
								}
							}
							else
							{
								if ($$ref_boundary1[2]>$$ref_boundary1[0])
								{
									my $left_bound_sr1 = $$ref_boundary1[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr1 = $$ref_boundary2[0] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[2] - ($position2-$position1+101) + $position3;
									my $right_bound_sr2 = $$ref_boundary2[3] - ($position2-$position1+101) + $position3;
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{0}});
										$reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $$ref_boundary1[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr1 = $$ref_boundary2[2] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[0] - ($position2-$position1+101) + $position3;
									my $right_bound_sr2 = $$ref_boundary2[1] - ($position2-$position1+101) + $position3;
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{1}});
										$reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
								}
							}
						}
						# single cluster
						else
						{
							my $left_bound_sr = $$ref_boundary1[1] - ($position2-$position1+101) + $position3 + $cut_sr;
							my $right_bound_sr = $$ref_boundary2[1] - ($position2-$position1+101) + $position3 + $cut_sr;
							$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
							$right_bound_sr = $data[5] unless ($$ref_boundary2[1]);
							my $size = $right_bound_sr - $left_bound_sr;
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
					}
				}
			}
		}
		if ($data[0] eq 'invers_r')
		{
			next unless ($cluster_region{$cl[0]});
			my ($trash, $position1, $position2, $trash, $position3, $position4, $orientation) = split (/__/, $cluster_region{$cl[0]});
			if (&covered($data[4], $data[4], $position1, $position2))
			{
				$bp_window{$data[4]} = 1;
			}
			if (&covered($data[4], $data[4], $position3, $position4))
			{
				$bp_window{$data[4]} = 2;
			}
			if (&covered($data[5], $data[5], $position1, $position2))
			{
				$bp_window{$data[5]} = 1;
			}
			if (&covered($data[5], $data[5], $position3, $position4))
			{
				$bp_window{$data[5]} = 2;
			}
			#print "$cluster_region{$cl[0]}\n$data[4]\t$bp_window{$data[4]}\n$data[5]\t$bp_window{$data[5]}\n";
			my @alignments;
			open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
			while ($newline1 = <SAM>)
			{
				chomp $newline1;
				push @alignments, $newline1;
			}
			close SAM;
			if ($bp_window{$data[4]} ne $bp_window{$data[5]})
			{
				if ($orientation == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my @src_return = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
						my $left_bound_sr = $position1 + $src_return[0];
						my $right_bound_sr = $position4 - ($src_return[3] - ($position2-$position1+101) + $cut_sr);
						$left_bound_sr = $data[4] unless ($src_return[0]);
						$right_bound_sr = $data[5] unless ($src_return[3]);
						my $size = $right_bound_sr - $left_bound_sr;
						if ($src_return[4] >= $support_reads)
						{
							my $reads = join("\t", @{$src_return[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
							$sri++;
						}
					}
				}
				else
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my @src_return = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
						my $left_bound_sr = $position1 + $src_return[0];
						my $right_bound_sr = $src_return[3] - ($position2-$position1+101) + $position3 + $cut_sr;
						$left_bound_sr = $data[4] unless ($src_return[0]);
						$right_bound_sr = $data[5] unless ($src_return[3]);
						my $size = $right_bound_sr - $left_bound_sr;
						if ($src_return[4] >= $support_reads)
						{
							my $reads = join("\t", @{$src_return[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
							$sri++;
						}
					}
				}
				
			}
			else
			{
				if ($bp_window{$data[4]} == 1 and $bp_window{$data[5]} == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start < ($position2-$position1) and $mstart < ($position2-$position1))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						# 2 cluster event
						if ($$ref_support_sr[1] >= $support_reads)
						{
							if (abs($$ref_boundary1[2]-$$ref_boundary1[1])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[2]-$$ref_boundary2[1])<($$ref_is{'rlu'}-$cut_sr) or abs($$ref_boundary1[0]-$$ref_boundary1[3])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[0]-$$ref_boundary2[3])<($$ref_is{'rlu'}-$cut_sr))
							{
								my $left_bound_sr = $position1 + $$ref_boundary1[2];
								my $right_bound_sr = $position1 + $$ref_boundary2[2];
								my $inv_size = $right_bound_sr - $left_bound_sr;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads0\t$reads1\n";
									$result_sr[$sri][0] = 'invers';
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$inv_size\n";
									$sri++;
								}
							}
							else
							{
								if ($$ref_boundary1[2]>$$ref_boundary1[0])
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr1 = $position1 + $$ref_boundary2[0] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[2];
									my $right_bound_sr2 = $position1 + $$ref_boundary2[3];
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{0}});
										$reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $position1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr1 = $position1 + $$ref_boundary2[2] + $cut_sr;
									my $left_bound_sr2 = $position1 + $$ref_boundary1[0];
									my $right_bound_sr2 = $position1 + $$ref_boundary2[1];
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{1}});
										$reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
								}
							}
						}
						# single cluster
						else
						{
							my $left_bound_sr = $position1 + $$ref_boundary1[0];
							my $right_bound_sr = $position1 + $$ref_boundary2[0];
							$left_bound_sr = $data[4] unless ($$ref_boundary1[0]);
							$right_bound_sr = $data[5] unless ($$ref_boundary2[0]);
							my $size = $right_bound_sr - $left_bound_sr;
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
					}
				}
				if ($bp_window{$data[4]} == 2 and $bp_window{$data[5]} == 2)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start > ($position2-$position1+100) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						# 2 cluster event
						if ($$ref_support_sr[1] >= $support_reads)
						{
							if (abs($$ref_boundary1[2]-$$ref_boundary1[1])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[2]-$$ref_boundary2[1])<($$ref_is{'rlu'}-$cut_sr) or abs($$ref_boundary1[0]-$$ref_boundary1[3])<($$ref_is{'rlu'}-$cut_sr) and abs($$ref_boundary2[0]-$$ref_boundary2[3])<($$ref_is{'rlu'}-$cut_sr))
							{
								my $left_bound_sr = $$ref_boundary1[2] - ($position2-$position1+101) + $position3;
								my $right_bound_sr = $$ref_boundary2[2] - ($position2-$position1+101) + $position3;
								my $inv_size = $right_bound_sr - $left_bound_sr;
								if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
								{
									my $reads0 = join("\t", @{$$ref_bpread{0}});
									my $reads1 = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads0\t$reads1\n";
									$result_sr[$sri][0] = 'invers';
									$result_sr[$sri][1] = $data[1];
									$result_sr[$sri][2] = $data[2];
									$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
									$result_sr[$sri][4] = $data[3];
									$result_sr[$sri][5] = $left_bound_sr;
									$result_sr[$sri][6] = $right_bound_sr;
									$result_sr[$sri][7] = "$inv_size\n";
									$sri++;
								}
							}
							else
							{
								if ($$ref_boundary1[2]>$$ref_boundary1[0])
								{
									my $left_bound_sr1 = $$ref_boundary1[1] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr1 = $$ref_boundary2[0] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[2] - ($position2-$position1+101) + $position3;
									my $right_bound_sr2 = $$ref_boundary2[3] - ($position2-$position1+101) + $position3;
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{0}});
										$reads1 = join("\t", @{$$ref_bpread{1}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
								}
								else
								{
									my $left_bound_sr1 = $$ref_boundary1[3] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $right_bound_sr1 = $$ref_boundary2[2] - ($position2-$position1+101) + $position3 + $cut_sr;
									my $left_bound_sr2 = $$ref_boundary1[0] - ($position2-$position1+101) + $position3;
									my $right_bound_sr2 = $$ref_boundary2[1] - ($position2-$position1+101) + $position3;
									my $del_size = $right_bound_sr2 - $left_bound_sr1 - 1;
									my $inv_size = $right_bound_sr1 - $left_bound_sr2;
									my $distance1 = $left_bound_sr2 - $left_bound_sr1;
									my $distance2 = $right_bound_sr2 - $right_bound_sr1;
									if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
									{
										my ($reads0, $reads1);
										$reads0 = join("\t", @{$$ref_bpread{1}});
										$reads1 = join("\t", @{$$ref_bpread{0}});
										print BPREAD "$data[3]\__$left_bound_sr1\__$data[3]\__$right_bound_sr1\n$reads0\n$data[3]\__$right_bound_sr2\__$data[3]\__$left_bound_sr2\n$reads1\n";
										$result_sr[$sri][0] = 'del_invers';
										$result_sr[$sri][1] = $data[1];
										$result_sr[$sri][2] = $data[2];
										$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
										$result_sr[$sri][4] = $data[3];
										$result_sr[$sri][5] = $left_bound_sr1;
										$result_sr[$sri][6] = $right_bound_sr2;
										$result_sr[$sri][7] = $del_size;
										$result_sr[$sri][8] = "$data[3]\t$left_bound_sr2\t$right_bound_sr1\t$inv_size\t$distance1\t$distance2\n";
										$sri++;
									}
								}
							}
						}
						# single cluster
						else
						{
							my $left_bound_sr = $$ref_boundary1[0] - ($position2-$position1+101) + $position3;
							my $right_bound_sr = $$ref_boundary2[0] - ($position2-$position1+101) + $position3;
							$left_bound_sr = $data[4] unless ($$ref_boundary1[0]);
							$right_bound_sr = $data[5] unless ($$ref_boundary2[0]);
							my $size = $right_bound_sr - $left_bound_sr;
							if ($$ref_support_sr[0] >= $support_reads)
							{
								my $reads = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr\__$data[3]\__$right_bound_sr\n$reads\n";
								$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t$right_bound_sr\t$size\n";
								$sri++;
							}
						}
					}
				}
			}
		}
		
	}
	close MPINTRASRD;
	
	$i = 0;
	open MPINTERSRD, "<$mpinter_outfile";
	while ($newline = <MPINTERSRD>)
	{
		chomp $newline;
		$i++;
		next if (defined($line) and $i < $line);
		last if (defined($line) and $i > $line);
		my @data = split (/\t/, $newline);
		my @cl = split (/\//, $data[1]);
		my (@discord_sr1, @discord_sr2);
		my $ref_discord_sr1 = \@discord_sr1;
		my $ref_discord_sr2 = \@discord_sr2;
		
		if ($data[0] =~ /inss/)
		{
			next unless ($cluster_region{$cl[0]} and $cluster_region{$cl[1]});
			my ($chr_bp_a1, $positiona1, $positiona2, $chr_bp_a2, $positiona3, $positiona4, $orientationa) = split (/__/, $cluster_region{$cl[0]});
			my ($chr_bp_b1, $positionb1, $positionb2, $chr_bp_b2, $positionb3, $positionb4, $orientationb) = split (/__/, $cluster_region{$cl[1]});
			if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
			{
				if ($cluster_region{$cl[0]} eq $cluster_region{$cl[1]})
				{
					print STDERR "inss $i\n";
				}
				#($positiona1, $positiona2, $positiona3, $positiona4) = ($positiona3, $positiona4, $positiona1, $positiona2);
			}
			#if ($chr_bp_b1 eq $data[7] and $chr_bp_b2 eq $data[3]){#($positionb1, $positionb2, $positionb3, $positionb4) = ($positionb3, $positionb4, $positionb1, $positionb2);}
			#print "$cluster_region{$cl[0]}\n$cluster_region{$cl[1]}\n$positiona1, $positiona2, $positiona3, $positiona4\n$positionb1, $positionb2, $positionb3, $positionb4\n";
			
			if ($cluster_region{$cl[0]} eq $cluster_region{$cl[1]})
			{
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				@discord_sr1 = undef;
				for my $a (@alignments)
				{
					my @data1 = split (/\t/, $a);
					my $start = $data1[3];
					my $strand = 1;
					$strand = -1 if ($data1[1] =~ /r/);
					my $mseqid = $data1[6];
					my $mstart = $data1[7];
					my $mstrand = 1;
					$mstrand = -1 if ($data1[1] =~ /R/);
					my $isize = $data1[8];
					if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
					{
						#print "$start $mstart\t$isize\t$strand $mstrand\n";
						push @discord_sr1, $a;
					}
				}
				if (@discord_sr1 >= $support_reads)
				{
					my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 1);
					#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
					if ($$ref_boundary1[2] > $$ref_boundary1[0])
					{
						my $left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
						my $right_bound_sr1 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
						my $left_bound_sr2 = $positiona1 + $$ref_boundary1[2];
						my $right_bound_sr2 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
						if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
						{
							$left_bound_sr1 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$right_bound_sr1 = $positiona1 + $$ref_boundary1[0];
							$left_bound_sr2 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
							$right_bound_sr2 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
						}
						if ($data[3] lt $data[7])
						{
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[1]);
							$right_bound_sr1 = $data[8] unless ($$ref_boundary2[0]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[2]);
							$right_bound_sr2 = $data[9] unless ($$ref_boundary2[3]);
						}
						else
						{
							$left_bound_sr1 = $data[5] unless ($$ref_boundary1[1]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[0]);
							$left_bound_sr2 = $data[4] unless ($$ref_boundary1[2]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[3]);
						}
						if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
						{
							my $reads0 = join("\t", @{$$ref_bpread{0}});
							my $reads1 = join("\t", @{$$ref_bpread{1}});
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr2 - $right_bound_sr1 + 1;
							print BPREAD "$data[3]\__$left_bound_sr1\__$data[7]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[7]\__$right_bound_sr2\n$reads1\n";
							$result_sr[$sri][0] = $data[0];
							$result_sr[$sri][1] = $data[1];
							$result_sr[$sri][2] = $data[2];
							$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
							$result_sr[$sri][4] = $data[3];
							$result_sr[$sri][5] = $left_bound_sr1;
							$result_sr[$sri][6] = $left_bound_sr2;
							$result_sr[$sri][7] = $del_size;
							$result_sr[$sri][8] = "$data[7]\t$right_bound_sr1\t$right_bound_sr2\t$ins_size\n";
							$sri++;
						}
						else
						{
							my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
							my ($mpd1, $mpd2) = split (/\//, $data[2]);
							if ($$ref_support_sr[0] >= $support_reads)
							{
								if ($data[3] lt $data[7])
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
										$right_bound_sr = $positiona1 + $$ref_boundary1[0];
									}
									$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
									$right_bound_sr = $data[8] unless ($$ref_boundary2[0]);
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t1\t$data[7]\t$right_bound_sr\t-1\n";
									$sri++;
								}
								else
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[1];
									my $right_bound_sr = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3;
										$right_bound_sr = $positiona1 + $$ref_boundary1[0] + $cut_sr;
									}
									$left_bound_sr = $data[5] unless ($$ref_boundary1[1]);
									$right_bound_sr = $data[9] unless ($$ref_boundary2[0]);
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[7]\t$right_bound_sr\t1\t$data[3]\t$left_bound_sr\t-1\n";
									$sri++;
								}
							}
							if ($$ref_support_sr[1] >= $support_reads)
							{
								if ($data[3] lt $data[7])
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[2];
									my $right_bound_sr = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
										$right_bound_sr = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									}
									$left_bound_sr = $data[5] unless ($$ref_boundary1[2]);
									$right_bound_sr = $data[9] unless ($$ref_boundary2[3]);
									my $reads = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t-1\t$data[7]\t$right_bound_sr\t1\n";
									$sri++;
								}
								else
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[2] + $cut_sr;
									my $right_bound_sr = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
										$right_bound_sr = $positiona1 + $$ref_boundary1[3];
									}
									$left_bound_sr = $data[4] unless ($$ref_boundary1[2]);
									$right_bound_sr = $data[8] unless ($$ref_boundary2[3]);
									my $reads = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[7]\t$right_bound_sr\t-1\t$data[3]\t$left_bound_sr\t1\n";
									$sri++;
								}
							}
						}
					}
					else
					{
						my $left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
						my $right_bound_sr1 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
						my $left_bound_sr2 = $positiona1 + $$ref_boundary1[0];
						my $right_bound_sr2 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
						if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
						{
							$left_bound_sr1 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$right_bound_sr1 = $positiona1 + $$ref_boundary1[2];
							$left_bound_sr2 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
							$right_bound_sr2 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
						}
						if ($data[3] lt $data[7])
						{
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[3]);
							$right_bound_sr1 = $data[8] unless ($$ref_boundary2[2]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[0]);
							$right_bound_sr2 = $data[9] unless ($$ref_boundary2[1]);
						}
						else
						{
							$left_bound_sr1 = $data[5] unless ($$ref_boundary1[3]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[2]);
							$left_bound_sr2 = $data[4] unless ($$ref_boundary1[0]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[1]);
						}
						if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
						{
							my $reads0 = join("\t", @{$$ref_bpread{1}});
							my $reads1 = join("\t", @{$$ref_bpread{0}});
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr2 - $right_bound_sr1 + 1;
							print BPREAD "$data[3]\__$left_bound_sr1\__$data[7]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[7]\__$right_bound_sr2\n$reads1\n";
							$result_sr[$sri][0] = $data[0];
							$result_sr[$sri][1] = $data[1];
							$result_sr[$sri][2] = $data[2];
							$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
							$result_sr[$sri][4] = $data[3];
							$result_sr[$sri][5] = $left_bound_sr1;
							$result_sr[$sri][6] = $left_bound_sr2;
							$result_sr[$sri][7] = $del_size;
							$result_sr[$sri][8] = "$data[7]\t$right_bound_sr1\t$right_bound_sr2\t$ins_size\n";
							$sri++;
						}
						else
						{
							my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
							my ($mpd1, $mpd2) = split (/\//, $data[2]);
							if ($$ref_support_sr[1] >= $support_reads)
							{
								if ($data[3] lt $data[7])
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
										$right_bound_sr = $positiona1 + $$ref_boundary1[2];
									}
									$left_bound_sr = $data[4] unless ($$ref_boundary1[3]);
									$right_bound_sr = $data[8] unless ($$ref_boundary2[2]);
									my $reads = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t1\t$data[7]\t$right_bound_sr\t-1\n";
									$sri++;
								}
								else
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[3];
									my $right_bound_sr = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3;
										$right_bound_sr = $positiona1 + $$ref_boundary1[2] + $cut_sr;
									}
									$left_bound_sr = $data[5] unless ($$ref_boundary1[3]);
									$right_bound_sr = $data[9] unless ($$ref_boundary2[2]);
									my $reads = join("\t", @{$$ref_bpread{1}});
									print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[7]\t$right_bound_sr\t1\t$data[3]\t$left_bound_sr\t-1\n";
									$sri++;
								}
							}
							if ($$ref_support_sr[0] >= $support_reads)
							{
								if ($data[3] lt $data[7])
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[0];
									my $right_bound_sr = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
										$right_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									}
									$left_bound_sr = $data[5] unless ($$ref_boundary1[0]);
									$right_bound_sr = $data[9] unless ($$ref_boundary2[1]);
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t-1\t$data[7]\t$right_bound_sr\t1\n";
									$sri++;
								}
								else
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[0] + $cut_sr;
									my $right_bound_sr = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
										$right_bound_sr = $positiona1 + $$ref_boundary1[1];
									}
									$left_bound_sr = $data[4] unless ($$ref_boundary1[0]);
									$right_bound_sr = $data[8] unless ($$ref_boundary2[1]);
									my $reads = join("\t", @{$$ref_bpread{0}});
									print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
									$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[7]\t$right_bound_sr\t-1\t$data[3]\t$left_bound_sr\t1\n";
									$sri++;
								}
							}
						}
					}
				}
			}
			else
			{
				my (@src_return1, $left_bound_sr1, $right_bound_sr1, @src_return2, $left_bound_sr2, $right_bound_sr2);
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				@discord_sr1 = undef;
				for my $a (@alignments)
				{
					my @data1 = split (/\t/, $a);
					my $start = $data1[3];
					my $strand = 1;
					$strand = -1 if ($data1[1] =~ /r/);
					my $mseqid = $data1[6];
					my $mstart = $data1[7];
					my $mstrand = 1;
					$mstrand = -1 if ($data1[1] =~ /R/);
					my $isize = $data1[8];
					if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
					{
						#print "1 $start $mstart\t$isize\t$strand $mstrand\n";
						push @discord_sr1, $a;
					}
				}
				if (@discord_sr1 >= $support_reads)
				{
					@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
					#print "@src_return1\n";
					if ($data[3] lt $data[7])
					{
						$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
						$right_bound_sr1 = $src_return1[2] - ($positiona2-$positiona1+101) + $positiona3;
						if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
						{
							$left_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$right_bound_sr1 = $positiona1 + $src_return1[0];
						}
						$left_bound_sr1 = $data[4] unless ($src_return1[1]);
						$right_bound_sr1 = $data[8] unless ($src_return1[2]);
					}
					else
					{
						$left_bound_sr2 = $positiona1 + $src_return1[0];
						$right_bound_sr2 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
						if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
						{
							$left_bound_sr2 = $src_return1[2] - ($positiona2-$positiona1+101) + $positiona3;
							$right_bound_sr2 = $positiona1 + $src_return1[1] + $cut_sr;
						}
						$left_bound_sr2 = $data[5] unless ($src_return1[0]);
						$right_bound_sr2 = $data[9] unless ($src_return1[3]);
					}
				}
				
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[1]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				@discord_sr2 = undef;
				for my $a (@alignments)
				{
					my @data1 = split (/\t/, $a);
					my $start = $data1[3];
					my $strand = 1;
					$strand = -1 if ($data1[1] =~ /r/);
					my $mseqid = $data1[6];
					my $mstart = $data1[7];
					my $mstrand = 1;
					$mstrand = -1 if ($data1[1] =~ /R/);
					my $isize = $data1[8];
					if ($isize>100 and $strand == $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
					{
						#print "2 $start $mstart\t$isize\t$strand $mstrand\n";
						push @discord_sr2, $a;
					}
				}
				if (@discord_sr2 >= $support_reads)
				{
					@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
					if ($data[3] lt $data[7])
					{
						$left_bound_sr2 = $positionb1 + $src_return2[0];
						$right_bound_sr2 = $src_return2[3] - ($positionb2-$positionb1+101) + $positionb3 + $cut_sr;
						if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
						{
							$left_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
							$right_bound_sr2 = $positionb1 + $src_return2[1] + $cut_sr;
						}
						$left_bound_sr2 = $data[5] unless ($src_return2[0]);
						$right_bound_sr2 = $data[9] unless ($src_return2[3]);
					}
					else
					{
						$left_bound_sr1 = $positionb1 + $src_return2[1] + $cut_sr;
						$right_bound_sr1 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
						if ($chr_bp_b1 eq $data[7] and $chr_bp_b2 eq $data[3])
						{
							$left_bound_sr1 = $src_return2[3] - ($positionb2-$positionb1+101) + $positionb3 + $cut_sr;
							$right_bound_sr1 = $positionb1 + $src_return2[0];
						}
						$left_bound_sr1 = $data[4] unless ($src_return2[1]);
						$right_bound_sr1 = $data[8] unless ($src_return2[2]);
					}
				}
				
				if ($src_return1[4] >= $support_reads and $src_return2[4] >= $support_reads)
				{
					my $reads0 = join("\t", @{$src_return1[5]{0}});
					my $reads1 = join("\t", @{$src_return2[5]{0}});
					if ($data[3] lt $data[7])
					{
						my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
						my $ins_size = $right_bound_sr2 - $right_bound_sr1 + 1;
						print BPREAD "$data[3]\__$left_bound_sr1\__$data[7]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[7]\__$right_bound_sr2\n$reads1\n";
						$result_sr[$sri][0] = $data[0];
						$result_sr[$sri][1] = $data[1];
						$result_sr[$sri][2] = $data[2];
						$result_sr[$sri][3] = "$src_return1[4]/$src_return2[4]";
						$result_sr[$sri][4] = $data[3];
						$result_sr[$sri][5] = $left_bound_sr1;
						$result_sr[$sri][6] = $left_bound_sr2;
						$result_sr[$sri][7] = $del_size;
						$result_sr[$sri][8] = "$data[7]\t$right_bound_sr1\t$right_bound_sr2\t$ins_size\n";
						$sri++;
					}
					else
					{
						($reads0, $reads1) = ($reads1, $reads0);
						my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
						my $ins_size = $right_bound_sr2 - $right_bound_sr1 + 1;
						print BPREAD "$data[3]\__$left_bound_sr1\__$data[7]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[7]\__$right_bound_sr2\n$reads1\n";
						$result_sr[$sri][0] = $data[0];
						$result_sr[$sri][1] = $data[1];
						$result_sr[$sri][2] = $data[2];
						$result_sr[$sri][3] = "$src_return2[4]/$src_return1[4]";
						$result_sr[$sri][4] = $data[3];
						$result_sr[$sri][5] = $left_bound_sr1;
						$result_sr[$sri][6] = $left_bound_sr2;
						$result_sr[$sri][7] = $del_size;
						$result_sr[$sri][8] = "$data[7]\t$right_bound_sr1\t$right_bound_sr2\t$ins_size\n";
						$sri++;
					}
				}
				else
				{
					my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
					my ($mpd1, $mpd2) = split (/\//, $data[2]);
					if ($src_return1[4] >= $support_reads)
					{
						if ($data[3] lt $data[7])
						{
							my $left_bound_sr = $positiona1 + $src_return1[1] + $cut_sr;
							my $right_bound_sr = $src_return1[2] - ($positiona2-$positiona1+101) + $positiona3;
							if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
							{
								$left_bound_sr = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
								$right_bound_sr = $positiona1 + $src_return1[0];
							}
							$left_bound_sr = $data[4] unless ($src_return1[1]);
							$right_bound_sr = $data[8] unless ($src_return1[2]);
							my $reads = join("\t", @{$src_return1[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[3]\t$left_bound_sr\t1\t$data[7]\t$right_bound_sr\t-1\n";
							$sri++;
						}
						else
						{
							my $left_bound_sr = $positiona1 + $src_return1[0];
							my $right_bound_sr = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
							{
								$left_bound_sr = $src_return1[2] - ($positiona2-$positiona1+101) + $positiona3;
								$right_bound_sr = $positiona1 + $src_return1[1] + $cut_sr;
							}
							$left_bound_sr = $data[5] unless ($src_return1[0]);
							$right_bound_sr = $data[9] unless ($src_return1[3]);
							my $reads = join("\t", @{$src_return1[5]{0}});
							print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[7]\t$right_bound_sr\t1\t$data[3]\t$left_bound_sr\t-1\n";
							$sri++;
						}
					}
					if ($src_return2[4] >= $support_reads)
					{
						if ($data[3] lt $data[7])
						{
							my $left_bound_sr = $positionb1 + $src_return2[0];
							my $right_bound_sr = $src_return2[3] - ($positionb2-$positionb1+101) + $positionb3 + $cut_sr;
							if ($chr_bp_b1 eq $data[7] and $chr_bp_b2 eq $data[3])
							{
								$left_bound_sr = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
								$right_bound_sr = $positionb1 + $src_return2[1] + $cut_sr;
							}
							$left_bound_sr = $data[5] unless ($src_return2[0]);
							$right_bound_sr = $data[9] unless ($src_return2[3]);
							my $reads = join("\t", @{$src_return2[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[3]\t$left_bound_sr\t-1\t$data[7]\t$right_bound_sr\t1\n";
							$sri++;
						}
						else
						{
							my $left_bound_sr = $positionb1 + $src_return2[1] + $cut_sr;
							my $right_bound_sr = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
							if ($chr_bp_b1 eq $data[7] and $chr_bp_b2 eq $data[3])
							{
								$left_bound_sr = $src_return2[3] - ($positionb2-$positionb1+101) + $positionb3 + $cut_sr;
								$right_bound_sr = $positionb1 + $src_return2[0];
							}
							$left_bound_sr = $data[4] unless ($src_return2[1]);
							$right_bound_sr = $data[8] unless ($src_return2[2]);
							my $reads = join("\t", @{$src_return2[5]{0}});
							print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[7]\t$right_bound_sr\t-1\t$data[3]\t$left_bound_sr\t1\n";
							$sri++;
						}
					}
				}
			}
		}
		if ($data[0] =~ /inso/)
		{
			next unless ($cluster_region{$cl[0]} and $cluster_region{$cl[1]});
			my ($chr_bp_a1, $positiona1, $positiona2, $chr_bp_a2, $positiona3, $positiona4, $orientationa) = split (/__/, $cluster_region{$cl[0]});
			my ($chr_bp_b1, $positionb1, $positionb2, $chr_bp_b2, $positionb3, $positionb4, $orientationb) = split (/__/, $cluster_region{$cl[1]});
			if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
			{
				if ($cluster_region{$cl[0]} eq $cluster_region{$cl[1]})
				{
					print STDERR "inso $i\n";
				}
				#($positiona1, $positiona2, $positiona3, $positiona4) = ($positiona3, $positiona4, $positiona1, $positiona2);
			}
				
			#if ($chr_bp_b1 eq $data[7] and $chr_bp_b2 eq $data[3]){#($positionb1, $positionb2, $positionb3, $positionb4) = ($positionb3, $positionb4, $positionb1, $positionb2);}
			#print "$cluster_region{$cl[0]}\n$cluster_region{$cl[1]}\n$positiona1, $positiona2, $positiona3, $positiona4\n$positionb1, $positionb2, $positionb3, $positionb4\n";
			
			if ($cluster_region{$cl[0]} eq $cluster_region{$cl[1]})
			{
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($orientationa == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
						{
							#print "$readname\t$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 1, 1);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						if ($$ref_boundary1[2] > $$ref_boundary1[0])
						{
							my $left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
							my $right_bound_sr1 = $positiona4 - ($$ref_boundary2[0] - ($positiona2-$positiona1+101));
							my $left_bound_sr2 = $positiona1 + $$ref_boundary1[2];
							my $right_bound_sr2 = $positiona4 - ($$ref_boundary2[3] - ($positiona2-$positiona1+101) + $cut_sr);
							if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
							{
								$left_bound_sr1 = $positiona4 - ($$ref_boundary2[0] - ($positiona2-$positiona1+101));
								$right_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
								$left_bound_sr2 = $positiona4 - ($$ref_boundary2[3] - ($positiona2-$positiona1+101) + $cut_sr);
								$right_bound_sr2 = $positiona1 + $$ref_boundary1[2];
							}
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[1]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[0]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[2]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[3]);
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
							{
								my $reads0 = join("\t", @{$$ref_bpread{0}});
								my $reads1 = join("\t", @{$$ref_bpread{1}});
								print BPREAD "$data[3]\__$left_bound_sr1\__$data[7]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[7]\__$right_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr1;
								$result_sr[$sri][6] = $left_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[7]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[0] >= $support_reads)
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr = $positiona4 - ($$ref_boundary2[0] - ($positiona2-$positiona1+101));
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $positiona4 - ($$ref_boundary2[0] - ($positiona2-$positiona1+101));
										$right_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									}
									$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
									$right_bound_sr = $data[9] unless ($$ref_boundary2[0]);
									my $reads = join("\t", @{$$ref_bpread{0}});
									if ($data[3] lt $data[7])
									{
										print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t1\t$data[7]\t$right_bound_sr\t1\n";
										$sri++;
									}
									else
									{
										print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[7]\t$right_bound_sr\t1\t$data[3]\t$left_bound_sr\t1\n";
										$sri++;
									}
								}
								if ($$ref_support_sr[1] >= $support_reads)
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[2];
									my $right_bound_sr = $positiona4 - ($$ref_boundary2[3] - ($positiona2-$positiona1+101) + $cut_sr);
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $positiona4 - ($$ref_boundary2[3] - ($positiona2-$positiona1+101) + $cut_sr);
										$right_bound_sr = $positiona1 + $$ref_boundary1[2];
									}
									$left_bound_sr = $data[5] unless ($$ref_boundary1[2]);
									$right_bound_sr = $data[8] unless ($$ref_boundary2[3]);
									my $reads = join("\t", @{$$ref_bpread{1}});
									if ($data[3] lt $data[7])
									{
										print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t-1\t$data[7]\t$right_bound_sr\t-1\n";
										$sri++;
									}
									else
									{
										print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[7]\t$right_bound_sr\t-1\t$data[3]\t$left_bound_sr\t-1\n";
										$sri++;
									}
								}
							}
						}
						else
						{
							my $left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
							my $right_bound_sr1 = $positiona4 - ($$ref_boundary2[2] - ($positiona2-$positiona1+101));
							my $left_bound_sr2 = $positiona1 + $$ref_boundary1[0];
							my $right_bound_sr2 = $positiona4 - ($$ref_boundary2[1] - ($positiona2-$positiona1+101) + $cut_sr);
							if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
							{
								$left_bound_sr1 = $positiona4 - ($$ref_boundary2[2] - ($positiona2-$positiona1+101));
								$right_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
								$left_bound_sr2 = $positiona4 - ($$ref_boundary2[1] - ($positiona2-$positiona1+101) + $cut_sr);
								$right_bound_sr2 = $positiona1 + $$ref_boundary1[0];
							}
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[3]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[2]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[0]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[1]);
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
							{
								my $reads0 = join("\t", @{$$ref_bpread{1}});
								my $reads1 = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr1\__$data[7]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[7]\__$right_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr1;
								$result_sr[$sri][6] = $left_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[7]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[1] >= $support_reads)
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr = $positiona4 - ($$ref_boundary2[2] - ($positiona2-$positiona1+101));
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $positiona4 - ($$ref_boundary2[2] - ($positiona2-$positiona1+101));
										$right_bound_sr = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									}
									$left_bound_sr = $data[4] unless ($$ref_boundary1[3]);
									$right_bound_sr = $data[9] unless ($$ref_boundary2[2]);
									my $reads = join("\t", @{$$ref_bpread{1}});
									if ($data[3] lt $data[7])
									{
										print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t1\t$data[7]\t$right_bound_sr\t1\n";
										$sri++;
									}
									else
									{
										print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[7]\t$right_bound_sr\t1\t$data[3]\t$left_bound_sr\t1\n";
										$sri++;
									}
								}
								if ($$ref_support_sr[0] >= $support_reads)
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[0];
									my $right_bound_sr = $positiona4 - ($$ref_boundary2[1] - ($positiona2-$positiona1+101) + $cut_sr);
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $positiona4 - ($$ref_boundary2[1] - ($positiona2-$positiona1+101) + $cut_sr);
										$right_bound_sr = $positiona1 + $$ref_boundary1[0];
									}
									$left_bound_sr = $data[5] unless ($$ref_boundary1[0]);
									$right_bound_sr = $data[8] unless ($$ref_boundary2[1]);
									my $reads = join("\t", @{$$ref_bpread{0}});
									if ($data[3] lt $data[7])
									{
										print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t-1\t$data[7]\t$right_bound_sr\t-1\n";
										$sri++;
									}
									else
									{
										print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[7]\t$right_bound_sr\t-1\t$data[3]\t$left_bound_sr\t-1\n";
										$sri++;
									}
								}
							}
						}
					}
				}
				else
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
						{
							#print "$readname\t$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my ($ref_boundary1, $ref_boundary2, $ref_support_sr, $ref_bpread) = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 1);
						#print "@$ref_boundary1\n@$ref_boundary2\n@$ref_support_sr\n";
						if ($$ref_boundary1[2] > $$ref_boundary1[0])
						{
							my $left_bound_sr1 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
							my $right_bound_sr1 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
							my $left_bound_sr2 = $positiona1 + $$ref_boundary1[2];
							my $right_bound_sr2 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
							{
								$left_bound_sr1 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
								$right_bound_sr1 = $positiona1 + $$ref_boundary1[0];
								$left_bound_sr2 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
								$right_bound_sr2 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
							}
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[1]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[0]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[2]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[3]);
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
							{
								my $reads0 = join("\t", @{$$ref_bpread{0}});
								my $reads1 = join("\t", @{$$ref_bpread{1}});
								print BPREAD "$data[3]\__$left_bound_sr1\__$data[7]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[7]\__$right_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[0]/$$ref_support_sr[1]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr1;
								$result_sr[$sri][6] = $left_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[7]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[0] >= $support_reads)
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									my $right_bound_sr = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
										$right_bound_sr = $positiona1 + $$ref_boundary1[0];
									}
									$left_bound_sr = $data[4] unless ($$ref_boundary1[1]);
									$right_bound_sr = $data[9] unless ($$ref_boundary2[0]);
									my $reads = join("\t", @{$$ref_bpread{0}});
									if ($data[3] lt $data[7])
									{
										print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t1\t$data[7]\t$right_bound_sr\t1\n";
										$sri++;
									}
									else
									{
										print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[0]\t$data[7]\t$right_bound_sr\t1\t$data[3]\t$left_bound_sr\t1\n";
										$sri++;
									}
								}
								if ($$ref_support_sr[1] >= $support_reads)
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[2];
									my $right_bound_sr = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
										$right_bound_sr = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									}
									$left_bound_sr = $data[5] unless ($$ref_boundary1[2]);
									$right_bound_sr = $data[8] unless ($$ref_boundary2[3]);
									my $reads = join("\t", @{$$ref_bpread{1}});
									if ($data[3] lt $data[7])
									{
										print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t-1\t$data[7]\t$right_bound_sr\t-1\n";
										$sri++;
									}
									else
									{
										print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[1]\t$data[7]\t$right_bound_sr\t-1\t$data[3]\t$left_bound_sr\t-1\n";
										$sri++;
									}
								}
							}
						}
						else
						{
							my $left_bound_sr1 = $positiona1 + $$ref_boundary1[3] + $cut_sr;
							my $right_bound_sr1 = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
							my $left_bound_sr2 = $positiona1 + $$ref_boundary1[0];
							my $right_bound_sr2 = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
							{
								$left_bound_sr1 = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
								$right_bound_sr1 = $positiona1 + $$ref_boundary1[2];
								$left_bound_sr2 = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
								$right_bound_sr2 = $positiona1 + $$ref_boundary1[1] + $cut_sr;
							}
							$left_bound_sr1 = $data[4] unless ($$ref_boundary1[3]);
							$right_bound_sr1 = $data[9] unless ($$ref_boundary2[2]);
							$left_bound_sr2 = $data[5] unless ($$ref_boundary1[0]);
							$right_bound_sr2 = $data[8] unless ($$ref_boundary2[1]);
							my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
							my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
							if ($$ref_support_sr[0] >= $support_reads and $$ref_support_sr[1] >= $support_reads)
							{
								my $reads0 = join("\t", @{$$ref_bpread{1}});
								my $reads1 = join("\t", @{$$ref_bpread{0}});
								print BPREAD "$data[3]\__$left_bound_sr1\__$data[7]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[7]\__$right_bound_sr2\n$reads1\n";
								$result_sr[$sri][0] = $data[0];
								$result_sr[$sri][1] = $data[1];
								$result_sr[$sri][2] = $data[2];
								$result_sr[$sri][3] = "$$ref_support_sr[1]/$$ref_support_sr[0]";
								$result_sr[$sri][4] = $data[3];
								$result_sr[$sri][5] = $left_bound_sr1;
								$result_sr[$sri][6] = $left_bound_sr2;
								$result_sr[$sri][7] = $del_size;
								$result_sr[$sri][8] = "$data[7]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
								$sri++;
							}
							else
							{
								my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
								my ($mpd1, $mpd2) = split (/\//, $data[2]);
								if ($$ref_support_sr[1] >= $support_reads)
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[3] + $cut_sr;
									my $right_bound_sr = $$ref_boundary2[2] - ($positiona2-$positiona1+101) + $positiona3;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
										$right_bound_sr = $positiona1 + $$ref_boundary1[2];
									}
									$left_bound_sr = $data[4] unless ($$ref_boundary1[3]);
									$right_bound_sr = $data[9] unless ($$ref_boundary2[2]);
									my $reads = join("\t", @{$$ref_bpread{1}});
									if ($data[3] lt $data[7])
									{
										print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[3]\t$left_bound_sr\t1\t$data[7]\t$right_bound_sr\t1\n";
										$sri++;
									}
									else
									{
										print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$$ref_support_sr[1]\t$data[7]\t$right_bound_sr\t1\t$data[3]\t$left_bound_sr\t1\n";
										$sri++;
									}
								}
								if ($$ref_support_sr[0] >= $support_reads)
								{
									my $left_bound_sr = $positiona1 + $$ref_boundary1[0];
									my $right_bound_sr = $$ref_boundary2[1] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
									if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
									{
										$left_bound_sr = $$ref_boundary2[0] - ($positiona2-$positiona1+101) + $positiona3;
										$right_bound_sr = $positiona1 + $$ref_boundary1[1] + $cut_sr;
									}
									$left_bound_sr = $data[5] unless ($$ref_boundary1[0]);
									$right_bound_sr = $data[8] unless ($$ref_boundary2[1]);
									my $reads = join("\t", @{$$ref_bpread{0}});
									if ($data[3] lt $data[7])
									{
										print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[3]\t$left_bound_sr\t-1\t$data[7]\t$right_bound_sr\t-1\n";
										$sri++;
									}
									else
									{
										print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
										$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$$ref_support_sr[0]\t$data[7]\t$right_bound_sr\t-1\t$data[3]\t$left_bound_sr\t-1\n";
										$sri++;
									}
								}
							}
						}
					}
				}
			}
			else
			{
				my (@src_return1, $left_bound_sr1, $right_bound_sr1, @src_return2, $left_bound_sr2, $right_bound_sr2);
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($orientationa == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						@src_return1 = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
						#print "@src_return1\n";
						$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
						$right_bound_sr1 = $positiona4 - ($src_return1[2] - ($positiona2-$positiona1+101));
						if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
						{
							$left_bound_sr1 = $positiona4 - ($src_return1[2] - ($positiona2-$positiona1+101));
							$right_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
						}
						$left_bound_sr1 = $data[4] unless ($src_return1[1]);
						$right_bound_sr1 = $data[9] unless ($src_return1[2]);
					}
				}
				else
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start < ($positiona2-$positiona1) and $mstart > ($positiona2-$positiona1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						@src_return1 = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
						#print "@src_return1\n";
						$left_bound_sr1 = $positiona1 + $src_return1[1] + $cut_sr;
						$right_bound_sr1 = $src_return1[2] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
						if ($chr_bp_a1 eq $data[7] and $chr_bp_a2 eq $data[3])
						{
							$left_bound_sr1 = $src_return1[3] - ($positiona2-$positiona1+101) + $positiona3 + $cut_sr;
							$right_bound_sr1 = $positiona1 + $src_return1[0] + $cut_sr;
						}
						$left_bound_sr1 = $data[4] unless ($src_return1[1]);
						$right_bound_sr1 = $data[9] unless ($src_return1[2]);
					}
				}
				
				
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[1]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($orientationb == 1)
				{
					@discord_sr2 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr2, $a;
						}
					}
					if (@discord_sr2 >= $support_reads)
					{
						@src_return2 = &sr_cluster($ref_discord_sr2, $ref_is, $cut_sr, 0, 0);
						$left_bound_sr2 = $positionb1 + $src_return2[0];
						$right_bound_sr2 = $positionb4 - ($src_return2[3] - ($positionb2-$positionb1+101) + $cut_sr);
						if ($chr_bp_b1 eq $data[7] and $chr_bp_b2 eq $data[3])
						{
							$left_bound_sr2 = $positionb4 - ($src_return2[3] - ($positionb2-$positionb1+101) + $cut_sr);
							$right_bound_sr2 = $positionb1 + $src_return2[0];
						}
						$left_bound_sr2 = $data[5] unless ($src_return2[0]);
						$right_bound_sr2 = $data[8] unless ($src_return2[3]);
					}
				}
				else
				{
					@discord_sr2 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start < ($positionb2-$positionb1) and $mstart > ($positionb2-$positionb1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr2, $a;
						}
					}
					if (@discord_sr2 >= $support_reads)
					{
						@src_return2 = &sr_cluster_1($ref_discord_sr2, $ref_is, $cut_sr, 0);
						$left_bound_sr2 = $positionb1 + $src_return2[0];
						$right_bound_sr2 = $src_return2[3] - ($positionb2-$positionb1+101) + $positionb3;
						if ($chr_bp_b1 eq $data[7] and $chr_bp_b2 eq $data[3])
						{
							$left_bound_sr2 = $src_return2[2] - ($positionb2-$positionb1+101) + $positionb3;
							$right_bound_sr2 = $positionb1 + $src_return2[1];
						}
						$left_bound_sr2 = $data[5] unless ($src_return2[0]);
						$right_bound_sr2 = $data[8] unless ($src_return2[3]);
					}
				}
				
				
				my $del_size = $left_bound_sr2 - $left_bound_sr1 - 1;
				my $ins_size = $right_bound_sr1 - $right_bound_sr2 + 1;
				if ($src_return1[4] >= $support_reads and $src_return2[4] >= $support_reads)
				{
					my $reads0 = join("\t", @{$src_return1[5]{0}});
					my $reads1 = join("\t", @{$src_return2[5]{0}});
					print BPREAD "$data[3]\__$left_bound_sr1\__$data[7]\__$right_bound_sr1\n$reads0\n$data[3]\__$left_bound_sr2\__$data[7]\__$right_bound_sr2\n$reads1\n";
					$result_sr[$sri][0] = $data[0];
					$result_sr[$sri][1] = $data[1];
					$result_sr[$sri][2] = $data[2];
					$result_sr[$sri][3] = "$src_return1[4]/$src_return2[4]";
					$result_sr[$sri][4] = $data[3];
					$result_sr[$sri][5] = $left_bound_sr1;
					$result_sr[$sri][6] = $left_bound_sr2;
					$result_sr[$sri][7] = $del_size;
					$result_sr[$sri][8] = "$data[7]\t$right_bound_sr2\t$right_bound_sr1\t$ins_size\n";
					$sri++;
				}
				else
				{
					my ($cluster_id1, $cluster_id2) = split (/\//, $data[1]);
					my ($mpd1, $mpd2) = split (/\//, $data[2]);
					if ($src_return1[4] >= $support_reads)
					{
						my $left_bound_sr = $left_bound_sr1;
						my $right_bound_sr = $right_bound_sr1;
						my $reads = join("\t", @{$src_return1[5]{0}});
						if ($data[3] lt $data[7])
						{
							print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[3]\t$left_bound_sr\t1\t$data[7]\t$right_bound_sr\t1\n";
							$sri++;
						}
						else
						{
							print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "transl_inter\t$cluster_id1\t$mpd1\t$src_return1[4]\t$data[7]\t$right_bound_sr\t1\t$data[3]\t$left_bound_sr\t1\n";
							$sri++;
						}
					}
					if ($src_return2[4] >= $support_reads)
					{
						my $left_bound_sr = $left_bound_sr2;
						my $right_bound_sr = $right_bound_sr2;
						my $reads = join("\t", @{$src_return2[5]{0}});
						if ($data[3] lt $data[7])
						{
							print BPREAD "$data[3]\__$left_bound_sr\__$data[7]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[3]\t$left_bound_sr\t-1\t$data[7]\t$right_bound_sr\t-1\n";
							$sri++;
						}
						else
						{
							print BPREAD "$data[7]\__$right_bound_sr\__$data[3]\__$left_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "transl_inter\t$cluster_id2\t$mpd2\t$src_return2[4]\t$data[7]\t$right_bound_sr\t-1\t$data[3]\t$left_bound_sr\t-1\n";
							$sri++;
						}
					}
				}
			}
		}
		if ($data[0] eq 'transl_inter')
		{
			if ($data[5] == 1 and $data[8] == -1)
			{
				next unless ($cluster_region{$cl[0]});
				my ($chr_bp_1, $position1, $position2, $chr_bp_2, $position3, $position4, $orientation) = split (/__/, $cluster_region{$cl[0]});
				if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
				{
					($position1, $position2, $position3, $position4) = ($position3, $position4, $position1, $position2);
				}
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				@discord_sr1 = undef;
				for my $a (@alignments)
				{
					my @data1 = split (/\t/, $a);
					my $start = $data1[3];
					my $strand = 1;
					$strand = -1 if ($data1[1] =~ /r/);
					my $mseqid = $data1[6];
					my $mstart = $data1[7];
					my $mstrand = 1;
					$mstrand = -1 if ($data1[1] =~ /R/);
					my $isize = $data1[8];
					if ($isize>100 and $strand == $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
					{
						#print "$start $mstart\t$isize\t$strand $mstrand\n";
						push @discord_sr1, $a;
					}
				}
				if (@discord_sr1 >= $support_reads)
				{
					my @src_return = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
					my $left_bound_sr = $position1 + $src_return[1] + $cut_sr;
					my $right_bound_sr = $src_return[2] - ($position2-$position1+101) + $position3;
					if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
					{
						$left_bound_sr = $position1 + $src_return[1];
						$right_bound_sr = $src_return[2] - ($position2-$position1+101) + $position3 + $cut_sr;
					}
					$left_bound_sr = $data[4] unless ($src_return[1]);
					$right_bound_sr = $data[7] unless ($src_return[2]);
					if ($src_return[4] >= $support_reads)
					{
						my $reads = join("\t", @{$src_return[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[6]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$data[5]\t$data[6]\t$right_bound_sr\t$data[8]\n";
						$sri++;
					}
				}
			}
			if ($data[5] == -1 and $data[8] == 1)
			{
				next unless ($cluster_region{$cl[0]});
				my ($chr_bp_1, $position1, $position2, $chr_bp_2, $position3, $position4, $orientation) = split (/__/, $cluster_region{$cl[0]});
				if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
				{
					($position1, $position2, $position3, $position4) = ($position3, $position4, $position1, $position2);
				}
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				@discord_sr1 = undef;
				for my $a (@alignments)
				{
					my @data1 = split (/\t/, $a);
					my $start = $data1[3];
					my $strand = 1;
					$strand = -1 if ($data1[1] =~ /r/);
					my $mseqid = $data1[6];
					my $mstart = $data1[7];
					my $mstrand = 1;
					$mstrand = -1 if ($data1[1] =~ /R/);
					my $isize = $data1[8];
					if ($isize>100 and $strand == $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
					{
						#print "$start $mstart\t$isize\t$strand $mstrand\n";
						push @discord_sr1, $a;
					}
				}
				if (@discord_sr1 >= $support_reads)
				{
					my @src_return = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
					my $left_bound_sr = $position1 + $src_return[0];
					my $right_bound_sr = $src_return[3] - ($position2-$position1+101) + $cut_sr + $position3;
					if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
					{
						$left_bound_sr = $position1 + $src_return[0] + $cut_sr;
						$right_bound_sr = $src_return[3] - ($position2-$position1+101) + $position3;
					}
					$left_bound_sr = $data[4] unless ($src_return[0]);
					$right_bound_sr = $data[7] unless ($src_return[3]);
					if ($src_return[4] >= $support_reads)
					{
						my $reads = join("\t", @{$src_return[5]{0}});
						print BPREAD "$data[3]\__$left_bound_sr\__$data[6]\__$right_bound_sr\n$reads\n";
						$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$data[5]\t$data[6]\t$right_bound_sr\t$data[8]\n";
						$sri++;
					}
				}
			}
			if ($data[5] == 1 and $data[8] == 1)
			{
				next unless ($cluster_region{$cl[0]});
				my ($chr_bp_1, $position1, $position2, $chr_bp_2, $position3, $position4, $orientation) = split (/__/, $cluster_region{$cl[0]});
				if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
				{
					($position1, $position2, $position3, $position4) = ($position3, $position4, $position1, $position2);
				}
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($orientation == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my @src_return = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
						my $left_bound_sr = $position1 + $src_return[1] + $cut_sr;
						my $right_bound_sr = $position4 - ($src_return[2] - ($position2-$position1+101));
						if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
						{
							$left_bound_sr = $position1 + $src_return[1];
							$right_bound_sr = $position4 - ($src_return[2] - ($position2-$position1+101) + $cut_sr);
						}
						$left_bound_sr = $data[4] unless ($src_return[1]);
						$right_bound_sr = $data[7] unless ($src_return[2]);
						if ($src_return[4] >= $support_reads)
						{
							my $reads = join("\t", @{$src_return[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[6]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$data[5]\t$data[6]\t$right_bound_sr\t$data[8]\n";
							$sri++;
						}
					}
				}
				else
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my @src_return = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
						my $left_bound_sr = $position1 + $src_return[1] + $cut_sr;
						my $right_bound_sr = $src_return[2] - ($position2-$position1+101) + $position3;
						if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
						{
							$left_bound_sr = $position1 + $src_return[1];
							$right_bound_sr = $src_return[2] - ($position2-$position1+101) + $position3 + $cut_sr;
						}
						$left_bound_sr = $data[4] unless ($src_return[1]);
						$right_bound_sr = $data[7] unless ($src_return[2]);
						if ($src_return[4] >= $support_reads)
						{
							my $reads = join("\t", @{$src_return[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[6]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$data[5]\t$data[6]\t$right_bound_sr\t$data[8]\n";
							$sri++;
						}
					}
				}
				
			}
			if ($data[5] == -1 and $data[8] == -1)
			{
				next unless ($cluster_region{$cl[0]});
				my ($chr_bp_1, $position1, $position2, $chr_bp_2, $position3, $position4, $orientation) = split (/__/, $cluster_region{$cl[0]});
				if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
				{
					($position1, $position2, $position3, $position4) = ($position3, $position4, $position1, $position2);
				}
				my @alignments;
				open (SAM, "$samtools_command view -X $sr_sortbam $cluster_region{$cl[0]}|");
				while ($newline1 = <SAM>)
				{
					chomp $newline1;
					push @alignments, $newline1;
				}
				close SAM;
				if ($orientation == 1)
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>100 and $strand == $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my @src_return = &sr_cluster($ref_discord_sr1, $ref_is, $cut_sr, 0, 0);
						my $left_bound_sr = $position1 + $src_return[0];
						my $right_bound_sr = $position4 - ($src_return[3] - ($position2-$position1+101) + $cut_sr);
						if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
						{
							$left_bound_sr = $position1 + $src_return[0] + $cut_sr;
							$right_bound_sr = $position4 - ($src_return[3] - ($position2-$position1+101));
						}
						$left_bound_sr = $data[4] unless ($src_return[0]);
						$right_bound_sr = $data[7] unless ($src_return[3]);
						if ($src_return[4] >= $support_reads)
						{
							my $reads = join("\t", @{$src_return[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[6]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$data[5]\t$data[6]\t$right_bound_sr\t$data[8]\n";
							$sri++;
						}
					}
				}
				else
				{
					@discord_sr1 = undef;
					for my $a (@alignments)
					{
						my @data1 = split (/\t/, $a);
						my $start = $data1[3];
						my $strand = 1;
						$strand = -1 if ($data1[1] =~ /r/);
						my $mseqid = $data1[6];
						my $mstart = $data1[7];
						my $mstrand = 1;
						$mstrand = -1 if ($data1[1] =~ /R/);
						my $isize = $data1[8];
						if ($isize>0 and $strand != $mstrand and $start < ($position2-$position1) and $mstart > ($position2-$position1+100))
						{
							#print "$start $mstart\t$isize\t$strand $mstrand\n";
							push @discord_sr1, $a;
						}
					}
					if (@discord_sr1 >= $support_reads)
					{
						my @src_return = &sr_cluster_1($ref_discord_sr1, $ref_is, $cut_sr, 0);
						my $left_bound_sr = $position1 + $src_return[0];
						my $right_bound_sr = $src_return[3] - ($position2-$position1+101) + $position3 + $cut_sr;
						if ($chr_bp_1 eq $data[6] and $chr_bp_2 eq $data[3])
						{
							$left_bound_sr = $position1 + $src_return[0] + $cut_sr;
							$right_bound_sr = $src_return[3] - ($position2-$position1+101) + $position3;
						}
						$left_bound_sr = $data[4] unless ($src_return[0]);
						$right_bound_sr = $data[7] unless ($src_return[3]);
						if ($src_return[4] >= $support_reads)
						{
							my $reads = join("\t", @{$src_return[5]{0}});
							print BPREAD "$data[3]\__$left_bound_sr\__$data[6]\__$right_bound_sr\n$reads\n";
							$result_sr[$sri][0] = "$data[0]\t$data[1]\t$data[2]\t$src_return[4]\t$data[3]\t$left_bound_sr\t$data[5]\t$data[6]\t$right_bound_sr\t$data[8]\n";
							$sri++;
						}
					}
				}
				
			}
		}
	}
	close MPINTERSRD;
	
	close BPREAD;
	
	# add deletion clusters without split read support into other events if supporting complex deletion events
	foreach (@unsupport_del)
	{
		my @del = @{$_};
		#my $toprint = join ("\t", @del); print "$toprint";
		foreach (@result_sr)
		{
			if ($$_[0] =~ /del/)
			{
				if ($del[3] eq $$_[4] and abs($del[4]-$$_[5])<$$ref_is{'rlu'} and abs($del[5]-$$_[6])<$$ref_is{'rlu'})
				{
					$$_[1] .= '/'.$del[1];
					$$_[2] .= '/'.$del[2];
				}
			}
		}
	}
	
	open SRINTEROUT, ">$srinter_outfile";
	open SRINTRAOUT, ">$srintra_outfile";
	foreach (@result_sr)
	{
		my $toprint = join ("\t", @{$_});
		if ($$_[0] =~ /transl_inter/ or $$_[0] eq 'del_inss' or $$_[0] eq 'del_inso' or $$_[0] eq 'inss' or $$_[0] eq 'inso')
		{
			print SRINTEROUT "$toprint";
		}
		else
		{
			print SRINTRAOUT "$toprint";
		}
	}
	close SRINTRAOUT;
	close SRINTEROUT;
}


#	construct search space and align split reads
sub alg
{
	my $line = shift;
	my $prefix = shift;
	my $bwa_command = shift;
	my $db = shift;
	my $ref_is = shift;
	my $cut_sr = shift;
	my $threads_bwa = shift;
	my $mpintra_outfile = $prefix.'.mp.intra.out';
	my $mpinter_outfile = $prefix.'.mp.inter.out';
	my $fafile = $prefix.'.bp.fasta';
	my $bpinfofile = $prefix.'.bp.info';
	my $sr_insertsize = 2*$cut_sr + 100;
	my $faifile = $prefix.'.bp.fasta.fai';
	my $sr_fq1 = $prefix.'.sr.1.fq.gz';
	my $sr_fq2 = $prefix.'.sr.2.fq.gz';
	my $fq1_pipe = $prefix.'.sr.1.fq.pipe';
	my $fq2_pipe = $prefix.'.sr.2.fq.pipe';
	my $sr_sai1 = $prefix.'.sr.1.sai';
	my $sr_sai2 = $prefix.'.sr.2.sai';
	my $sr_rawsam = $prefix.'.sr.raw.sam';
	my $sr_sam = $prefix.'.sr.sam';
	my $sr_bam = $prefix.'.sr.bam';
	my $sr_sort = $prefix.'.sr.sorted';
	my $sr_sortbam = $prefix.'.sr.sorted.bam';
	goto ALGS1 if ($algs == 2);
	
	my $n;
	for (my $i=0;$i<100;$i++)
	{
		$n .= 'N';
	}
	my ($newline, $i, %bp_weight, %regions, %cluster_exist);
#	%bp_weight: where is the real break point more likely to be located
#	$bp_weight{regionname}[0]: left break point
#	$bp_weight{regionname}[1]: right break point
#	%regions, a list of existing regions
#	$regions{$cluster_id}: region name, chr_p_p_chr_p_p_0/1, 0 same orientation, 1 opposite orientation
#	%cluster_exist: the cluster with key id exists
#	$cluster_exist{$cluster_id}: 1 exist, 0 not
	open MPINTRAALG, "<$mpintra_outfile";
	while ($newline = <MPINTRAALG>)
	{
		chomp $newline;
		$i++;
		next if ($i < $line);
		my @data = split (/\t/, $newline);
		#print "$data[1]\t$data[7]\n";
		my @cl = split (/\//, $data[1]);
		if ($data[0] eq 'del')
		{
			my ($position1, $position2, $position3, $position4) = split (/:/, $data[7]);
			$position1 = $position1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$position2 = $position2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$position3 = $position3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$position4 = $position4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $name = $data[3].'__'.$position1.'__'.$position2.'__'.$data[3].'__'.$position3.'__'.$position4.'__0';
			$bp_weight{$name}[0] = 1;
			$bp_weight{$name}[1] = 4;
			unless ($cluster_exist{$cl[0]})
			{
				$regions{a}{$data[3]}{$cl[0]} = $name;
				$cluster_exist{$cl[0]} = 1;
			}
		}
		if ($data[0] =~ 'inss')
		{
			my ($positiona1, $positiona2, $positiona3, $positiona4) = split (/:/, $data[12]);
			$positiona1 = $positiona1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona2 = $positiona2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona3 = $positiona3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona4 = $positiona4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $namea = $data[3].'__'.$positiona1.'__'.$positiona2.'__'.$data[3].'__'.$positiona3.'__'.$positiona4.'__0';
			$bp_weight{$namea}[0] = 2;
			$bp_weight{$namea}[1] = 3;
			my ($positionb1, $positionb2, $positionb3, $positionb4) = split (/:/, $data[13]);
			$positionb1 = $positionb1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb2 = $positionb2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb3 = $positionb3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb4 = $positionb4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $nameb = $data[3].'__'.$positionb1.'__'.$positionb2.'__'.$data[3].'__'.$positionb3.'__'.$positionb4.'__0';
			$bp_weight{$nameb}[0] = 1;
			$bp_weight{$nameb}[1] = 4;
			unless ($cluster_exist{$cl[0]})
			{
				$regions{a}{$data[3]}{$cl[0]} = $namea;
				$cluster_exist{$cl[0]} = 1;
			}
			unless ($cluster_exist{$cl[1]})
			{
				$regions{a}{$data[3]}{$cl[1]} = $nameb;
				$cluster_exist{$cl[1]} = 1;
			}
		}
		if ($data[0] =~ 'inso')
		{
			my ($positiona1, $positiona2, $positiona3, $positiona4) = split (/:/, $data[12]);
			$positiona1 = $positiona1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona2 = $positiona2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona3 = $positiona3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona4 = $positiona4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $namea = $data[3].'__'.$positiona1.'__'.$positiona2.'__'.$data[3].'__'.$positiona3.'__'.$positiona4.'__1';
			my ($positionb1, $positionb2, $positionb3, $positionb4) = split (/:/, $data[13]);
			$positionb1 = $positionb1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb2 = $positionb2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb3 = $positionb3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb4 = $positionb4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $nameb = $data[3].'__'.$positionb1.'__'.$positionb2.'__'.$data[3].'__'.$positionb3.'__'.$positionb4.'__1';
			if ($data[0] =~ 'insod')
			{
				$bp_weight{$namea}[0] = 1;
				$bp_weight{$namea}[1] = 3;
				$bp_weight{$nameb}[0] = 2;
				$bp_weight{$nameb}[1] = 4;
			}
			else
			{
				$bp_weight{$namea}[0] = 2;
				$bp_weight{$namea}[1] = 4;
				$bp_weight{$nameb}[0] = 1;
				$bp_weight{$nameb}[1] = 3;
			}
			unless ($cluster_exist{$cl[0]})
			{
				$regions{a}{$data[3]}{$cl[0]} = $namea;
				$cluster_exist{$cl[0]} = 1;
			}
			unless ($cluster_exist{$cl[1]})
			{
				$regions{a}{$data[3]}{$cl[1]} = $nameb;
				$cluster_exist{$cl[1]} = 1;
			}
		}
		if ($data[0] eq 'invers' or $data[0] eq 'del_invers')
		{
			my ($positiona1, $positiona2, $positiona3, $positiona4) = split (/:/, $data[10]);
			$positiona1 = $positiona1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona2 = $positiona2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona3 = $positiona3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona4 = $positiona4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $namea = $data[3].'__'.$positiona1.'__'.$positiona2.'__'.$data[3].'__'.$positiona3.'__'.$positiona4.'__1';
			my ($positionb1, $positionb2, $positionb3, $positionb4) = split (/:/, $data[11]);
			$positionb1 = $positionb1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb2 = $positionb2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb3 = $positionb3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb4 = $positionb4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $nameb = $data[3].'__'.$positionb1.'__'.$positionb2.'__'.$data[3].'__'.$positionb3.'__'.$positionb4.'__1';
			unless (&covered($positiona1, $positiona2, $positionb1, $positionb2) and &covered($positiona3, $positiona4, $positionb3, $positionb4))
			{
				$bp_weight{$namea}[0] = 1;
				$bp_weight{$namea}[1] = 3;
				$bp_weight{$nameb}[0] = 2;
				$bp_weight{$nameb}[1] = 4;
			}
			unless ($cluster_exist{$cl[0]})
			{
				$regions{a}{$data[3]}{$cl[0]} = $namea;
				$cluster_exist{$cl[0]} = 1;
			}
			unless ($cluster_exist{$cl[1]})
			{
				$regions{a}{$data[3]}{$cl[1]} = $nameb;
				$cluster_exist{$cl[1]} = 1;
			}
		}
		if ($data[0] eq 'tandem_dup')
		{
			my ($position1, $position2, $position3, $position4) = split (/:/, $data[7]);
			$position1 = $position1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$position2 = $position2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$position3 = $position3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$position4 = $position4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $name = $data[3].'__'.$position1.'__'.$position2.'__'.$data[3].'__'.$position3.'__'.$position4.'__0';
			$bp_weight{$name}[0] = 2;
			$bp_weight{$name}[1] = 3;
			unless ($cluster_exist{$cl[0]})
			{
				$regions{a}{$data[3]}{$cl[0]} = $name;
				$cluster_exist{$cl[0]} = 1;
			}
		}
		if ($data[0] eq 'invers_f')
		{
			my ($position1, $position2, $position3, $position4) = split (/:/, $data[7]);
			$position1 = $position1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$position2 = $position2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$position3 = $position3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$position4 = $position4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $name = $data[3].'__'.$position1.'__'.$position2.'__'.$data[3].'__'.$position3.'__'.$position4.'__1';
			$bp_weight{$name}[0] = 1;
			$bp_weight{$name}[1] = 3;
			unless ($cluster_exist{$cl[0]})
			{
				$regions{a}{$data[3]}{$cl[0]} = $name;
				$cluster_exist{$cl[0]} = 1;
			}
		}
		if ($data[0] eq 'invers_r')
		{
			my ($position1, $position2, $position3, $position4) = split (/:/, $data[7]);
			$position1 = $position1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$position2 = $position2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$position3 = $position3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$position4 = $position4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $name = $data[3].'__'.$position1.'__'.$position2.'__'.$data[3].'__'.$position3.'__'.$position4.'__1';
			$bp_weight{$name}[0] = 2;
			$bp_weight{$name}[1] = 4;
			unless ($cluster_exist{$cl[0]})
			{
				$regions{a}{$data[3]}{$cl[0]} = $name;
				$cluster_exist{$cl[0]} = 1;
			}
		}
	}
	close MPINTRAALG;
	
	open MPINTERALG, "<$mpinter_outfile";
	while ($newline = <MPINTERALG>)
	{
		chomp $newline;
		$i++;
		next if ($i < $line);
		my @data = split (/\t/, $newline);
		my @cl = split (/\//, $data[1]);
		if ($data[0] =~ /inss/)
		{
			my ($positiona1, $positiona2, $positiona3, $positiona4) = split (/:/, $data[11]);
			$positiona1 = $positiona1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona2 = $positiona2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona3 = $positiona3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona4 = $positiona4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $namea = $data[3].'__'.$positiona1.'__'.$positiona2.'__'.$data[7].'__'.$positiona3.'__'.$positiona4.'__0';
			$bp_weight{$namea}[0] = 2;
			$bp_weight{$namea}[1] = 3;
			my ($positionb1, $positionb2, $positionb3, $positionb4) = split (/:/, $data[12]);
			$positionb1 = $positionb1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb2 = $positionb2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb3 = $positionb3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb4 = $positionb4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $nameb = $data[3].'__'.$positionb1.'__'.$positionb2.'__'.$data[7].'__'.$positionb3.'__'.$positionb4.'__0';
			$bp_weight{$nameb}[0] = 1;
			$bp_weight{$nameb}[1] = 4;
			unless ($cluster_exist{$cl[0]})
			{
				$regions{e}{$data[3]}{$data[7]}{$cl[0]} = $namea;
				$cluster_exist{$cl[0]} = 1;
			}
			unless ($cluster_exist{$cl[1]})
			{
				$regions{e}{$data[3]}{$data[7]}{$cl[1]} = $nameb;
				$cluster_exist{$cl[1]} = 1;
			}
		}
		if ($data[0] =~ /inso/)
		{
			my ($positiona1, $positiona2, $positiona3, $positiona4) = split (/:/, $data[11]);
			$positiona1 = $positiona1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona2 = $positiona2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona3 = $positiona3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positiona4 = $positiona4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $namea = $data[3].'__'.$positiona1.'__'.$positiona2.'__'.$data[7].'__'.$positiona3.'__'.$positiona4.'__1';
			$bp_weight{$namea}[0] = 2;
			$bp_weight{$namea}[1] = 3;
			my ($positionb1, $positionb2, $positionb3, $positionb4) = split (/:/, $data[12]);
			$positionb1 = $positionb1 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb2 = $positionb2 + ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb3 = $positionb3 - ($$ref_is{'rlu'} - $cut_sr + 10);
			$positionb4 = $positionb4 + ($$ref_is{'rlu'} - $cut_sr + 10);
			my $nameb = $data[3].'__'.$positionb1.'__'.$positionb2.'__'.$data[7].'__'.$positionb3.'__'.$positionb4.'__1';
			$bp_weight{$nameb}[0] = 1;
			$bp_weight{$nameb}[1] = 4;
			unless ($cluster_exist{$cl[0]})
			{
				$regions{e}{$data[3]}{$data[7]}{$cl[0]} = $namea;
				$cluster_exist{$cl[0]} = 1;
			}
			unless ($cluster_exist{$cl[1]})
			{
				$regions{e}{$data[3]}{$data[7]}{$cl[1]} = $nameb;
				$cluster_exist{$cl[1]} = 1;
			}
		}
		if ($data[0] eq 'transl_inter')
		{
			if ($data[5] == 1 and $data[8] == -1)
			{
				my ($position1, $position2, $position3, $position4) = split (/:/, $data[9]);
				$position1 = $position1 - ($$ref_is{'rlu'} - $cut_sr + 10);
				$position2 = $position2 + ($$ref_is{'rlu'} - $cut_sr + 10);
				$position3 = $position3 - ($$ref_is{'rlu'} - $cut_sr + 10);
				$position4 = $position4 + ($$ref_is{'rlu'} - $cut_sr + 10);
				my $name = $data[3].'__'.$position1.'__'.$position2.'__'.$data[6].'__'.$position3.'__'.$position4.'__0';
				$bp_weight{$name}[0] = 1;
				$bp_weight{$name}[1] = 4;
				unless ($cluster_exist{$cl[0]})
				{
					$regions{e}{$data[3]}{$data[6]}{$cl[0]} = $name;
					$cluster_exist{$cl[0]} = 1;
				}
			}
			if ($data[5] == -1 and $data[8] == 1)
			{
				my ($position1, $position2, $position3, $position4) = split (/:/, $data[9]);
				$position1 = $position1 - ($$ref_is{'rlu'} - $cut_sr + 10);
				$position2 = $position2 + ($$ref_is{'rlu'} - $cut_sr + 10);
				$position3 = $position3 - ($$ref_is{'rlu'} - $cut_sr + 10);
				$position4 = $position4 + ($$ref_is{'rlu'} - $cut_sr + 10);
				my $name = $data[3].'__'.$position1.'__'.$position2.'__'.$data[6].'__'.$position3.'__'.$position4.'__0';
				$bp_weight{$name}[0] = 2;
				$bp_weight{$name}[1] = 3;
				unless ($cluster_exist{$cl[0]})
				{
					$regions{e}{$data[3]}{$data[6]}{$cl[0]} = $name;
					$cluster_exist{$cl[0]} = 1;
				}
			}
			if ($data[5] == 1 and $data[8] == 1)
			{
				my ($position1, $position2, $position3, $position4) = split (/:/, $data[9]);
				$position1 = $position1 - ($$ref_is{'rlu'} - $cut_sr + 10);
				$position2 = $position2 + ($$ref_is{'rlu'} - $cut_sr + 10);
				$position3 = $position3 - ($$ref_is{'rlu'} - $cut_sr + 10);
				$position4 = $position4 + ($$ref_is{'rlu'} - $cut_sr + 10);
				my $name = $data[3].'__'.$position1.'__'.$position2.'__'.$data[6].'__'.$position3.'__'.$position4.'__1';
				$bp_weight{$name}[0] = 1;
				$bp_weight{$name}[1] = 3;
				unless ($cluster_exist{$cl[0]})
				{
					$regions{e}{$data[3]}{$data[6]}{$cl[0]} = $name;
					$cluster_exist{$cl[0]} = 1;
				}
			}
			if ($data[5] == -1 and $data[8] == -1)
			{
				my ($position1, $position2, $position3, $position4) = split (/:/, $data[9]);
				$position1 = $position1 - ($$ref_is{'rlu'} - $cut_sr + 10);
				$position2 = $position2 + ($$ref_is{'rlu'} - $cut_sr + 10);
				$position3 = $position3 - ($$ref_is{'rlu'} - $cut_sr + 10);
				$position4 = $position4 + ($$ref_is{'rlu'} - $cut_sr + 10);
				my $name = $data[3].'__'.$position1.'__'.$position2.'__'.$data[6].'__'.$position3.'__'.$position4.'__1';
				$bp_weight{$name}[0] = 2;
				$bp_weight{$name}[1] = 4;
				unless ($cluster_exist{$cl[0]})
				{
					$regions{e}{$data[3]}{$data[6]}{$cl[0]} = $name;
					$cluster_exist{$cl[0]} = 1;
				}
			}
		}
	}
	close MPINTERALG;
	
	# merge overlapping regions
	open FA, ">$fafile";
	open BPDT, ">$bpinfofile";
	foreach my $chr (keys(%{$regions{a}}))
	{
		my $merge = 1;
		while ($merge)
		{
			$merge = 0;
			foreach my $key (keys(%{$regions{a}{$chr}}))
			{
				my @data = split (/__/, $regions{a}{$chr}{$key});
				if (&covered($data[1], $data[2], $data[4], $data[5]))
				{
					$merge = 1;
					my @temp = ($data[1], $data[2], $data[4], $data[5]);
					@temp = sort {$a <=> $b} @temp;
					$regions{a}{$chr}{$key} = $data[0].'__'.$temp[0].'__'.$temp[-1];
					#print STDERR "@data\n$regions{a}{$chr}{$key}\n";
				}
			}
MA:			foreach my $key1 (keys(%{$regions{a}{$chr}}))
			{
				my @data1 = split (/__/, $regions{a}{$chr}{$key1});
				foreach my $key2 (keys(%{$regions{a}{$chr}}))
				{
					next if ($key1 eq $key2);
					my @data2 = split (/__/, $regions{a}{$chr}{$key2});
					if ($data1[4] and $data2[4])
					{
						if (&covered($data1[1], $data1[2], $data2[1], $data2[2]) and &covered($data1[4], $data1[5], $data2[4], $data2[5]))
						{
							$merge = 1;
							my @temp1 = ($data1[1], $data1[2], $data2[1], $data2[2]);
							@temp1 = sort {$a <=> $b} @temp1;
							my @temp2 = ($data1[4], $data1[5], $data2[4], $data2[5]);
							@temp2 = sort {$a <=> $b} @temp2;
							delete $regions{a}{$chr}{$key1};
							delete $regions{a}{$chr}{$key2};
							my $newkey = $key1.'/'.$key2;
							if ($data1[6] eq $data2[6])
							{
								$regions{a}{$chr}{$newkey} = $data1[0].'__'.$temp1[0].'__'.$temp1[-1].'__'.$data1[3].'__'.$temp2[0].'__'.$temp2[-1].'__'.$data1[6];
							}
							else
							{
								$regions{a}{$chr}{$newkey} = $data1[0].'__'.$temp1[0].'__'.$temp1[-1].'__'.$data1[3].'__'.$temp2[0].'__'.$temp2[-1].'__10';
							}
							#print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
							next MA;
						}
					}
					if ($data1[4] and !($data2[4]))
					{
						if (&covered($data1[1], $data1[2], $data2[1], $data2[2]))
						{
							$merge = 1;
							my @temp = ($data1[1], $data1[2], $data2[1], $data2[2]);
							@temp = sort {$a <=> $b} @temp;
							delete $regions{a}{$chr}{$key1};
							delete $regions{a}{$chr}{$key2};
							my $newkey = $key1.'/'.$key2;
							$regions{a}{$chr}{$newkey} = $data1[0].'__'.$temp[0].'__'.$temp[-1].'__'.$data1[3].'__'.$data1[4].'__'.$data1[5].'__'.$data1[6];
							#print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
							next MA;
						}
						if (&covered($data1[4], $data1[5], $data2[1], $data2[2]))
						{
							$merge = 1;
							my @temp = ($data1[4], $data1[5], $data2[1], $data2[2]);
							@temp = sort {$a <=> $b} @temp;
							delete $regions{a}{$chr}{$key1};
							delete $regions{a}{$chr}{$key2};
							my $newkey = $key1.'/'.$key2;
							$regions{a}{$chr}{$newkey} = $data1[0].'__'.$data1[1].'__'.$data1[2].'__'.$data1[3].'__'.$temp[0].'__'.$temp[-1].'__'.$data1[6];
							#print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
							next MA;
						}
					}
					if (!($data1[4]) and $data2[4])
					{
						if (&covered($data1[1], $data1[2], $data2[1], $data2[2]))
						{
							$merge = 1;
							my @temp = ($data1[1], $data1[2], $data2[1], $data2[2]);
							@temp = sort {$a <=> $b} @temp;
							delete $regions{a}{$chr}{$key1};
							delete $regions{a}{$chr}{$key2};
							my $newkey = $key1.'/'.$key2;
							$regions{a}{$chr}{$newkey} = $data2[0].'__'.$temp[0].'__'.$temp[-1].'__'.$data2[3].'__'.$data2[4].'__'.$data2[5].'__'.$data2[6];
							#print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
							next MA;
						}
						if (&covered($data1[4], $data1[5], $data2[1], $data2[2]))
						{
							$merge = 1;
							my @temp = ($data1[4], $data1[5], $data2[1], $data2[2]);
							@temp = sort {$a <=> $b} @temp;
							delete $regions{a}{$chr}{$key1};
							delete $regions{a}{$chr}{$key2};
							my $newkey = $key1.'/'.$key2;
							$regions{a}{$chr}{$newkey} = $data2[0].'__'.$data2[1].'__'.$data2[2].'__'.$data2[3].'__'.$temp[0].'__'.$temp[-1].'__'.$data2[6];
							#print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
							next MA;
						}
					}
					if (!($data1[4]) and !($data2[4]))
					{
						if (&covered($data1[1], $data1[2], $data2[1], $data2[2]))
						{
							$merge = 1;
							my @temp = ($data1[1], $data1[2], $data2[1], $data2[2]);
							@temp = sort {$a <=> $b} @temp;
							delete $regions{a}{$chr}{$key1};
							delete $regions{a}{$chr}{$key2};
							my $newkey = $key1.'/'.$key2;
							$regions{a}{$chr}{$newkey} = $data1[0].'__'.$temp[0].'__'.$temp[-1];
							#print STDERR "@data1\n@data2\n$regions{a}{$chr}{$newkey}\n";
							next MA;
						}
					}
				}
			}
		}
		foreach my $key (keys(%{$regions{a}{$chr}}))
		{
			my @data = split (/__/, $regions{a}{$chr}{$key});
			next if ($data[1]<0 or $data[2]<0 or $data[4]<0 or $data[5]<0);
			my $seq = $db->seq($data[0], $data[1] => $data[2]);
			if ($data[3])
			{
				$seq .= $n;
				if ($data[6] eq '1')
				{
					$seq .= $db->seq($data[3], $data[5] => $data[4]);
				}
				else
				{
					$seq .= $db->seq($data[3], $data[4] => $data[5]);
				}
			}
			print BPDT "$key\t$regions{a}{$chr}{$key}\n";
			print FA ">$regions{a}{$chr}{$key}\n$seq\n" if ($seq);
		}
	}
	foreach my $chr1 (keys(%{$regions{e}}))
	{
		foreach my $chr2 (keys(%{$regions{e}{$chr1}}))
		{
			my $merge = 1;
			while ($merge)
			{
				$merge = 0;
MI:				foreach my $key1 (keys(%{$regions{e}{$chr1}{$chr2}}))
				{
					my @data1 = split (/__/, $regions{e}{$chr1}{$chr2}{$key1});
					foreach my $key2 (keys(%{$regions{e}{$chr1}{$chr2}}))
					{
						next if ($key1 eq $key2);
						my @data2 = split (/__/, $regions{e}{$chr1}{$chr2}{$key2});
						if (&covered($data1[1], $data1[2], $data2[1], $data2[2]) and &covered($data1[4], $data1[5], $data2[4], $data2[5]))
						{
							$merge = 1;
							my @temp1 = ($data1[1], $data1[2], $data2[1], $data2[2]);
							@temp1 = sort {$a <=> $b} @temp1;
							my @temp2 = ($data1[4], $data1[5], $data2[4], $data2[5]);
							@temp2 = sort {$a <=> $b} @temp2;
							delete $regions{e}{$chr1}{$chr2}{$key1};
							delete $regions{e}{$chr1}{$chr2}{$key2};
							my $newkey = $key1.'/'.$key2;
							if ($data1[6] eq $data2[6])
							{
								$regions{e}{$chr1}{$chr2}{$newkey} = $data1[0].'__'.$temp1[0].'__'.$temp1[-1].'__'.$data1[3].'__'.$temp2[0].'__'.$temp2[-1].'__'.$data1[6];
							}
							else
							{
								$regions{e}{$chr1}{$chr2}{$newkey} = $data1[0].'__'.$temp1[0].'__'.$temp1[-1].'__'.$data1[3].'__'.$temp2[0].'__'.$temp2[-1].'__10';
							}
							#print STDERR "@data1\n@data2\n$regions{e}{$chr1}{$chr2}{$newkey}\n";
							next MI;
						}
					}
				}
			}
			foreach my $key (keys(%{$regions{e}{$chr1}{$chr2}}))
			{
				my @data = split (/__/, $regions{e}{$chr1}{$chr2}{$key});
				next if ($data[1]<0 or $data[2]<0 or $data[4]<0 or $data[5]<0);
				my $seq = $db->seq($data[0], $data[1] => $data[2]);
				$seq .= $n;
				if ($data[6] eq '1')
				{
					$seq .= $db->seq($data[3], $data[5] => $data[4]);
				}
				else
				{
					$seq .= $db->seq($data[3], $data[4] => $data[5]);
				}
				print BPDT "$key\t$regions{e}{$chr1}{$chr2}{$key}\n";
				print FA ">$regions{e}{$chr1}{$chr2}{$key}\n$seq\n" if ($seq);
			}
		}
	}
	close BPDT;
	close FA;
	
	
	system "$bwa_command index $fafile";
	system "$bwa_command aln $fafile -l $cut_sr -k 1 -t $threads_bwa $sr_fq1 > $sr_sai1  2>>bwa.err";
	system "$bwa_command aln $fafile -l $cut_sr -k 1 -t $threads_bwa $sr_fq2 > $sr_sai2  2>>bwa.err";
	return if ($algs == 1);
ALGS1:	system "$bwa_command sampe -a $sr_insertsize -P -N 20 $fafile $sr_sai1 $sr_sai2 $sr_fq1 $sr_fq2 2>>bwa.err >$sr_rawsam";
	
#	adjust pair end mapping for mis-aligned reads
	my %alg_mis;
#	%alg_mis: mis-aligned reads
#	$alg_mis{readname}[0]: read id mis-aligned
#	$alg_mis{readname}[1]: read to be placed
	open RAWSAM, "<$sr_rawsam";
	while ($newline = <RAWSAM>)
	{
		chomp $newline;
		my @data = split (/\t/, $newline);
		$data[0] =~ /_(\d{2,3})$/;
		my $readlength = $1;
		if ($data[6] ne '=' and $data[11] and $newline =~ /XA:Z:/)
		{
			my $sr_alt = $';#'
			my @sr_alt = split (/;/, $sr_alt);
			my ($select_chr, $select_p, $trash);
			foreach (@sr_alt)
			{
				if ($_ =~ /$data[6]\S{0,100}/)
				{
					#print "$newline\n";
					my $current_match = $&;
					unless ($alg_mis{$data[0]}[1])
					{
						$alg_mis{$data[0]}[1] = $current_match;
						if ($data[1] & 0x40) # first in mate
						{
							$alg_mis{$data[0]}[0] = 1;
						}
						if ($data[1] & 0x80) # second in mate
						{
							$alg_mis{$data[0]}[0] = 2;
						}
						($select_chr, $select_p, $trash) = split (/,/, $current_match);
						$select_p = substr ($select_p, 1);
						next;
					}
					my @p = split (/__/, $select_chr);
					if ($p[5])
					{
						$p[3] = $p[4];
						$p[4] = $p[5];
						$p[5] = undef;
					}
					my $left_size = $p[2] - $p[1] + 1;
					$p[4] = $p[4] - $p[3] + $left_size + 100;
					$p[3] = $left_size + 100;
					$p[2] = $p[2] - $p[1];
					$p[1] = 1;
					my ($trash, $current_p, $trash) = split (/,/, $current_match);
					$current_p = substr ($current_p, 1);
					if (abs(abs($current_p-$data[7])-($readlength-$cut_sr))<2)
					{
						if ($data[1] & 0x40) # first in mate
						{
							$alg_mis{$data[0]}[0] = 1;
						}
						if ($data[1] & 0x80) # second in mate
						{
							$alg_mis{$data[0]}[0] = 2;
						}
						$alg_mis{$data[0]}[1] = $current_match;
						$select_p = $current_p;
						last;
					}
					if ($select_p < $left_size and $current_p < $left_size)
					{
						if (abs($current_p-$p[$bp_weight{$data[6]}[0]])<abs($select_p-$p[$bp_weight{$data[6]}[0]]))
						{
							if ($data[1] & 0x40) # first in mate
							{
								$alg_mis{$data[0]}[0] = 1;
							}
							if ($data[1] & 0x80) # second in mate
							{
								$alg_mis{$data[0]}[0] = 2;
							}
							$alg_mis{$data[0]}[1] = $current_match;
							$select_p = $current_p;
						}
					}
					if ($select_p > ($left_size + 100) and $current_p > ($left_size + 100))
					{
						if (abs($current_p-$p[$bp_weight{$data[6]}[1]])<abs($select_p-$p[$bp_weight{$data[6]}[1]]))
						{
							if ($data[1] & 0x40) # first in mate
							{
								$alg_mis{$data[0]}[0] = 1;
							}
							if ($data[1] & 0x80) # second in mate
							{
								$alg_mis{$data[0]}[0] = 2;
							}
							$alg_mis{$data[0]}[1] = $current_match;
							$select_p = $current_p;
						}
					}
				}
			}
			#print "$data[0]\t$data[1]\t$alg_mis{$data[0]}[0]\t$alg_mis{$data[0]}[1]\n";
		}
		elsif ($data[6] eq '=' and $data[11] and $newline =~ /XA:Z:/)
		{
			my $sr_alt = $';#'
			my @sr_alt = split (/;/, $sr_alt);
			my $select_p = $data[3];
			foreach (@sr_alt)
			{
				if ($_ =~ /$data[2]\S{0,100}/)
				{
					#print "$newline\n";
					my $current_match = $&;
					my @p = split (/__/, $data[2]);
					if ($p[5])
					{
						$p[3] = $p[4];
						$p[4] = $p[5];
						$p[5] = undef;
					}
					my $left_size = $p[2] - $p[1] + 1;
					$p[4] = $p[4] - $p[3] + $left_size + 100;
					$p[3] = $left_size + 100;
					$p[2] = $p[2] - $p[1];
					$p[1] = 1;
					my ($trash, $current_p, $trash) = split (/,/, $current_match);
					$current_p = substr ($current_p, 1);
					if (abs(abs($current_p-$data[7])-($readlength-$cut_sr))<2)
					{
						if ($data[1] & 0x40) # first in mate
						{
							$alg_mis{$data[0]}[0] = 1;
						}
						if ($data[1] & 0x80) # second in mate
						{
							$alg_mis{$data[0]}[0] = 2;
						}
						$alg_mis{$data[0]}[1] = $current_match;
						$select_p = $current_p;
						last;
					}
					if ($select_p < $left_size and $current_p < $left_size)
					{
						if (abs($current_p-$p[$bp_weight{$data[2]}[0]])<abs($select_p-$p[$bp_weight{$data[2]}[0]]))
						{
							if ($data[1] & 0x40) # first in mate
							{
								$alg_mis{$data[0]}[0] = 1;
							}
							if ($data[1] & 0x80) # second in mate
							{
								$alg_mis{$data[0]}[0] = 2;
							}
							$alg_mis{$data[0]}[1] = $current_match;
							$select_p = $current_p;
						}
					}
					if ($select_p > ($left_size + 100) and $current_p > ($left_size + 100))
					{
						if (abs($current_p-$p[$bp_weight{$data[2]}[1]])<abs($select_p-$p[$bp_weight{$data[2]}[1]]))
						{
							#print "$data[2]\t$data[0]\t$select_p\t$current_p\t@p\t$p[$bp_weight{$data[2]}[1]]\n";
							if ($data[1] & 0x40) # first in mate
							{
								$alg_mis{$data[0]}[0] = 1;
							}
							if ($data[1] & 0x80) # second in mate
							{
								$alg_mis{$data[0]}[0] = 2;
							}
							$alg_mis{$data[0]}[1] = $current_match;
							$select_p = $current_p;
						}
					}
				}
			}
			#print "$data[0]\t$data[1]\t$alg_mis{$data[0]}[0]\t$alg_mis{$data[0]}[1]\n";
		}
	}
	close RAWSAM;
	
	open RAWSAM, "<$sr_rawsam";
	open MDSAM, ">$sr_sam";
	while ($newline = <RAWSAM>)
	{
		chomp $newline;
		my @data = split (/\t/, $newline);
		my @data1 = @data;
		if ($alg_mis{$data[0]})
		{
			#print "$newline\n";
			# modify read
			if (($data[1] & 0x40 and $alg_mis{$data[0]}[0]==1) or ($data[1] & 0x80 and $alg_mis{$data[0]}[0]==2))
			{
				#print "$data[0]\t1\t$data[1]\t$alg_mis{$data[0]}[0]\n";
				if ($data[1] & 0x10 and $alg_mis{$data[0]}[1] =~ /\+/)
				{
					$data1[1] -= 16;
				}
				if (!($data[1] & 0x10) and $alg_mis{$data[0]}[1] =~ /-/)
				{
					$data1[1] += 16;
				}
				if ($alg_mis{$data[0]}[1] =~ /(\S{0,100}?),(\S{0,100}?),(\S{0,10}?),(\S{0,10})/)
				{
					#print "$1\t$2\t$3\t$4\n";
					$data1[6] = '=';
					$data1[2] = $1;
					$data1[3] = substr ($2, 1);
					$data1[8] = $data1[7] - $data1[3];
				}
				$data1[1] += 2 unless ($data[1] & 0x2);
				#print "$data1[0]\t$data1[1]\t$data1[6]\t$data1[7]\t$data1[8]\n";
			}
			# modify mate
			else
			{
				#print "$data[0]\t2\t$data[1]\t$alg_mis{$data[0]}[0]\n";
				if ($data[1] & 0x20 and $alg_mis{$data[0]}[1] =~ /\+/)
				{
					$data1[1] -= 32;
				}
				if (!($data[1] & 0x20) and $alg_mis{$data[0]}[1] =~ /-/)
				{
					$data1[1] += 32;
				}
				if ($alg_mis{$data[0]}[1] =~ /(\S{0,100}?),(\S{0,100}?),(\S{0,10}?),(\S{0,10})/)
				{
					#print "$data[0]\t$1\t$2\t$3\t$4\n";
					$data1[6] = '=';
					$data1[7] = substr ($2, 1);
					$data1[8] = $data1[7] - $data1[3];
				}
				$data1[1] -= 2 if ($data[1] & 0x2);
				#print "$data1[0]\t$data1[1]\t$data1[6]\t$data1[7]\t$data1[8]\n";
			}
			my $data = join ("\t", @data1);
			print MDSAM "$data\n";
		}
		else
		{
			print MDSAM "$newline\n";
		}
	}
	close RAWSAM;
	close MDSAM;
	
	system "$samtools_command faidx $fafile";
	system "$samtools_command view -bt $faifile $sr_sam -o $sr_bam";
	system "$samtools_command sort $sr_bam $sr_sort";
	system "$samtools_command index $sr_sortbam";
	#system "rm $sr_sai1";
	#system "rm $sr_sai2";
	#system "rm $sr_sam";
	system "rm $sr_rawsam";
	#system "rm $sr_bam";
	system "rm $fafile.*";
}

#	extract discordant read pairs
sub discord
{
	my $bamfile = shift;
	my $prefix = shift;
	my $blacklistfile = shift;
	my $clip = shift;
	my $dre_command = shift;
	my $sd_cutoff_disc = shift;
	my $isinfofile = shift;
	my $samtools_command = shift;
	my $remove_dup = shift;
	my $black_rg = shift;
	my $drelog = $prefix.'.dre.log';
	my $cl_bam = $prefix.'.cl.sorted.bam';
	my $disc_bam = $prefix.'.disc.bam';
	my $dup_file = $prefix.'.dup.bam';
	my $disc_cl_bam = $prefix.'.cl.disc.bam';
	my $dup_cl_file = $prefix.'.cl.dup.bam';
	
	die "bam file doesn't exist\n" unless (-e $bamfile);
	
	if ($blacklistfile)
	{
		if ($remove_dup)
		{
			system "$dre_command -f -v -P -s $sd_cutoff_disc -b $blacklistfile -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file $bamfile >$drelog";
			my $time = time;
			my $local = localtime($time);
			print STDERR "$local extracted discordant read pairs\n";
			if ($clip)
			{
				system "$dre_command -f -v -P -s $sd_cutoff_disc -b $blacklistfile -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file $cl_bam >>$drelog";
				my $time = time;
				my $local = localtime($time);
				print STDERR "$local extracted discordant read pairs from remapped bam file\n";
			}
		}
		else
		{
			system "$dre_command -f -v -m -s $sd_cutoff_disc -b $blacklistfile -r $black_rg -i $isinfofile -o $disc_bam $bamfile >$drelog";
			my $time = time;
			my $local = localtime($time);
			print STDERR "$local extracted discordant read pairs\n";
			if ($clip)
			{
				system "$dre_command -f -v -P -s $sd_cutoff_disc -b $blacklistfile -r $black_rg -i $isinfofile -o $disc_cl_bam $cl_bam >>$drelog";
				my $time = time;
				my $local = localtime($time);
				print STDERR "$local extracted discordant read pairs from remapped bam file\n";
			}
		}
	}
	else
	{
		if ($remove_dup)
		{
			system "$dre_command -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam -d $dup_file $bamfile >$drelog";
			my $time = time;
			my $local = localtime($time);
			print STDERR "$local extracted discordant read pairs\n";
			if ($clip)
			{
				system "$dre_command -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam -d $dup_cl_file $cl_bam >>$drelog";
				my $time = time;
				my $local = localtime($time);
				print STDERR "$local extracted discordant read pairs from remapped bam file\n";
			}
		}
		else
		{
			system "$dre_command -f -v -m -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_bam $bamfile >$drelog";
			my $time = time;
			my $local = localtime($time);
			print STDERR "$local extracted discordant read pairs\n";
			if ($clip)
			{
				system "$dre_command -f -v -P -s $sd_cutoff_disc -r $black_rg -i $isinfofile -o $disc_cl_bam $cl_bam >>$drelog";
				my $time = time;
				my $local = localtime($time);
				print STDERR "$local extracted discordant read pairs from remapped bam file\n";
			}
		}
	}
	my $disc_sort = $prefix.'.disc.sorted';
	my $disc_cl_sort = $prefix.'.cl.disc.sorted';
	
	system "$samtools_command sort -n $disc_bam $disc_sort";
	system "rm $disc_bam";
	my $time = time;
	my $local = localtime($time);
	print STDERR "$local sorted discordant bam by name\n";
	if ($clip)
	{
		system "$samtools_command sort -n $disc_cl_bam $disc_cl_sort";
		system "rm $disc_cl_bam";
		my $time = time;
		my $local = localtime($time);
		print STDERR "$local sorted discordant bam by name from remapped bam file\n";
	}
}

#	cluster discordant read pairs
sub cluster
{
#	prefix.discord file format: cluster id, number of read pairs, chr, position, strand, mchr, mposition, mstrand, size
#	prefix.clusters file format: cluster id primary, cluster id secondary, old cluster id, weight, mismatch, readname, read group, chr, strand, start, length, mchr, mstrand, mstart, mlength, isize

# hash keys
# n: readname
# c: seqid
# t: start
# e: end
# d: strand
# m: nm
# l: len
# C: mseqid
# T: mstart
# E: mend
# D: mstrad
# M: mnm
# L: mlen
# i: isize
# r: rg
# w: weight
# a: alternative mapping
	
	my $prefix = shift;
	my $blacklistfile = shift;
	my $ref_is = shift;
	my $cut_sr = shift;
	my $ad_align = shift;
	my $alt_map_max = shift;
	my $alt_map_max_clip = shift;
	my $clip = shift;
	my $samtools_command = shift;
	my $cl = shift;
	my $scluster_command = shift;
	my $disc_bam = $prefix.'.disc.sorted.bam';
	my $disc_cl_bam = $prefix.'.cl.disc.sorted.bam';
	my $raw_alt_map_file = $prefix.'.alt.map.raw';
	my $sort_alt_map_file = $prefix.'.alt.map.sort';
	my $raw_mapping = $prefix.'.mapping.raw';
	my $raw_cluster = $prefix.'.clusters.raw'; # for -a 0 only
	my $clusterfile = $prefix.'.clusters'; # final output
	goto AS2 if ($cl == 2 or $cl == 3);
	
	# parse positions to be ignored
	my %blacklist;
	if (-e $blacklistfile)
	{
		open(FILE, "gunzip -c $blacklistfile |");
		my $newline;
		while ($newline = <FILE>)
		{
			chomp $newline;
			my @data = split (/\t/, $newline);
			$blacklist{$data[0]}{$data[1]} = 1;
			#print "$data[0]\t$data[1]\n";
		}
		close FILE;
	}
	
	open(my $rawalt_fh, ">", "$raw_alt_map_file");
	my (@disc_alg, $last_readname);
	open (SAM, "$samtools_command view -X $disc_bam|");
	while ($newline = <SAM>)
	{
		chomp $newline;
		my @data = split (/\t/, $newline);
		if ($data[0] ne $last_readname)
		{
			&alt_map(\@disc_alg, $rawalt_fh, \%blacklist, $ref_is);
			@disc_alg = undef;
			$disc_alg[0] = $newline;
			$last_readname = $data[0];
		}
		else
		{
			push @disc_alg, $newline;
		}
	}
	close SAM;
	&alt_map(\@disc_alg, $rawalt_fh, \%blacklist, $ref_is);
	
	if ($clip)
	{
		open (CLSAM, "$samtools_command view -X $disc_cl_bam|");
		while ($newline = <CLSAM>)
		{
			chomp $newline;
			my @data = split (/\t/, $newline);
			if ($data[0] ne $last_readname)
			{
				&alt_map(\@disc_alg, $rawalt_fh, \%blacklist, $ref_is);
				@disc_alg = undef;
				$disc_alg[0] = $newline;
				$last_readname = $data[0];
			}
			else
			{
				push @disc_alg, $newline;
			}
		}
		close CLSAM;
		&alt_map(\@disc_alg, $rawalt_fh, \%blacklist, $ref_is);
	}
	close $rawalt_fh;
	
	my $time = time;
	my $local = localtime($time);
	print STDERR "$local constructed all mappings\n";
		
	system "sort -k 2,2 -k 3,3n $raw_alt_map_file > $sort_alt_map_file";
	system "rm $raw_alt_map_file";
	return if ($cl == 1);
	
#	construct clusters for all alternative mappings include uniq mapped reads
AS2:	goto AS3 if ($cl == 3);
	my %cluster_alt;
	my $i = 0;
#	%cluster_alt: cluster of all alternative mappings, same format as %cluster
	
	my $count_alt = 0;
	my $lastseqid;
	open MAPPING, ">$raw_mapping";
	open SORTALT, "<$sort_alt_map_file";
CLA4:	while (<SORTALT>)
	{
		chomp;
		my @data = split (/\t/, $_);
		my $readname = $data[0]; my $seqid = $data[1]; my $start = $data[2]; my $strand = $data[3]; my $mseqid = $data[4]; my $mstart = $data[5]; my $mstrand = $data[6]; my $isize = $data[7]; my $nm = $data[8]; my $mnm = $data[9]; my $rg = $data[10]; my $len = $data[11]; my $mlen = $data[12]; my $weight = $data[13]; 
		%cluster_alt = undef if ($seqid ne $lastseqid);
		$lastseqid = $seqid;
		my ($incluster, $topush);
CLA5:		for (my $i=$count_alt;$i>=0;$i--)
		{
			my $outofrange;
			foreach (@{$cluster_alt{$i}})
			{
				my $readnamet = $$_[0]; my $seqidt = $$_[1]; my $startt = $$_[2]; my $strandt = $$_[3]; my $mseqidt = $$_[4]; my $mstartt = $$_[5]; my $mstrandt = $$_[6]; my $isizet = $$_[7]; my $rgt = $$_[10]; 
				my $isize_cutoff_u = $$ref_is{$rg}{'median'}-$$ref_is{$rgt}{'median'}+$sd_cutoff_cl*sqrt($$ref_is{$rg}{'sd'}**2+$$ref_is{$rgt}{'sd'}**2);
				my $isize_cutoff_d = $$ref_is{$rg}{'median'}-$$ref_is{$rgt}{'median'}-$sd_cutoff_cl*sqrt($$ref_is{$rg}{'sd'}**2+$$ref_is{$rgt}{'sd'}**2);
				my $map_cutoff_u = ($$ref_is{$rgt}{'isu'} > $$ref_is{$rg}{'isu'})?$$ref_is{$rgt}{'isu'}:$$ref_is{$rg}{'isu'};
				if ($seqid eq $seqidt)
				{
					if ($readnamet eq $readname and $mseqid eq $mseqidt and $startt == $start and $mstartt == $mstart and $strand == $strandt and $mstrand == $mstrandt)
					{
						next CLA4;
					}
					elsif ($mseqid eq $mseqidt and $strand == $strandt and $mstrand == $mstrandt and abs($start-$startt)<=$map_cutoff_u and abs($mstart-$mstartt)<=$map_cutoff_u)
					{
						if ($strand == $mstrand)
						{
							$topush = $i;
							$incluster = 1;
							last CLA5;
						}
						else
						{
							if ($isize-$isizet<=$isize_cutoff_u and $isize-$isizet>=$isize_cutoff_d)
							{
								$topush = $i;
								$incluster = 1;
								last CLA5;
							}
						}
					}
					if (abs($start-$startt)>$$ref_is{'isu'})
					{
						$outofrange = 1;
					}
				}
				else
				{
					last CLA5;
				}
			}
			last if ($outofrange);
		}
		unless ($incluster)
		{
			$count_alt++;
			$cluster_alt{$count_alt}[0][0] = $readname;
			$cluster_alt{$count_alt}[0][1] = $seqid;
			$cluster_alt{$count_alt}[0][2] = $start;
			$cluster_alt{$count_alt}[0][3] = $strand;
			$cluster_alt{$count_alt}[0][4] = $mseqid;
			$cluster_alt{$count_alt}[0][5] = $mstart;
			$cluster_alt{$count_alt}[0][6] = $mstrand;
			$cluster_alt{$count_alt}[0][7] = $isize;
			$cluster_alt{$count_alt}[0][8] = $nm;
			$cluster_alt{$count_alt}[0][9] = $mnm;
			$cluster_alt{$count_alt}[0][10] = $rg;
			$cluster_alt{$count_alt}[0][11] = $len;
			$cluster_alt{$count_alt}[0][12] = $mlen;
			$cluster_alt{$count_alt}[0][13] = $weight;
			print MAPPING "$count_alt\t$readname\t$rg\t$seqid\t$strand\t$start\t$len\t$mseqid\t$mstrand\t$mstart\t$mlen\t$isize\t$nm\t$mnm\t$weight\n";
		}
		if (defined($topush))
		{
			my $k = @{$cluster_alt{$topush}};
			$cluster_alt{$topush}[$k][0] = $readname;
			$cluster_alt{$topush}[$k][1] = $seqid;
			$cluster_alt{$topush}[$k][2] = $start;
			$cluster_alt{$topush}[$k][3] = $strand;
			$cluster_alt{$topush}[$k][4] = $mseqid;
			$cluster_alt{$topush}[$k][5] = $mstart;
			$cluster_alt{$topush}[$k][6] = $mstrand;
			$cluster_alt{$topush}[$k][7] = $isize;
			$cluster_alt{$topush}[$k][8] = $nm;
			$cluster_alt{$topush}[$k][9] = $mnm;
			$cluster_alt{$topush}[$k][10] = $rg;
			$cluster_alt{$topush}[$k][11] = $len;
			$cluster_alt{$topush}[$k][12] = $mlen;
			$cluster_alt{$topush}[$k][13] = $weight;
			print MAPPING "$topush\t$readname\t$rg\t$seqid\t$strand\t$start\t$len\t$mseqid\t$mstrand\t$mstart\t$mlen\t$isize\t$nm\t$mnm\t$weight\n";
		}
	}
	close MAPPING;
	close SORTALT;
	system "rm $sort_alt_map_file";
	my $time = time;
	my $local = localtime($time);
	print STDERR "$local constructed discordant clusters for all alternative mappings\n";
	return if ($cl == 2);
	
AS3:	if ($ad_align)
	{
		system "$scluster_command $raw_mapping > $clusterfile";
		#system "rm $raw_mapping";
		my $time = time;
		my $local = localtime($time);
		print STDERR "$local called discordant clusters\n";
	}
	else
	{
		system "sort -k 1,1n $raw_mapping > $raw_cluster";
		open CLFILE, ">$clusterfile";
		open RAWCL, "<$raw_cluster";
		while ($newline = <RAWCL>)
		{
			chomp $newline;
			my @data = split (/\t/, $newline);
			print CLFILE "$data[0]\t0\t\t\t\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t$data[7]\t$data[8]\t$data[9]\t$data[10]\t$data[11]\n";
		}
		close CLFILE;
		system "rm $raw_cluster";
		system "rm $raw_mapping";
		my $time = time;
		my $local = localtime($time);
		print STDERR "$local called discordant clusters\n";
	}
}

sub alt_map
{
	my $ref_disc_alg = shift;
	my $rawalt_fh = shift;
	my $ref_blacklist = shift;
	my $ref_is = shift;
	
	my @data1 = split (/\t/, $$ref_disc_alg[0]);
	my @data2 = split (/\t/, $$ref_disc_alg[1]);
	my @data3 = split (/\t/, $$ref_disc_alg[2]);
	my @data4 = split (/\t/, $$ref_disc_alg[3]);
	
	if ($data3[0])
	{
		# 4 entries per pair
		if ($data4[0])
		{
			#print "$$ref_disc_alg[0]\n$$ref_disc_alg[1]\n$$ref_disc_alg[2]\n$$ref_disc_alg[3]\n\n";
			return 0;
			if (substr($data1[1], -1) != substr($data2[1], -1) and substr($data1[1], -1) != substr($data3[1], -1) and substr($data1[1], -1) != substr($data4[1], -1))
			{
				
			}
			if (substr($data2[1], -1) != substr($data1[1], -1) and substr($data2[1], -1) != substr($data3[1], -1) and substr($data2[1], -1) != substr($data4[1], -1))
			{
				# no such pair
			}
			if (substr($data3[1], -1) != substr($data1[1], -1) and substr($data3[1], -1) != substr($data2[1], -1) and substr($data3[1], -1) != substr($data4[1], -1))
			{
				# no such pair
			}
			if (substr($data4[1], -1) != substr($data1[1], -1) and substr($data4[1], -1) != substr($data2[1], -1) and substr($data4[1], -1) != substr($data3[1], -1))
			{
			#print "$data1[0]\t$data1[1]\t$data2[1]\t$data3[1]\t$data4[1]\n";
				
			}
		}
		# 3 entries per pair
		else
		{
			# entry 2 and 3 from the same read
			if (substr($data1[1], -1) != substr($data2[1], -1) and substr($data1[1], -1) != substr($data3[1], -1))
			{
				#print "$$ref_disc_alg[0]\n$$ref_disc_alg[1]\n$$ref_disc_alg[2]\n\n";
				my (@cigar2, @cigar3, $maxmatch_value2, $maxmatch_i2, $maxmatch_value3, $maxmatch_i3);
				my $i = 0;
				while ($data2[5] =~ /(\d{1,5})(\D)/g)
				{
					$cigar2[$i][0] = $1;
					$cigar2[$i][1] = $2;
					$i++;
					#print "$data2[5]\t$i\t$cigar2[$i][0]\t$cigar2[$i][1]\n";
				}
				for ($i=0;$cigar2[$i][1];$i++)
				{
					#print "$data2[5]\t$i\t$cigar2[$i][0]\t$cigar2[$i][1]\n";
					if ($cigar2[$i][1] eq 'M' and $cigar2[$i][0] > $maxmatch_value2)
					{
						$maxmatch_value2 = $cigar2[$i][0];
						$maxmatch_i2 = $i;
					}
				}
				$i = 0;
				while ($data3[5] =~ /(\d{1,5})(\D)/g)
				{
					$cigar3[$i][0] = $1;
					$cigar3[$i][1] = $2;
					$i++;
					#print "$data3[5]\t$i\t$cigar3[$i][0]\t$cigar3[$i][1]\n";
				}
				for ($i=0;$cigar3[$i][1];$i++)
				{
					#print "$data3[5]\t$i\t$cigar3[$i][0]\t$cigar3[$i][1]\n";
					if ($cigar3[$i][1] eq 'M' and $cigar3[$i][0] > $maxmatch_value3)
					{
						$maxmatch_value3 = $cigar3[$i][0];
						$maxmatch_i3 = $i;
					}
				}
				
				my ($matchup2, $matchup3, $matchdown2, $matchdown3);
				for ($i=0;$cigar2[$i][1];$i++)
				{
					if ($cigar2[$i][1] eq 'M' and $cigar2[$i][0] == $maxmatch_value2)
					{
						last;
					}
					$matchup2 += $cigar2[$i][0];
				}
				for ($i=0;$cigar3[$i][1];$i++)
				{
					if ($cigar3[$i][1] eq 'M' and $cigar3[$i][0] == $maxmatch_value3)
					{
						last;
					}
					$matchup3 += $cigar3[$i][0];
				}
				my $tocount = 0;
				for ($i=0;$cigar2[$i][1];$i++)
				{
					if ($cigar2[$i][1] eq 'M' and $cigar2[$i][0] == $maxmatch_value2)
					{
						$tocount = 1;
						next;
					}
					$matchdown2 += $cigar2[$i][0] if ($tocount == 1);
				}
				$tocount = 0;
				for ($i=0;$cigar3[$i][1];$i++)
				{
					if ($cigar3[$i][1] eq 'M' and $cigar3[$i][0] == $maxmatch_value3)
					{
						$tocount = 1;
						next;
					}
					$matchdown3 += $cigar3[$i][0] if ($tocount == 1);
				}
				#print "$data1[0]\t$data2[1]\t$data2[5]\t$maxmatch_value2\t$maxmatch_i2\n$data1[0]\t$data3[1]\t$data3[5]\t$maxmatch_value3\t$maxmatch_i3\n";
					
				if ($data2[1] =~ /r/ and $data3[1] =~ /r/)
				{
					#print "$matchup2\t$matchup3\n";
					if ($matchup2 < $matchup3)
					{
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					goto ALTMAP;
				}
				if ($data2[1] =~ /r/ and !($data3[1] =~ /r/))
				{
					#print "$matchup2\t$matchdown3\n";
					if ($matchup2 < $matchdown3)
					{
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					goto ALTMAP;
				}
				if (!($data2[1] =~ /r/) and $data3[1] =~ /r/)
				{
					#print "$matchdown2\t$matchup3\n";
					if ($matchdown2 < $matchup3)
					{
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					goto ALTMAP;
				}
				if (!($data2[1] =~ /r/) and !($data3[1] =~ /r/))
				{
					#print "$matchdown2\t$matchdown3\n";
					if ($matchdown2 < $matchdown3)
					{
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					goto ALTMAP;
				}
			}
			if (substr($data2[1], -1) != substr($data1[1], -1) and substr($data2[1], -1) != substr($data3[1], -1))
			{
				# no such pair
			}
			# entry 1 and 2 from the same read
			if (substr($data3[1], -1) != substr($data1[1], -1) and substr($data3[1], -1) != substr($data2[1], -1))
			{
				#print "$$ref_disc_alg[0]\n$$ref_disc_alg[1]\n$$ref_disc_alg[2]\n\n";
				my (@cigar1, @cigar2, $maxmatch_value1, $maxmatch_i1, $maxmatch_value2, $maxmatch_i2);
				my $i = 0;
				while ($data1[5] =~ /(\d{1,5})(\D)/g)
				{
					$cigar1[$i][0] = $1;
					$cigar1[$i][1] = $2;
					$i++;
					#print "$data1[5]\t$i\t$cigar1[$i][0]\t$cigar1[$i][1]\n";
				}
				for ($i=0;$cigar1[$i][1];$i++)
				{
					#print "$data1[5]\t$i\t$cigar1[$i][0]\t$cigar1[$i][1]\n";
					if ($cigar1[$i][1] eq 'M' and $cigar1[$i][0] > $maxmatch_value1)
					{
						$maxmatch_value1 = $cigar1[$i][0];
						$maxmatch_i1 = $i;
					}
				}
				$i = 0;
				while ($data2[5] =~ /(\d{1,5})(\D)/g)
				{
					$cigar2[$i][0] = $1;
					$cigar2[$i][1] = $2;
					$i++;
					#print "$data2[5]\t$i\t$cigar2[$i][0]\t$cigar2[$i][1]\n";
				}
				for ($i=0;$cigar2[$i][1];$i++)
				{
					#print "$data2[5]\t$i\t$cigar2[$i][0]\t$cigar2[$i][1]\n";
					if ($cigar2[$i][1] eq 'M' and $cigar2[$i][0] > $maxmatch_value2)
					{
						$maxmatch_value2 = $cigar2[$i][0];
						$maxmatch_i2 = $i;
					}
				}
				
				my ($matchup1, $matchup2, $matchdown1, $matchdown2);
				for ($i=0;$cigar1[$i][1];$i++)
				{
					if ($cigar1[$i][1] eq 'M' and $cigar1[$i][0] == $maxmatch_value1)
					{
						last;
					}
					$matchup1 += $cigar1[$i][0];
				}
				for ($i=0;$cigar2[$i][1];$i++)
				{
					if ($cigar2[$i][1] eq 'M' and $cigar2[$i][0] == $maxmatch_value2)
					{
						last;
					}
					$matchup2 += $cigar2[$i][0];
				}
				my $tocount = 0;
				for ($i=0;$cigar1[$i][1];$i++)
				{
					if ($cigar1[$i][1] eq 'M' and $cigar1[$i][0] == $maxmatch_value1)
					{
						$tocount = 1;
						next;
					}
					$matchdown1 += $cigar1[$i][0] if ($tocount == 1);
				}
				$tocount = 0;
				for ($i=0;$cigar2[$i][1];$i++)
				{
					if ($cigar2[$i][1] eq 'M' and $cigar2[$i][0] == $maxmatch_value2)
					{
						$tocount = 1;
						next;
					}
					$matchdown2 += $cigar2[$i][0] if ($tocount == 1);
				}
				#print "$data1[0]\t$data1[1]\t$data1[5]\t$maxmatch_value1\t$maxmatch_i1\n$data1[0]\t$data2[1]\t$data2[5]\t$maxmatch_value2\t$maxmatch_i2\n";
					
				if ($data1[1] =~ /r/ and $data2[1] =~ /r/)
				{
					#print "$matchup1\t$matchup2\n";
					if ($matchup1 < $matchup2)
					{
						@data1 = @data2;
						$$ref_disc_alg[0] = $$ref_disc_alg[1];
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					else
					{
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					goto ALTMAP;
				}
				if ($data1[1] =~ /r/ and !($data2[1] =~ /r/))
				{
					#print "$matchup1\t$matchdown2\n";
					if ($matchup1 < $matchdown2)
					{
						@data1 = @data2;
						$$ref_disc_alg[0] = $$ref_disc_alg[1];
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					else
					{
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					goto ALTMAP;
				}
				if (!($data1[1] =~ /r/) and $data2[1] =~ /r/)
				{
					#print "$matchdown1\t$matchup2\n";
					if ($matchdown1 < $matchup2)
					{
						@data1 = @data2;
						$$ref_disc_alg[0] = $$ref_disc_alg[1];
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					else
					{
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					goto ALTMAP;
				}
				if (!($data1[1] =~ /r/) and !($data2[1] =~ /r/))
				{
					#print "$matchdown1\t$matchdown2\n";
					if ($matchdown1 < $matchdown2)
					{
						@data1 = @data2;
						$$ref_disc_alg[0] = $$ref_disc_alg[1];
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					else
					{
						@data2 = @data3;
						$$ref_disc_alg[1] = $$ref_disc_alg[2];
					}
					goto ALTMAP;
				}
			}
ALTMAP:			$data1[0] .= 'sc'; $data2[0] .= 'sc';
			#print "$$ref_disc_alg[0]\n$$ref_disc_alg[1]\n$$ref_disc_alg[2]\n\n\n";
		}
	}
	
	return 0 unless (defined($data2[0])); # skip unpaired reads
	if (($data1[2] gt $data2[2]) or ($data1[2] eq $data2[2] and $data1[3] > $data2[3]))
	{
		($$ref_disc_alg[0], $$ref_disc_alg[1]) = ($$ref_disc_alg[1], $$ref_disc_alg[0]);
		my @temp = @data2;
		@data2 = @data1;
		@data1 = @temp;
	}
	if ($min_mapq)
	{
		return 0 if ($data1[4] < $min_mapq or $data2[4] < $min_mapq);
	}
	my $readname = $data1[0];
	# a read is non-uniq mapped if XT tag with R, or has XA tag or mapping qual is 0
	my $xt = 'U';
	my $mxt = 'U';
	if ($$ref_disc_alg[0] =~ /XT:A:(\w)/)
	{
		$xt = $1;
	}
	if ($$ref_disc_alg[0] =~ /XA:Z:(\w)/)
	{
		$xt = 'R';
	}
	if ($$ref_disc_alg[1] =~ /XT:A:(\w)/)
	{
		$mxt = $1;
	}
	if ($$ref_disc_alg[1] =~ /XA:Z:(\w)/)
	{
		$mxt = 'R';
	}
	$xt = 'R' unless ($data1[4]);
	$mxt = 'R' unless ($data2[4]);
	if ($ad_align)
	{
		return 0 if ($xt eq 'R' and $mxt eq 'R');
		return 0 if ($xt eq 'M' or $mxt eq 'M');
	}
	else
	{
		unless ($use_all_align)
		{
			return 0 unless ($xt eq 'U' and $mxt eq 'U');
		}
	}
	my $seqid = $data1[2];
	my $mseqid = $data2[2];
	my $start = $data1[3];
	my $mstart = $data2[3];
	return 0 if ($$ref_blacklist{$seqid}{$start} or $$ref_blacklist{$mseqid}{$mstart} or $$ref_blacklist{$seqid}{$start+2*$cut_sr} or $$ref_blacklist{$mseqid}{$mstart+2*$cut_sr});
	my $strand = 1;
	$strand = -1 if ($data1[1] =~ /r/);
	my $mstrand = 1;
	$mstrand = -1 if ($data2[1] =~ /r/);
	my ($len, $mlen);
	while ($data1[5] =~ /(\d{1,5})(\D)/g)
	{
		if ($2 eq 'M' and $1 > $len)
		{
			$len = $1;
		}
	}
	while ($data2[5] =~ /(\d{1,5})(\D)/g)
	{
		if ($2 eq 'M' and $1 > $mlen)
		{
			$mlen = $1;
		}
	}
	
	my $isize;
	if ($data1[2] eq $data2[2])
	{
		if ($strand == $mstrand)
		{
			$isize = abs($data1[3] - $data2[3]) + 1;
		}
		elsif ($strand == 1 and $mstrand == -1)
		{
			$isize = abs($data1[3] - $data2[3]) + $mlen;
		}
		elsif ($strand == -1 and $mstrand == 1)
		{
			$isize = abs($data1[3] - $data2[3]) - $mlen + 2;
		}
	}
	my $rg = 'none';
	if ($$ref_disc_alg[0] =~ /RG:Z:(\S{1,500})/)
	{
		$rg = $1;
	}
	return 0 if ($seqid eq $mseqid and $strand == 1 and $mstrand == -1 and $isize <= $$ref_is{$rg}{'isu'});
	
	my ($xa, $mxa, $nm, $mnm, $tr, $mtr, @alt_mapping, @malt_mapping, $alt_mapping, $malt_mapping);
	if ($ad_align)
	{
		$$ref_disc_alg[0] =~ /NM:i:(\d{1,3})/;
		$nm = $1;
		$$ref_disc_alg[1] =~ /NM:i:(\d{1,3})/;
		$mnm = $1;
		if ($xt eq 'R')
		{
			if ($$ref_disc_alg[0] =~ /XA:Z:(\S{1,10000})/)
			{
				$xa = $1;
			}
			$tr = 1 unless ($xa);
			@alt_mapping = split (/;/, $xa);
			my $n_alt_mapping = @alt_mapping;
			if (defined($alt_map_max) and $n_alt_mapping>$alt_map_max)
			{
				$tr = 1;
			}
		}
		if ($mxt eq 'R')
		{
			if ($$ref_disc_alg[1] =~ /XA:Z:(\S{1,10000})/)
			{
				$mxa = $1;
			}
			$mtr = 1 unless ($mxa);
			@malt_mapping = split (/;/, $mxa);
			my $n_malt_mapping = @malt_mapping;
			if (defined($alt_map_max) and $n_malt_mapping>$alt_map_max)
			{
				$tr = 1;
			}
		}
		if ($xt eq 'R' or $mxt eq 'R')
		{
			return 0 if ($tr == 1 or $mtr == 1);
		}
		$alt_mapping = @alt_mapping;
		$malt_mapping = @malt_mapping;
		my $weight = 1;
		if ($xt eq 'R' or $mxt eq 'R')
		{
			$weight = ($alt_mapping>$malt_mapping)?$alt_mapping+1:$malt_mapping+1;
			$weight = sprintf("%.9f", 1/$weight);
		}
		my @toprint;
		push @toprint, "$readname\t$seqid\t$start\t$strand\t$mseqid\t$mstart\t$mstrand\t$isize\t$nm\t$mnm\t$rg\t$len\t$mlen\t$weight\n";
		foreach (@alt_mapping)
		{
			my ($seqid_alt, $start_alt, $mseqid_alt, $mstart_alt, $isize_alt, $strand_alt, $mstrand_alt, $nm_alt, $mnm_alt, $len_alt, $mlen_alt, $mapping, $mmapping);
			($seqid_alt, $start_alt, $mapping, $nm_alt) = split (/,/, $_);
			my $nosign_start_alt = substr ($start_alt, 1);
			return 0 if ($$ref_blacklist{$seqid_alt}{$nosign_start_alt} or $$ref_blacklist{$seqid_alt}{$nosign_start_alt+2*$cut_sr});
			($mseqid_alt, $mstart_alt, $mstrand_alt, $mnm_alt, $len_alt, $mlen_alt) = ($mseqid, $mstart, $mstrand, $mnm, $len, $mlen);
			next if ($seqid eq $seqid_alt and $start == $start_alt);
			$strand_alt = 1;
			$strand_alt = -1 if ($start_alt =~ /-/);
			$start_alt = $nosign_start_alt;
			$isize_alt = 0;
			if ($seqid_alt eq $mseqid_alt)
			{
				if ($start_alt <= $mstart_alt)
				{
					if ($strand_alt == $mstrand_alt)
					{
						$isize_alt = abs($mstart_alt - $start_alt) + 1;
					}
					else
					{
						if ($strand_alt == -1 and $mstrand_alt == 1)
						{
							$isize_alt = abs($mstart_alt - $start_alt) - $mlen_alt + 2;
						}
						else
						{
							$isize_alt = abs($mstart_alt - $start_alt) + $mlen_alt;
						}
					}
				}
				else
				{
					if ($strand_alt == $mstrand_alt)
					{
						$isize_alt = abs($start_alt - $mstart_alt) + 1;
					}
					else
					{
						if ($strand_alt == -1 and $mstrand_alt == 1)
						{
							$isize_alt = abs($start_alt - $mstart_alt) + $mlen_alt;
						}
						else
						{
							$isize_alt = abs($start_alt - $mstart_alt) - $mlen_alt + 2;
						}
					}
					($seqid_alt, $start_alt, $strand_alt, $nm_alt, $len_alt, $mseqid_alt, $mstart_alt, $mstrand_alt, $mnm_alt, $mlen_alt) = ($mseqid_alt, $mstart_alt, $mstrand_alt, $mnm_alt, $mlen_alt, $seqid_alt, $start_alt, $strand_alt, $nm_alt, $len_alt);
				}
			}
			if ($seqid_alt gt $mseqid_alt)
			{
				($seqid_alt, $start_alt, $strand_alt, $nm_alt, $len_alt, $mseqid_alt, $mstart_alt, $mstrand_alt, $mnm_alt, $mlen_alt) = ($mseqid_alt, $mstart_alt, $mstrand_alt, $mnm_alt, $mlen_alt, $seqid_alt, $start_alt, $strand_alt, $nm_alt, $len_alt);
			}
			if ($seqid_alt eq $mseqid_alt and $strand_alt == 1 and $mstrand_alt == -1 and $isize_alt <= $$ref_is{$rg}{'isu'})
			{
				return 0;
			}
			push @toprint, "$readname\t$seqid_alt\t$start_alt\t$strand_alt\t$mseqid_alt\t$mstart_alt\t$mstrand_alt\t$isize_alt\t$nm_alt\t$mnm_alt\t$rg\t$len_alt\t$mlen_alt\t$weight\n";
		}
		foreach (@malt_mapping)
		{
			my ($seqid_alt, $start_alt, $mseqid_alt, $mstart_alt, $isize_alt, $strand_alt, $mstrand_alt, $nm_alt, $mnm_alt, $len_alt, $mlen_alt, $mapping, $mmapping);
			($mseqid_alt, $mstart_alt, $mmapping, $mnm_alt) = split (/,/, $_);
			my $nosign_mstart_alt = substr ($mstart_alt, 1);
			return 0 if ($$ref_blacklist{$mseqid_alt}{$nosign_mstart_alt} or $$ref_blacklist{$mseqid_alt}{$nosign_mstart_alt+2*$cut_sr});
			($seqid_alt, $start_alt, $strand_alt, $nm_alt, $len_alt, $mlen_alt) = ($seqid, $start, $strand, $nm, $len, $mlen);
			next if ($mseqid eq $mseqid_alt and $mstart == $mstart_alt);
			$mstrand_alt = 1;
			$mstrand_alt = -1 if ($mstart_alt =~ /-/);
			$mstart_alt = $nosign_mstart_alt;
			$isize_alt = 0;
			if ($seqid_alt eq $mseqid_alt)
			{
				if ($start_alt <= $mstart_alt)
				{
					if ($strand_alt == $mstrand_alt)
					{
						$isize_alt = abs($mstart_alt - $start_alt) + 1;
					}
					else
					{
						if ($strand_alt == -1 and $mstrand_alt == 1)
						{
							$isize_alt = abs($mstart_alt - $start_alt) - $mlen_alt + 2;
						}
						else
						{
							$isize_alt = abs($mstart_alt - $start_alt) + $mlen_alt;
						}
					}
				}
				else
				{
					if ($strand_alt == $mstrand_alt)
					{
						$isize_alt = abs($start_alt - $mstart_alt) + 1;
					}
					else
					{
						if ($strand_alt == -1 and $mstrand_alt == 1)
						{
							$isize_alt = abs($start_alt - $mstart_alt) + $mlen_alt;
						}
						else
						{
							$isize_alt = abs($start_alt - $mstart_alt) - $mlen_alt + 2;
						}
					}
					($seqid_alt, $start_alt, $strand_alt, $nm_alt, $len_alt, $mseqid_alt, $mstart_alt, $mstrand_alt, $mnm_alt, $mlen_alt) = ($mseqid_alt, $mstart_alt, $mstrand_alt, $mnm_alt, $mlen_alt, $seqid_alt, $start_alt, $strand_alt, $nm_alt, $len_alt);
				}
			}
			if ($seqid_alt gt $mseqid_alt)
			{
				($seqid_alt, $start_alt, $strand_alt, $nm_alt, $len_alt, $mseqid_alt, $mstart_alt, $mstrand_alt, $mnm_alt, $mlen_alt) = ($mseqid_alt, $mstart_alt, $mstrand_alt, $mnm_alt, $mlen_alt, $seqid_alt, $start_alt, $strand_alt, $nm_alt, $len_alt);
			}
			if ($seqid_alt eq $mseqid_alt and $strand_alt == 1 and $mstrand_alt == -1 and $isize_alt <= $$ref_is{$rg}{'isu'})
			{
				return 0;
			}
			push @toprint, "$readname\t$seqid_alt\t$start_alt\t$strand_alt\t$mseqid_alt\t$mstart_alt\t$mstrand_alt\t$isize_alt\t$nm_alt\t$mnm_alt\t$rg\t$len_alt\t$mlen_alt\t$weight\n";
		}
		foreach (@toprint)
		{
			print $rawalt_fh "$_";
		}
	}
	else
	{
		print $rawalt_fh "$readname\t$seqid\t$start\t$strand\t$mseqid\t$mstart\t$mstrand\t$isize\t$nm\t$mnm\t$rg\t$len\t$mlen\n";
	}
}

#	identify candidates from read pairs
sub mpd
{
#	file format of prefix.mp.intra.out
#	del: simple deletion
#	del_inssd: deletion with insertion in the same orientation, insertion come from downstream of deletion
#	del_inssu: deletion with insertion in the same orientation, insertion come from upstream of deletion
#	del_insod: deletion with insertion in the opposite orientation, insertion come from downstream of deletion
#	del_insou: deletion with insertion in the opposite orientation, insertion come from upstream of deletion
#	invers: inversion supported by 2 reciprocal clusters, include simple inversion and deletion with inversion
#	tandem_dup: tandem duplication
#	invers_f: inversion supported by one cluster on forward strand
#	invers_r: inversion supported by one cluster on reverse strand
#	del, cluster id, number of supporting read pairs, chr, range of deletion (2 col), deletion size, region of break point
#	del_ins*, cluster id, number of supporting read pairs, chr, range of deletion (2 col), deletion size, chr (donor), range of insertion (2 col), insert size, distance of deletion and insertion, region of break point (2 col)
#	*invers, cluster id, number of supporting read pairs, chr, inversion left boundary (2 col), inversion size, inversion right boundary (2 col), region of break point (2 col)
#	invers_*, cluster id, number of supporting read pairs, chr, inversion left boundary, inversion right boundary, inversion size, region of break point
#	tandem_dup, cluster id, number of supporting read pairs, chr, tandem duplication boundary 1, tandem duplication boundary 2, tandem duplication size, region of break point

#	file format of prefix.mp.inter.out
#	del_inss: deletion with insertion in the same orientation
#	del_inso: deletion with insertion in the opposite orientation
#	transl_inter: inter-chromosomal translocation
#	del_ins*, cluster id, number of supporting read pairs, chr of deletion, range of deletion (2 col), deletion size, chr of insertion donor, range of insertion (2 col), insert size, region of break point (2 col)
#	transl_inter, cluster id, number of supporting read pairs, chr of 1st cluster, range (2 col), orientation of 1st cluster, chr of 2nd cluster, range (2 col), orientation of 2nd cluster, region of break point
	my $ref_is = shift;
	my $prefix = shift;
	my $support_mps = shift;
	my $support_mpf = shift;
	my $clusterfile = $prefix.'.clusters';
	my $discclfile = $prefix.'.discord';
	my $mpintra_outfile = $prefix.'.mp.intra.out';
	my $mpinter_outfile = $prefix.'.mp.inter.out';
	
	my (@cluster, %cluster_type, $count, %readname, %support, %support_f, $lastpid, $lastsid);
	my (%start, %mstart, %orientation, %cbp);
#	@cluster, detail of each cluster
#	$cluster[k][j]
#	%cluster_type: index of cluster type
#	$cluster_type{chr1}{chr2}{0-3}: array of cluster id, primary_secondary
#	0: 1 -1; 1: 1 1; 2: -1 -1; 3: -1 1
	my $k = 0; # number of supporting entries
	my $l = 0; # number of supporting read pairs
	my $j = 0; # number of supporting full length read pairs
	open DISCCL, ">$discclfile";
	open CLMPD, "<$clusterfile";
	while ($newline = <CLMPD>)
	{
		chomp $newline;
		my @data = split (/\t/, $newline);
		if ($lastpid ne $data[0] or $lastsid ne $data[1])
		{
			# last cluster
			$support{$lastpid} = $l;
			$support_f{$lastpid} = $j;
			#if ($support{$lastpid} >= $support_mps)
			{
				my (@start, @mstart, @isize, @cbp);
				foreach (@cluster)
				{
					my $readname = $$_[0]; my $rg = $$_[1]; my $seqid = $$_[2]; my $strand = $$_[3]; my $start = $$_[4]; my $end = $$_[5]; my $len = $$_[6]; my $mseqid = $$_[7]; my $mstrand = $$_[8]; my $mstart = $$_[9]; my $mend = $$_[10]; my $mlen = $$_[11]; my $isize = $$_[12];
					next unless ($seqid);
					push @start, $start if ($start);
					push @start, $end if ($end);
					push @mstart, $mstart if ($mstart);
					push @mstart, $mend if ($mend);
					push @isize, $isize;
					$end = $start unless ($end);
					$mend = $mstart unless ($mend);
					if ($strand == 1)
					{
						($cbp[0], $cbp[1]) = &overlap($end-10, $start+$$ref_is{$rg}{'isu'}-$mlen, $cbp[0], $cbp[1]);
						# -10 or +10: expand the window to avoid one read too close to a break point
					}
					else
					{
						($cbp[0], $cbp[1]) = &overlap($end-$$ref_is{$rg}{'isu'}+$mlen, $start+10, $cbp[0], $cbp[1]);
					}
					if ($mstrand == 1)
					{
						($cbp[2], $cbp[3]) = &overlap($mend-10, $mstart+$$ref_is{$rg}{'isu'}-$len, $cbp[2], $cbp[3]);
					}
					else
					{
						($cbp[2], $cbp[3]) = &overlap($mend-$$ref_is{$rg}{'isu'}+$len, $mstart+10, $cbp[2], $cbp[3]);
					}
					#print STDERR "$lastpid $readname\t$end-$$ref_is{$rg}{'isu'}+$mlen, $start+10, $cbp[0]\t$cbp[1]\t$cbp[2]\t$cbp[3]\n";
					$orientation{$lastpid}{$lastsid}[0] = $strand;
					$orientation{$lastpid}{$lastsid}[1] = $mstrand;
					$orientation{$lastpid}{$lastsid}[2] = $seqid;
					$orientation{$lastpid}{$lastsid}[3] = $mseqid;
				}
				@start = sort {$a <=> $b} @start;
				@mstart = sort {$a <=> $b} @mstart;
				@isize = sort {$a <=> $b} @isize;
				#print "$lastpid\t$$ref_is{'isu'}\n$start[0]\t$start[-1]\t$mstart[0]\t$mstart[-1]\n" if ($i == 79);
				if ($cluster[0])
				{
					$start{$lastpid}{$lastsid}[0] = $start[0];
					$start{$lastpid}{$lastsid}[1] = $start[-1];
					$mstart{$lastpid}{$lastsid}[0] = $mstart[0];
					$mstart{$lastpid}{$lastsid}[1] = $mstart[-1];
					$cbp{$lastpid}{$lastsid}[0] = $cbp[0];
					$cbp{$lastpid}{$lastsid}[1] = $cbp[1];
					$cbp{$lastpid}{$lastsid}[2] = $cbp[2];
					$cbp{$lastpid}{$lastsid}[3] = $cbp[3];
				}
				if ($orientation{$lastpid}{$lastsid}[0] == 1 and $orientation{$lastpid}{$lastsid}[1] == -1)
				{
					my $size = $mstart{$lastpid}{$lastsid}[0] - $start{$lastpid}{$lastsid}[1] if ($orientation{$lastpid}{$lastsid}[2] eq $orientation{$lastpid}{$lastsid}[3]);
					my $p_s_id = $lastpid.'_'.$lastsid;
					push @{$cluster_type{$orientation{$lastpid}{$lastsid}[2]}{$orientation{$lastpid}{$lastsid}[3]}{0}}, $p_s_id;
					print DISCCL "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n" if ($size > 0 or $orientation{$lastpid}{$lastsid}[2] ne $orientation{$lastpid}{$lastsid}[3]);
				}
				if ($orientation{$lastpid}{$lastsid}[0] == 1 and $orientation{$lastpid}{$lastsid}[1] == 1)
				{
					my $size = $mstart{$lastpid}{$lastsid}[1] - $start{$lastpid}{$lastsid}[1] if ($orientation{$lastpid}{$lastsid}[2] eq $orientation{$lastpid}{$lastsid}[3]);
					my $p_s_id = $lastpid.'_'.$lastsid;
					push @{$cluster_type{$orientation{$lastpid}{$lastsid}[2]}{$orientation{$lastpid}{$lastsid}[3]}{1}}, $p_s_id;
					print DISCCL "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
				}
				if ($orientation{$lastpid}{$lastsid}[0] == -1 and $orientation{$lastpid}{$lastsid}[1] == -1)
				{
					my $size = $mstart{$lastpid}{$lastsid}[0] - $start{$lastpid}{$lastsid}[0] if ($orientation{$lastpid}{$lastsid}[2] eq $orientation{$lastpid}{$lastsid}[3]);
					my $p_s_id = $lastpid.'_'.$lastsid;
					push @{$cluster_type{$orientation{$lastpid}{$lastsid}[2]}{$orientation{$lastpid}{$lastsid}[3]}{2}}, $p_s_id;
					print DISCCL "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
				}
				if ($orientation{$lastpid}{$lastsid}[0] == -1 and $orientation{$lastpid}{$lastsid}[1] == 1)
				{
					my $size = $mstart{$lastpid}{$lastsid}[1] - $start{$lastpid}{$lastsid}[0] if ($orientation{$lastpid}{$lastsid}[2] eq $orientation{$lastpid}{$lastsid}[3]);
					my $p_s_id = $lastpid.'_'.$lastsid;
					push @{$cluster_type{$orientation{$lastpid}{$lastsid}[2]}{$orientation{$lastpid}{$lastsid}[3]}{3}}, $p_s_id;
					print DISCCL "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
				}
				if ($l < $support_mps or $j < $support_mpf)
				{
					$start{$lastpid} = undef;
					$mstart{$lastpid} = undef;
					$cbp{$lastpid} = undef;
					$orientation{$lastpid} = undef;
					$support{$lastpid} = undef;
					$support_f{$lastpid} = undef;
				}
			}
			@cluster = undef;
			
			# next cluster
			$lastpid = $data[0];
			$lastsid = $data[1];
			$k = 0;
			$l = 0;
			$j = 0;
			%readname = undef;
		}
		$cluster[$k][0] = $data[5];
		$cluster[$k][1] = $data[6];
		$cluster[$k][2] = $data[7];
		$cluster[$k][3] = $data[8];
		$cluster[$k][4] = $data[9];
		$cluster[$k][5] = $data[9]+$data[10]-1;
		$cluster[$k][6] = $data[10];
		$cluster[$k][7] = $data[11];
		$cluster[$k][8] = $data[12];
		$cluster[$k][9] = $data[13];
		$cluster[$k][10] = $data[13]+$data[14]-1;
		$cluster[$k][11] = $data[14];
		$cluster[$k][12] = $data[15];
		$count = $data[0];
		#print STDERR "$cluster[$k][10]\n" if ($data[0] == 125);
		$k++ unless ($readname{$data[5]});
		if ($data[5] =~ /mu1|mu2|sc$/)
		{
			$l++ unless ($readname{$`});
			$readname{$`} = 1;
		}
		else
		{
			$j++;
			$l++ unless ($readname{$data[5]});
			$readname{$data[5]} = 1;
		}
	}
	close CLMPD;
	# parse last cluster
	$support{$lastpid} = $l;
	$support_f{$lastpid} = $j;
	#if ($support{$lastpid} >= $support_mps)
	{
		my (@start, @mstart, @isize, @cbp);
		foreach (@cluster)
		{
			my $readname = $$_[0]; my $rg = $$_[1]; my $seqid = $$_[2]; my $strand = $$_[3]; my $start = $$_[4]; my $end = $$_[5]; my $len = $$_[6]; my $mseqid = $$_[7]; my $mstrand = $$_[8]; my $mstart = $$_[9]; my $mend = $$_[10]; my $mlen = $$_[11]; my $isize = $$_[12];
			next unless ($seqid);
			push @start, $start if ($start);
			push @start, $end if ($end);
			push @mstart, $mstart if ($mstart);
			push @mstart, $mend if ($mend);
			push @isize, $isize;
			$end = $start unless ($end);
			$mend = $mstart unless ($mend);
			if ($strand == 1)
			{
				($cbp[0], $cbp[1]) = &overlap($end-10, $start+$$ref_is{$rg}{'isu'}-$mlen, $cbp[0], $cbp[1]);
				# -10 or +10: expand the window to avoid one read too close to a break point
			}
			else
			{
				($cbp[0], $cbp[1]) = &overlap($end-$$ref_is{$rg}{'isu'}+$mlen, $start+10, $cbp[0], $cbp[1]);
			}
			if ($mstrand == 1)
			{
				($cbp[2], $cbp[3]) = &overlap($mend-10, $mstart+$$ref_is{$rg}{'isu'}-$len, $cbp[2], $cbp[3]);
			}
			else
			{
				($cbp[2], $cbp[3]) = &overlap($mend-$$ref_is{$rg}{'isu'}+$len, $mstart+10, $cbp[2], $cbp[3]);
			}
			#print "$lastpid $readname\t$len $mlen $cbp[0]\t$cbp[1]\t$cbp[2]\t$cbp[3]\n";
			$orientation{$lastpid}{$lastsid}[0] = $strand;
			$orientation{$lastpid}{$lastsid}[1] = $mstrand;
			$orientation{$lastpid}{$lastsid}[2] = $seqid;
			$orientation{$lastpid}{$lastsid}[3] = $mseqid;
		}
		@start = sort {$a <=> $b} @start;
		@mstart = sort {$a <=> $b} @mstart;
		@isize = sort {$a <=> $b} @isize;
		#print "$lastpid\t$$ref_is{'isu'}\n$start[0]\t$start[-1]\t$mstart[0]\t$mstart[-1]\n" if ($i == 79);
		if ($cluster[0])
		{
			$start{$lastpid}{$lastsid}[0] = $start[0];
			$start{$lastpid}{$lastsid}[1] = $start[-1];
			$mstart{$lastpid}{$lastsid}[0] = $mstart[0];
			$mstart{$lastpid}{$lastsid}[1] = $mstart[-1];
			$cbp{$lastpid}{$lastsid}[0] = $cbp[0];
			$cbp{$lastpid}{$lastsid}[1] = $cbp[1];
			$cbp{$lastpid}{$lastsid}[2] = $cbp[2];
			$cbp{$lastpid}{$lastsid}[3] = $cbp[3];
		}
		if ($orientation{$lastpid}{$lastsid}[0] == 1 and $orientation{$lastpid}{$lastsid}[1] == -1)
		{
			my $size = $mstart{$lastpid}{$lastsid}[0] - $start{$lastpid}{$lastsid}[1] if ($orientation{$lastpid}{$lastsid}[2] eq $orientation{$lastpid}{$lastsid}[3]);
			my $p_s_id = $lastpid.'_'.$lastsid;
			push @{$cluster_type{$orientation{$lastpid}{$lastsid}[2]}{$orientation{$lastpid}{$lastsid}[3]}{0}}, $p_s_id;
			print DISCCL "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n" if ($size > 0 or $orientation{$lastpid}{$lastsid}[2] ne $orientation{$lastpid}{$lastsid}[3]);
		}
		if ($orientation{$lastpid}{$lastsid}[0] == 1 and $orientation{$lastpid}{$lastsid}[1] == 1)
		{
			my $size = $mstart{$lastpid}{$lastsid}[1] - $start{$lastpid}{$lastsid}[1] if ($orientation{$lastpid}{$lastsid}[2] eq $orientation{$lastpid}{$lastsid}[3]);
			my $p_s_id = $lastpid.'_'.$lastsid;
			push @{$cluster_type{$orientation{$lastpid}{$lastsid}[2]}{$orientation{$lastpid}{$lastsid}[3]}{1}}, $p_s_id;
			print DISCCL "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
		}
		if ($orientation{$lastpid}{$lastsid}[0] == -1 and $orientation{$lastpid}{$lastsid}[1] == -1)
		{
			my $size = $mstart{$lastpid}{$lastsid}[0] - $start{$lastpid}{$lastsid}[0] if ($orientation{$lastpid}{$lastsid}[2] eq $orientation{$lastpid}{$lastsid}[3]);
			my $p_s_id = $lastpid.'_'.$lastsid;
			push @{$cluster_type{$orientation{$lastpid}{$lastsid}[2]}{$orientation{$lastpid}{$lastsid}[3]}{2}}, $p_s_id;
			print DISCCL "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
		}
		if ($orientation{$lastpid}{$lastsid}[0] == -1 and $orientation{$lastpid}{$lastsid}[1] == 1)
		{
			my $size = $mstart{$lastpid}{$lastsid}[1] - $start{$lastpid}{$lastsid}[0] if ($orientation{$lastpid}{$lastsid}[2] eq $orientation{$lastpid}{$lastsid}[3]);
			my $p_s_id = $lastpid.'_'.$lastsid;
			push @{$cluster_type{$orientation{$lastpid}{$lastsid}[2]}{$orientation{$lastpid}{$lastsid}[3]}{3}}, $p_s_id;
			print DISCCL "$lastpid\t$lastsid\t$support{$lastpid}\t$orientation{$lastpid}{$lastsid}[2]\t$start{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[0]\t$orientation{$lastpid}{$lastsid}[3]\t$mstart{$lastpid}{$lastsid}[1]\t$orientation{$lastpid}{$lastsid}[1]\t$size\n";
		}
		if ($k < $support_mps or $j < $support_mpf)
		{
			$start{$lastpid} = undef;
			$mstart{$lastpid} = undef;
			$cbp{$lastpid} = undef;
			$orientation{$lastpid} = undef;
			$support{$lastpid} = undef;
			$support_f{$lastpid} = undef;
		}
	}
	@cluster = undef;
	close DISCCL;	
	#print "@{$cluster_type{chr21}{chr21}{0}}\n";
	
	my @used_cluster;
	my @del_cluster;
	my @result_mp; # store all results
	my $mpi = 0; # index in @result_mp
	
	# intra-chr events
	foreach my $chr (keys %cluster_type)
	{
		if ($cluster_type{$chr}{$chr}{0} and $cluster_type{$chr}{$chr}{3})
		{
			foreach (@{$cluster_type{$chr}{$chr}{3}})
			{
				my $id3 = $_;
				my ($i, $ii) = split (/_/, $id3);
				next if ($support{$i} < $support_mps or $support_f{$i} < $support_mpf);
				#print "$chr\t$id3\n";
				next if (($cbp{$i}{$ii}[0] == 0 and $cbp{$i}{$ii}[1] == 0) or ($cbp{$i}{$ii}[2] == 0 and $cbp{$i}{$ii}[3] == 0));
				
				my (@current_intra_event, %ca_event_size);
				my $cei = 0; # index for current event
				# @current_intra_event: all possible events from cluster i
				# %ca_event_size: all event sizes for @current_intra_event
				foreach (@{$cluster_type{$chr}{$chr}{0}})
				{
					my $id0 = $_;
					my ($k, $kk) = split (/_/, $id0);
					next if ($support{$k} < $support_mps or $support_f{$k} < $support_mpf);
					next if (($cbp{$k}{$kk}[0] == 0 and $cbp{$k}{$kk}[1] == 0) or ($cbp{$k}{$kk}[2] == 0 and $cbp{$k}{$kk}[3] == 0));
					next if ($i == $k);
					
					# deletion with insertion in the same orientation
					if ($orientation{$i}{$ii}[2] eq $orientation{$i}{$ii}[3] and $orientation{$k}{$kk}[2] eq $orientation{$k}{$kk}[3] and $orientation{$i}{$ii}[0] == -1 and $orientation{$i}{$ii}[1] == 1 and $orientation{$k}{$kk}[0] == 1 and $orientation{$k}{$kk}[1] == -1)
					{
						if (&covered($start{$i}{$ii}[0], $start{$i}{$ii}[1], $start{$k}{$kk}[1], $mstart{$k}{$kk}[0]) and &covered($mstart{$k}{$kk}[0], $mstart{$k}{$kk}[1], $start{$i}{$ii}[0], $mstart{$i}{$ii}[1]) and $start{$k}{$kk}[0] < $start{$i}{$ii}[0])
						{
							my $del_size = $start{$i}{$ii}[0] - $start{$k}{$kk}[1] - 1;
							my $ins_size = $mstart{$i}{$ii}[1] - $mstart{$k}{$kk}[0] + 1;
							my $distance = $mstart{$k}{$kk}[0] - $start{$i}{$ii}[0];
							if ($del_size < $sv_size_cutoff and $ins_size < $sv_size_cutoff)
							{
								if ($del_size < $$ref_is{'rlu'})
								{
									$current_intra_event[$cei][0] = 'inssd';
								}
								else
								{
									$current_intra_event[$cei][0] = 'del_inssd';
								}
								$current_intra_event[$cei][1] = "$i\_$ii/$k\_$kk";
								$current_intra_event[$cei][2] = "$support{$i}/$support{$k}";
								$current_intra_event[$cei][3] = $orientation{$i}{$ii}[2];
								$current_intra_event[$cei][4] = $start{$k}{$kk}[1];
								$current_intra_event[$cei][5] = $start{$i}{$ii}[0];
								$current_intra_event[$cei][6] = $del_size;
								$current_intra_event[$cei][7] = "$orientation{$i}{$ii}[2]\t$mstart{$k}{$kk}[0]\t$mstart{$i}{$ii}[1]\t$ins_size\t$distance\t$cbp{$i}{$ii}[0]:$cbp{$i}{$ii}[1]:$cbp{$i}{$ii}[2]:$cbp{$i}{$ii}[3]\t$cbp{$k}{$kk}[0]:$cbp{$k}{$kk}[1]:$cbp{$k}{$kk}[2]:$cbp{$k}{$kk}[3]";
								$ca_event_size{$cei} = abs($del_size) + abs($ins_size);
								$cei++;
							}
						}
						if (&covered($start{$k}{$kk}[0], $start{$k}{$kk}[1], $start{$i}{$ii}[0], $mstart{$i}{$ii}[1]) and &covered($mstart{$i}{$ii}[0], $mstart{$i}{$ii}[1], $start{$k}{$kk}[1], $mstart{$k}{$kk}[0]) and $mstart{$i}{$ii}[1] < $mstart{$k}{$kk}[1])
						{
							my $del_size = $mstart{$k}{$kk}[0] - $mstart{$i}{$ii}[1] - 1;
							my $ins_size = $start{$k}{$kk}[1] - $start{$i}{$ii}[0] + 1;
							my $distance = $mstart{$i}{$ii}[1] - $start{$k}{$kk}[1];
							if ($del_size < $sv_size_cutoff and $ins_size < $sv_size_cutoff)
							{
								if ($del_size < $$ref_is{'rlu'})
								{
									$current_intra_event[$cei][0] = 'inssu';
								}
								else
								{
									$current_intra_event[$cei][0] = 'del_inssu';
								}
								$current_intra_event[$cei][1] = "$i\_$ii/$k\_$kk";
								$current_intra_event[$cei][2] = "$support{$i}/$support{$k}";
								$current_intra_event[$cei][3] = $orientation{$i}{$ii}[2];
								$current_intra_event[$cei][4] = $mstart{$i}{$ii}[1];
								$current_intra_event[$cei][5] = $mstart{$k}{$kk}[0];
								$current_intra_event[$cei][6] = $del_size;
								$current_intra_event[$cei][7] = "$orientation{$i}{$ii}[2]\t$start{$i}{$ii}[0]\t$start{$k}{$kk}[1]\t$ins_size\t$distance\t$cbp{$i}{$ii}[0]:$cbp{$i}{$ii}[1]:$cbp{$i}{$ii}[2]:$cbp{$i}{$ii}[3]\t$cbp{$k}{$kk}[0]:$cbp{$k}{$kk}[1]:$cbp{$k}{$kk}[2]:$cbp{$k}{$kk}[3]";
								$ca_event_size{$cei} = abs($del_size) + abs($ins_size);
								$cei++;
							}
						}
					}
				}
				
				# call smallest event for cluster i
				my ($top1, $top2, $top3, $select_cei);
				my $chi = 0;
				foreach my $cei (sort {$ca_event_size{$a} <=> $ca_event_size{$b}} (keys(%ca_event_size)))
				{
					$top1 = $cei if ($chi == 0);
					$top2 = $cei if ($chi == 1);
					$top3 = $cei if ($chi == 2);
					$chi++;
					last if ($chi == 2);
				}
				if (($current_intra_event[$top1][0] =~ /insod/ or $current_intra_event[$top1][0] =~ /insou/) and ($current_intra_event[$top2][0] eq 'invers' or $current_intra_event[$top3][0] eq 'invers'))
				{
					if ($current_intra_event[$top1][6] < 20 and $current_intra_event[$top1][9] < 20 and $current_intra_event[$top1][10] > 20)
					{
						if ($current_intra_event[$top2][0] eq 'invers')
						{
							@{$result_mp[$mpi]} = @{$current_intra_event[$top2]};
							my @support_cl = split (/\//, $current_intra_event[$top2][1]);
							push @used_cluster, $support_cl[0];
							push @used_cluster, $support_cl[1];
							$mpi++;
						}
						elsif ($current_intra_event[$top3][0] eq 'invers')
						{
							@{$result_mp[$mpi]} = @{$current_intra_event[$top3]};
							my @support_cl = split (/\//, $current_intra_event[$top3][1]);
							push @used_cluster, $support_cl[0];
							push @used_cluster, $support_cl[1];
							$mpi++;
						}
					}
					else
					{
						@{$result_mp[$mpi]} = @{$current_intra_event[$top1]};
						my @support_cl = split (/\//, $current_intra_event[$top1][1]);
						push @used_cluster, $support_cl[0];
						push @used_cluster, $support_cl[1];
						$mpi++;
					}
				}
				else
				{
					if ($current_intra_event[$top1][0])
					{
						@{$result_mp[$mpi]} = @{$current_intra_event[$top1]};
						my @support_cl = split (/\//, $current_intra_event[$top1][1]);
						push @used_cluster, $support_cl[0];
						push @used_cluster, $support_cl[1];
						$mpi++;
					}
				}
			}
		}
		
		if ($cluster_type{$chr}{$chr}{1} and $cluster_type{$chr}{$chr}{2})
		{
			foreach (@{$cluster_type{$chr}{$chr}{1}})
			{
				my $id1 = $_;
				my ($i, $ii) = split (/_/, $id1);
				next if ($support{$i} < $support_mps or $support_f{$i} < $support_mpf);
				next if (($cbp{$i}{$ii}[0] == 0 and $cbp{$i}{$ii}[1] == 0) or ($cbp{$i}{$ii}[2] == 0 and $cbp{$i}{$ii}[3] == 0));
				
				my (@current_intra_event, %ca_event_size);
				my $cei = 0; # index for current event
				# @current_intra_event: all possible events from cluster i
				# %ca_event_size: all event sizes for @current_intra_event
				foreach (@{$cluster_type{$chr}{$chr}{2}})
				{
					my $id2 = $_;
					#print "id1:$id1 id2:$id2 ";
					my ($k, $kk) = split (/_/, $id2);
					next if ($support{$k} < $support_mps or $support_f{$k} < $support_mpf);
					next if (($cbp{$k}{$kk}[0] == 0 and $cbp{$k}{$kk}[1] == 0) or ($cbp{$k}{$kk}[2] == 0 and $cbp{$k}{$kk}[3] == 0));
					next if ($i == $k);
					
					# inversions				
					if (&covered($start{$k}{$kk}[0], $start{$k}{$kk}[1], $start{$i}{$ii}[1], $mstart{$i}{$ii}[1]) and &covered($mstart{$i}{$ii}[0], $mstart{$i}{$ii}[1], $start{$k}{$kk}[0], $mstart{$k}{$kk}[0]) and $start{$i}{$ii}[0] < $start{$k}{$kk}[0] and $start{$i}{$ii}[1] < $start{$k}{$kk}[1] and $mstart{$i}{$ii}[0] < $mstart{$k}{$kk}[0] and $mstart{$i}{$ii}[1] < $mstart{$k}{$kk}[1])
					{
						my $inv_size = $mstart{$i}{$ii}[1] - $start{$k}{$kk}[0];
						if ($inv_size < $sv_size_cutoff)
						{
							if (&covered($start{$i}{$ii}[1] - 10, $start{$i}{$ii}[1] + $$ref_is{'rlu'}, $start{$k}{$kk}[0] - $$ref_is{'rlu'}, $start{$k}{$kk}[0] + 10) and &covered($mstart{$i}{$ii}[1] - 10, $mstart{$i}{$ii}[1] + $$ref_is{'rlu'}, $mstart{$k}{$kk}[0] - $$ref_is{'rlu'}, $mstart{$k}{$kk}[0] + 10))
							{
								$current_intra_event[$cei][0] = 'invers';
							}
							else
							{
								$current_intra_event[$cei][0] = 'del_invers';
							}
							$current_intra_event[$cei][1] = "$i\_$ii/$k\_$kk";
							$current_intra_event[$cei][2] = "$support{$i}/$support{$k}";
							$current_intra_event[$cei][3] = $orientation{$i}{$ii}[2];
							$current_intra_event[$cei][4] = $start{$i}{$ii}[1];
							$current_intra_event[$cei][5] = $start{$k}{$kk}[0];
							$current_intra_event[$cei][6] = $inv_size;
							$current_intra_event[$cei][7] = $orientation{$i}{$ii}[2];
							$current_intra_event[$cei][8] = $mstart{$i}{$ii}[1];
							$current_intra_event[$cei][9] = $mstart{$k}{$kk}[0];
							$current_intra_event[$cei][10] = "$cbp{$i}{$ii}[0]:$cbp{$i}{$ii}[1]:$cbp{$i}{$ii}[2]:$cbp{$i}{$ii}[3]\t$cbp{$k}{$kk}[0]:$cbp{$k}{$kk}[1]:$cbp{$k}{$kk}[2]:$cbp{$k}{$kk}[3]";
							$ca_event_size{$cei} = abs($inv_size);
							$cei++;
						}
					}
					
					# deletion with insertion in opposite orientation
					if (&covered($start{$k}{$kk}[0], $start{$k}{$kk}[1], $start{$i}{$ii}[1], $mstart{$i}{$ii}[1]) and &covered($mstart{$k}{$kk}[0], $mstart{$k}{$kk}[1], $start{$i}{$ii}[1], $mstart{$i}{$ii}[1]) and $start{$i}{$ii}[0] < $start{$k}{$kk}[0])
					{
						my $del_size = $start{$k}{$kk}[0] - $start{$i}{$ii}[1] - 1;
						my $ins_size = $mstart{$i}{$ii}[1] - $mstart{$k}{$kk}[0] + 1;
						my $distance = $mstart{$k}{$kk}[0] - $start{$k}{$kk}[0];
						if ($del_size < $sv_size_cutoff and $ins_size < $sv_size_cutoff)
						{
							if ($del_size < $$ref_is{'rlu'})
							{
								$current_intra_event[$cei][0] = 'insod';
							}
							else
							{
								$current_intra_event[$cei][0] = 'del_insod';
							}
							$current_intra_event[$cei][1] = "$i\_$ii/$k\_$kk";
							$current_intra_event[$cei][2] = "$support{$i}/$support{$k}";
							$current_intra_event[$cei][3] = $orientation{$i}{$ii}[2];
							$current_intra_event[$cei][4] = $start{$i}{$ii}[1];
							$current_intra_event[$cei][5] = $start{$k}{$kk}[0];
							$current_intra_event[$cei][6] = $del_size;
							$current_intra_event[$cei][7] = $orientation{$i}{$ii}[2];
							$current_intra_event[$cei][8] = $mstart{$k}{$kk}[0];
							$current_intra_event[$cei][9] = $mstart{$i}{$ii}[1];
							$current_intra_event[$cei][10] = $ins_size;
							$current_intra_event[$cei][11] = $distance;
							$current_intra_event[$cei][12] = "$cbp{$i}{$ii}[0]:$cbp{$i}{$ii}[1]:$cbp{$i}{$ii}[2]:$cbp{$i}{$ii}[3]\t$cbp{$k}{$kk}[0]:$cbp{$k}{$kk}[1]:$cbp{$k}{$kk}[2]:$cbp{$k}{$kk}[3]";
							$ca_event_size{$cei} = abs($del_size) + abs($ins_size);
							$cei++;
						}
					}
					
					# deletion with insertion in opposite orientation
					if (&covered($start{$i}{$ii}[0], $start{$i}{$ii}[1], $start{$k}{$kk}[0], $mstart{$k}{$kk}[0]) and &covered($mstart{$i}{$ii}[0], $mstart{$i}{$ii}[1], $start{$k}{$kk}[0], $mstart{$k}{$kk}[0]) and $mstart{$i}{$ii}[1] < $mstart{$k}{$kk}[1])
					{
						my $del_size = $mstart{$k}{$kk}[0] - $mstart{$i}{$ii}[1] - 1;
						my $ins_size = $start{$i}{$ii}[1] - $start{$k}{$kk}[0] + 1;
						my $distance = $mstart{$i}{$ii}[1] - $start{$i}{$ii}[1];
						if ($del_size < $sv_size_cutoff and $ins_size < $sv_size_cutoff)
						{
							if ($del_size < $$ref_is{'rlu'})
							{
								$current_intra_event[$cei][0] = 'insou';
							}
							else
							{
								$current_intra_event[$cei][0] = 'del_insou';
							}
							$current_intra_event[$cei][1] = "$i\_$ii/$k\_$kk";
							$current_intra_event[$cei][2] = "$support{$i}/$support{$k}";
							$current_intra_event[$cei][3] = $orientation{$i}{$ii}[2];
							$current_intra_event[$cei][4] = $mstart{$i}{$ii}[1];
							$current_intra_event[$cei][5] = $mstart{$k}{$kk}[0];
							$current_intra_event[$cei][6] = $del_size;
							$current_intra_event[$cei][7] = $orientation{$i}{$ii}[2];
							$current_intra_event[$cei][8] = $start{$k}{$kk}[0];
							$current_intra_event[$cei][9] = $start{$i}{$ii}[1];
							$current_intra_event[$cei][10] = $ins_size;
							$current_intra_event[$cei][11] = $distance;
							$current_intra_event[$cei][12] = "$cbp{$i}{$ii}[0]:$cbp{$i}{$ii}[1]:$cbp{$i}{$ii}[2]:$cbp{$i}{$ii}[3]\t$cbp{$k}{$kk}[0]:$cbp{$k}{$kk}[1]:$cbp{$k}{$kk}[2]:$cbp{$k}{$kk}[3]";
							$ca_event_size{$cei} = abs($del_size) + abs($ins_size);
							$cei++;
						}
					}
				}
				
				# call smallest event for cluster i
				my ($top1, $top2, $top3, $select_cei);
				my $chi = 0;
				foreach my $cei (sort {$ca_event_size{$a} <=> $ca_event_size{$b}} (keys(%ca_event_size)))
				{
					$top1 = $cei if ($chi == 0);
					$top2 = $cei if ($chi == 1);
					$top3 = $cei if ($chi == 2);
					$chi++;
					last if ($chi == 2);
				}
				if (($current_intra_event[$top1][0] =~ /insod/ or $current_intra_event[$top1][0] =~ /insou/) and ($current_intra_event[$top2][0] eq 'invers' or $current_intra_event[$top3][0] eq 'invers'))
				{
					if ($current_intra_event[$top1][6] < 20 and $current_intra_event[$top1][10] < 20 and $current_intra_event[$top1][11] > 20)
					{
						if ($current_intra_event[$top2][0] eq 'invers')
						{
							@{$result_mp[$mpi]} = @{$current_intra_event[$top2]};
							my @support_cl = split (/\//, $current_intra_event[$top2][1]);
							push @used_cluster, $support_cl[0];
							push @used_cluster, $support_cl[1];
							$mpi++;
						}
						elsif ($current_intra_event[$top3][0] eq 'invers')
						{
							@{$result_mp[$mpi]} = @{$current_intra_event[$top3]};
							my @support_cl = split (/\//, $current_intra_event[$top3][1]);
							push @used_cluster, $support_cl[0];
							push @used_cluster, $support_cl[1];
							$mpi++;
						}
					}
					else
					{
						@{$result_mp[$mpi]} = @{$current_intra_event[$top1]};
						my @support_cl = split (/\//, $current_intra_event[$top1][1]);
						push @used_cluster, $support_cl[0];
						push @used_cluster, $support_cl[1];
						$mpi++;
					}
				}
				else
				{
					if ($current_intra_event[$top1][0])
					{
						@{$result_mp[$mpi]} = @{$current_intra_event[$top1]};
						my @support_cl = split (/\//, $current_intra_event[$top1][1]);
						push @used_cluster, $support_cl[0];
						push @used_cluster, $support_cl[1];
						$mpi++;
					}
				}
			}
		}
		
	}
	
	# inter-chr events
	my %used_inter;
	foreach my $chr1 (keys %cluster_type)
	{
		foreach my $chr2 (keys %{$cluster_type{$chr1}})
		{
			next if ($chr1 eq $chr2);
			next unless ($cluster_type{$chr1}{$chr2}{0} or $cluster_type{$chr1}{$chr2}{1} or $cluster_type{$chr1}{$chr2}{2} or $cluster_type{$chr1}{$chr2}{3});
			my @all_inter;
			push @all_inter, @{$cluster_type{$chr1}{$chr2}{0}} if ($cluster_type{$chr1}{$chr2}{0});
			push @all_inter, @{$cluster_type{$chr1}{$chr2}{1}} if ($cluster_type{$chr1}{$chr2}{1});
			push @all_inter, @{$cluster_type{$chr1}{$chr2}{2}} if ($cluster_type{$chr1}{$chr2}{2});
			push @all_inter, @{$cluster_type{$chr1}{$chr2}{3}} if ($cluster_type{$chr1}{$chr2}{3});
			foreach (@all_inter)
			{
				my $id = $_;
				my ($i, $ii) = split (/_/, $id);
				next if ($support{$i} < $support_mps or $support_f{$i} < $support_mpf);
				next if ($ii); # start from 0 secondary id
				my (%distance1, %distance2, $m, $n);
				# %distance1, %distance2: event size of possible paired clusters to query cluster
				# $m, $n: cluster id of all possible paired clusters
				for (my $ii=0;1;$ii++)
				{
					last unless ($start{$i}{$ii}[0]);
					next if (($cbp{$i}{$ii}[0] == 0 and $cbp{$i}{$ii}[1] == 0) or ($cbp{$i}{$ii}[2] == 0 and $cbp{$i}{$ii}[3] == 0));
					
					if ($orientation{$i}{$ii}[0] == 1 and $orientation{$i}{$ii}[1] == -1)
					{
						if ($cluster_type{$orientation{$i}{$ii}[2]}{$orientation{$i}{$ii}[3]}{3})
						{
							foreach (@{$cluster_type{$orientation{$i}{$ii}[2]}{$orientation{$i}{$ii}[3]}{3}})
							{
								my $id3 = $_;
								my ($k, $kk) = split (/_/, $id3);
								next if ($support{$k} < $support_mps or $support_f{$k} < $support_mpf);
								next if (($cbp{$k}{$kk}[0] == 0 and $cbp{$k}{$kk}[1] == 0) or ($cbp{$k}{$kk}[2] == 0 and $cbp{$k}{$kk}[3] == 0));
								next if ($i == $k);
								
								my $dis_id = $i.'_'.$ii.'_'.$k.'_'.$kk;
								$distance1{$dis_id} = abs($start{$k}{$kk}[0] - $start{$i}{$ii}[0]);
								$distance2{$dis_id} = abs($mstart{$k}{$kk}[0] - $mstart{$i}{$ii}[0]);
								#print "$dis_id\t$distance1{$dis_id}\t$distance2{$dis_id}\n";
							}
						}
					}
					elsif ($orientation{$i}{$ii}[0] == 1 and $orientation{$i}{$ii}[1] == 1)
					{
						if ($cluster_type{$orientation{$i}{$ii}[2]}{$orientation{$i}{$ii}[3]}{2})
						{
							foreach (@{$cluster_type{$orientation{$i}{$ii}[2]}{$orientation{$i}{$ii}[3]}{2}})
							{
								my $id2 = $_;
								my ($k, $kk) = split (/_/, $id2);
								next if ($support{$k} < $support_mps or $support_f{$k} < $support_mpf);
								next if (($cbp{$k}{$kk}[0] == 0 and $cbp{$k}{$kk}[1] == 0) or ($cbp{$k}{$kk}[2] == 0 and $cbp{$k}{$kk}[3] == 0));
								next if ($i == $k);
								
								my $dis_id = $i.'_'.$ii.'_'.$k.'_'.$kk;
								$distance1{$dis_id} = abs($start{$k}{$kk}[0] - $start{$i}{$ii}[0]);
								$distance2{$dis_id} = abs($mstart{$k}{$kk}[0] - $mstart{$i}{$ii}[0]);
#								print "$dis_id\t$distance1{$dis_id}\t$distance2{$dis_id}\n";
							}
						}
					}
					
				}
				my ($smallestm, $smallestn, @m, @n);
				my $sei = 0;
				foreach my $key (sort { $distance1 {$a} <=> $distance1 {$b}} keys %distance1)
				{
					if ($sei == 0)
					{
						push @m, $key;
						$smallestm = $distance1{$key};
						$sei++;
						next;
					}
					if ($smallestm == $distance1{$key})
					{
						push @m, $key;
					}
					last if ($distance1{$key} > $smallestm);
				}
				my $sei = 0;
				foreach my $key (sort { $distance2 {$a} <=> $distance2 {$b}} keys %distance2)
				{
					if ($sei == 0)
					{
						push @n, $key;
						$smallestn = $distance2{$key};
						$sei++;
						next;
					}
					if ($smallestn == $distance2{$key})
					{
						push @n, $key;
					}
					last if ($distance2{$key} > $smallestn);
				}
				#print "$i\t@m\t@n\n";
SEI:				foreach (@m)
				{
					$m = $_;
					foreach (@n)
					{
						$n = $_;
						if ($m eq $n)
						{
							last SEI;
						}
					}
				}
				#print "$i\t$m\t$n\t$distance1{$m}\t$distance2{$m}\n";
				if ($m eq $n and defined($m))
				{
					my ($trash, $ii, $k, $kk) = split (/_/, $m);
					if ($k != $i and $orientation{$i}{$ii}[2] lt $orientation{$i}{$ii}[3])
					{
						if ($orientation{$i}{$ii}[0] == $orientation{$i}{$ii}[1])
						{
							# acceptor on smaller chr
							if ($start{$k}{$kk}[0] - $start{$i}{$ii}[1] - 1 > -$$ref_is{'rlu'} and $mstart{$i}{$ii}[1] > $mstart{$k}{$kk}[0])
							{
								my $del_size = $start{$k}{$kk}[0] - $start{$i}{$ii}[1] - 1;
								my $ins_size = $mstart{$i}{$ii}[1] - $mstart{$k}{$kk}[0] + 1;
								if ($del_size < $sv_size_cutoff and $ins_size < $sv_size_cutoff)
								{
									push @used_cluster, $i;
									push @used_cluster, $k;
									$used_inter{$i} = 1;
									$used_inter{$k} = 1;
									if ($del_size < $$ref_is{'rlu'})
									{
										$result_mp[$mpi][0] = 'inso';
									}
									else
									{
										$result_mp[$mpi][0] = 'del_inso';
									}
									$result_mp[$mpi][1] = "$i\_$ii/$k\_$kk";
									$result_mp[$mpi][2] = "$support{$i}/$support{$k}";
									$result_mp[$mpi][3] = $orientation{$i}{$ii}[2];
									$result_mp[$mpi][4] = $start{$i}{$ii}[1];
									$result_mp[$mpi][5] = $start{$k}{$kk}[0];
									$result_mp[$mpi][6] = $del_size;
									$result_mp[$mpi][7] = "$orientation{$i}{$ii}[3]\t$mstart{$k}{$kk}[0]\t$mstart{$i}{$ii}[1]\t$ins_size\t$cbp{$i}{$ii}[0]:$cbp{$i}{$ii}[1]:$cbp{$i}{$ii}[2]:$cbp{$i}{$ii}[3]\t$cbp{$k}{$kk}[0]:$cbp{$k}{$kk}[1]:$cbp{$k}{$kk}[2]:$cbp{$k}{$kk}[3]";
									$mpi++;
								}
							}
							# donor on smaller chr
							elsif ($mstart{$k}{$kk}[0] - $mstart{$i}{$ii}[1] - 1 > -$$ref_is{'rlu'} and $start{$i}{$ii}[1] > $start{$k}{$kk}[0])
							{
								my $del_size = $mstart{$k}{$kk}[0] - $mstart{$i}{$ii}[1] - 1;
								my $ins_size = $start{$i}{$ii}[1] - $start{$k}{$kk}[0] + 1;
								if ($del_size < $sv_size_cutoff and $ins_size < $sv_size_cutoff)
								{
									push @used_cluster, $i;
									push @used_cluster, $k;
									$used_inter{$i} = 1;
									$used_inter{$k} = 1;
									if ($del_size < $$ref_is{'rlu'})
									{
										$result_mp[$mpi][0] = 'inso';
									}
									else
									{
										$result_mp[$mpi][0] = 'del_inso';
									}
									$result_mp[$mpi][1] = "$i\_$ii/$k\_$kk";
									$result_mp[$mpi][2] = "$support{$i}/$support{$k}";
									$result_mp[$mpi][3] = $orientation{$i}{$ii}[3];
									$result_mp[$mpi][4] = $mstart{$i}{$ii}[1];
									$result_mp[$mpi][5] = $mstart{$k}{$kk}[0];
									$result_mp[$mpi][6] = $del_size;
									$result_mp[$mpi][7] = "$orientation{$i}{$ii}[2]\t$start{$k}{$kk}[0]\t$start{$i}{$ii}[1]\t$ins_size\t$cbp{$i}{$ii}[2]:$cbp{$i}{$ii}[3]:$cbp{$i}{$ii}[0]:$cbp{$i}{$ii}[1]\t$cbp{$k}{$kk}[2]:$cbp{$k}{$kk}[3]:$cbp{$k}{$kk}[0]:$cbp{$k}{$kk}[1]";
									$mpi++;
								}
							}
						}
						else
						{
							# acceptor on smaller chr
							if ($start{$k}{$kk}[0] - $start{$i}{$ii}[1] - 1 > -$$ref_is{'rlu'} and $mstart{$k}{$kk}[1] > $mstart{$i}{$ii}[0])
							{
								my $del_size = $start{$k}{$kk}[0] - $start{$i}{$ii}[1] - 1;
								my $ins_size = $mstart{$k}{$kk}[1] - $mstart{$i}{$ii}[0] + 1;
								if ($del_size < $sv_size_cutoff and $ins_size < $sv_size_cutoff)
								{
									push @used_cluster, $i;
									push @used_cluster, $k;
									$used_inter{$i} = 1;
									$used_inter{$k} = 1;
									if ($del_size < $$ref_is{'rlu'})
									{
										$result_mp[$mpi][0] = 'inss';
									}
									else
									{
										$result_mp[$mpi][0] = 'del_inss';
									}
									$result_mp[$mpi][1] = "$i\_$ii/$k\_$kk";
									$result_mp[$mpi][2] = "$support{$i}/$support{$k}";
									$result_mp[$mpi][3] = $orientation{$i}{$ii}[2];
									$result_mp[$mpi][4] = $start{$i}{$ii}[1];
									$result_mp[$mpi][5] = $start{$k}{$kk}[0];
									$result_mp[$mpi][6] = $del_size;
									$result_mp[$mpi][7] = "$orientation{$i}{$ii}[3]\t$mstart{$i}{$ii}[0]\t$mstart{$k}{$kk}[1]\t$ins_size\t$cbp{$i}{$ii}[0]:$cbp{$i}{$ii}[1]:$cbp{$i}{$ii}[2]:$cbp{$i}{$ii}[3]\t$cbp{$k}{$kk}[0]:$cbp{$k}{$kk}[1]:$cbp{$k}{$kk}[2]:$cbp{$k}{$kk}[3]";
									$mpi++;
								}
							}
							# donor on smaller chr
							elsif ($mstart{$i}{$ii}[0] - $mstart{$k}{$kk}[1] - 1 > -$$ref_is{'rlu'} and $start{$i}{$ii}[1] > $start{$k}{$kk}[0])
							{
								my $del_size = $mstart{$i}{$ii}[0] - $mstart{$k}{$kk}[1] - 1;
								my $ins_size = $start{$i}{$ii}[1] - $start{$k}{$kk}[0] + 1;
								if ($del_size < $sv_size_cutoff and $ins_size < $sv_size_cutoff)
								{
									push @used_cluster, $i;
									push @used_cluster, $k;
									$used_inter{$i} = 1;
									$used_inter{$k} = 1;
									if ($del_size < $$ref_is{'rlu'})
									{
										$result_mp[$mpi][0] = 'inss';
									}
									else
									{
										$result_mp[$mpi][0] = 'del_inss';
									}
									$result_mp[$mpi][1] = "$i\_$ii/$k\_$kk";
									$result_mp[$mpi][2] = "$support{$i}/$support{$k}";
									$result_mp[$mpi][3] = $orientation{$i}{$ii}[3];
									$result_mp[$mpi][4] = $mstart{$k}{$kk}[1];
									$result_mp[$mpi][5] = $mstart{$i}{$ii}[0];
									$result_mp[$mpi][6] = $del_size;
									$result_mp[$mpi][7] = "$orientation{$i}{$ii}[2]\t$start{$k}{$kk}[0]\t$start{$i}{$ii}[1]\t$ins_size\t$cbp{$i}{$ii}[2]:$cbp{$i}{$ii}[3]:$cbp{$i}{$ii}[0]:$cbp{$i}{$ii}[1]\t$cbp{$k}{$kk}[2]:$cbp{$k}{$kk}[3]:$cbp{$k}{$kk}[0]:$cbp{$k}{$kk}[1]";
									$mpi++;
								}
							}
						}
					}
				}
			}
		}
	}
	
#	left over clusters
RESD:	for (my $i=0;$i<=$count;$i++)
	{
		next if ($support{$i} < $support_mps or $support_f{$i} < $support_mpf);
		next if (($cbp{$i}{0}[0] == 0 and $cbp{$i}{0}[1] == 0) or ($cbp{$i}{0}[2] == 0 and $cbp{$i}{0}[3] == 0));
		foreach (@used_cluster)
		{
			next RESD if ($i == $_);
		}
		#print "$i\t$start{$i}{0}[0]\t$start{$i}{0}[1]\t$orientation{$i}{0}[2]\t$orientation{$i}{0}[3]\n";
		if ($support{$i})
		{
			if ($orientation{$i}{0}[2] eq $orientation{$i}{0}[3])
			{
				# simple deletion
				if ($orientation{$i}{0}[0] == 1 and $orientation{$i}{0}[1] == -1)
				{
					my $printed;
					foreach (@del_cluster)
					{
						$printed = 1 if ($_ == $i);
					}
					unless ($printed)
					{
						my $del_size = $mstart{$i}{0}[0] - $start{$i}{0}[1] - 1;
						if ($del_size > 0 and $del_size < $sv_size_cutoff)
						{
							push @del_cluster, $i;
							$result_mp[$mpi][0] = 'del';
							$result_mp[$mpi][1] = $i;
							$result_mp[$mpi][2] = $support{$i};
							$result_mp[$mpi][3] = $orientation{$i}{0}[2];
							$result_mp[$mpi][4] = $start{$i}{0}[1];
							$result_mp[$mpi][5] = $mstart{$i}{0}[0];
							$result_mp[$mpi][6] = $del_size;
							$result_mp[$mpi][7] = "$cbp{$i}{0}[0]:$cbp{$i}{0}[1]:$cbp{$i}{0}[2]:$cbp{$i}{0}[3]";
							$mpi++;
						}
					}
				}
				if ($orientation{$i}{0}[0] == 1 and $orientation{$i}{0}[1] == 1)
				{
					my $invers_size = $mstart{$i}{0}[1] - $start{$i}{0}[1];
					if ($invers_size < $sv_size_cutoff)
					{
						$result_mp[$mpi][0] = 'invers_f';
						$result_mp[$mpi][1] = $i;
						$result_mp[$mpi][2] = $support{$i};
						$result_mp[$mpi][3] = $orientation{$i}{0}[2];
						$result_mp[$mpi][4] = $start{$i}{0}[1];
						$result_mp[$mpi][5] = $mstart{$i}{0}[1];
						$result_mp[$mpi][6] = $invers_size;
						$result_mp[$mpi][7] = "$cbp{$i}{0}[0]:$cbp{$i}{0}[1]:$cbp{$i}{0}[2]:$cbp{$i}{0}[3]";
						$mpi++;
					}
				}
				if ($orientation{$i}{0}[0] == -1 and $orientation{$i}{0}[1] == -1)
				{
					my $invers_size = $mstart{$i}{0}[0] - $start{$i}{0}[0];
					if ($invers_size < $sv_size_cutoff)
					{
						$result_mp[$mpi][0] = 'invers_r';
						$result_mp[$mpi][1] = $i;
						$result_mp[$mpi][2] = $support{$i};
						$result_mp[$mpi][3] = $orientation{$i}{0}[2];
						$result_mp[$mpi][4] = $start{$i}{0}[0];
						$result_mp[$mpi][5] = $mstart{$i}{0}[0];
						$result_mp[$mpi][6] = $invers_size;
						$result_mp[$mpi][7] = "$cbp{$i}{0}[0]:$cbp{$i}{0}[1]:$cbp{$i}{0}[2]:$cbp{$i}{0}[3]";
						$mpi++;
					}
				}
				if ($orientation{$i}{0}[0] == -1 and $orientation{$i}{0}[1] == 1)
				{
					my $transl_size = $mstart{$i}{0}[0] - $start{$i}{0}[0];
					if ($transl_size < $sv_size_cutoff)
					{
						$result_mp[$mpi][0] = 'tandem_dup';
						$result_mp[$mpi][1] = $i;
						$result_mp[$mpi][2] = $support{$i};
						$result_mp[$mpi][3] = $orientation{$i}{0}[2];
						$result_mp[$mpi][4] = $start{$i}{0}[0];
						$result_mp[$mpi][5] = $mstart{$i}{0}[1];
						$result_mp[$mpi][6] = $transl_size;
						$result_mp[$mpi][7] = "$cbp{$i}{0}[0]:$cbp{$i}{0}[1]:$cbp{$i}{0}[2]:$cbp{$i}{0}[3]";
						$mpi++;
					}
				}
			}
			
			# inter chromosomal events
			else
			{
				if ($orientation{$i}{0}[0] == 1 and $orientation{$i}{0}[1] == -1)
				{
					$result_mp[$mpi][0] = 'transl_inter';
					$result_mp[$mpi][1] = $i;
					$result_mp[$mpi][2] = $support{$i};
					$result_mp[$mpi][3] = $orientation{$i}{0}[2];
					$result_mp[$mpi][4] = $start{$i}{0}[1];
					$result_mp[$mpi][5] = $orientation{$i}{0}[0];
					$result_mp[$mpi][6] = $orientation{$i}{0}[3];
					$result_mp[$mpi][7] = $mstart{$i}{0}[0];
					$result_mp[$mpi][8] = $orientation{$i}{0}[1];
					$result_mp[$mpi][9] = "$cbp{$i}{0}[0]:$cbp{$i}{0}[1]:$cbp{$i}{0}[2]:$cbp{$i}{0}[3]";
					$mpi++;
				}
				if ($orientation{$i}{0}[0] == -1 and $orientation{$i}{0}[1] == 1)
				{
					$result_mp[$mpi][0] = 'transl_inter';
					$result_mp[$mpi][1] = $i;
					$result_mp[$mpi][2] = $support{$i};
					$result_mp[$mpi][3] = $orientation{$i}{0}[2];
					$result_mp[$mpi][4] = $start{$i}{0}[0];
					$result_mp[$mpi][5] = $orientation{$i}{0}[0];
					$result_mp[$mpi][6] = $orientation{$i}{0}[3];
					$result_mp[$mpi][7] = $mstart{$i}{0}[1];
					$result_mp[$mpi][8] = $orientation{$i}{0}[1];
					$result_mp[$mpi][9] = "$cbp{$i}{0}[0]:$cbp{$i}{0}[1]:$cbp{$i}{0}[2]:$cbp{$i}{0}[3]";
					$mpi++;
				}
				if ($orientation{$i}{0}[0] == 1 and $orientation{$i}{0}[1] == 1)
				{
					$result_mp[$mpi][0] = 'transl_inter';
					$result_mp[$mpi][1] = $i;
					$result_mp[$mpi][2] = $support{$i};
					$result_mp[$mpi][3] = $orientation{$i}{0}[2];
					$result_mp[$mpi][4] = $start{$i}{0}[1];
					$result_mp[$mpi][5] = $orientation{$i}{0}[0];
					$result_mp[$mpi][6] = $orientation{$i}{0}[3];
					$result_mp[$mpi][7] = $mstart{$i}{0}[1];
					$result_mp[$mpi][8] = $orientation{$i}{0}[1];
					$result_mp[$mpi][9] = "$cbp{$i}{0}[0]:$cbp{$i}{0}[1]:$cbp{$i}{0}[2]:$cbp{$i}{0}[3]";
					$mpi++;
				}
				if ($orientation{$i}{0}[0] == -1 and $orientation{$i}{0}[1] == -1)
				{
					$result_mp[$mpi][0] = 'transl_inter';
					$result_mp[$mpi][1] = $i;
					$result_mp[$mpi][2] = $support{$i};
					$result_mp[$mpi][3] = $orientation{$i}{0}[2];
					$result_mp[$mpi][4] = $start{$i}{0}[0];
					$result_mp[$mpi][5] = $orientation{$i}{0}[0];
					$result_mp[$mpi][6] = $orientation{$i}{0}[3];
					$result_mp[$mpi][7] = $mstart{$i}{0}[0];
					$result_mp[$mpi][8] = $orientation{$i}{0}[1];
					$result_mp[$mpi][9] = "$cbp{$i}{0}[0]:$cbp{$i}{0}[1]:$cbp{$i}{0}[2]:$cbp{$i}{0}[3]";
					$mpi++;
				}
			}
		}
	}
	#foreach (@result_mp){my $toprint = join ("\t", @{$_});unless ($$_[0] eq 'transl_inter' or $$_[0] eq 'del_inss' or $$_[0] eq 'del_inso' or $$_[0] eq 'inss' or $$_[0] eq 'inso'){print "$toprint\n";}}
	
	# call smaller event with same primary cluster id for intra-chr events
	my (%data, %cluster_event_map, %eventsize);
	# $data{event id} = complete event string
	# $cluster_event_map{$cluster_id} = array of event ids
	# $eventsize{event id} = event size
	my $mpai=0;
	foreach (@result_mp)
	{
		my @data = @$_;
		unless ($$_[0] eq 'transl_inter' or $$_[0] eq 'del_inss' or $$_[0] eq 'del_inso' or $$_[0] eq 'inss' or $$_[0] eq 'inso')
		{
			if ($data[1] =~ /\//)
			{
				my ($a, $b) = ($`, $');#'
				if ($a =~ /_/)
				{
					$a = $`;
				}
				if ($b =~ /_/)
				{
					$b = $`;
				}
				push @{$cluster_event_map{$a}}, $mpai;
				push @{$cluster_event_map{$b}}, $mpai;
				$eventsize{$mpai} += abs($data[6]);
				$eventsize{$mpai} += abs($data[10]);
			}
			else
			{
				my $a = $data[1];
				if ($a =~ /_/)
				{
					$a = $`;
				}
				push @{$cluster_event_map{$a}}, $mpai;
				$eventsize{$mpai} += abs($data[6]);
			}
			$data{$mpai} = join ("\t", @data);
			$mpai++;
		}
	}
	foreach my $cluster_id (%cluster_event_map)
	{
		if ($cluster_event_map{$cluster_id})
		{
			if (@{$cluster_event_map{$cluster_id}} > 1)
			{
				#print "$cluster_id\t@{$cluster_event_map{$cluster_id}}\n";
				my ($select_id, @delete_id);
				$select_id = $cluster_event_map{$cluster_id}[0];
				foreach (@{$cluster_event_map{$cluster_id}})
				{
					if ($eventsize{$_}<$eventsize{$select_id})
					{
						$select_id = $_;
					}
				}
				foreach (@{$cluster_event_map{$cluster_id}})
				{
					 if ($_ ne $select_id)
					 {
					 	push @delete_id, $_;
					 }
				}
				foreach (@delete_id)
				{
					#print "$cluster_id\t$select_id\t$_\n";
					my @current_event = split (/\t/, $data{$_});
					if ($current_event[1] =~ /\//)
					{
						if ($current_event[0] =~ /inssd/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my @temp_event;
							if ($cl_id1 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'del';
								$temp_event[1] = $cl_id2;
								$temp_event[2] = $mpd2;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[4];
								$temp_event[5] = $current_event[8];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[13];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "inssd\t$cluster_id\t$_\t$modified_event\n";
							}
							elsif ($cl_id2 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'tandem_dup';
								$temp_event[1] = $cl_id1;
								$temp_event[2] = $mpd1;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[5];
								$temp_event[5] = $current_event[9];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[12];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "inssd\t$cluster_id\t$_\t$modified_event\n";
							}
						}
						if ($current_event[0] =~ /inssu/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my @temp_event;
							if ($cl_id1 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'del';
								$temp_event[1] = $cl_id2;
								$temp_event[2] = $mpd2;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[9];
								$temp_event[5] = $current_event[5];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[13];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "inssu\t$cluster_id\t$_\t$modified_event\n";
							}
							elsif ($cl_id2 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'tandem_dup';
								$temp_event[1] = $cl_id1;
								$temp_event[2] = $mpd1;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[8];
								$temp_event[5] = $current_event[4];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[12];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "inssu\t$cluster_id\t$_\t$modified_event\n";
							}
							
						}
						if ($current_event[0] =~ /insod/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my @temp_event;
							if ($cl_id1 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'invers_r';
								$temp_event[1] = $cl_id2;
								$temp_event[2] = $mpd2;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[5];
								$temp_event[5] = $current_event[8];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[13];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "insod\t$cluster_id\t$_\t$modified_event\n";
							}
							elsif ($cl_id2 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'invers_f';
								$temp_event[1] = $cl_id1;
								$temp_event[2] = $mpd1;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[4];
								$temp_event[5] = $current_event[9];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[12];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "insod\t$cluster_id\t$_\t$modified_event\n";
							}
						}
						if ($current_event[0] =~ /insou/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my @temp_event;
							if ($cl_id1 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'invers_r';
								$temp_event[1] = $cl_id2;
								$temp_event[2] = $mpd2;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[8];
								$temp_event[5] = $current_event[5];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[13];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "insou\t$cluster_id\t$_\t$modified_event\n";
							}
							elsif ($cl_id2 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'invers_f';
								$temp_event[1] = $cl_id1;
								$temp_event[2] = $mpd1;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[9];
								$temp_event[5] = $current_event[4];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[12];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "insou\t$cluster_id\t$_\t$modified_event\n";
							}
						}
						if ($current_event[0] eq 'del_invers')
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my @temp_event;
							if ($cl_id1 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'invers_r';
								$temp_event[1] = $cl_id2;
								$temp_event[2] = $mpd2;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[5];
								$temp_event[5] = $current_event[9];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[11];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "del_invers\t$cluster_id\t$_\t$modified_event\n";
							}
							elsif ($cl_id2 =~ /$cluster_id\_/)
							{
								$temp_event[0] = 'invers_f';
								$temp_event[1] = $cl_id1;
								$temp_event[2] = $mpd1;
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[4];
								$temp_event[5] = $current_event[8];
								$temp_event[6] = $temp_event[5] - $temp_event[4];
								$temp_event[7] = $current_event[10];
								my $modified_event = join("\t", @temp_event);
								$data{$_} = $modified_event;
								#print STDERR "del_invers\t$cluster_id\t$_\t$modified_event\n";
							}
						}
					}
					else
					{
						delete $data{$_};
					}
				}
			}
		}
	}
	
	open MPINTRA, ">$mpintra_outfile";
	for (my $i=0;$i<=$mpai;$i++)
	{
		print MPINTRA "$data{$i}\n" if ($data{$i});
	}
	close MPINTRA;
	
	# call smaller event if same cluster is used for inter-chr events
	my (%data, %cluster_event_map, %eventsize);
	# $data{event id} = complete event string
	# $cluster_event_map{$cluster_id} = array of event ids
	# $eventsize{event id} = event size
	my $mpei=0;
	foreach (@result_mp)
	{
		my @data = @$_;
		if ($$_[0] eq 'transl_inter' or $$_[0] eq 'del_inss' or $$_[0] eq 'del_inso' or $$_[0] eq 'inss' or $$_[0] eq 'inso')
		{
			if ($data[1] =~ /\//)
			{
				my ($a, $b) = ($`, $');#'
				if ($a =~ /_/)
				{
					$a = $`;
				}
				if ($b =~ /_/)
				{
					$b = $`;
				}
				push @{$cluster_event_map{$a}}, $mpei;
				push @{$cluster_event_map{$b}}, $mpei;
				$eventsize{$mpei} += abs($data[6]);
				$eventsize{$mpei} += abs($data[10]);
			}
			else
			{
				my $a = $data[1];
				if ($a =~ /_/)
				{
					$a = $`;
				}
				push @{$cluster_event_map{$a}}, $mpei;
				$eventsize{$mpei} = 1000000000;
			}
			$data{$mpei} = join ("\t", @data);
			$mpei++;
		}
	}
	foreach my $cluster_id (%cluster_event_map)
	{
		if ($cluster_event_map{$cluster_id})
		{
			if (@{$cluster_event_map{$cluster_id}} > 1)
			{
				#print "$cluster_id\t@{$cluster_event_map{$cluster_id}}\n";
				my ($select_id, @delete_id);
				$select_id = $cluster_event_map{$cluster_id}[0];
				foreach (@{$cluster_event_map{$cluster_id}})
				{
					if ($eventsize{$_}<$eventsize{$select_id})
					{
						$select_id = $_;
					}
				}
				foreach (@{$cluster_event_map{$cluster_id}})
				{
					 if ($_ ne $select_id)
					 {
					 	push @delete_id, $_;
					 }
				}
				foreach (@delete_id)
				{
					my @current_event = split (/\t/, $data{$_});
					if ($current_event[1] =~ /\//)
					{
						if ($current_event[0] =~ /inss/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my @temp_event;
							$temp_event[0] = 'transl_inter';
							$temp_event[1] = $cl_id1;
							$temp_event[2] = $mpd1;
							if ($current_event[3] lt $current_event[7])
							{
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[4];
								$temp_event[5] = 1;
								$temp_event[6] = $current_event[7];
								$temp_event[7] = $current_event[8];
								$temp_event[8] = -1;
								$temp_event[9] = $current_event[11];
							}
							else
							{
								$temp_event[3] = $current_event[7];
								$temp_event[4] = $current_event[9];
								$temp_event[5] = 1;
								$temp_event[6] = $current_event[3];
								$temp_event[7] = $current_event[5];
								$temp_event[8] = -1;
								my @tempcbp = split (":", $current_event[11]);
								$temp_event[9] = $tempcbp[2].':'.$tempcbp[3].':'.$tempcbp[0].':'.$tempcbp[1];
							}
							my $modified_event = join("\t", @temp_event);
							$data{$_} = $modified_event;
							#print "inss $modified_event\n";
						}
						if ($current_event[0] =~ /inso/)
						{
							my ($cl_id1, $cl_id2) = split (/\//, $current_event[1]);
							my ($mpd1, $mpd2) = split (/\//, $current_event[2]);
							my @temp_event;
							$temp_event[0] = 'transl_inter';
							$temp_event[1] = $cl_id1;
							$temp_event[2] = $mpd1;
							if ($current_event[3] lt $current_event[7])
							{
								$temp_event[3] = $current_event[3];
								$temp_event[4] = $current_event[4];
								$temp_event[5] = 1;
								$temp_event[6] = $current_event[7];
								$temp_event[7] = $current_event[9];
								$temp_event[8] = 1;
								$temp_event[9] = $current_event[11];
							}
							else
							{
								$temp_event[3] = $current_event[7];
								$temp_event[4] = $current_event[9];
								$temp_event[5] = 1;
								$temp_event[6] = $current_event[3];
								$temp_event[7] = $current_event[4];
								$temp_event[8] = 1;
								my @tempcbp = split (":", $current_event[11]);
								$temp_event[9] = $tempcbp[2].':'.$tempcbp[3].':'.$tempcbp[0].':'.$tempcbp[1];
							}
							my $modified_event = join("\t", @temp_event);
							$data{$_} = $modified_event;
							#print "inso $modified_event\n";
						}
					}
					else
					{
						delete $data{$_};
					}
				}
			}
		}
	}
		
	open MPINTER, ">$mpinter_outfile";
	for (my $i=0;$i<=$mpei;$i++)
	{
		print MPINTER "$data{$i}\n" if ($data{$i});
	}
	close MPINTER;
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

sub overlap
{
	my $a1 = shift;
	my $a2 = shift;
	my $b1 = shift;
	my $b2 = shift;
	return ($a1, $a2) unless (defined($b1) and defined($b2));
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

sub print_usage
{
	die "Usage:\nperl ./scripts/meerkat.pl [options]
	-b FILE	sorted and indexed bam file, required
	-k INT	[0/1], use black list generated in pre_process.pl, default 1
	-d FLT	standard deviation cutoff to call discordant read pairs, default 3
	-c FLT	standard deviation cutoff to cluster discordant read pairs, default equal to -d, it is recommended to use same -c as -d if -d<=5, -c 5 if -d > 5
	-p INT	number of supporting read pairs required for an event to be called, default 2
	-o INT	number of supporting full length read pairs, default 0, specify this option will decrease sensitivity on small complex events
	-q INT	number of supporting split reads required, default 1
	-z INT	event size cutoff, default 1,000,000,000
	-s INT	bp to be cut off from beginning and end of unmapped reads, must be same as -s in pre_process.pl, default 20
	-m INT	[0/1], if set to 1, use Meerkat to remove duplicates; if set to 0, use flag 'd' marked by Picard or other tools to remove duplicates. If the bam file is aligned by bwa mem, it has to be processed by Picard to mark duplicates and use 0 option. The bwa mem aligned bam file won't work with option 1. Default 1
	-a INT	[0/1], adjust non-uniq mapped reads, default 1
	-u INT  [0/1], use all alignments in the BAM file, turn this option on if the BAM file is not generated by BWA, turn on this option will force turning off option a, default 0
	-Q INT	minimum mapping quality for reads to be used, default 0
	-g INT	number of alternative mappings to consider in main bam file, number of alternative mappings printed out in XA tag by bwa is controlled by -N, default use all in bam file
	-f INT	number of alternative mappings to consider in clipped alignments, default use all in bam file
	-l INT	[0/1], consider clipped alignments, default 1
	-t INT	number of threads used in bwa alignment, default 1
	-R STR	file name of read group to be ignored, one read group ID per line
	-F STR	/path/to/reference/fasta/files, path only, not the files, required
	-S STR	/path/to/samtools, path only, not the command, no need to specify if samtools is in PATH
	-W STR	/path/to/bwa, path only, not the command, no need to specify if bwa is in PATH
	-B STR	/path/to/blastall and formatdb, path only, not the command, no need to specify if blastall and formatdb is in PATH
	-P STR	specify step to run, dc|cl|mpd|alg|srd|rf|all, default 'all', each step require results from previous steps
		dc: extract discordant read pairs
		cl: construct clusters of discordant read pairs
		mpd: call events based on read pairs
		alg: align split reads to candidate break point regions
		srd: confirm events based on split reads and filter results
		rf:  refine break points by local alignments
		all: run all above steps
	-h help\n";
}
