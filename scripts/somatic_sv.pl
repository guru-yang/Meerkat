#!/usr/bin/perl
# somatic_sv.pl
# cross check tumor with all normal discordant clusters
# filter SR
# filter by homology level
# filter complex events from single cluster
# filter by event size for single cluster events
# split complex event by insertion size or distance
# filter none smallest single events
# filter by matched normal bam file, number of discordant pairs (uniq mapped, ignore too many mismatch reads) at SV break points based on orientation of tumor call in normal should be less than specified
# filter by matched normal bam file, ratio of non-uniq mapped read should be less than specified
# filter by matched normal bam file, number of soft-clipped reads in defined window must be less than specified, count q15 base only
# filter by tumor bam file, there must be at least one read pair with both end uniquely mappable in cluster
# filter by number of supporting read pairs (p), number of supporting split reads (q) or the sum of of those two (P)
# filter by blat results
# call mechanism for leftover events
# more TEI events after recall mechenisms, need to manually remove
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, USA
# Email: lixing_yang@hms.harvard.edu

use strict;
use FindBin '$Bin';
use Getopt::Std;
my $version = 'v.0.189';

my %opts = (F=>'null', D=>3, l=>500, x=>2, s=>100, E=>1, d=>100, t=>100, m=>1, n=>0, y=>0.1, u=>0, v=>100, r=>0.25, f=>0, g=>10, j=>3, e=>0, k=>800, a=>0, z=>0, p=>3, q=>1, P=>6, M=>0, N=>2000, Q=>0, R=>'/db/hg18/rmsk-hg18.txt');
getopts("i:o:D:I:K:F:l:x:s:E:d:t:m:n:y:u:v:r:f:g:j:e:k:B:a:c:z:p:q:P:M:N:C:Q:b:S:R:L:V:T:A:h", \%opts);

my $inputfile = $opts{i};
my $outputfile = $opts{o};
my $sd_cutoff = $opts{D};
my $isinfofile = $opts{I};
my $blackrgfile = $opts{K};
my $tmpdir = 'tmp_'.int(rand()*1000000);


# master normal
my $normal_disc_folder = $opts{F};
my $distance = $opts{l}; # distance of tumor and normal break point cutoff
my $number_discord = $opts{x}; # number of discordant read pairs in a cluster to be filtered
$normal_disc_folder .= '/' unless ($normal_disc_folder =~ /\/$/);

# SV size
my $sv_min = $opts{s}; # minimum SV size for simple events
my $di_min = 50; # minimum distance or insertion size for complex events

# TE
my $filter_te = $opts{E};

# homology
my $del_homology = $opts{d}; # max homology allowed for deletion and intra-chr events, 0 to disable
my $inter_transl_homology = $opts{t}; # max homology allowed for inter chromosomal translocation events, 0 to disable

# SR
my $filter_sr = $opts{m}; # filter satellite, simple repeats

# discordant in matched normal
my $bamfile = $opts{B};
my $filter_discord = $opts{n};
my $match_normal_disc = $opts{y}; # number of discordant pairs in normal should be no more than this value

# non-uniq mapped reads in matched normal
my $filter_non_uniq = $opts{u};
my $nu_window = $opts{v};
my $nu_ratio = $opts{r}; # ratio cutoff of fraction of non-uniq mapped reads out of mappable reads

# soft-clipped reads in matched normal
my $filter_soft_clip = $opts{f};
my $sc_window = $opts{g};
my $sc_cutoff = $opts{j};
my $clip_size_cutoff = 5; # soft clipped has to be at least this number of bp
my $qual_code = 33; # sanger 33-73
my $qual_cutoff = 15; # q15 base

# discordant in tumor
my $disc_tumor = $opts{e};
my $dt_window = $opts{k};

# tie on non-uniq mapped clusters
my $no_tie = $opts{a};
my $complex_sv = $no_tie; # filter complex SV as well
my $discordfile = $opts{c};

# filter out low confidence calls
my $no_lc = $opts{z};
my $mpd = $opts{p};
my $srd = $opts{q};
my $mpsrd = $opts{P};

# filter by mate of sr read
my $filter_mate = $opts{M};
my $filter_mate_window = $opts{N};
my $bp_readsfile = $opts{C};
my $min_mapq = $opts{Q};
my %readslist;

# filter by minor allele freq
my $alle_freq_cutoff = $opts{b};

&print_usage if (defined($opts{h}));
die "Please specify inputfile\n" unless (defined($inputfile));
die "Please specify outputfile\n" unless (defined($outputfile));
die "The inputfile must be .variants file\n" unless ($inputfile =~ /.variants$/);
die "The outputfile must be .variants file\n" unless ($outputfile =~ /.variants$/);
if ($opts{n} or $opts{e})
{
	die "Please specify isinfo file\n" unless (defined($isinfofile));
	die "$isinfofile file not exist\n" unless (-e $isinfofile);
}

my $time0 = time;
my $local0 = localtime($time0);
print STDERR "$local0 somatic_sv.pl $version started\n";

my ($newline, %is);
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
		#print "$rg\t$is{$rg}{'isu'}\n";
	}
}
close FILE;
my $window_size = $is{'isu'}; # upstream and downstream of break points

my %blackrg;
open FILE, "<$blackrgfile";
while ($newline = <FILE>)
{
	chomp $newline;
	$blackrg{$newline} = 1;
}
close FILE;

if ($filter_mate and $bamfile)
{
	open BPREADSRF, "<$bp_readsfile";
	while ($newline = <BPREADSRF>)
	{
		chomp $newline;
		my $bpname = $newline;
		$newline = <BPREADSRF>;
		chomp $newline;
		my @reads_name = split (/\t/, $newline);
		@{$readslist{$bpname}} = @reads_name;
	}
	close BPREADSRF;
}
	
my $samtools_path = $opts{S};
my $rmsk = $opts{R};
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

my $blat_path = $opts{L};
my $gfClient_command;
if (defined($opts{L}))
{
	$blat_path .= '/' unless ($blat_path =~ /\/$/);
	$gfClient_command = $blat_path.'gfClient';
	system "mkdir $tmpdir" unless (-e $tmpdir);
}
else
{
	$gfClient_command = 'gfClient';
}
my $blat_server = $opts{V};
my $blat_port = $opts{T};
my $srout_dir = $opts{A};
if (defined($opts{A}))
{
	$srout_dir .= '/' unless ($srout_dir =~ /\/$/);
}

$bamfile =~ s/ /\\ / unless ($bamfile =~ /\\ /);

opendir(DIR, $normal_disc_folder);
my @samplelist=readdir(DIR);
closedir(DIR);

my %normal_cluster;
# $normal_cluster{sample}{chr1}{chr2}{clustertype} is an array of cluster
foreach (@samplelist)
{
	next if ($_ =~ m/^\./);
	my $normal_variant = $normal_disc_folder.$_;
	#print STDERR "Reading clusters in $normal_variant\n";
	open FILE, "<$normal_variant";
	while ($newline = <FILE>)
	{
		chomp $newline;
		my @data = split ('\t', $newline);
		#print "@data\n";
		next if ($data[2] < $number_discord);
		if ($data[5] == 1 and $data[8] == -1)
		{
			push @{$normal_cluster{$normal_variant}{$data[3]}{$data[6]}{0}}, \@data;
		}
		if ($data[5] == -1 and $data[8] == 1)
		{
			push @{$normal_cluster{$normal_variant}{$data[3]}{$data[6]}{1}}, \@data;
		}
		if ($data[5] == 1 and $data[8] == 1)
		{
			push @{$normal_cluster{$normal_variant}{$data[3]}{$data[6]}{2}}, \@data;
		}
		if ($data[5] == -1 and $data[8] == -1)
		{
			push @{$normal_cluster{$normal_variant}{$data[3]}{$data[6]}{3}}, \@data;
		}
	}
	close FILE;
}

my %discord;
if ($no_tie)
{
	open FILE, "cat $discordfile|cut -f1,2|";
	while ($newline = <FILE>)
	{
		chomp $newline;
		my @data = split ('\t', $newline);
		$discord{$data[0]} = $data[1];
	}
	close FILE;
}

my $tumor_out = $outputfile.'.temp';
open FILE, "<$inputfile";
open TUMOR, ">$tumor_out";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split ('\t', $newline);
	my @data1 = @data;
	my @rp;
	if ($data[-1] =~ /RP:/)
	{
		@rp = split ('\/', $');#'
		delete $data[-1];
		
	}
	if ($data[-1] =~ /BP:/)
	{
		delete $data[-1];
	}
	#print "@rp\n@data\n";
	my $found;
	if ($data[0] eq 'del' or $data[0] eq 'del_ins')
	{
		my $hm;
		$hm = $data[9] if ($data[0] eq 'del');
		$hm = $data[12] if ($data[0] eq 'del_ins');
		if ($filter_sr)
		{
			$found = 1 if ($newline =~ /SR_SR/);
		}
		my $cluster_id = $data[2];
		$found = 1 if ($newline =~ /TEI/ and $filter_te);
		$found = 1 if ($del_homology and ($data[9]>$del_homology or $data[12]>$del_homology));
		$found = 1 if ($data[8] < $sv_min);
		$found = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
		$found = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[5], $data[7], -1, $hm));
		$found = 1 if ($no_lc and $data[3] < $mpd);
		$found = 1 if ($no_lc and $data[4] < $srd);
		$found = 1 if ($no_lc and ($data[3] + $data[4] < $mpsrd));
		if (defined($blat_server) and defined($blat_port))
		{
			$found = 1 unless (&local_align ($data[2], $data[5], $data[6], $data[5], $data[7]));
		}
		if ($filter_mate)
		{
			$found = 1 unless (&mate($data[2], $data[5], $data[6], $data[5], $data[7]));
		}
		if ($data[2] =~ /_/)
		{
			my $smallest = $';#'
			$cluster_id = $`;
			$found = 1 if ($smallest > 0);
		}
		$found = 1 if ($no_tie and $discord{$cluster_id});
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord = $rp1[2];
			$max_concord = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq = $rp1[0]/($rp1[0]+$max_concord);
			$found = 1 if ($alle_freq < $alle_freq_cutoff);
		}
LOOP1:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP1 if ($found);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{0}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[6] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
				{
					$found = 1;
					last LOOP1;
				}
			}
		}
		unless ($found)
		{
			my $eventtype = $data[0];
			shift(@data); shift(@data);
			unshift(@data, $eventtype);
			my $toprint = join ("\t", @data);
			$toprint .= "\tRP:".$rp[0] if ($rp[0]);
			print TUMOR "$toprint\n";
		}
	}
	if ($data[0] eq 'tandem_dup')
	{
		my $hm = $data[9];
		if ($filter_sr)
		{
			$found = 1 if ($newline =~ /SR_SR/);
		}
		$found = 1 if ($del_homology and abs($data[9])>$del_homology);
		$found = 1 if ($data[8] < $sv_min);
		$found = 1 if ($filter_discord and &dcn($data[5], $data[6], -1));
		$found = 1 if ($filter_discord and &dcn($data[5], $data[7], 1));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found = 1 if ($disc_tumor and &dt($data[5], $data[6], -1, $data[5], $data[7], 1, $hm));
		$found = 1 if ($no_lc and $data[3] < $mpd);
		$found = 1 if ($no_lc and $data[4] < $srd);
		$found = 1 if ($no_lc and ($data[3] + $data[4] < $mpsrd));
		my $cluster_id = $data[2];
		if (defined($blat_server) and defined($blat_port))
		{
			$found = 1 unless (&local_align ($data[2], $data[5], $data[6], $data[5], $data[7]));
		}
		if ($filter_mate)
		{
			$found = 1 unless (&mate($data[2], $data[5], $data[6], $data[5], $data[7]));
		}
		if ($data[2] =~ /_/)
		{
			my $smallest = $';#'
			$cluster_id = $`;
			$found = 1 if ($smallest > 0);
		}
		$found = 1 if ($no_tie and $discord{$cluster_id});
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord = $rp1[2];
			$max_concord = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq = $rp1[0]/($rp1[0]+$max_concord);
			$found = 1 if ($alle_freq < $alle_freq_cutoff);
		}
LOOP2:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP2 if ($found);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{1}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[6] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
				{
					$found = 1;
					last LOOP2;
				}
			}
		}
		unless ($found)
		{
			my $eventtype = $data[0];
			shift(@data); shift(@data);
			unshift(@data, $eventtype);
			my $toprint = join ("\t", @data);
			$toprint .= "\tRP:".$rp[0] if ($rp[0]);
			print TUMOR "$toprint\n";
		}
	}
	if ($data[0] eq 'invers')
	{
		my @hm = split ("\/", $data[9]);
		if ($filter_sr)
		{
			$found = 1 if ($newline =~ /SR_SR/);
		}
		my @hm = split ("\/", $data[9]);
		$found = 1 if ($del_homology and (abs($hm[0])>$del_homology or abs($hm[1])>$del_homology));
		$found = 1 if ($data[8] < $sv_min);
		$found = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
		$found = 1 if ($filter_discord and &dcn($data[5], $data[7], 1));
		$found = 1 if ($filter_discord and &dcn($data[5], $data[6], -1));
		$found = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[5], $data[7], 1, $hm[0]));
		$found = 1 if ($disc_tumor and &dt($data[5], $data[6], -1, $data[5], $data[7], -1, $hm[0]));
		my @mpd = split ("\/", $data[3]);
		my @srd = split ("\/", $data[4]);
		if (defined($blat_server) and defined($blat_port))
		{
			$found = 1 unless (&local_align ($mpd[0], $data[5], $data[6], $data[5], $data[7]));
		}
		if ($filter_mate)
		{
			$found = 1 unless (&mate($data[2], $data[5], $data[6], $data[5], $data[7]));
		}
		$found = 1 if ($no_lc and $mpd[0] < $mpd);
		$found = 1 if ($no_lc and $srd[0] < $srd);
		$found = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
		unless ($data[2] =~ /\//)
		{
			$found = 1;
		}
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord1 = $rp1[2];
			$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
			$found = 1 if ($alle_freq1 < $alle_freq_cutoff);
			if ($rp[1])
			{
				my @rp2 = split ('_', $rp[1]);
				my $max_concord2 = $rp2[2];
				$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
				my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
				$found = 1 if ($alle_freq2 < $alle_freq_cutoff);
			}
		}
LOOP27:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP27 if ($found);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{2}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[6] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
				{
					$found = 1;
					last LOOP27;
				}
			}
		}
LOOP28:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP28 if ($found);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{3}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[6] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
				{
					$found = 1;
					last LOOP28;
				}
			}
		}
		unless ($found)
		{
			my $eventtype = $data[0];
			shift(@data); shift(@data);
			unshift(@data, $eventtype);
			my $toprint = join ("\t", @data);
			$toprint .= "\tRP:".$rp[0] if ($rp[0]);
			if ($rp[1])
			{
				$toprint .= "/".$rp[1];
			}
			print TUMOR "$toprint\n";
		}
	}
	if ($data[0] eq 'invers_f')
	{
		my $hm = $data[9];
		if ($filter_sr)
		{
			$found = 1 if ($newline =~ /SR_SR/);
		}
		$found = 1 if ($del_homology and abs($data[9])>$del_homology);
		$found = 1 if ($data[8] < $sv_min);
		$found = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
		$found = 1 if ($filter_discord and &dcn($data[5], $data[7], 1));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[5], $data[7], 1, $hm));
		$found = 1 if ($no_lc and $data[3] < $mpd);
		$found = 1 if ($no_lc and $data[4] < $srd);
		$found = 1 if ($no_lc and ($data[3] + $data[4] < $mpsrd));
		my $cluster_id = $data[2];
		if (defined($blat_server) and defined($blat_port))
		{
			$found = 1 unless (&local_align ($data[2], $data[5], $data[6], $data[5], $data[7]));
		}
		if ($filter_mate)
		{
			$found = 1 unless (&mate($data[2], $data[5], $data[6], $data[5], $data[7]));
		}
		if ($data[2] =~ /_/)
		{
			my $smallest = $';#'
			$cluster_id = $`;
			$found = 1 if ($smallest > 0);
		}
		$found = 1 if ($no_tie and $discord{$cluster_id});
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord = $rp1[2];
			$max_concord = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq = $rp1[0]/($rp1[0]+$max_concord);
			$found = 1 if ($alle_freq < $alle_freq_cutoff);
		}
LOOP3:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP3 if ($found);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{2}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[6] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
				{
					$found = 1;
					last LOOP3;
				}
			}
		}
		unless ($found)
		{
			my $eventtype = $data[0];
			shift(@data); shift(@data);
			unshift(@data, $eventtype);
			my $toprint = join ("\t", @data);
			$toprint .= "\tRP:".$rp[0] if ($rp[0]);
			print TUMOR "$toprint\n";
		}
	}
	if ($data[0] eq 'invers_r')
	{
		my $hm = $data[9];
		if ($filter_sr)
		{
			$found = 1 if ($newline =~ /SR_SR/);
		}
		$found = 1 if ($del_homology and abs($data[9])>$del_homology);
		$found = 1 if ($data[8] < $sv_min);
		$found = 1 if ($filter_discord and &dcn($data[5], $data[6], -1));
		$found = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found = 1 if ($disc_tumor and &dt($data[5], $data[6], -1, $data[5], $data[7], -1, $hm));
		$found = 1 if ($no_lc and $data[3] < $mpd);
		$found = 1 if ($no_lc and $data[4] < $srd);
		$found = 1 if ($no_lc and ($data[3] + $data[4] < $mpsrd));
		my $cluster_id = $data[2];
		if (defined($blat_server) and defined($blat_port))
		{
			$found = 1 unless (&local_align ($data[2], $data[5], $data[6], $data[5], $data[7]));
		}
		if ($filter_mate)
		{
			$found = 1 unless (&mate($data[2], $data[5], $data[6], $data[5], $data[7]));
		}
		if ($data[2] =~ /_/)
		{
			my $smallest = $';#'
			$cluster_id = $`;
			$found = 1 if ($smallest > 0);
		}
		$found = 1 if ($no_tie and $discord{$cluster_id});
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord = $rp1[2];
			$max_concord = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq = $rp1[0]/($rp1[0]+$max_concord);
			$found = 1 if ($alle_freq < $alle_freq_cutoff);
		}
LOOP4:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP4 if ($found);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{3}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[6] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
				{
					$found = 1;
					last LOOP4;
				}
			}
		}
		unless ($found)
		{
			my $eventtype = $data[0];
			shift(@data); shift(@data);
			unshift(@data, $eventtype);
			my $toprint = join ("\t", @data);
			$toprint .= "\tRP:".$rp[0] if ($rp[0]);
			print TUMOR "$toprint\n";
		}
	}
	my ($found1, $found2);
	if ($data[0] =~ /inssd/)
	{
		my @hm = split ("\/", $data[14]);
		my $sr = substr ($data1[15], 3);
		my @sr = split ("_", $sr);
		if ($filter_sr)
		{
			$found1 = 1 if ($sr[1] eq 'SR' and $sr[3] eq 'SR');
			$found2 = 1 if ($sr[0] eq 'SR' and $sr[2] eq 'SR');
		}
		my @hm = split ("\/", $data[14]);
		$found1 = 1 if ($del_homology and abs($hm[0])>$del_homology);
		$found2 = 1 if ($del_homology and abs($hm[1])>$del_homology);
		$found1 = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
		$found1 = 1 if ($filter_discord and &dcn($data[9], $data[11], 1));
		$found2 = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
		$found2 = 1 if ($filter_discord and &dcn($data[9], $data[10], -1));
		$found1 = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found1 = 1 if ($filter_non_uniq and &nud($data[9], $data[11]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[9], $data[10]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[9], $data[11]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[9], $data[10]));
		$found1 = 1 if ($disc_tumor and &dt($data[5], $data[7], -1, $data[9], $data[11], 1, $hm[0]));
		$found2 = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[9], $data[10], -1, $hm[1]));
		my @mpd = split ("\/", $data[3]);
		my @srd = split ("\/", $data[4]);
		if (defined($blat_server) and defined($blat_port))
		{
			$found1 = 1 unless (&local_align ($mpd[0], $data[5], $data[7], $data[9], $data[11]));
			$found2 = 1 unless (&local_align ($mpd[1], $data[5], $data[6], $data[9], $data[10]));
		}
		if ($filter_mate)
		{
			$found1 = 1 unless (&mate($mpd[0], $data[5], $data[7], $data[9], $data[11]));
			$found2 = 1 unless (&mate($mpd[1], $data[5], $data[6], $data[9], $data[10]));
		}
		$found1 = 1 if ($no_lc and $mpd[0] < $mpd);
		$found1 = 1 if ($no_lc and $srd[0] < $srd);
		$found1 = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
		$found2 = 1 if ($no_lc and $mpd[1] < $mpd);
		$found2 = 1 if ($no_lc and $srd[1] < $srd);
		$found2 = 1 if ($no_lc and ($mpd[1] + $srd[1] < $mpsrd));
		if ($complex_sv)
		{
			my @cluster_id = split ("\/", $data[2]);
			if ($cluster_id[0] =~ /_/)
			{
				$cluster_id[0] = $`;
			}
			$found1 = 1 if ($no_tie and $discord{$cluster_id[0]});
			if ($cluster_id[1] =~ /_/)
			{
				$cluster_id[1] = $`;
			}
			$found2 = 1 if ($no_tie and $discord{$cluster_id[1]});
		}
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord1 = $rp1[2];
			$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
			$found1 = 1 if ($alle_freq1 < $alle_freq_cutoff);
			my @rp2 = split ('_', $rp[1]);
			my $max_concord2 = $rp2[2];
			$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
			my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
			$found2 = 1 if ($alle_freq2 < $alle_freq_cutoff);
		}
LOOP5:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP5 if ($found1);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{1}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[7] - $$cluster[4]) < $distance and abs($data[11] - $$cluster[7]) < $distance)
				{
					$found1 = 1;
					last LOOP5;
				}
			}
		}
LOOP6:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP6 if ($found2);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{0}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[6] - $$cluster[4]) < $distance and abs($data[10] - $$cluster[7]) < $distance)
				{
					$found2 = 1;
					last LOOP6;
				}
			}
		}
		if ($sv_min > 0)
		{
			if ($data[12] < $di_min or $data[13] < $di_min)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found2)
						{
							my @tempdata;
							$tempdata[0] = 'del';
							$tempdata[1] = $id2;
							$tempdata[2] = $mpd2;
							$tempdata[3] = $srd2;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[6];
							$tempdata[6] = $data[10];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							if ($hm2 >= 0)
							{
								$tempdata[8] = $hm2;
							}
							else
							{
								$tempdata[0] = 'del_ins';
								$tempdata[8] = '-';
								$tempdata[9] = '-';
								$tempdata[10] = '-';
								$tempdata[11] = -$hm2;
							}
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[1] if ($rp[1]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found1)
						{
							my @tempdata;
							$tempdata[0] = 'tandem_dup';
							$tempdata[1] = $id1;
							$tempdata[2] = $mpd1;
							$tempdata[3] = $srd1;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[7];
							$tempdata[6] = $data[11];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							$tempdata[8] = $hm1;
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[0] if ($rp[0]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				next;
			}
		}
		unless ($data[2] =~ /\//)
		{
			$found1 = 1;
			$found2 = 1;
		}
		if ($found1)
		{
			# both found, ignore whole event
			if ($found2)
			{
				
			}
			# -+ found
			else
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'del';
				$tempdata[1] = $id2;
				$tempdata[2] = $mpd2;
				$tempdata[3] = $srd2;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[6];
				$tempdata[6] = $data[10];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				if ($hm2 >= 0)
				{
					$tempdata[8] = $hm2;
				}
				else
				{
					$tempdata[0] = 'del_ins';
					$tempdata[8] = '-';
					$tempdata[9] = '-';
					$tempdata[10] = '-';
					$tempdata[11] = -$hm2;
				}
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[1] if ($rp[1]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
		}
		else
		{
			# +- found
			if ($found2)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'tandem_dup';
				$tempdata[1] = $id1;
				$tempdata[2] = $mpd1;
				$tempdata[3] = $srd1;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[7];
				$tempdata[6] = $data[11];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				$tempdata[8] = $hm1;
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[0] if ($rp[0]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
			# both not found, keep whole event
			else
			{
				my $eventtype = $data[0];
				shift(@data); shift(@data);
				unshift(@data, $eventtype);
				my $toprint = join ("\t", @data);
				$toprint .= "\tRP:".$rp[0].'/'.$rp[1] if ($rp[0] and $rp[1]);
				print TUMOR "$toprint\n";
			}
		}
	}
	if ($data[0] =~ /inssu/)
	{
		my @hm = split ("\/", $data[14]);
		my $sr = substr ($data1[15], 3);
		my @sr = split ("_", $sr);
		if ($filter_sr)
		{
			$found1 = 1 if ($sr[0] eq 'SR' and $sr[2] eq 'SR');
			$found2 = 1 if ($sr[1] eq 'SR' and $sr[3] eq 'SR');
		}
		my @hm = split ("\/", $data[14]);
		$found1 = 1 if ($del_homology and abs($hm[0])>$del_homology);
		$found2 = 1 if ($del_homology and abs($hm[1])>$del_homology);
		$found1 = 1 if ($filter_discord and &dcn($data[9], $data[10], -1));
		$found1 = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
		$found2 = 1 if ($filter_discord and &dcn($data[9], $data[11], 1));
		$found2 = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
		$found1 = 1 if ($filter_non_uniq and &nud($data[9], $data[10]));
		$found1 = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[9], $data[11]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[9], $data[10]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[9], $data[11]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found1 = 1 if ($disc_tumor and &dt($data[9], $data[10], -1, $data[5], $data[6], 1, $hm[0]));
		$found2 = 1 if ($disc_tumor and &dt($data[9], $data[11], 1, $data[5], $data[7], -1, $hm[1]));
		my @mpd = split ("\/", $data[3]);
		my @srd = split ("\/", $data[4]);
		if (defined($blat_server) and defined($blat_port))
		{
			$found1 = 1 unless (&local_align ($mpd[0], $data[5], $data[6], $data[9], $data[10]));
			$found2 = 1 unless (&local_align ($mpd[1], $data[5], $data[7], $data[9], $data[11]));
		}
		if ($filter_mate)
		{
			$found1 = 1 unless (&mate($mpd[0], $data[5], $data[6], $data[9], $data[10]));
			$found2 = 1 unless (&mate($mpd[1], $data[5], $data[7], $data[9], $data[11]));
		}
		$found1 = 1 if ($no_lc and $mpd[0] < $mpd);
		$found1 = 1 if ($no_lc and $srd[0] < $srd);
		$found1 = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
		$found2 = 1 if ($no_lc and $mpd[1] < $mpd);
		$found2 = 1 if ($no_lc and $srd[1] < $srd);
		$found2 = 1 if ($no_lc and ($mpd[1] + $srd[1] < $mpsrd));
		if ($complex_sv)
		{
			my @cluster_id = split ("\/", $data[2]);
			if ($cluster_id[0] =~ /_/)
			{
				$cluster_id[0] = $`;
			}
			$found1 = 1 if ($no_tie and $discord{$cluster_id[0]});
			if ($cluster_id[1] =~ /_/)
			{
				$cluster_id[1] = $`;
			}
			$found2 = 1 if ($no_tie and $discord{$cluster_id[1]});
		}
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord1 = $rp1[2];
			$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
			$found1 = 1 if ($alle_freq1 < $alle_freq_cutoff);
			my @rp2 = split ('_', $rp[1]);
			my $max_concord2 = $rp2[2];
			$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
			my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
			$found2 = 1 if ($alle_freq2 < $alle_freq_cutoff);
		}
LOOP7:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP7 if ($found1);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{1}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[10] - $$cluster[4]) < $distance and abs($data[6] - $$cluster[7]) < $distance)
				{
					$found1 = 1;
					last LOOP7;
				}
			}
		}
LOOP8:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP8 if ($found2);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{0}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[11] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
				{
					$found2 = 1;
					last LOOP8;
				}
			}
		}
		unless ($data[2] =~ /\//)
		{
			$found1 = 1;
			$found2 = 1;
		}
		if ($sv_min > 0)
		{
			if ($data[12] < $di_min or $data[13] < $di_min)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found2)
						{
							my @tempdata;
							$tempdata[0] = 'del';
							$tempdata[1] = $id2;
							$tempdata[2] = $mpd2;
							$tempdata[3] = $srd2;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[11];
							$tempdata[6] = $data[7];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							if ($hm2 >= 0)
							{
								$tempdata[8] = $hm2;
							}
							else
							{
								$tempdata[0] = 'del_ins';
								$tempdata[8] = '-';
								$tempdata[9] = '-';
								$tempdata[10] = '-';
								$tempdata[11] = -$hm2;
							}
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[1] if ($rp[1]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found1)
						{
							my @tempdata;
							$tempdata[0] = 'tandem_dup';
							$tempdata[1] = $id1;
							$tempdata[2] = $mpd1;
							$tempdata[3] = $srd1;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[10];
							$tempdata[6] = $data[6];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							$tempdata[8] = $hm1;
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[0] if ($rp[0]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				next;
			}
		}
		if ($found1)
		{
			# both found, ignore whole event
			if ($found2)
			{
				
			}
			# -+ found
			else
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'del';
				$tempdata[1] = $id2;
				$tempdata[2] = $mpd2;
				$tempdata[3] = $srd2;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[11];
				$tempdata[6] = $data[7];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				if ($hm2 >= 0)
				{
					$tempdata[8] = $hm2;
				}
				else
				{
					$tempdata[0] = 'del_ins';
					$tempdata[8] = '-';
					$tempdata[9] = '-';
					$tempdata[10] = '-';
					$tempdata[11] = -$hm2;
				}
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[1] if ($rp[1]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
		}
		else
		{
			# +- found
			if ($found2)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'tandem_dup';
				$tempdata[1] = $id1;
				$tempdata[2] = $mpd1;
				$tempdata[3] = $srd1;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[10];
				$tempdata[6] = $data[6];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				$tempdata[8] = $hm1;
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[0] if ($rp[0]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
			# both not found, keep whole event
			else
			{
				my $eventtype = $data[0];
				shift(@data); shift(@data);
				unshift(@data, $eventtype);
				my $toprint = join ("\t", @data);
				$toprint .= "\tRP:".$rp[0].'/'.$rp[1] if ($rp[0] and $rp[1]);
				print TUMOR "$toprint\n";
			}
		}
	}
	if ($data[0] =~ /insod/)
	{
		my @hm = split ("\/", $data[14]);
		my $sr = substr ($data1[15], 3);
		my @sr = split ("_", $sr);
		if ($filter_sr)
		{
			$found1 = 1 if ($sr[0] eq 'SR' and $sr[3] eq 'SR');
			$found2 = 1 if ($sr[1] eq 'SR' and $sr[2] eq 'SR');
		}
		my @hm = split ("\/", $data[14]);
		$found1 = 1 if ($del_homology and abs($hm[0])>$del_homology);
		$found2 = 1 if ($del_homology and abs($hm[1])>$del_homology);
		$found1 = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
		$found1 = 1 if ($filter_discord and &dcn($data[9], $data[11], 1));
		$found2 = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
		$found2 = 1 if ($filter_discord and &dcn($data[9], $data[10], -1));
		$found1 = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found1 = 1 if ($filter_non_uniq and &nud($data[9], $data[11]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[9], $data[10]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[9], $data[11]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[9], $data[10]));
		$found1 = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[9], $data[11], 1, $hm[0]));
		$found2 = 1 if ($disc_tumor and &dt($data[5], $data[7], -1, $data[9], $data[10], -1, $hm[1]));
		my @mpd = split ("\/", $data[3]);
		my @srd = split ("\/", $data[4]);
		if (defined($blat_server) and defined($blat_port))
		{
			$found1 = 1 unless (&local_align ($mpd[0], $data[5], $data[6], $data[9], $data[11]));
			$found2 = 1 unless (&local_align ($mpd[1], $data[5], $data[7], $data[9], $data[10]));
		}
		if ($filter_mate)
		{
			$found1 = 1 unless (&mate($mpd[0], $data[5], $data[6], $data[9], $data[11]));
			$found2 = 1 unless (&mate($mpd[1], $data[5], $data[7], $data[9], $data[10]));
		}
		$found1 = 1 if ($no_lc and $mpd[0] < $mpd);
		$found1 = 1 if ($no_lc and $srd[0] < $srd);
		$found1 = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
		$found2 = 1 if ($no_lc and $mpd[1] < $mpd);
		$found2 = 1 if ($no_lc and $srd[1] < $srd);
		$found2 = 1 if ($no_lc and ($mpd[1] + $srd[1] < $mpsrd));
		if ($complex_sv)
		{
			my @cluster_id = split ("\/", $data[2]);
			if ($cluster_id[0] =~ /_/)
			{
				$cluster_id[0] = $`;
			}
			$found1 = 1 if ($no_tie and $discord{$cluster_id[0]});
			if ($cluster_id[1] =~ /_/)
			{
				$cluster_id[1] = $`;
			}
			$found2 = 1 if ($no_tie and $discord{$cluster_id[1]});
		}
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord1 = $rp1[2];
			$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
			$found1 = 1 if ($alle_freq1 < $alle_freq_cutoff);
			my @rp2 = split ('_', $rp[1]);
			my $max_concord2 = $rp2[2];
			$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
			my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
			$found2 = 1 if ($alle_freq2 < $alle_freq_cutoff);
		}
LOOP9:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP9 if ($found1);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{2}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[6] - $$cluster[4]) < $distance and abs($data[11] - $$cluster[7]) < $distance)
				{
					$found1 = 1;
					last LOOP9;
				}
			}
		}
LOOP10:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP10 if ($found2);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{3}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[7] - $$cluster[4]) < $distance and abs($data[10] - $$cluster[7]) < $distance)
				{
					$found2 = 1;
					last LOOP10;
				}
			}
		}
		unless ($data[2] =~ /\//)
		{
			$found1 = 1;
			$found2 = 1;
		}
		if ($sv_min > 0)
		{
			if ($data[12] < $di_min or $data[13] < $di_min)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found2)
						{
							my @tempdata;
							$tempdata[0] = 'invers_r';
							$tempdata[1] = $id2;
							$tempdata[2] = $mpd2;
							$tempdata[3] = $srd2;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[7];
							$tempdata[6] = $data[10];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							$tempdata[8] = $hm2;
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[1] if ($rp[1]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found1)
						{
							my @tempdata;
							$tempdata[0] = 'invers_f';
							$tempdata[1] = $id1;
							$tempdata[2] = $mpd1;
							$tempdata[3] = $srd1;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[6];
							$tempdata[6] = $data[11];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							$tempdata[8] = $hm1;
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[0] if ($rp[0]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				next;
			}
		}
		if ($found1)
		{
			# both found, ignore whole event
			if ($found2)
			{
				
			}
			# ++ found
			else
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'invers_r';
				$tempdata[1] = $id2;
				$tempdata[2] = $mpd2;
				$tempdata[3] = $srd2;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[7];
				$tempdata[6] = $data[10];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				$tempdata[8] = $hm2;
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[1] if ($rp[1]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
		}
		else
		{
			# -- found
			if ($found2)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'invers_f';
				$tempdata[1] = $id1;
				$tempdata[2] = $mpd1;
				$tempdata[3] = $srd1;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[6];
				$tempdata[6] = $data[11];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				$tempdata[8] = $hm1;
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[0] if ($rp[0]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
			# both not found, keep whole event
			else
			{
				my $eventtype = $data[0];
				shift(@data); shift(@data);
				unshift(@data, $eventtype);
				my $toprint = join ("\t", @data);
				$toprint .= "\tRP:".$rp[0].'/'.$rp[1] if ($rp[0] and $rp[1]);
				print TUMOR "$toprint\n";
			}
		}
	}
	if ($data[0] =~ /insou/)
	{
		my @hm = split ("\/", $data[14]);
		my $sr = substr ($data1[15], 3);
		my @sr = split ("_", $sr);
		if ($filter_sr)
		{
			$found1 = 1 if ($sr[0] eq 'SR' and $sr[3] eq 'SR');
			$found2 = 1 if ($sr[1] eq 'SR' and $sr[2] eq 'SR');
		}
		my @hm = split ("\/", $data[14]);
		$found1 = 1 if ($del_homology and abs($hm[0])>$del_homology);
		$found2 = 1 if ($del_homology and abs($hm[1])>$del_homology);
		$found1 = 1 if ($filter_discord and &dcn($data[9], $data[11], 1));
		$found1 = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
		$found2 = 1 if ($filter_discord and &dcn($data[9], $data[10], -1));
		$found2 = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
		$found1 = 1 if ($filter_non_uniq and &nud($data[9], $data[11]));
		$found1 = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[9], $data[10]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[9], $data[11]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[9], $data[10]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found1 = 1 if ($disc_tumor and &dt($data[9], $data[11], 1, $data[5], $data[6], 1, $hm[0]));
		$found2 = 1 if ($disc_tumor and &dt($data[9], $data[10], -1, $data[5], $data[7], -1, $hm[1]));
		my @mpd = split ("\/", $data[3]);
		my @srd = split ("\/", $data[4]);
		if (defined($blat_server) and defined($blat_port))
		{
			$found1 = 1 unless (&local_align ($mpd[0], $data[5], $data[6], $data[9], $data[11]));
			$found2 = 1 unless (&local_align ($mpd[1], $data[5], $data[7], $data[9], $data[10]));
		}
		if ($filter_mate)
		{
			$found1 = 1 unless (&mate($mpd[0], $data[5], $data[6], $data[9], $data[11]));
			$found2 = 1 unless (&mate($mpd[1], $data[5], $data[7], $data[9], $data[10]));
		}
		$found1 = 1 if ($no_lc and $mpd[0] < $mpd);
		$found1 = 1 if ($no_lc and $srd[0] < $srd);
		$found1 = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
		$found2 = 1 if ($no_lc and $mpd[1] < $mpd);
		$found2 = 1 if ($no_lc and $srd[1] < $srd);
		$found2 = 1 if ($no_lc and ($mpd[1] + $srd[1] < $mpsrd));
		if ($complex_sv)
		{
			my @cluster_id = split ("\/", $data[2]);
			if ($cluster_id[0] =~ /_/)
			{
				$cluster_id[0] = $`;
			}
			$found1 = 1 if ($no_tie and $discord{$cluster_id[0]});
			if ($cluster_id[1] =~ /_/)
			{
				$cluster_id[1] = $`;
			}
			$found2 = 1 if ($no_tie and $discord{$cluster_id[1]});
		}
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord1 = $rp1[2];
			$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
			$found1 = 1 if ($alle_freq1 < $alle_freq_cutoff);
			my @rp2 = split ('_', $rp[1]);
			my $max_concord2 = $rp2[2];
			$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
			my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
			$found2 = 1 if ($alle_freq2 < $alle_freq_cutoff);
		}
LOOP11:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP11 if ($found1);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{2}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[11] - $$cluster[4]) < $distance and abs($data[6] - $$cluster[7]) < $distance)
				{
					$found1 = 1;
					last LOOP11;
				}
			}
		}
LOOP12:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP12 if ($found2);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{3}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[10] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
				{
					$found2 = 1;
					last LOOP12;
				}
			}
		}
		unless ($data[2] =~ /\//)
		{
			$found1 = 1;
			$found2 = 1;
		}
		if ($sv_min > 0)
		{
			if ($data[12] < $di_min or $data[13] < $di_min)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found2)
						{
							my @tempdata;
							$tempdata[0] = 'invers_r';
							$tempdata[1] = $id2;
							$tempdata[2] = $mpd2;
							$tempdata[3] = $srd2;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[10];
							$tempdata[6] = $data[7];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							$tempdata[8] = $hm2;
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[1] if ($rp[1]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found1)
						{
							my @tempdata;
							$tempdata[0] = 'invers_f';
							$tempdata[1] = $id1;
							$tempdata[2] = $mpd1;
							$tempdata[3] = $srd1;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[11];
							$tempdata[6] = $data[6];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							$tempdata[8] = $hm1;
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[0] if ($rp[0]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				next;
			}
		}
		if ($found1)
		{
			# both found, ignore whole event
			if ($found2)
			{
				
			}
			# ++ found
			else
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'invers_r';
				$tempdata[1] = $id2;
				$tempdata[2] = $mpd2;
				$tempdata[3] = $srd2;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[10];
				$tempdata[6] = $data[7];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				$tempdata[8] = $hm2;
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[1] if ($rp[1]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
		}
		else
		{
			# -- found
			if ($found2)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[14]);
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'invers_f';
				$tempdata[1] = $id1;
				$tempdata[2] = $mpd1;
				$tempdata[3] = $srd1;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[11];
				$tempdata[6] = $data[6];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				$tempdata[8] = $hm1;
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[0] if ($rp[0]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
			# both not found, keep whole event
			else
			{
				my $eventtype = $data[0];
				shift(@data); shift(@data);
				unshift(@data, $eventtype);
				my $toprint = join ("\t", @data);
				$toprint .= "\tRP:".$rp[0].'/'.$rp[1] if ($rp[0] and $rp[1]);
				print TUMOR "$toprint\n";
			}
		}
	}
	if ($data[0] eq 'del_invers')
	{
		my @hm = split ("\/", $data[15]);
		my $sr = substr ($data1[16], 3);
		my @sr = split ("_", $sr);
		if ($filter_sr)
		{
			$found1 = 1 if ($sr[0] eq 'SR' and $sr[3] eq 'SR');
			$found2 = 1 if ($sr[1] eq 'SR' and $sr[2] eq 'SR');
		}
		my @hm = split ("\/", $data[15]);
		$found1 = 1 if ($del_homology and abs($hm[0])>$del_homology);
		$found2 = 1 if ($del_homology and abs($hm[1])>$del_homology);
		$found1 = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
		$found1 = 1 if ($filter_discord and &dcn($data[9], $data[11], 1));
		$found2 = 1 if ($filter_discord and &dcn($data[9], $data[10], -1));
		$found2 = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
		$found1 = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found1 = 1 if ($filter_non_uniq and &nud($data[9], $data[11]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[9], $data[10]));
		$found2 = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[9], $data[11]));
		$found1 = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[9], $data[10]));
		$found2 = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
		$found1 = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[9], $data[11], 1, $hm[0]));
		$found2 = 1 if ($disc_tumor and &dt($data[9], $data[10], -1, $data[5], $data[7], -1, $hm[1]));
		my @mpd = split ("\/", $data[3]);
		my @srd = split ("\/", $data[4]);
		if (defined($blat_server) and defined($blat_port))
		{
			$found1 = 1 unless (&local_align ($mpd[0], $data[5], $data[6], $data[9], $data[11]));
			$found2 = 1 unless (&local_align ($mpd[1], $data[5], $data[7], $data[9], $data[10]));
		}
		if ($filter_mate)
		{
			$found1 = 1 unless (&mate($mpd[0], $data[5], $data[6], $data[9], $data[11]));
			$found2 = 1 unless (&mate($mpd[1], $data[5], $data[7], $data[9], $data[10]));
		}
		$found1 = 1 if ($no_lc and $mpd[0] < $mpd);
		$found1 = 1 if ($no_lc and $srd[0] < $srd);
		$found1 = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
		$found2 = 1 if ($no_lc and $mpd[1] < $mpd);
		$found2 = 1 if ($no_lc and $srd[1] < $srd);
		$found2 = 1 if ($no_lc and ($mpd[1] + $srd[1] < $mpsrd));
		if ($complex_sv)
		{
			my @cluster_id = split ("\/", $data[2]);
			if ($cluster_id[0] =~ /_/)
			{
				$cluster_id[0] = $`;
			}
			$found1 = 1 if ($no_tie and $discord{$cluster_id[0]});
			if ($cluster_id[1] =~ /_/)
			{
				$cluster_id[1] = $`;
			}
			$found2 = 1 if ($no_tie and $discord{$cluster_id[1]});
		}
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord1 = $rp1[2];
			$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
			$found1 = 1 if ($alle_freq1 < $alle_freq_cutoff);
			my @rp2 = split ('_', $rp[1]);
			my $max_concord2 = $rp2[2];
			$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
			my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
			$found2 = 1 if ($alle_freq2 < $alle_freq_cutoff);
		}
LOOP13:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP13 if ($found1);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{2}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[6] - $$cluster[4]) < $distance and abs($data[11] - $$cluster[7]) < $distance)
				{
					$found1 = 1;
					last LOOP13;
				}
			}
		}
LOOP14:		foreach my $sample (keys %normal_cluster)
		{
			last LOOP14 if ($found2);
			foreach (@{$normal_cluster{$sample}{$data[5]}{$data[5]}{3}})
			{
				my $cluster = $_;
				#print "@$_\n";
				if (abs($data[10] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
				{
					$found2 = 1;
					last LOOP14;
				}
			}
		}
		unless ($data[2] =~ /\//)
		{
			$found1 = 1;
			$found2 = 1;
		}
		unless ($data[15] =~ /\//)
		{
			$found1 = 1;
			$found2 = 1;
		}
		if ($sv_min > 0)
		{
			if ($data[12] < $di_min)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[15]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found2)
						{
							my @tempdata;
							$tempdata[0] = 'invers_r';
							$tempdata[1] = $id2;
							$tempdata[2] = $mpd2;
							$tempdata[3] = $srd2;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[10];
							$tempdata[6] = $data[7];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							$tempdata[8] = $hm2;
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[1] if ($rp[1]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					if ($smallest == 0)
					{
						unless ($found1)
						{
							my @tempdata;
							$tempdata[0] = 'invers_f';
							$tempdata[1] = $id1;
							$tempdata[2] = $mpd1;
							$tempdata[3] = $srd1;
							$tempdata[4] = $data[5];
							$tempdata[5] = $data[6];
							$tempdata[6] = $data[11];
							$tempdata[7] = $tempdata[6] - $tempdata[5];
							$tempdata[8] = $hm1;
							my $toprint = join ("\t", @tempdata);
							$toprint .= "\tRP:".$rp[0] if ($rp[0]);
							print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
						}
					}
				}
				next;
			}
		}
		if ($found1)
		{
			# both found, ignore whole event
			if ($found2)
			{
				
			}
			# ++ found
			else
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[15]);
				if ($id2 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'invers_r';
				$tempdata[1] = $id2;
				$tempdata[2] = $mpd2;
				$tempdata[3] = $srd2;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[10];
				$tempdata[6] = $data[7];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				$tempdata[8] = $hm2;
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[1] if ($rp[1]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
		}
		else
		{
			# -- found
			if ($found2)
			{
				my ($id1, $id2) = split (/\//, $data[2]);
				my ($mpd1, $mpd2) = split (/\//, $data[3]);
				my ($srd1, $srd2) = split (/\//, $data[4]);
				my ($hm1, $hm2) = split (/\//, $data[15]);
				if ($id1 =~ /_/)
				{
					my $smallest = $';#'
					next if ($smallest > 0);
				}
				my @tempdata;
				$tempdata[0] = 'invers_f';
				$tempdata[1] = $id1;
				$tempdata[2] = $mpd1;
				$tempdata[3] = $srd1;
				$tempdata[4] = $data[5];
				$tempdata[5] = $data[6];
				$tempdata[6] = $data[11];
				$tempdata[7] = $tempdata[6] - $tempdata[5];
				$tempdata[8] = $hm1;
				my $toprint = join ("\t", @tempdata);
				$toprint .= "\tRP:".$rp[0] if ($rp[0]);
				print TUMOR "$toprint\n" if ($tempdata[7] >= $sv_min);
			}
			# both not found, keep whole event
			else
			{
				my $eventtype = $data[0];
				shift(@data); shift(@data);
				unshift(@data, $eventtype);
				my $toprint = join ("\t", @data);
				$toprint .= "\tRP:".$rp[0].'/'.$rp[1] if ($rp[0] and $rp[1]);
				print TUMOR "$toprint\n";
			}
		}
	}
	if ($data[0] eq 'transl_inter')
	{
		if ($filter_sr)
		{
			$found = 1 if ($newline =~ /SR_SR/);
		}
		$found = 1 if ($inter_transl_homology and abs($data[11])>$inter_transl_homology);
		$found = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
		$found = 1 if ($filter_non_uniq and &nud($data[8], $data[9]));
		$found = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
		$found = 1 if ($filter_soft_clip and &scd($data[8], $data[9]));
		$found = 1 if ($disc_tumor and &dt($data[5], $data[6], $data[7], $data[8], $data[9], $data[10], $data[11]));
		$found = 1 if ($no_lc and $data[3] < $mpd);
		$found = 1 if ($no_lc and $data[4] < $srd);
		$found = 1 if ($no_lc and ($data[3] + $data[4] < $mpsrd));
		my $cluster_id = $data[2];
		if (defined($blat_server) and defined($blat_port))
		{
			$found = 1 unless (&local_align ($data[2], $data[5], $data[6], $data[8], $data[9]));
		}
		if ($filter_mate)
		{
			$found = 1 unless (&mate($data[2], $data[5], $data[6], $data[8], $data[9]));
		}
		if ($data[2] =~ /_/)
		{
			my $smallest = $';#'
			$cluster_id = $`;
			$found = 1 if ($smallest > 0);
		}
		$found = 1 if ($no_tie and $discord{$cluster_id});
		if ($alle_freq_cutoff)
		{
			my @rp1 = split ('_', $rp[0]);
			my $max_concord = $rp1[2];
			$max_concord = $rp1[3] if ($rp1[2] < $rp1[3]);
			my $alle_freq = $rp1[0]/($rp1[0]+$max_concord);
			$found = 1 if ($alle_freq < $alle_freq_cutoff);
		}
		if ($data[7] == 1 and $data[10] == -1)
		{
			$found = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
			$found = 1 if ($filter_discord and &dcn($data[8], $data[9], -1));
LOOP15:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP15 if ($found);
				foreach (@{$normal_cluster{$sample}{$data[5]}{$data[8]}{0}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[6] - $$cluster[4]) < $distance and abs($data[9] - $$cluster[7]) < $distance)
					{
						$found = 1;
						last LOOP15;
					}
				}
			}
			unless ($found)
			{
				my $eventtype = $data[0];
				shift(@data); shift(@data);
				unshift(@data, $eventtype);
				my $toprint = join ("\t", @data);
				$toprint .= "\tRP:".$rp[0] if ($rp[0]);
				print TUMOR "$toprint\n";
			}
		}
		if ($data[7] == -1 and $data[10] == 1)
		{
			$found = 1 if ($filter_discord and &dcn($data[5], $data[6], -1));
			$found = 1 if ($filter_discord and &dcn($data[8], $data[9], 1));
LOOP16:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP16 if ($found);
				foreach (@{$normal_cluster{$sample}{$data[5]}{$data[8]}{1}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[6] - $$cluster[4]) < $distance and abs($data[9] - $$cluster[7]) < $distance)
					{
						$found = 1;
						last LOOP16;
					}
				}
			}
			unless ($found)
			{
				my $eventtype = $data[0];
				shift(@data); shift(@data);
				unshift(@data, $eventtype);
				my $toprint = join ("\t", @data);
				$toprint .= "\tRP:".$rp[0] if ($rp[0]);
				print TUMOR "$toprint\n";
			}
		}
		if ($data[7] == 1 and $data[10] == 1)
		{
			$found = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
			$found = 1 if ($filter_discord and &dcn($data[8], $data[9], 1));
LOOP17:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP17 if ($found);
				foreach (@{$normal_cluster{$sample}{$data[5]}{$data[8]}{2}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[6] - $$cluster[4]) < $distance and abs($data[9] - $$cluster[7]) < $distance)
					{
						$found = 1;
						last LOOP17;
					}
				}
			}
			unless ($found)
			{
				my $eventtype = $data[0];
				shift(@data); shift(@data);
				unshift(@data, $eventtype);
				my $toprint = join ("\t", @data);
				$toprint .= "\tRP:".$rp[0] if ($rp[0]);
				print TUMOR "$toprint\n";
			}
		}
		if ($data[7] == -1 and $data[10] == -1)
		{
			$found = 1 if ($filter_discord and &dcn($data[5], $data[6], -1));
			$found = 1 if ($filter_discord and &dcn($data[8], $data[9], -1));
LOOP18:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP18 if ($found);
				foreach (@{$normal_cluster{$sample}{$data[5]}{$data[8]}{3}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[6] - $$cluster[4]) < $distance and abs($data[9] - $$cluster[7]) < $distance)
					{
						$found = 1;
						last LOOP18;
					}
				}
			}
			unless ($found)
			{
				my $eventtype = $data[0];
				shift(@data); shift(@data);
				unshift(@data, $eventtype);
				my $toprint = join ("\t", @data);
				$toprint .= "\tRP:".$rp[0] if ($rp[0]);
				print TUMOR "$toprint\n";
			}
		}
	}
	if ($data[0] eq 'inss' or $data[0] eq 'del_inss')
	{
		my @hm = split ("\/", $data[13]);
		if ($data[5] lt $data[9])
		{
			my $sr = substr ($data1[14], 3);
			my @sr = split ("_", $sr);
			if ($filter_sr)
			{
				$found1 = 1 if ($sr[0] eq 'SR' and $sr[2] eq 'SR');
				$found2 = 1 if ($sr[1] eq 'SR' and $sr[3] eq 'SR');
			}
			my @hm = split ("\/", $data[13]);
			$found1 = 1 if ($inter_transl_homology and abs($hm[0])>$inter_transl_homology);
			$found2 = 1 if ($inter_transl_homology and abs($hm[1])>$inter_transl_homology);
			$found1 = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
			$found1 = 1 if ($filter_discord and &dcn($data[9], $data[10], -1));
			$found2 = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
			$found2 = 1 if ($filter_discord and &dcn($data[9], $data[11], 1));
			$found1 = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
			$found1 = 1 if ($filter_non_uniq and &nud($data[9], $data[10]));
			$found2 = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
			$found2 = 1 if ($filter_non_uniq and &nud($data[9], $data[11]));
			$found1 = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
			$found1 = 1 if ($filter_soft_clip and &scd($data[9], $data[10]));
			$found2 = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
			$found2 = 1 if ($filter_soft_clip and &scd($data[9], $data[11]));
			$found1 = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[9], $data[10], -1, $hm[0]));
			$found2 = 1 if ($disc_tumor and &dt($data[5], $data[7], -1, $data[9], $data[11], 1, $hm[1]));
			my @mpd = split ("\/", $data[3]);
			my @srd = split ("\/", $data[4]);
			if (defined($blat_server) and defined($blat_port))
			{
				$found1 = 1 unless (&local_align ($mpd[0], $data[5], $data[6], $data[9], $data[10]));
				$found2 = 1 unless (&local_align ($mpd[1], $data[5], $data[7], $data[9], $data[11]));
			}
			if ($filter_mate)
			{
				$found1 = 1 unless (&mate($mpd[0], $data[5], $data[6], $data[9], $data[10]));
				$found2 = 1 unless (&mate($mpd[1], $data[5], $data[7], $data[9], $data[11]));
			}
			$found1 = 1 if ($no_lc and $mpd[0] < $mpd);
			$found1 = 1 if ($no_lc and $srd[0] < $srd);
			$found1 = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
			$found2 = 1 if ($no_lc and $mpd[1] < $mpd);
			$found2 = 1 if ($no_lc and $srd[1] < $srd);
			$found2 = 1 if ($no_lc and ($mpd[1] + $srd[1] < $mpsrd));
			if ($complex_sv)
			{
				my @cluster_id = split ("\/", $data[2]);
				if ($cluster_id[0] =~ /_/)
				{
					$cluster_id[0] = $`;
				}
				$found1 = 1 if ($no_tie and $discord{$cluster_id[0]});
				if ($cluster_id[1] =~ /_/)
				{
					$cluster_id[1] = $`;
				}
				$found2 = 1 if ($no_tie and $discord{$cluster_id[1]});
			}
			if ($alle_freq_cutoff)
			{
				my @rp1 = split ('_', $rp[0]);
				my $max_concord1 = $rp1[2];
				$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
				my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
				$found1 = 1 if ($alle_freq1 < $alle_freq_cutoff);
				my @rp2 = split ('_', $rp[1]);
				my $max_concord2 = $rp2[2];
				$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
				my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
				$found2 = 1 if ($alle_freq2 < $alle_freq_cutoff);
			}
LOOP19:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP19 if ($found1);
				foreach (@{$normal_cluster{$sample}{$data[5]}{$data[9]}{0}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[6] - $$cluster[4]) < $distance and abs($data[10] - $$cluster[7]) < $distance)
					{
						$found1 = 1;
						last LOOP19;
					}
				}
			}
LOOP20:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP20 if ($found2);
				foreach (@{$normal_cluster{$sample}{$data[5]}{$data[9]}{1}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[7] - $$cluster[4]) < $distance and abs($data[11] - $$cluster[7]) < $distance)
					{
						$found2 = 1;
						last LOOP20;
					}
				}
			}
			unless ($data[2] =~ /\//)
			{
				$found1 = 1;
				$found2 = 1;
			}
			if ($sv_min > 0)
			{
				if ($data[12] < $di_min)
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id2 =~ /_/)
					{
						my $smallest = $';#'
						if ($smallest == 0)
						{
							unless ($found2)
							{
								my @tempdata;
								$tempdata[0] = 'transl_inter';
								$tempdata[1] = $id2;
								$tempdata[2] = $mpd2;
								$tempdata[3] = $srd2;
								$tempdata[4] = $data[5];
								$tempdata[5] = $data[7];
								$tempdata[6] = '-1';
								$tempdata[7] = $data[9];
								$tempdata[8] = $data[11];
								$tempdata[9] = '1';
								$tempdata[10] = $hm2;
								my $toprint = join ("\t", @tempdata);
								$toprint .= "\tRP:".$rp[1] if ($rp[1]);
								print TUMOR "$toprint\n";
							}
						}
					}
					
					if ($id1 =~ /_/)
					{
						my $smallest = $';#'
						if ($smallest == 0)
						{
							unless ($found1)
							{
								my @tempdata;
								$tempdata[0] = 'transl_inter';
								$tempdata[1] = $id1;
								$tempdata[2] = $mpd1;
								$tempdata[3] = $srd1;
								$tempdata[4] = $data[5];
								$tempdata[5] = $data[6];
								$tempdata[6] = '1';
								$tempdata[7] = $data[9];
								$tempdata[8] = $data[10];
								$tempdata[9] = '-1';
								$tempdata[10] = $hm1;
								my $toprint = join ("\t", @tempdata);
								$toprint .= "\tRP:".$rp[0] if ($rp[0]);
								print TUMOR "$toprint\n";
							}
						}
					}
					next;
				}
			}
			if ($found1)
			{
				# both found, ignore whole event
				if ($found2)
				{
					
				}
				# +- found
				else
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id2 =~ /_/)
					{
						my $smallest = $';#'
						next if ($smallest > 0);
					}
					my @tempdata;
					$tempdata[0] = 'transl_inter';
					$tempdata[1] = $id2;
					$tempdata[2] = $mpd2;
					$tempdata[3] = $srd2;
					$tempdata[4] = $data[5];
					$tempdata[5] = $data[7];
					$tempdata[6] = '-1';
					$tempdata[7] = $data[9];
					$tempdata[8] = $data[11];
					$tempdata[9] = '1';
					$tempdata[10] = $hm2;
					my $toprint = join ("\t", @tempdata);
					$toprint .= "\tRP:".$rp[1] if ($rp[1]);
					print TUMOR "$toprint\n";
				}
			}
			else
			{
				# -+ found
				if ($found2)
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id1 =~ /_/)
					{
						my $smallest = $';#'
						next if ($smallest > 0);
					}
					my @tempdata;
					$tempdata[0] = 'transl_inter';
					$tempdata[1] = $id1;
					$tempdata[2] = $mpd1;
					$tempdata[3] = $srd1;
					$tempdata[4] = $data[5];
					$tempdata[5] = $data[6];
					$tempdata[6] = '1';
					$tempdata[7] = $data[9];
					$tempdata[8] = $data[10];
					$tempdata[9] = '-1';
					$tempdata[10] = $hm1;
					my $toprint = join ("\t", @tempdata);
					$toprint .= "\tRP:".$rp[0] if ($rp[0]);
					print TUMOR "$toprint\n";
				}
				# both not found, keep whole event
				else
				{
					my $eventtype = $data[0];
					shift(@data); shift(@data);
					unshift(@data, $eventtype);
					my $toprint = join ("\t", @data);
					$toprint .= "\tRP:".$rp[0].'/'.$rp[1] if ($rp[0] and $rp[1]);
					print TUMOR "$toprint\n";
				}
			}
		}
		else
		{
			my $sr = substr ($data1[14], 3);
			my @sr = split ("_", $sr);
			if ($filter_sr)
			{
				$found1 = 1 if ($sr[1] eq 'SR' and $sr[3] eq 'SR');
				$found2 = 1 if ($sr[0] eq 'SR' and $sr[2] eq 'SR');
			}
			my @hm = split ("\/", $data[13]);
			$found2 = 1 if ($inter_transl_homology and abs($hm[0])>$inter_transl_homology);
			$found1 = 1 if ($inter_transl_homology and abs($hm[1])>$inter_transl_homology); # inss events homology alway first for 6-10, second for 7-11, not in order of cluster id
			$found1 = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
			$found1 = 1 if ($filter_discord and &dcn($data[9], $data[11], 1));
			$found2 = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
			$found2 = 1 if ($filter_discord and &dcn($data[9], $data[10], -1));
			$found1 = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
			$found1 = 1 if ($filter_non_uniq and &nud($data[9], $data[11]));
			$found2 = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
			$found2 = 1 if ($filter_non_uniq and &nud($data[9], $data[10]));
			$found1 = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
			$found1 = 1 if ($filter_soft_clip and &scd($data[9], $data[11]));
			$found2 = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
			$found2 = 1 if ($filter_soft_clip and &scd($data[9], $data[10]));
			$found1 = 1 if ($disc_tumor and &dt($data[5], $data[7], -1, $data[9], $data[11], 1, $hm[0]));
			$found2 = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[9], $data[10], -1, $hm[1]));
			my @mpd = split ("\/", $data[3]);
			my @srd = split ("\/", $data[4]);
			if (defined($blat_server) and defined($blat_port))
			{
				$found1 = 1 unless (&local_align ($mpd[0], $data[5], $data[7], $data[9], $data[11]));
				$found2 = 1 unless (&local_align ($mpd[1], $data[5], $data[6], $data[9], $data[10]));
			}
			if ($filter_mate)
			{
				$found1 = 1 unless (&mate($mpd[0], $data[5], $data[7], $data[9], $data[11]));
				$found2 = 1 unless (&mate($mpd[1], $data[5], $data[6], $data[9], $data[10]));
			}
			$found1 = 1 if ($no_lc and $mpd[0] < $mpd);
			$found1 = 1 if ($no_lc and $srd[0] < $srd);
			$found1 = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
			$found2 = 1 if ($no_lc and $mpd[1] < $mpd);
			$found2 = 1 if ($no_lc and $srd[1] < $srd);
			$found2 = 1 if ($no_lc and ($mpd[1] + $srd[1] < $mpsrd));
			if ($complex_sv)
			{
				my @cluster_id = split ("\/", $data[2]);
				if ($cluster_id[0] =~ /_/)
				{
					$cluster_id[0] = $`;
				}
				$found1 = 1 if ($no_tie and $discord{$cluster_id[0]});
				if ($cluster_id[1] =~ /_/)
				{
					$cluster_id[1] = $`;
				}
				$found2 = 1 if ($no_tie and $discord{$cluster_id[1]});
			}
			if ($alle_freq_cutoff)
			{
				my @rp1 = split ('_', $rp[0]);
				my $max_concord1 = $rp1[2];
				$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
				my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
				$found1 = 1 if ($alle_freq1 < $alle_freq_cutoff);
				my @rp2 = split ('_', $rp[1]);
				my $max_concord2 = $rp2[2];
				$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
				my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
				$found2 = 1 if ($alle_freq2 < $alle_freq_cutoff);
			}
LOOP21:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP21 if ($found1);
				foreach (@{$normal_cluster{$sample}{$data[9]}{$data[5]}{0}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[11] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
					{
						$found1 = 1;
						last LOOP21;
					}
				}
			}
LOOP22:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP22 if ($found2);
				foreach (@{$normal_cluster{$sample}{$data[9]}{$data[5]}{1}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[10] - $$cluster[4]) < $distance and abs($data[6] - $$cluster[7]) < $distance)
					{
						$found2 = 1;
						last LOOP22;
					}
				}
			}
			unless ($data[2] =~ /\//)
			{
				$found1 = 1;
				$found2 = 1;
			}
			if ($sv_min > 0)
			{
				if ($data[12] < $di_min)
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id2 =~ /_/)
					{
						my $smallest = $';#'
						if ($smallest == 0)
						{
							unless ($found2)
							{
								my @tempdata;
								$tempdata[0] = 'transl_inter';
								$tempdata[1] = $id2;
								$tempdata[2] = $mpd2;
								$tempdata[3] = $srd2;
								$tempdata[4] = $data[9];
								$tempdata[5] = $data[10];
								$tempdata[6] = '-1';
								$tempdata[7] = $data[5];
								$tempdata[8] = $data[6];
								$tempdata[9] = '1';
								$tempdata[10] = $hm1;
								my $toprint = join ("\t", @tempdata);
								$toprint .= "\tRP:".$rp[1] if ($rp[1]);
								print TUMOR "$toprint\n";
							}
						}
					}
					
					if ($id1 =~ /_/)
					{
						my $smallest = $';#'
						if ($smallest == 0)
						{
							unless ($found1)
							{
								my @tempdata;
								$tempdata[0] = 'transl_inter';
								$tempdata[1] = $id1;
								$tempdata[2] = $mpd1;
								$tempdata[3] = $srd1;
								$tempdata[4] = $data[9];
								$tempdata[5] = $data[11];
								$tempdata[6] = '1';
								$tempdata[7] = $data[5];
								$tempdata[8] = $data[7];
								$tempdata[9] = '-1';
								$tempdata[10] = $hm2;
								my $toprint = join ("\t", @tempdata);
								$toprint .= "\tRP:".$rp[0] if ($rp[0]);
								print TUMOR "$toprint\n";
							}
						}
					}
					next;
				}
			}
			if ($found1)
			{
				# both found, ignore whole event
				if ($found2)
				{
					
				}
				# +- found
				else
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id2 =~ /_/)
					{
						my $smallest = $';#'
						next if ($smallest > 0);
					}
					my @tempdata;
					$tempdata[0] = 'transl_inter';
					$tempdata[1] = $id2;
					$tempdata[2] = $mpd2;
					$tempdata[3] = $srd2;
					$tempdata[4] = $data[9];
					$tempdata[5] = $data[10];
					$tempdata[6] = '-1';
					$tempdata[7] = $data[5];
					$tempdata[8] = $data[6];
					$tempdata[9] = '1';
					$tempdata[10] = $hm1;
					my $toprint = join ("\t", @tempdata);
					$toprint .= "\tRP:".$rp[1] if ($rp[1]);
					print TUMOR "$toprint\n";
				}
			}
			else
			{
				# -+ found
				if ($found2)
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id1 =~ /_/)
					{
						my $smallest = $';#'
						next if ($smallest > 0);
					}
					my @tempdata;
					$tempdata[0] = 'transl_inter';
					$tempdata[1] = $id1;
					$tempdata[2] = $mpd1;
					$tempdata[3] = $srd1;
					$tempdata[4] = $data[9];
					$tempdata[5] = $data[11];
					$tempdata[6] = '1';
					$tempdata[7] = $data[5];
					$tempdata[8] = $data[7];
					$tempdata[9] = '-1';
					$tempdata[10] = $hm2;
					my $toprint = join ("\t", @tempdata);
					$toprint .= "\tRP:".$rp[0] if ($rp[0]);
					print TUMOR "$toprint\n";
				}
				# both not found, keep whole event
				else
				{
					my $eventtype = $data[0];
					shift(@data); shift(@data);
					unshift(@data, $eventtype);
					my $toprint = join ("\t", @data);
					$toprint .= "\tRP:".$rp[0].'/'.$rp[1] if ($rp[0] and $rp[1]);
					print TUMOR "$toprint\n";
				}
			}
		}
	}
	if ($data[0] eq 'inso' or $data[0] eq 'del_inso')
	{
		my @hm = split ("\/", $data[13]);
		if ($data[5] lt $data[9])
		{
			my $sr = substr ($data1[14], 3);
			my @sr = split ("_", $sr);
			if ($filter_sr)
			{
				$found1 = 1 if ($sr[0] eq 'SR' and $sr[3] eq 'SR');
				$found2 = 1 if ($sr[1] eq 'SR' and $sr[2] eq 'SR');
			}
			my @hm = split ("\/", $data[13]);
			$found1 = 1 if ($inter_transl_homology and abs($hm[0])>$inter_transl_homology);
			$found2 = 1 if ($inter_transl_homology and abs($hm[1])>$inter_transl_homology);
			$found1 = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
			$found1 = 1 if ($filter_discord and &dcn($data[9], $data[11], 1));
			$found2 = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
			$found2 = 1 if ($filter_discord and &dcn($data[9], $data[10], -1));
			$found1 = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
			$found1 = 1 if ($filter_non_uniq and &nud($data[9], $data[11]));
			$found2 = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
			$found2 = 1 if ($filter_non_uniq and &nud($data[9], $data[10]));
			$found1 = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
			$found1 = 1 if ($filter_soft_clip and &scd($data[9], $data[11]));
			$found2 = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
			$found2 = 1 if ($filter_soft_clip and &scd($data[9], $data[10]));
			$found1 = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[9], $data[11], 1, $hm[0]));
			$found2 = 1 if ($disc_tumor and &dt($data[5], $data[7], -1, $data[9], $data[10], -1, $hm[1]));
			my @mpd = split ("\/", $data[3]);
			my @srd = split ("\/", $data[4]);
			if (defined($blat_server) and defined($blat_port))
			{
				$found1 = 1 unless (&local_align ($mpd[0], $data[5], $data[6], $data[9], $data[11]));
				$found2 = 1 unless (&local_align ($mpd[1], $data[5], $data[7], $data[9], $data[10]));
			}
			if ($filter_mate)
			{
				$found1 = 1 unless (&mate($mpd[0], $data[5], $data[6], $data[9], $data[11]));
				$found2 = 1 unless (&mate($mpd[1], $data[5], $data[7], $data[9], $data[10]));
			}
			$found1 = 1 if ($no_lc and $mpd[0] < $mpd);
			$found1 = 1 if ($no_lc and $srd[0] < $srd);
			$found1 = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
			$found2 = 1 if ($no_lc and $mpd[1] < $mpd);
			$found2 = 1 if ($no_lc and $srd[1] < $srd);
			$found2 = 1 if ($no_lc and ($mpd[1] + $srd[1] < $mpsrd));
			if ($complex_sv)
			{
				my @cluster_id = split ("\/", $data[2]);
				if ($cluster_id[0] =~ /_/)
				{
					$cluster_id[0] = $`;
				}
				$found1 = 1 if ($no_tie and $discord{$cluster_id[0]});
				if ($cluster_id[1] =~ /_/)
				{
					$cluster_id[1] = $`;
				}
				$found2 = 1 if ($no_tie and $discord{$cluster_id[1]});
			}
			if ($alle_freq_cutoff)
			{
				my @rp1 = split ('_', $rp[0]);
				my $max_concord1 = $rp1[2];
				$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
				my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
				$found1 = 1 if ($alle_freq1 < $alle_freq_cutoff);
				my @rp2 = split ('_', $rp[1]);
				my $max_concord2 = $rp2[2];
				$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
				my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
				$found2 = 1 if ($alle_freq2 < $alle_freq_cutoff);
			}
LOOP23:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP23 if ($found1);
				foreach (@{$normal_cluster{$sample}{$data[5]}{$data[9]}{2}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[6] - $$cluster[4]) < $distance and abs($data[11] - $$cluster[7]) < $distance)
					{
						$found1 = 1;
						last LOOP23;
					}
				}
			}
LOOP24:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP24 if ($found2);
				foreach (@{$normal_cluster{$sample}{$data[5]}{$data[9]}{3}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[7] - $$cluster[4]) < $distance and abs($data[10] - $$cluster[7]) < $distance)
					{
						$found2 = 1;
						last LOOP24;
					}
				}
			}
			unless ($data[2] =~ /\//)
			{
				$found1 = 1;
				$found2 = 1;
			}
			if ($sv_min > 0)
			{
				if ($data[12] < $di_min)
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id2 =~ /_/)
					{
						my $smallest = $';#'
						if ($smallest == 0)
						{
							unless ($found2)
							{
								my @tempdata;
								$tempdata[0] = 'transl_inter';
								$tempdata[1] = $id2;
								$tempdata[2] = $mpd2;
								$tempdata[3] = $srd2;
								$tempdata[4] = $data[5];
								$tempdata[5] = $data[7];
								$tempdata[6] = '-1';
								$tempdata[7] = $data[9];
								$tempdata[8] = $data[10];
								$tempdata[9] = '-1';
								$tempdata[10] = $hm2;
								my $toprint = join ("\t", @tempdata);
								$toprint .= "\tRP:".$rp[1] if ($rp[1]);
								print TUMOR "$toprint\n";
							}
						}
					}
					
					if ($id1 =~ /_/)
					{
						my $smallest = $';#'
						if ($smallest == 0)
						{
							unless ($found1)
							{
								my @tempdata;
								$tempdata[0] = 'transl_inter';
								$tempdata[1] = $id1;
								$tempdata[2] = $mpd1;
								$tempdata[3] = $srd1;
								$tempdata[4] = $data[5];
								$tempdata[5] = $data[6];
								$tempdata[6] = '1';
								$tempdata[7] = $data[9];
								$tempdata[8] = $data[11];
								$tempdata[9] = '1';
								$tempdata[10] = $hm1;
								my $toprint = join ("\t", @tempdata);
								$toprint .= "\tRP:".$rp[0] if ($rp[0]);
								print TUMOR "$toprint\n";
							}
						}
					}
					next;
				}
			}
			if ($found1)
			{
				# both found, ignore whole event
				if ($found2)
				{
					
				}
				# ++ found
				else
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id2 =~ /_/)
					{
						my $smallest = $';#'
						next if ($smallest > 0);
					}
					my @tempdata;
					$tempdata[0] = 'transl_inter';
					$tempdata[1] = $id2;
					$tempdata[2] = $mpd2;
					$tempdata[3] = $srd2;
					$tempdata[4] = $data[5];
					$tempdata[5] = $data[7];
					$tempdata[6] = '-1';
					$tempdata[7] = $data[9];
					$tempdata[8] = $data[10];
					$tempdata[9] = '-1';
					$tempdata[10] = $hm2;
					my $toprint = join ("\t", @tempdata);
					$toprint .= "\tRP:".$rp[1] if ($rp[1]);
					print TUMOR "$toprint\n";
				}
			}
			else
			{
				# -- found
				if ($found2)
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id1 =~ /_/)
					{
						my $smallest = $';#'
						next if ($smallest > 0);
					}
					my @tempdata;
					$tempdata[0] = 'transl_inter';
					$tempdata[1] = $id1;
					$tempdata[2] = $mpd1;
					$tempdata[3] = $srd1;
					$tempdata[4] = $data[5];
					$tempdata[5] = $data[6];
					$tempdata[6] = '1';
					$tempdata[7] = $data[9];
					$tempdata[8] = $data[11];
					$tempdata[9] = '1';
					$tempdata[10] = $hm1;
					my $toprint = join ("\t", @tempdata);
					$toprint .= "\tRP:".$rp[0] if ($rp[0]);
					print TUMOR "$toprint\n";
				}
				# both not found, keep whole event
				else
				{
					my $eventtype = $data[0];
					shift(@data); shift(@data);
					unshift(@data, $eventtype);
					my $toprint = join ("\t", @data);
					$toprint .= "\tRP:".$rp[0].'/'.$rp[1] if ($rp[0] and $rp[1]);
					print TUMOR "$toprint\n";
				}
			}
		}
		else
		{
			my $sr = substr ($data1[14], 3);
			my @sr = split ("_", $sr);
			if ($filter_sr)
			{
				$found1 = 1 if ($sr[0] eq 'SR' and $sr[3] eq 'SR');
				$found2 = 1 if ($sr[1] eq 'SR' and $sr[2] eq 'SR');
			}
			my @hm = split ("\/", $data[13]);
			$found1 = 1 if ($inter_transl_homology and abs($hm[0])>$inter_transl_homology);
			$found2 = 1 if ($inter_transl_homology and abs($hm[1])>$inter_transl_homology);
			$found1 = 1 if ($filter_discord and &dcn($data[5], $data[6], 1));
			$found1 = 1 if ($filter_discord and &dcn($data[9], $data[11], 1));
			$found2 = 1 if ($filter_discord and &dcn($data[5], $data[7], -1));
			$found2 = 1 if ($filter_discord and &dcn($data[9], $data[10], -1));
			$found1 = 1 if ($filter_non_uniq and &nud($data[5], $data[6]));
			$found1 = 1 if ($filter_non_uniq and &nud($data[9], $data[11]));
			$found2 = 1 if ($filter_non_uniq and &nud($data[5], $data[7]));
			$found2 = 1 if ($filter_non_uniq and &nud($data[9], $data[10]));
			$found1 = 1 if ($filter_soft_clip and &scd($data[5], $data[6]));
			$found1 = 1 if ($filter_soft_clip and &scd($data[9], $data[11]));
			$found2 = 1 if ($filter_soft_clip and &scd($data[5], $data[7]));
			$found2 = 1 if ($filter_soft_clip and &scd($data[9], $data[10]));
			$found1 = 1 if ($disc_tumor and &dt($data[5], $data[6], 1, $data[9], $data[11], 1, $hm[0]));
			$found2 = 1 if ($disc_tumor and &dt($data[5], $data[7], -1, $data[9], $data[10], -1, $hm[1]));
			my @mpd = split ("\/", $data[3]);
			my @srd = split ("\/", $data[4]);
			if (defined($blat_server) and defined($blat_port))
			{
				$found1 = 1 unless (&local_align ($mpd[0], $data[5], $data[6], $data[9], $data[11]));
				$found2 = 1 unless (&local_align ($mpd[1], $data[5], $data[7], $data[9], $data[10]));
			}
			if ($filter_mate)
			{
				$found1 = 1 unless (&mate($mpd[0], $data[5], $data[6], $data[9], $data[11]));
				$found2 = 1 unless (&mate($mpd[1], $data[5], $data[7], $data[9], $data[10]));
			}
			$found1 = 1 if ($no_lc and $mpd[0] < $mpd);
			$found1 = 1 if ($no_lc and $srd[0] < $srd);
			$found1 = 1 if ($no_lc and ($mpd[0] + $srd[0] < $mpsrd));
			$found2 = 1 if ($no_lc and $mpd[1] < $mpd);
			$found2 = 1 if ($no_lc and $srd[1] < $srd);
			$found2 = 1 if ($no_lc and ($mpd[1] + $srd[1] < $mpsrd));
			if ($complex_sv)
			{
				my @cluster_id = split ("\/", $data[2]);
				if ($cluster_id[0] =~ /_/)
				{
					$cluster_id[0] = $`;
				}
				$found1 = 1 if ($no_tie and $discord{$cluster_id[0]});
				if ($cluster_id[1] =~ /_/)
				{
					$cluster_id[1] = $`;
				}
				$found2 = 1 if ($no_tie and $discord{$cluster_id[1]});
			}
			if ($alle_freq_cutoff)
			{
				my @rp1 = split ('_', $rp[0]);
				my $max_concord1 = $rp1[2];
				$max_concord1 = $rp1[3] if ($rp1[2] < $rp1[3]);
				my $alle_freq1 = $rp1[0]/($rp1[0]+$max_concord1);
				$found1 = 1 if ($alle_freq1 < $alle_freq_cutoff);
				my @rp2 = split ('_', $rp[1]);
				my $max_concord2 = $rp2[2];
				$max_concord2 = $rp2[3] if ($rp2[2] < $rp2[3]);
				my $alle_freq2 = $rp2[0]/($rp2[0]+$max_concord2);
				$found2 = 1 if ($alle_freq2 < $alle_freq_cutoff);
			}
LOOP25:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP25 if ($found1);
				foreach (@{$normal_cluster{$sample}{$data[9]}{$data[5]}{2}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[11] - $$cluster[4]) < $distance and abs($data[6] - $$cluster[7]) < $distance)
					{
						$found1 = 1;
						last LOOP25;
					}
				}
			}
LOOP26:			foreach my $sample (keys %normal_cluster)
			{
				last LOOP26 if ($found2);
				foreach (@{$normal_cluster{$sample}{$data[9]}{$data[5]}{3}})
				{
					my $cluster = $_;
					#print "@$_\n";
					if (abs($data[10] - $$cluster[4]) < $distance and abs($data[7] - $$cluster[7]) < $distance)
					{
						$found2 = 1;
						last LOOP26;
					}
				}
			}
			unless ($data[2] =~ /\//)
			{
				$found1 = 1;
				$found2 = 1;
			}
			if ($sv_min > 0)
			{
				if ($data[12] < $di_min)
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id2 =~ /_/)
					{
						my $smallest = $';#'
						if ($smallest == 0)
						{
							unless ($found2)
							{
								my @tempdata;
								$tempdata[0] = 'transl_inter';
								$tempdata[1] = $id2;
								$tempdata[2] = $mpd2;
								$tempdata[3] = $srd2;
								$tempdata[4] = $data[9];
								$tempdata[5] = $data[10];
								$tempdata[6] = '-1';
								$tempdata[7] = $data[5];
								$tempdata[8] = $data[7];
								$tempdata[9] = '-1';
								$tempdata[10] = $hm2;
								my $toprint = join ("\t", @tempdata);
								$toprint .= "\tRP:".$rp[1] if ($rp[1]);
								print TUMOR "$toprint\n";
							}
						}
					}
					
					if ($id1 =~ /_/)
					{
						my $smallest = $';#'
						if ($smallest == 0)
						{
							unless ($found1)
							{
								my @tempdata;
								$tempdata[0] = 'transl_inter';
								$tempdata[1] = $id1;
								$tempdata[2] = $mpd1;
								$tempdata[3] = $srd1;
								$tempdata[4] = $data[9];
								$tempdata[5] = $data[11];
								$tempdata[6] = '1';
								$tempdata[7] = $data[5];
								$tempdata[8] = $data[6];
								$tempdata[9] = '1';
								$tempdata[10] = $hm1;
								my $toprint = join ("\t", @tempdata);
								$toprint .= "\tRP:".$rp[0] if ($rp[0]);
								print TUMOR "$toprint\n";
							}
						}
					}
					next;
				}
			}
			if ($found1)
			{
				# both found, ignore whole event
				if ($found2)
				{
					
				}
				# ++ found
				else
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id2 =~ /_/)
					{
						my $smallest = $';#'
						next if ($smallest > 0);
					}
					my @tempdata;
					$tempdata[0] = 'transl_inter';
					$tempdata[1] = $id2;
					$tempdata[2] = $mpd2;
					$tempdata[3] = $srd2;
					$tempdata[4] = $data[9];
					$tempdata[5] = $data[10];
					$tempdata[6] = '-1';
					$tempdata[7] = $data[5];
					$tempdata[8] = $data[7];
					$tempdata[9] = '-1';
					$tempdata[10] = $hm2;
					my $toprint = join ("\t", @tempdata);
					$toprint .= "\tRP:".$rp[1] if ($rp[1]);
					print TUMOR "$toprint\n";
				}
			}
			else
			{
				# -- found
				if ($found2)
				{
					my ($id1, $id2) = split (/\//, $data[2]);
					my ($mpd1, $mpd2) = split (/\//, $data[3]);
					my ($srd1, $srd2) = split (/\//, $data[4]);
					my ($hm1, $hm2) = split (/\//, $data[13]);
					if ($id1 =~ /_/)
					{
						my $smallest = $';#'
						next if ($smallest > 0);
					}
					my @tempdata;
					$tempdata[0] = 'transl_inter';
					$tempdata[1] = $id1;
					$tempdata[2] = $mpd1;
					$tempdata[3] = $srd1;
					$tempdata[4] = $data[9];
					$tempdata[5] = $data[11];
					$tempdata[6] = '1';
					$tempdata[7] = $data[5];
					$tempdata[8] = $data[6];
					$tempdata[9] = '1';
					$tempdata[10] = $hm1;
					my $toprint = join ("\t", @tempdata);
					$toprint .= "\tRP:".$rp[0] if ($rp[0]);
					print TUMOR "$toprint\n";
				}
				# both not found, keep whole event
				else
				{
					my $eventtype = $data[0];
					shift(@data); shift(@data);
					unshift(@data, $eventtype);
					my $toprint = join ("\t", @data);
					$toprint .= "\tRP:".$rp[0].'/'.$rp[1] if ($rp[0] and $rp[1]);
					print TUMOR "$toprint\n";
				}
			}
		}
	}
}
close FILE;
close TUMOR;

system "perl $Bin/mechanism2.pl -i $tumor_out -R $rmsk";
system "rm $tumor_out";
system "rm -rf $tmpdir" if (-e $tmpdir);

my $time1 = time;
my $local1 = localtime($time1);
my $difference = $time1 - $time0;
my $seconds    =  $difference % 60;
$difference = ($difference - $seconds) / 60;
my $minutes    =  $difference % 60;
$difference = ($difference - $minutes) / 60;
print STDERR "$local1 somatic_sv.pl Finished\n";
print STDERR "Time used: $difference:$minutes:$seconds\n";

sub mate
{
	my $cluster_id = shift;
	my $chra = shift;
	my $posa = shift;
	my $chrb = shift;
	my $posb = shift;
	my ($region1, $region2);
	foreach my $key (keys %readslist)
	{
		my @coordinates = split (/__/, $key);
		if
		(
		($coordinates[0] eq $chra and abs($coordinates[1] - $posa) < $filter_mate_window and $coordinates[2] eq $chrb and abs($coordinates[3] - $posb) < $filter_mate_window)
		or
		($coordinates[0] eq $chrb and abs($coordinates[1] - $posb) < $filter_mate_window and $coordinates[2] eq $chra and abs($coordinates[3] - $posa) < $filter_mate_window)
		)
		{
			$region1 = $chra.':'.($posa-$filter_mate_window).'-'.($posa+$filter_mate_window);
			$region2 = $chrb.':'.($posb-$filter_mate_window).'-'.($posb+$filter_mate_window);
			#print "$key\t$region1\t$region2\n@{$readslist{$key}}\n";
			foreach (@{$readslist{$key}})
			{
				my @readname = split (/_/, $_); 
				my $newline;
				open TEMP, "$samtools_command view -X $bamfile $region1|grep $readname[0]|";
				$newline = <TEMP>;
				chomp $newline;
				close TEMP;
				if ($newline)
				{
					#print "$newline\n";
					return 1;
				}
				open TEMP, "$samtools_command view -X $bamfile $region2|grep $readname[0]|";
				$newline = <TEMP>;
				chomp $newline;
				close TEMP;
				if ($newline)
				{
					#print "$newline\n";
					return 1;
				}
			}
		}
	}
}

sub local_align
{
	my $cluster_id = shift;
	my $chra = shift;
	my $posa = shift;
	my $chrb = shift;
	my $posb = shift;
	my $al_cutoff = 500;
	my $overlapping = 0.5;
	my $align_file = $srout_dir.$cluster_id.'.srout';
	my $reads_file = $tmpdir.'/'.$cluster_id.'.fa';
	my $outfile = $tmpdir.'/'.$cluster_id.'.psl';
	#print "$cluster_id\t$align_file\n";
	
	my ($newline, @reads);
	open TEMP, "<$align_file";
	$newline = <TEMP>; chomp $newline;
	return 0 unless ($newline);
	$newline = <TEMP>;
	while ($newline = <TEMP>)
	{
		chomp $newline;
		last unless ($newline);
		$newline =~ s/ //g;
		push @reads, $newline;
	}
	close TEMP;
	
	open READS, ">$reads_file";
	my $toprint;
	for (my $i=0; $reads[$i]; $i++)
	{
		$toprint .= '>'.$i."\n".$reads[$i]."\n";
	}
	print READS "$toprint";
	close READS;
	
	system "$gfClient_command $blat_server $blat_port / $reads_file $outfile -nohead -minScore=10 -maxIntron=10 >/dev/null";
	
	my ($align_line, %alg, $last_readname);
	my $i=0;
	open ALIGN, "<$outfile";
	while ($align_line = <ALIGN>)
	{
		chomp $align_line;
		my @data = split (/\t/, $align_line);
		if ($last_readname ne $data[9])
		{
			$i = 0;
		}
		$last_readname = $data[9];
		my @size = split (",", $data[18]);
		my @query_start = split (",", $data[19]);
		my @target_start = split (",", $data[20]);
		if (@size > 1)
		{
			my $total_size;
			foreach (@size)
			{
				$total_size += $_;
			}
			for (my $j=0; defined($size[$j]); $j++)
			{
				next if ($size[$j] < 10);
				@{$alg{$data[9]}{$i}} = @data;
				$alg{$data[9]}{$i}[1] = int($size[$j]/$total_size*$data[1]+0.5);
				$alg{$data[9]}{$i}[0] = $size[$j] - $alg{$data[9]}{$i}[1];
				$alg{$data[9]}{$i}[17] = 1;
				$alg{$data[9]}{$i}[18] = $size[$j];
				$alg{$data[9]}{$i}[19] = $query_start[$j];
				$alg{$data[9]}{$i}[20] = $target_start[$j];
				if ($alg{$data[9]}{$i}[8] eq '+')
				{
					$alg{$data[9]}{$i}[11] = $alg{$data[9]}{$i}[19];
					$alg{$data[9]}{$i}[12] = $alg{$data[9]}{$i}[19] + $alg{$data[9]}{$i}[18];
				}
				else
				{
					$alg{$data[9]}{$i}[11] = $alg{$data[9]}{$i}[10] - ($alg{$data[9]}{$i}[19] + $alg{$data[9]}{$i}[18]);
					$alg{$data[9]}{$i}[12] = $alg{$data[9]}{$i}[10] - $alg{$data[9]}{$i}[19];
				}
				$alg{$data[9]}{$i}[15] = $alg{$data[9]}{$i}[20];
				$alg{$data[9]}{$i}[16] = $alg{$data[9]}{$i}[20] + $alg{$data[9]}{$i}[18];
				$alg{$data[9]}{$i}[11] = $alg{$data[9]}{$i}[11] + 1;
				$alg{$data[9]}{$i}[15] = $alg{$data[9]}{$i}[15] + 1;
				#print "@{$alg{$data[9]}{$i}}\n";
				$i++;
			}
		}
		else
		{
			@{$alg{$data[9]}{$i}} = @data;
			$alg{$data[9]}{$i}[11] = $alg{$data[9]}{$i}[11] + 1;
			$alg{$data[9]}{$i}[15] = $alg{$data[9]}{$i}[15] + 1;
			#print "@{$alg{$data[9]}{$i}}\n";
			$i++;
		}
	}
	close ALIGN;
	
	for (my $read=0; $alg{$read}{0}; $read++)
	{
		#print "@{$alg{$read}{0}}\n";
		my %alg_top2;
		my $maxk = keys %{$alg{$read}};
		if ($maxk > 1)
		{
			# initialize first entry
			@{$alg_top2{0}} = @{$alg{$read}{0}};
			delete $alg{$read}{0};
			# get best hit for first entry
			for (my $k=1;$k<$maxk;$k++)
			{
				my ($overlap1, $overlap2) = &overlap($alg{$read}{$k}[11], $alg{$read}{$k}[12], $alg_top2{0}[11], $alg_top2{0}[12]);
				#print "1 @{$alg_top2{0}}\n2 @{$alg{$read}{$k}}\n";
				#print "$alg{$read}{$k}[12] - $alg{$read}{$k}[11]\t$alg_top2{0}[12] - $alg_top2{0}[11]\n";
				if (($overlap2-$overlap1)/($alg{$read}{$k}[12] - $alg{$read}{$k}[11]) >= $overlapping or ($overlap2-$overlap1)/($alg_top2{0}[12] - $alg_top2{0}[11]) >= $overlapping)
				{
					if ($alg{$read}{$k}[0] > $alg_top2{0}[0])
					{
						@{$alg_top2{0}} = @{$alg{$read}{$k}};
					}
					delete $alg{$read}{$k};
				}
			}
			for (my $k=1;$k<$maxk;$k++)
			{
				next unless ($alg{$read}{$k}[0]);
				my ($overlap1, $overlap2) = &overlap($alg{$read}{$k}[11], $alg{$read}{$k}[12], $alg_top2{0}[11], $alg_top2{0}[12]);
				if (($overlap2-$overlap1)/($alg{$read}{$k}[12] - $alg{$read}{$k}[11]) >= $overlapping or ($overlap2-$overlap1)/($alg_top2{0}[12] - $alg_top2{0}[11]) >= $overlapping)
				{
					#print "3 @{$alg_top2{0}}\n4 @{$alg{$read}{$k}}\n";
					if ($alg{$read}{$k}[0] > $alg_top2{0}[0])
					{
						@{$alg_top2{0}} = @{$alg{$read}{$k}};
					}
					delete $alg{$read}{$k};
				}
			}
			# initialize second entry
			for (my $k=1;$k<$maxk;$k++)
			{
				if (defined($alg{$read}{$k}[0]))
				{
					@{$alg_top2{1}} = @{$alg{$read}{$k}};
					delete $alg{$read}{$k};
					last;
				}
			}
			for (my $k=2;$k<$maxk;$k++)
			{
				next unless ($alg{$read}{$k}[0]);
				my ($overlap1, $overlap2) = &overlap($alg{$read}{$k}[11], $alg{$read}{$k}[12], $alg_top2{1}[11], $alg_top2{1}[12]);
				if (($overlap2-$overlap1)/($alg{$read}{$k}[12] - $alg{$read}{$k}[11]) >= $overlapping or ($overlap2-$overlap1)/($alg_top2{1}[12] - $alg_top2{1}[11]) >= $overlapping)
				{
					if ($alg{$read}{$k}[0] > $alg_top2{1}[0])
					{
						@{$alg_top2{1}} = @{$alg{$read}{$k}};
					}
					delete $alg{$read}{$k};
				}
			}
			for (my $k=2;$k<$maxk;$k++)
			{
				next unless ($alg{$read}{$k}[0]);
				my ($overlap1, $overlap2) = &overlap($alg{$read}{$k}[11], $alg{$read}{$k}[12], $alg_top2{1}[11], $alg_top2{1}[12]);
				if (($overlap2-$overlap1)/($alg{$read}{$k}[12] - $alg{$read}{$k}[11]) >= $overlapping or ($overlap2-$overlap1)/($alg_top2{1}[12] - $alg_top2{1}[11]) >= $overlapping)
				{
					if ($alg{$read}{$k}[0] > $alg_top2{1}[0])
					{
						@{$alg_top2{1}} = @{$alg{$read}{$k}};
					}
					delete $alg{$read}{$k};
				}
			}
			#print "$cluster_id\n@{$alg_top2{0}}\n" if (defined($alg_top2{0})); print "@{$alg_top2{1}}\n" if (defined($alg_top2{1}));
		}
		if 
		(
			($alg_top2{0}[13] eq $chra and abs($alg_top2{0}[15] - $posa) < $al_cutoff and $alg_top2{1}[13] eq $chrb and abs($alg_top2{1}[15] - $posb) < $al_cutoff)
			or
			($alg_top2{0}[13] eq $chrb and abs($alg_top2{0}[15] - $posb) < $al_cutoff and $alg_top2{1}[13] eq $chra and abs($alg_top2{1}[15] - $posa) < $al_cutoff)
		)
		{
			return 1;
		}
	}
}

sub dcn
{
	my $chr = shift;
	my $pos = shift;
	my $orient = shift;
	
	my ($p1, $p2);
	if ($orient == 1)
	{
		$p1 = $pos - $window_size;
		$p2 = $pos;
	}
	elsif ($orient == -1)
	{
		$p1 = $pos;
		$p2 = $pos + $window_size;
	}
	else
	{
		$p1 = $pos - $window_size;
		$p2 = $pos + $window_size;
	}
	my $region = $chr.':'.$p1.'-'.$p2;
	
	my ($discord, $total_pair);
	my $newline1;
	open TEMP, "$samtools_command view -X $bamfile $region|";
	while ($newline1 = <TEMP>)
	{
		chomp $newline1;
		my $rg = 'none';
		if ($newline1 =~ /RG:Z:(\S{1,500})/)
		{
			$rg = $1;
			next if ($blackrg{$rg});
		}
		my @tempdata = split ('\t', $newline1);
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
		next unless ($tempdata[1] =~ /p/);
		
		my $strand = 1;
		$strand = -1 if ($tempdata[1] =~ /r/);
		my $mstrand = 1;
		$mstrand = -1 if ($tempdata[1] =~ /R/);
		#next if ($orient != $strand);
		
		$total_pair++;
		if ($tempdata[6] eq '=')
		{
			if ($strand == $mstrand and $tempdata[8] > 0)
			{
				#print "$rg\t$is{$rg}{'isu'}\t$newline1\n";
				$discord++;
			}
			else
			{
				if ($tempdata[8] > $is{$rg}{'isu'})
				{
					#print "$rg\t$is{$rg}{'isu'}\t$newline1\n";
					$discord++;
				}
			}
		}
		else
		{
			#print "$rg\t$is{$rg}{'isu'}\t$newline1\n";
			$discord++;
		}
	}
	close TEMP;
	$total_pair = 0.01 unless ($total_pair);
	
	#print "$region\t$discord\t$total_pair\n";
	return 1 if ($discord/$total_pair > $match_normal_disc);
	return 0;
}

sub nud
{
	my $chr = shift;
	my $pos = shift;
	
	my $p1 = $pos - $nu_window;
	my $p2 = $pos + $nu_window;
	my $region = $chr.':'.$p1.'-'.$p2;
	
	#system "samtools view -X $bamfile $region|grep XT:A:U|wc";
	
	my ($ratio, $uniq, $nonuniq);
	my $newline1;
	if ($min_mapq)
	{
		open TEMP, "$samtools_command view -X $bamfile $region|";
		while ($newline1 = <TEMP>)
		{
			chomp $newline1;
			my $rg = 'none';
			if ($newline1 =~ /RG:Z:(\S{1,500})/)
			{
				$rg = $1;
				next if ($blackrg{$rg});
			}
			my @query = split ('\t', $newline1);
			if ($query[4] > $min_mapq)
			{
				$uniq++;
			}
			else
			{
				$nonuniq++;
			}
		}
		close TEMP;
	}
	else
	{
		open TEMP, "$samtools_command view -X $bamfile $region|";
		while ($newline1 = <TEMP>)
		{
			chomp $newline1;
			my $rg = 'none';
			if ($newline1 =~ /RG:Z:(\S{1,500})/)
			{
				$rg = $1;
				next if ($blackrg{$rg});
			}
			if ($newline1 =~ /XT:A:U/)
			{
				$uniq++;
			}
			if ($newline1 =~ /XT:A:R/)
			{
				$nonuniq++;
			}
		}
		close TEMP;
	}
	
	return 0 if (($uniq + $nonuniq) == 0);
	$ratio = $nonuniq/($uniq + $nonuniq);
	#print "$uniq\t$nonuniq\t$ratio\n";
	return 1 if ($ratio > $nu_ratio);
	return 0;
}

sub scd
{
	my $chr = shift;
	my $pos = shift;
	
	my $p1 = $pos - $sc_window;
	my $p2 = $pos + $sc_window;
	my $region = $chr.':'.$p1.'-'.$p2;
	
	
	my ($newline1, $sc_number);
	open TEMP, "$samtools_command view -X $bamfile $region|";
	while ($newline1 = <TEMP>)
	{
		chomp $newline1;
		my $rg = 'none';
		if ($newline1 =~ /RG:Z:(\S{1,500})/)
		{
			$rg = $1;
			next if ($blackrg{$rg});
		}
		my @query = split ('\t', $newline1);
		my @clipped;
		while ($query[5] =~ /(\d{1,4})S/g)
		{
			push @clipped, $1;
			
		}
		next unless ($clipped[0]);
		my $clip_side; 
		# 0: both side, 1: left, 2: right
		if ($query[5] =~ /^(\d{1,4})S.*(\d{1,4})S$/)
		{
			$clip_side = 0;
		}
		elsif ($query[5] =~ /^(\d{1,4})S/)
		{
			$clip_side = 1;
		}
		elsif ($query[5] =~ /(\d{1,4})S$/)
		{
			$clip_side = 2;
		}
		if ($clip_side == 0) 
		{
			my $position = $query[3] + length($query[9]) - $clipped[0] - $clipped[1];
			if ($position > $p1 and $position < $p2 and $clipped[1] >= $clip_size_cutoff)
			{
				my $qual = substr ($query[10], -$clipped[1]);
				my $high_qual = &fq($qual);
				if ($high_qual >= $clip_size_cutoff)
				{
					$sc_number++ ;
					next;
				}
			}
			if ($query[3] > $p1 and $query[3] < $p2 and $clipped[0] >= $clip_size_cutoff)
			{
				my $qual = substr ($query[10], 0, $clipped[0]);
				my $high_qual = &fq($qual);
				if ($high_qual >= $clip_size_cutoff)
				{
					$sc_number++;
					next;
				}
			}
			#print "0 $query[1]\t$query[3]\t$query[5]\t@clipped\t$position\n";
		}
		if ($clip_side == 1)
		{
			my $position = $query[3];
			if ($position > $p1 and $position < $p2 and $clipped[0] >= $clip_size_cutoff)
			{
				my $qual = substr ($query[10], 0, $clipped[0]);
				my $high_qual = &fq($qual);
				if ($high_qual >= $clip_size_cutoff)
				{
					$sc_number++;
				}
			}
			#print "1 $query[1]\t$query[3]\t$query[5]\t@clipped\n";
		}
		if ($clip_side == 2)
		{
			my $position = $query[3] + length($query[9]) - $clipped[0];
			if ($position > $p1 and $position < $p2 and $clipped[0] >= $clip_size_cutoff)
			{
				my $qual = substr ($query[10], -$clipped[0]);
				my $high_qual = &fq($qual);
				if ($high_qual >= $clip_size_cutoff)
				{
					$sc_number++;
				}
			}
			#print "2 $query[1]\t$query[3]\t$query[5]\t@clipped\t$position\n";
		}
	}
	close TEMP;
		
	#print "$sc_number\n";
	return 1 if ($sc_number > $sc_cutoff);
	return 0;
}

sub fq
{
	my $qual = shift;
	
	my @qual = split(//,$qual);
	my $high_qual;
	foreach (@qual)
	{
		$high_qual++ if (ord($_)-$qual_code >= $qual_cutoff);
	}
	#print "$high_qual\n";
	return $high_qual;
}

sub dt
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
	#print STDERR "$region1\n$region2\n";
	
	my $discord = 0;
	my $newline1;
	open TEMP, "$samtools_command view -X $bamfile $region1|";
	while ($newline1 = <TEMP>)
	{
		chomp $newline1;
		my $rg = 'none';
		if ($newline1 =~ /RG:Z:(\S{1,500})/)
		{
			$rg = $1;
			next if ($blackrg{$rg});
		}
		my @tempdata = split ('\t', $newline1);
		#print "@tempdata\n";
		next unless ($tempdata[1] =~ /p/);
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
					#print "@tempdata\n$ori1\t$ori2\t$pos1\t$pos2\t$tempdata[3]\t$tempdata[7]\t$disc1 and $disc2\n";
					if ($disc1 and $disc2)
					{
						$discord++;
					}
				}
			}
		}
	}
	close TEMP;
	#print "$discord\n";
	return 0 if ($discord);
	return 1;
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
	die "somatic_sv.pl [options]
	-i FILE	input file, required
	-o FILE	output file, required
	-D FLT	standard deviation cutoff to call discordant read pairs, use the same value as used in d option of meerkat.pl, default 3
	-I FILE	isinfo file from Meerkat run, required if n or e option is turned on
	-K FILE	file name of read group to be ignored, one read group ID per line, same as R option in pre_process.pl and meerkat.pl. If all the read groups are of high quality, you don't need to specify this option. If there is any blacklist rg and any of the options n, u, f, e is enabled, K option is required.
	-F STR	name of folder contains all *.discord files from normal genomes to filter germline events, we recommend to filter against all normal genomes from one tumor type
	-x INT	number of discordant read pairs supporting the same event to be used as filter from *.discord file, default 2
	-l INT	distance of breakpoints to filter germline events from *.discord file, default 500
	-s INT	minimum size of deletions, default 100
	-E INT	filter TEI for deletions, these are typically germline events, default 1
	-d INT	max homology allowed for deletion and intra-chr events, 0 to disable, default 100
	-t INT	max homology allowed for inter chromosomal translocation events, 0 to disable, default 100
	-m INT	filter events that both breakpoints fall in satellite or simple repeats, default 1
	-n INT	filter by total number of discordant read pairs in matched normal genome, if certain number of discordant read pairs are observed in the given genome, discard the event, if enable this option, B and I options are required, default 0
	-y INT	maximum fraction of discordant pairs in normal bam, default 0.1
	-u INT	filter by non-uniq mapped reads in matched normal genome, determined by XT tag or mapping quality, if too many non-uniq mapped reads are observed in given genome, discard the event, if enable this option, B option is required, default 0
	-v INT	window size to look for  non-uniq mapped reads in normal bam file, default 100
	-r FLT	cutoff of ratio of non-uniq mapped reads to all mappable reads, default 0.25
	-f INT	filter by soft-clipped reads in matched normal genome, if certain number of soft-clipped reads are observed in the given genome, discard the event, if enable this option, B option is required, default 0
	-g INT window size to look for soft-clipped reads in normal bam file, default 10
	-j INT	cutoff for number of  soft-clipped reads in normal bam file, default 3
	-e INT	filter by discordant read pairs in tumor genome, if certain breakpoint has no discordant read pair, it is an artifact of split read mapping, discard the event, if enable this option, B and I options are required, default 0
	-B FILE	bam file
	-k INT	window size to look for  discordant read pairs in tumor bam file, default 800
	-z INT	to enable parameter p, q and P, default 0
	-p INT	filter by number of supporting discordant read pairs, default 3
	-q INT	filter by number of supporting split reads, default 1
	-P INT	filter by sum of supporting discordant read pairs and supporting split reads, default 6
	-M INT	filter by mate position of split reads, for reads covering breakpoint junctions, their mate should map near breakpoints, if enable this option, B option is required, default 0
	-N INT	window size to look for mate in tumor bam file, default 2000
	-C FILE	bp_reads file from Meerkat output, required if M option is enabled
	-Q INT	minimum mapping quality for reads to be used, default 0
	-b INT	allele frequency cutoff
	-V STR	blat server, if enable option V and T, will filter by blat the split reads against whole genome, and require both side of breakpoints to be best hit. One needs to set up a blat server before using this filter (i.e. gfServer start 10.11.240.76 17777 /reference/hg18/hg18.2bit -stepSize=5)
	-T STR	blat port
	-L STR	/path/to/blat, path only, not the command, no need to specify if blat is in PATH
	-A STR	location of srout folder from Meerkat run.
	-S STR	/path/to/samtools, path only, not the command, no need to specify if samtools is in PATH
	-h help\n";
}