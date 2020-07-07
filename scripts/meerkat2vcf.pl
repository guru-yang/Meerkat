#!/usr/bin/perl
# meerkat2vcf.pl
# convert meerkat output to vcf format
# Author: Lixing Yang, The Center for Biomedical Informatics, Harvard Medical School, Boston, MA, 02115, USA
# Email: lixing_yang@hms.harvard.edu

use strict;
use Getopt::Std;
use IO::File;
use Bio::DB::Fasta;

my %opts = ();
getopts("i:o:H:F:h", \%opts);
my $infile = $opts{i};
my $outfile = $opts{o};
my $vcfheader = $opts{H};
my $reference_path = $opts{F};
if (defined($opts{h}))
{
	&print_usage;
	die;
}
die "Input file not specified or not Meerkat variant file\n" unless ($infile =~ /variants$/);
die "Header file not specified\n" unless (defined($vcfheader));
if (!($outfile) and $infile =~ /.variants$/)
{
	$outfile = $`.'.vcf';
}
my $reference_db = Bio::DB::Fasta->new($reference_path) if (defined($reference_path));

system "cat $vcfheader > $outfile";

my $fh_out = IO::File->new(">>$outfile");
my $newline;
my $mateid = 1;
open FILE, "<$infile";
while ($newline = <FILE>)
{
	chomp $newline;
	my @data = split (/\t/, $newline);
	my $event = $data[0].'_'.$data[2];
	my @ids;
	if ($data[3] =~ /\//)
	{
		@ids = split (/\//, $data[3]);
	}
	if ($data[0] eq 'del')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[11], $data[3], $event);
	}
	elsif ($data[0] eq 'del_ins')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[14], $data[3], $event);
	}
	elsif ($data[0] eq 'invers_f')
	{
		$event = 'transl_intra_'.$data[2];
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 1;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[11], $data[3], $event);
	}
	elsif ($data[0] eq 'invers_r')
	{
		$event = 'transl_intra_'.$data[2];
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 2;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[11], $data[3], $event);
	}
	elsif ($data[0] eq 'tandem_dup')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 3;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7], $data[11], $data[3], $event);
	}
	elsif ($data[0] eq 'invers')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 1;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6]-1, $data[5], $data[7], $data[11], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 2;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[5], $data[7]+1, $data[11], $ids[1], $event);
	}
	elsif ($data[0] =~ /inssd/)
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 3;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[11], $data[16], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[10], $data[17], $ids[1], $event);
	}
	elsif ($data[0] =~ /inssu/)
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 3;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[10], $data[16], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[11], $data[17], $ids[1], $event);
	}
	elsif ($data[0] =~ /inso/)
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 1;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[11], $data[16], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 2;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[10], $data[17], $ids[1], $event);
	}
	elsif ($data[0] eq 'del_invers')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 1;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[11], $data[17], $ids[0], $event);
		
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 2;
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[10], $data[18], $ids[1], $event);
	}
	elsif ($data[0] =~ /inss/)
	{
		if ($data[5] lt $data[9])
		{
			my $mateida = 'bnd_'.$mateid;
			$mateid++;
			my $mateidb = 'bnd_'.$mateid;
			$mateid++;
			my $mate_type = 3;
			&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[10], $data[16], $ids[0], $event);
			
			my $mateida = 'bnd_'.$mateid;
			$mateid++;
			my $mateidb = 'bnd_'.$mateid;
			$mateid++;
			my $mate_type = 0;
			&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[11], $data[17], $ids[1], $event);
		}
		else
		{
			my $mateida = 'bnd_'.$mateid;
			$mateid++;
			my $mateidb = 'bnd_'.$mateid;
			$mateid++;
			my $mate_type = 3;
			&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[7], $data[9], $data[11], $data[17], $ids[0], $event);
			
			my $mateida = 'bnd_'.$mateid;
			$mateid++;
			my $mateidb = 'bnd_'.$mateid;
			$mateid++;
			my $mate_type = 0;
			&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[9], $data[10], $data[16], $ids[1], $event);
		}
	}
	elsif ($data[0] eq 'transl_inter')
	{
		my $mateida = 'bnd_'.$mateid;
		$mateid++;
		my $mateidb = 'bnd_'.$mateid;
		$mateid++;
		my $mate_type = 0;
		$mate_type = 1 if ($data[7] == 1 and $data[10] == 1);
		$mate_type = 2 if ($data[7] == -1 and $data[10] == -1);
		$mate_type = 3 if ($data[7] == -1 and $data[10] == 1);
		&bnd ($data[0], $mate_type, $mateida, $mateidb, $data[5], $data[6], $data[8], $data[9], $data[13], $data[3], $event);
	}
}
close FILE;
$fh_out->close;

sub bnd
{
	my $event_type = shift;
	my $mate_type = shift;
	my $mateida = shift;
	my $mateidb = shift;
	my $chra = shift;
	my $posa = shift;
	my $chrb = shift;
	my $posb = shift;
	my $gene_string = shift;
	my $ad = shift;
	my $event = shift;
	
	my ($genea, $geneb);
	if ($gene_string =~ /GN:(\S{0,20}),.*;(\S{0,20}),/)
	{
		($genea, $geneb) = ($1, $2);
		$genea = substr ($genea, 0, -1);
		$geneb = substr ($geneb, 0, -1);
	}
	#print "$gene_string\t$genea\t$geneb\n";
	
	my ($alta, $altb);
	my $nta = '.';
	my $ntb = '.';
	$nta = $reference_db->seq($chra, $posa => $posa) if (defined($reference_path));
	$ntb = $reference_db->seq($chrb, $posb => $posb) if (defined($reference_path));
	if ($mate_type == 0)
	{
		$altb = $nta.'['.$chrb.':'.$posb.'[';
		$alta = ']'.$chra.':'.$posa.']'.$ntb;
	}
	if ($mate_type == 1)
	{
		$altb = $nta.']'.$chrb.':'.$posb.']';
		$alta = $ntb.']'.$chra.':'.$posa.']';
	}
	if ($mate_type == 2)
	{
		$altb = '['.$chrb.':'.$posb.'['.$nta;
		$alta = '['.$chra.':'.$posa.'['.$ntb;
	}
	if ($mate_type == 3)
	{
		$altb = ']'.$chrb.':'.$posb.']'.$nta;
		$alta = $ntb.'['.$chra.':'.$posa.'[';
	}
	
	my $toprint = "$chra\t$posa\t$mateida\t$nta\t$altb\t100\tPASS\tSVTYPE=BND;MATEID=$mateidb;EVENT=$event;GENE=$genea;\tGT:AD:DP:SS:SSC:BQ\t.:.:.:.:.:.\t.:$ad:.:2:100:.\n";
	$toprint .= "$chrb\t$posb\t$mateidb\t$ntb\t$alta\t100\tPASS\tSVTYPE=BND;MATEID=$mateida;EVENT=$event;GENE=$geneb;\tGT:AD:DP:SS:SSC:BQ\t.:.:.:.:.:.\t.:$ad:.:2:100:.\n";
	print $fh_out $toprint;
	return;
}

sub print_usage
{
	die "Usage:\nperl ./scripts/meerkat2vcf.pl [options]
	-i FILE	input file name, required
	-o FILE	output file name
	-H FILE	header file name
	-F STR	/path/to/reference/fasta/files, path only, not the files
	-h help\n";
}
