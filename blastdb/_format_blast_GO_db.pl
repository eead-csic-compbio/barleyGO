#!/usr/bin/perl -w
use strict;

# Creates a FASTA file with sequences and their GO annotations in this order:
# barley (B2G), plaza monocots, Aegilops tauschii

# from http://bioinfo.cipf.es/node/839
# please edit paths
my $B2GSEQS  = 'b2gfar/barley_b2gfar.faa';
my $B2GANNOT = 'b2gfar/4513.annot';

# from http://bioinformatics.psb.ugent.be/plaza/versions/plaza_v3_monocots/download/index
# please edit paths
my $PLAZASEQS  = 'plaza/monocots_plaza3.faa';
my $PLAZAANNOT = 'plaza/monocots_plaza3.go.csv';

# from http://plants.ensembl.org/index.html
# please edit paths 
my $ATAUSCHIISEQS  = 'ensembl_plants/Aegilops_tauschii.GCA_000347335.1.21.pep.all.fa';
my $ATAUSCHIIANNOT = 'ensembl_plants/go_ensembl_aegilops_tauschii.gaf';


#####################################

my (%syn,%GO,%FAA,%annot,@IDs,$id,$annot);

open(GO,$B2GANNOT) || die "# cannot find GO file ($B2GANNOT)\n";
while(<GO>)
{
	#1092232	GO:0009627	gi|167048|gb|AAA32958.1|	genbank	167048	AAA32958.1	gb	75
	#57656429	GO:0006355	MLOC_6708.5	Hordeum vulgar
	$_ =~ s/"//g;
	my @data = split(/\t/,$_);
	$GO{$data[0]} .= "$data[1](IEA),"; 
	
	# add IBSC gene names, as these were used in _format_barley_B2G_annotations.pl
	if(/Hordeum vulgare/)
	{
		$syn{$data[0]} .= "$data[2]," if(!$syn{$data[0]} || $syn{$data[0]} !~ /$data[2]/); 
	}	
}
close(GO); 

open(SEQS,$B2GSEQS) || die "# cannot find file ($B2GSEQS)\n";
while(<SEQS>)
{
	#>1092232
	chomp;
	if(/^>(\S+)/)
	{
		$id = $1;
		if($syn{$id}){ $annot = 'B2G:'.$syn{$id}.'|'; }
		else{ $annot = '' }
		if($GO{$id}){ $annot .= $GO{$id} }
		$annot{$id} = $annot; 
		push(@IDs,$id);
	}
	else{ $FAA{$id} .= $_; }
}
close(SEQS); 

foreach $id (@IDs){ printf(">%s|%d|%s\n%s\n",$id,length($FAA{$id}),$annot{$id},$FAA{$id}) if($FAA{$id}) }

%GO = %FAA = %annot = @IDs = %syn = ();
open(GO,$PLAZAANNOT) || die "# cannot find GO file ($PLAZAANNOT)\n";
while(<GO>)
{
  #"id";"species";"gene_id";"go";"evidence";"go_source";"provider";"comment";"is_shown"
  #"1";"cpa";"CP00074G00650";"GO:0005634";"IEA";"primary";"from_annotation";"";"1"
	$_ =~ s/"//g;
	my @data = split(/;/,$_);
	$GO{$data[2]} .= "$data[3]($data[4]),";
}
close(GO); 

open(SEQS,$PLAZASEQS) || die "# cannot find file ($PLAZASEQS)\n";
while(<SEQS>)
{
	#>AT5G16970
	chomp;
	if(/^>(\S+)/)
	{
		$id = $1;
		$annot = $GO{$id} || '';
		$annot{$id} = $annot;
		push(@IDs,$id);
	}
	else{ $FAA{$id} .= $_; }
}
close(SEQS); 

foreach $id (@IDs){ printf(">%s|%d|%s\n%s\n",$id,length($FAA{$id}),$annot{$id},$FAA{$id}) if($FAA{$id}) }



%GO = %FAA = %annot = @IDs = ();
open(GO,$ATAUSCHIIANNOT) || die "# cannot find GO file ($ATAUSCHIIANNOT)\n";
while(<GO>)
{
	#ensembl_aegilops_tauschii	F775_11690			GO:0006812	GR_REF:8396	IEA					gene	37682	20140110	Gramene
	my @data = split(/\t/,$_);
	next if(scalar(@data)<7);
	$GO{$data[1]} .= "$data[4]($data[6]),"; 
}
close(GO); 

open(SEQS,$ATAUSCHIISEQS) || die "# cannot find file ($ATAUSCHIISEQS)\n";
while(<SEQS>)
{
	#>EMT30525 pep:novel supercontig:GCA_000347335.1:Scaffold5955:89268:90667:-1 gene:F775_24069 transc
	chomp;
	if(/^>(\S+).*?gene:(\S+)/)
	{
		$id = $1;
		$annot = $GO{$2} || '';
		$annot{$id} = $annot;
		push(@IDs,$id);
	}
	else{ $FAA{$id} .= $_; }
}
close(SEQS); 

foreach $id (@IDs){ printf(">%s|%d|%s\n%s\n",$id,length($FAA{$id}),$annot{$id},$FAA{$id}) if($FAA{$id}) }
