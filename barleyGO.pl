#!/usr/bin/perl 
# EEAD-CSIC 2014 Bruno Contreras Moreira
# http://www.eead.csic.es/compbio/material/bioperl/node52.html
# http://bioinfoperl.blogspot.com.es/2011/08/decodificando-la-gene-ontology-cc.html
# http://bioinfoperl.blogspot.com.es/2013/04/splitblastpl-real-multicore-blast.html

use strict;
use IO::File;
use Getopt::Std;
use FindBin '$Bin';
use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OboParser;
use GO::TermFinderReport::Text;
use GO::Utils::File    qw (GenesFromFile);
use GO::Utils::General qw (CategorizeGenes);

my $GOfile            = $Bin.'/gene_ontology.1_2.obo'; 

my @annot_files       = ( '_barleyGO.B2G.go','_barleyGO.annot.go' ); 
my @annot_files_descr = ( 'B2G-FAR annotation for barley IBSC genes PMID:21335611', 
									'original GO annotation from IBSC PMID:23075845' );
									
my $BLASTDBGOfile     = $Bin.'/_grassesAthaGOdb.faa'; 								
								
my $GETANNOTEXE       = $Bin.'/_get_go_annotation';
my $SPLITBLASTEXE     = $Bin.'/_split_blast.pl';
my $BLASTPATH         = '/home/contrera/soft/ncbi-blast-2.2.27+/bin/'; #paths $PATH by default/home/contrera/soft/ncbi-blast-2.2.27+/bin/
my $BLASTPEXE         = ($BLASTPATH && -s $BLASTPATH) ? $BLASTPATH.'/blastp' : 'blastp';
my $BLASTXEXE         = ($BLASTPATH && -s $BLASTPATH) ? $BLASTPATH.'/blastx' : 'blastx';

my $BLASTPARAMS       = '-max_target_seqs 1 -seg yes -soft_masking true';
my $BLASTFORMAT       = '6 qseqid sseqid qlen qstart qend pident evalue';

###########################################################

my ($INP_Pvalue,$INP_file,$INP_term,$INP_verb,$INP_threads) = (0,'','',0,2);  
my ($INP_seqfile,$INP_index,$INP_fdr,$INP_annot,$INP_nucl) = ('',0,0.01,0,0); 
my ($INP_cover,$INP_ident) = (50,70); 

my ($annot_file,%opts) = ($Bin.'/'.$annot_files[0]);

getopts('SnhavC:I:T:A:F:f:p:f:t:s:', \%opts);
if(($opts{'h'})||(scalar(keys(%opts))==0))
{
	print   "\n[options]:\n\n";
	print   "-h \t this message\n";
	print   "-t \t GO term to analyze, or file with 1 terms per line   (optional, example: -t GO:0003924)\n";
	print   "-f \t file with gene identifiers for GO enrichment test   (optional, choose either -f,-s or -t)\n";
	print   "-s \t FASTA file with sequences to be assigned GO terms   (optional, choose either -f,-s or -t)\n";
	print   "   \t and matched IBSC genes from annotations in $BLASTDBGOfile\n\n";
	
	print   "Options that affect BLAST searches conducted with -s:\n";
	print   "-n \t FASTA file contains nucleotide sequences            (optional, default: protein)\n";
	print   "-C \t minimum \%sequence coverage of best hit              (optional, default: -C $INP_cover)\n";
	print   "-I \t minimum \%sequence identity to best hit              (optional, default: -I $INP_ident)\n";
	print   "-T \t number of threads for BLAST searches                (optional, default: -T $INP_threads)\n\n";
	
	print   "Options that control GO enrichment tests conducted with -f:\n";
	print   "-F \t False Discovery Rate (FDR) for enrichment test      (optional, default: -F $INP_fdr)\n";
	print   "-p \t Bonferroni-corrected p-value for enrichment test    (optional, overrides FDR, example -p 0.01)\n";
	print   "-S \t show available GO annotation files and exit         (optional)\n";
	print   "-A \t use indicated GO annotation file                    (optional, default: -A 1\n";
	print   "-a \t show annotation of genes in -f file                 (optional, requires -f )\n";
	print   "-v \t increase verbosity                                  (optional)\n";
	exit(-1);
}

if(defined($opts{'S'}))
{ 
	foreach my $annot (0 .. $#annot_files)
	{
		printf(" %d:\t%s\t%s\n",$annot+1,$annot_files[$annot],$annot_files_descr[$annot]);
	}
	exit(0);
}

if(defined($opts{'f'}))
{
	if(-e $opts{'f'})
	{ 
		$INP_file = $opts{'f'};
		if(defined($opts{'a'})){ $INP_annot = 1 }
	}
	else{ die "# need a valid input file, cannot find $opts{'f'}\n" }
}
elsif(defined($opts{'s'}))
{
	if(-e $opts{'s'}){ $INP_seqfile = $opts{'s'} }
	else{ die "# need a valid file with input FASTA sequences, cannot find $opts{'s'}\n" }
	
	if(defined($opts{'n'})){ $INP_nucl = 1 }
	
	# check blast binary
	my $binary = $BLASTPEXE;
	if($INP_nucl){ $binary = $BLASTXEXE }
	
	my $blastStatus = `which $binary 2>&1`;
	if(!$blastStatus)
	{
		die "# WARNING: make sure that BLAST ($binary) is on your \$PATH (check README file)\n";
	}
}
elsif(defined($opts{'t'}))
{
	$INP_term = $opts{'t'};
	if($INP_term =~ /^\d+$/){ $INP_term = 'GO:'.$INP_term }
}

if(defined($opts{'p'}) && $opts{'p'} > 0)
{
	$INP_Pvalue = $opts{'p'};
}
elsif(defined($opts{'F'}) && $opts{'F'} > 0)
{
	$INP_fdr = $opts{'F'};
}

if(defined($opts{'T'}) && $opts{'T'} > 0)
{
	$INP_threads = $opts{'T'};
}

if(defined($opts{'C'}) && $opts{'C'} >= 0 && $opts{'C'} <= 100)
{
	$INP_cover = $opts{'C'};
}

if(defined($opts{'I'}) && $opts{'I'} >= 0 && $opts{'I'} <= 100)
{
	$INP_ident = $opts{'I'};
}

if(defined($opts{'v'})){ $INP_verb = 1 }

if(defined($opts{'A'}) && $opts{'A'} > 0 && $opts{'A'} <= scalar(@annot_files))
{ 
	$annot_file = $Bin.'/'.$annot_files[$opts{'A'}-1];
}

#############################################

if($INP_annot)
{
	my (%annot,$id,$ann);
	open(ANNOT,$annot_file);
	while(<ANNOT>)
	{
			#IBSC2012	AK248138.1	Cytochrome P450		GO:0005506	PMID:23075845	IEA		F	Cytochrome P450		
			my @data = split(/\t/,$_);
			$annot{$data[1]} = $data[2];
	}
	close(ANNOT);

	print "# annotations of genes in $INP_file:\n";
	open(LIST,$INP_file);
	while(<LIST>)
	{
		$id = (split)[0];
		$ann = $annot{$id} || '';
		print "$id\t$ann\n";
	}
	close(LIST);
	print "\n";
}

if($INP_seqfile)
{
	my ($cover,$besthit,$slen,$annot,$s,@syns);
	my $GOseqfile = $INP_seqfile . '.go';
	my $matchedISBCgenesfile = $INP_seqfile . '.ISBC.list';
	my ($nseqs,$nbesthits,$nIBSChits,$isProtein) = (0,0,0,0);
	
	open(SEQS,$INP_seqfile) || die "# cannot open file $INP_seqfile\n";
	while(<SEQS>)
	{
		chomp;
		if(/^>/){ $nseqs++; }
		elsif(/[^ACTGNX-]/){ $isProtein = 1 }
	}
	close(SEQS); 
	
	print "# input sequences: $nseqs (-s $INP_seqfile)\n";
	if($isProtein){ print "# input seems to be: amino acid sequences\n"; }
	else{ print "# input seems to be: nucleotide sequences\n"; }
	
	if(!$isProtein && !$INP_nucl)
	{
		die "# WARNING: please use option -n if your input are nucleotide sequences ($isProtein,$INP_nucl)\n"
	}
	elsif($isProtein && $INP_nucl)
	{
		die "# WARNING: please don't use option -n if your input are protein sequences ($isProtein,$INP_nucl)\n"
	}
	
	print "# GO sequence library: $BLASTDBGOfile\n\n";
	
	open(GOANNOT,">$GOseqfile") || die "# cannot create file $GOseqfile\n";
	print GOANNOT "# input sequence file: -s $INP_seqfile\n";
	print GOANNOT "# GO sequence library: $BLASTDBGOfile\n";
	print GOANNOT "# INP_ident=$INP_ident INP_cover=$INP_cover\n";
	print GOANNOT "# identifier\tGOterm\tbesthit\tcover(besthit)\t\%identity\tEvalue\n"; 
	
	open(IBSCMATCHES,">$matchedISBCgenesfile") || die "# cannot create file $matchedISBCgenesfile\n";
	print IBSCMATCHES "# input sequence file: -s $INP_seqfile\n";
	print IBSCMATCHES "# GO sequence library: $BLASTDBGOfile\n";
	print IBSCMATCHES "# INP_ident=$INP_ident INP_cover=$INP_cover\n";
	print IBSCMATCHES "# ISBC best hits of input sequences\n"; 
	
	my $blast_command = sprintf("%s -query %s -db %s -outfmt \'\"%s\"\' %s",
		$BLASTPEXE,$INP_seqfile,$BLASTDBGOfile,$BLASTFORMAT,$BLASTPARAMS);
		
	if($INP_nucl)
	{
		$blast_command = sprintf("%s -query %s -db %s -outfmt \'\"%s\"\' %s",
			$BLASTXEXE, $INP_seqfile,$BLASTDBGOfile,$BLASTFORMAT,$BLASTPARAMS);
	}
	
	print "# running BLAST jobs...\n\n";
	
	open(BLAST,"$SPLITBLASTEXE $INP_threads 100 $blast_command 2>&1 |") ||
		die "# cannot run $SPLITBLASTEXE $INP_threads 100 $blast_command, exit\n";
	while(<BLAST>)
	{
		#MLOC_31.1	SB05G021010|222|GO:0005840(IEA),GO:0006412(IEA),	66	22.73	7.8
		#MLOC_44.1	56334045|400|B2G:AK362988,|GO:0005794(IEA),GO:0003824(IEA),GO:0016020(IEA),	95	82.11	5e-36
		next if($_ =~ /^#/ || $_ =~ /^$/);
		chomp;
		my @data = split(/\t/,$_); #qseqid sseqid qlen qstart qend pident evalue
		($besthit,$slen,$annot) = split(/\|/,$data[1],3);
		$cover  = sprintf("%1.1f",100*(abs($data[4]-$data[3])/$data[2]));
		
		print "#BLASTHIT: $_ $cover\n" if($INP_verb);
		
		next if($INP_ident && $data[5] < $INP_ident);
		next if($INP_cover && $cover < $INP_cover);
		
		if($annot =~ /GO:\d+/)
		{ 
			while($annot =~ /(GO:\d+)/g)
			{ 
				print GOANNOT "$data[0]\t$1\t$besthit\t$cover\t$data[5]\t$data[6]\n"; 
			}
			$nbesthits++; 
		}
		
		if($annot =~ /B2G:(\S+?)\|/g)
		{ 
			@syns = split(/,/,$1);
			foreach $s (@syns)
			{
				print IBSCMATCHES  "$s\n"; 
			}
			$nIBSChits++;
		}
	}
	close(BLAST);
	
	close(IBSCMATCHES);
	close(GOANNOT);
		
	print "# output file with GO terms associated to input sequences: $GOseqfile ($nbesthits sequences)\n\n";
	print "# output file with IBSC genes matched in input sequences: $matchedISBCgenesfile ($nIBSChits sequences)\n\n";
}

if($INP_file)
{
	my ($fh,$num_assoc,@list,@notFound,@ambiguous);

	print "# parsing ontology file $GOfile ...\n";

	my $process   = GO::OntologyProvider::OboParser->new(ontologyFile => $GOfile, aspect => 'P');
	my $component = GO::OntologyProvider::OboParser->new(ontologyFile => $GOfile, aspect => 'C');
	my $function  = GO::OntologyProvider::OboParser->new(ontologyFile => $GOfile, aspect => 'F');

	print "# reading annotations from $annot_file ...\n";
	my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annot_file);

	my $termFinderP = GO::TermFinder->new(annotationProvider=> $annotation, ontologyProvider => $process, aspect => 'P');
	my $termFinderC = GO::TermFinder->new(annotationProvider=> $annotation, ontologyProvider => $component, aspect => 'C');
	my $termFinderF = GO::TermFinder->new(annotationProvider=> $annotation, ontologyProvider => $function, aspect => 'F');
	my $totalNumGenes = $termFinderF->totalNumGenes();

	my $report = GO::TermFinderReport::Text->new();
	my @genes = GenesFromFile($INP_file); 

	CategorizeGenes(annotation => $annotation, genes => \@genes, ambiguous => \@ambiguous,
		unambiguous => \@list, notFound => \@notFound);
		
	print "# parsing ".scalar(@genes)." gene identifiers of $INP_file ...\n";

	if(@list)
	{
		printf("# matched/recognized gene identifiers: %d\n",scalar(@list));
		if($INP_verb)
		{
			foreach my $gen (@list)
			{ 
				print $gen."\t",$annotation->standardNameByName($gen)."\n"; 
			}	
			print "\n";
		}	
	}
	else{ die "# WARNING: could not match any of the identifiers in $INP_file, please check them\n" }

	if(@ambiguous)
	{
		printf("# ambiguous gene identifiers: %d\n",scalar(@ambiguous)); 
		if($INP_verb){ print join("\n",@ambiguous)."\n\n"; }
	}

	if (@notFound)
	{
		printf("# unrecognized gene identifiers: %d\n",scalar(@notFound));
		if($INP_verb){ print join("\n",@notFound)."\n\n"; }
	}

	# check associations in all ontology branches: P,C,F
	foreach my $termFinder ($termFinderP, $termFinderC, $termFinderF) 
	{
		if(!$INP_Pvalue)
		{
			print "# associations for branch ".$termFinder->aspect()." (FDR<=$INP_fdr):\n";
		}
		else
		{
			print "# associations for branch ".$termFinder->aspect()." (corrected P-value<=$INP_Pvalue):\n";
		}	
		
		
		my @pvalues = $termFinder->findTerms(genes=>\@list,calculateFDR=>1); 
		if(@pvalues)
		{
			my @filtered_pvalues;
			foreach my $it (@pvalues)
			{ 
				if(!$INP_Pvalue)
				{
					if($it->{'FDR_RATE'} <= $INP_fdr){ push(@filtered_pvalues,$it) }
				}
				else
				{
					if($it->{'CORRECTED_PVALUE'} <= $INP_Pvalue){ push(@filtered_pvalues,$it) }
				} 
			}
		
			if(@filtered_pvalues)
			{
				$num_assoc = $report->print(
					pvalues  => \@filtered_pvalues,
					numGenes => scalar(@list),
					totalNum => $totalNumGenes,
					table    => 1 );
			}			
		
			if(!@filtered_pvalues || !$num_assoc)
			{
				print "# no significant associations found\n\n";
			}
		}
	}	
}
elsif($INP_term)
{
	if(-s $INP_term)
	{
		system("$GETANNOTEXE -f $INP_term");
	}
	else
	{
		system("$GETANNOTEXE -t $INP_term");
	}	
}
