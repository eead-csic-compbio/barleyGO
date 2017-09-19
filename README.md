# barleyGO

barleyGO is a Perl script to annotate barley sequences with Gene Ontology terms
inferred by homology. It uses the IBSC2012 barley GO annotation and supports both nucleotide
and peptide sequences. It should be trivial to modify it for any other species.

Authors Bruno Contreras Moreira (1,2), Carlos P Cantalapiedra (1), MJ Garcia-Pereira (1) 

    1. Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
    2. Fundacion ARAID, Zaragoza, Spain


## 1) software dependencies

+ Perl module 'GO::TermFinder', which in turns requires GraphViz. In Ubuntu run:

    $ sudo apt-get install graphviz
    $ sudo cpan -i GO::TermFinder 

+ _get_go_annotation, which reads a local copy of gene_ontology.1_2.obo, and can be compiled with g++ as follows:
    
    $ gcc -o _get_go_annotation get_go_annotation.cpp myGOparser.cpp -lstdc++

+ ncbi-blast+ blastp and blastx, expected to be in $PATH, otherwise edit variable $BLASTPATH accordingly on barleyGO.pl


## 2) data dependencies

This software requires several files which should be included in your downloaded bundle:

    gene_ontology.1_2.obo     from http://www.geneontology.org/GO.downloads.ontology.shtml
 
    _barleyGO.annot.go        derived from barley_HighConf_genes_MIPS_23Mar12_HumReadDesc.txt and
                              barley_HighConf_genes_MIPS_23Mar12_IPROScan_GOs.txt, originally from
                              ftp://ftpmips.helmholtz-muenchen.de/plants/barley/public_data/genes,
                              made with perl $ _format_barley_GO_annotations.pl > _barleyGO.annot.go
                              PMID:23075845 (PubMed identifier)

    _barleyGO.B2G.go          derived from the MIPS file above and 4513.annot from 
                              http://www.b2gfar.org/showspecies?species=4513, made with 
                              $ perl _format_barley_B2G_annotations.pl > _barleyGO.B2G.go 
                              PMID:21335611

    _grassesAthaGOdb.faa.p*   derived from barleyGO.B2G annotated sequences, grasses + A.thaliana in
                              http://bioinformatics.psb.ugent.be/plaza/ and Aegilops tauschii sequences
                              and annotations retrieved from http://plants.ensembl.org/info/data/ftp/index.html
                              and ftp://ftp.gramene.org/pub/gramene/CURRENT_RELEASE/data/ontology/go
                              made with $ perl _format_blast_GO_db.pl > _blast_GO_db.faa and then 
                              $ bin/makeblastdb -in _blast_GO_db.faa -dbtype prot
                              PMIDs:21335611,20040540,24163254,24217918

## 3) Examples of use:

    $ perl barleyGO.pl 

    [options]:
    
    -h    this message
    -t    GO term to analyze, or file with 1 terms per line   (optional, example: -t GO:0003924)
    -f    file with gene identifiers for GO enrichment test   (optional, choose either -f,-s or -t)
    -s    FASTA file with sequences to be assigned GO terms   (optional, choose either -f,-s or -t)
          and matched IBSC genes from annotations in /home/contrera/codigo/cvs/barleyGO/_grassesAthaGOdb.faa
    
    Options that affect BLAST searches conducted with -s:
    -n    FASTA file contains nucleotide sequences            (optional, default: protein)
    -C    minimum %sequence coverage of best hit              (optional, default: -C 50)
    -O    minimum %sequence identity to best hit              (optional, default: -I 70)
    -T    number of threads for BLAST searches                (optional, default: -T 2)
    
    Options that control GO enrichment tests conducted with -f:
    -F    False Discovery Rate (FDR) for enrichment test      (optional, default: -F 0.01)
    -p    Bonferroni-corrected p-value for enrichment test    (optional, overrides FDR, example -p 0.01)
    -S    show available GO annotation files and exit         (optional)
    -A    use indicated GO annotation file                    (optional, default: -A 1
    -a    show annotation of genes in -f file                 (optional)
    -v    increase verbosity                                  (optional)

$ perl barleyGO.pl -S

    1:  _barleyGO.B2G.go    B2G-FAR annotation for barley IBSC genes PMID:21335611
    2:  _barleyGO.annot.go  original GO annotation from IBSC PMID:23075845

$ perl barleyGO.pl -t GO:0003924

    # parsing GO flat file ... done (34858 records)
    >term: GO:0003924
    GTPase activity(GO:0003924), nucleoside-triphosphatase activity(GO:0017111), pyrophosphatase activity(GO:0016462), hydrolase activity, acting on acid anhydrides, in phosphorus-containing anhydrides(GO:0016818), hydrolase activity, acting on acid anhydrides(GO:0016817), hydrolase activity(GO:0016787), catalytic activity(GO:0003824), molecular_function(GO:0003674) |

$ perl barleyGO.pl -s test.fna -n

    # input sequences: 10 (-s test.fna)
    # input seems to be: nucleotide sequences
    # GO sequence library: /home/contrera/codigo/cvs/barleyGO/_grassesAthaGOdb.faa
    
    # running BLAST jobs...
    
    # output file with GO terms associated to input sequences: test.fna.go (5 sequences)
    
    # output file with IBSC genes matched in input sequences: test.fna.ISBC.list (4 sequences)


Note that output *test.fna.ISBC.list* can be further used to calculate GO enrichment, with a command such as:

$ perl barleyGO.pl -f test.fna.ISBC.list

    # parsing ontology file /home/contrera/codigo/github/barleyGO/gene_ontology.1_2.obo ...
    # reading annotations from /home/contrera/codigo/github/barleyGO/_barleyGO.B2G.go ...
    # parsing 8 gene identifiers of test.fna.ISBC.list ...
    # matched/recognized gene identifiers: 4
    # unrecognized gene identifiers: 4
    # associations for branch P (FDR<=0.01):
    # no significant associations found

    # associations for branch C (FDR<=0.01):
    # no significant associations found

    # associations for branch F (FDR<=0.01):
    # no significant associations found

