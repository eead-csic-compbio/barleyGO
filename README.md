# barleyGO

barleyGO is a Perl script to annotate barley sequences with Gene Ontology terms
inferred by homology.

Authors Bruno Contreras Moreira (1,2), Carlos P Cantalapiedra (1), MJ Garcia-Pereira (1) 

    1. Estacion Experimental de Aula Dei-CSIC, Zaragoza, Spain
    2. Fundacion ARAID, Zaragoza, Spain


## 1) software dependencies

+ Perl module 'GO::TermFinder', which in turns requires GraphViz. In Ubuntu run:

    $ sudo apt-get install graphviz
    $ sudo cpan -i GO::TermFinder 

+ _get_go_annotation, which reads a local copy of gene_ontology.1_2.obo, and can be compiled as follows:
    
    $ gcc -o _get_go_annotation get_go_annotation.cpp myGOparser.cpp -lstdc++ -static

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

$ perl barleyGO.pl -f test.list

    # parsing ontology file /home/contrera/codigo/cvs/barleyGO/gene_ontology.1_2.obo ...
    # reading annotations from /home/contrera/codigo/cvs/barleyGO/_barleyGO.B2G.go ...
    # parsing 212 gene identifiers of test.list ...
    # matched/recognized gene identifiers: 187
    # unrecognized gene identifiers: 25
    # associations for branch P (FDR<=0.01):
    GOID TERM CORRECTED_PVALUE UNCORR_PVALUE NUM_LIST_ANNOTS LIST_SIZE TOTAL_NUM_ANNOTS POPULATION_SIZEFDR_RATE EXPECTED_FALSE_POS ANNOTATED_GENES
    GO:0005975 carbohydrate metabolic process  1.48882927117949e-08  2.91927308074409e-11  36  187 1080  19808 0.00% 0.00  MLOC_61506.1, MLOC_18499.1,... 
    GO:0044262  cellular carbohydrate metabolic process 0.000793050470563112  1.55500092267277e-06  22  187 728 19808 0.00% 0.00  MLOC_18499.1, ...
    # associations for branch C (FDR<=0.01):
    GO:0005737 cytoplasm 0.000934996748518741  7.85711553377093e-06  113 187 8813  19808 0.00% 0.00  MLOC_5326.1, MLOC_8093.1, MLOC_66714.1, ...
    GO:0030312 external encapsulating structure  0.00295124696574557 2.48003946701309e-05  29  187 1344  19808 0.00% 0.00  MLOC_11692.1, MLOC_75285.1,...
    GO:0005618 cell wall 0.0124082751070936  0.000104271219387341  27  187 1309  19808 0.00% 0.00  MLOC_11692.1, MLOC_75285.1, MLOC_61506.1, ...
    GO:0048046 apoplast  0.0140914406818047  0.000118415467914325  17  187 637 19808 0.00% 0.00  MLOC_18499.1, MLOC_70866.1, MLOC_44256.2, ...
    # associations for branch F (FDR<=0.01):
    GO:0003824 catalytic activity  0.00246797600274457 1.54248500171536e-05  113 187 8921  19808 0.00% 0.00  MLOC_12033.1, MLOC_5326.1, ...


$ perl barleyGO.pl -f test.list -A 2

    # parsing ontology file /home/contrera/codigo/cvs/barleyGO/gene_ontology.1_2.obo ...
    # reading annotations from /home/contrera/codigo/cvs/barleyGO/_barleyGO.annot.go ...
    # parsing 212 gene identifiers of test.list ...
    # matched/recognized gene identifiers: 141
    # unrecognized gene identifiers: 71
    # associations for branch P (FDR<=0.01):
    GO:0005975 carbohydrate metabolic process  8.45021282146773e-06  9.60251456984969e-08  23  141 659 14450 0.00% 0.00  MLOC_61506.1, MLOC_18499.1,...
    # associations for branch C (FDR<=0.01):
    # no significant associations found
    # associations for branch F (FDR<=0.01):
    GO:0003824 catalytic activity  0.00238566080823399 2.11120425507433e-05  96  141 7335  14450 0.00% 0.00  MLOC_12033.1, MLOC_5326.1, MLOC_18499.1,...

$ perl barleyGO.pl -s test.fna -n

    # input sequences: 10 (-s test.fna)
    # input seems to be: nucleotide sequences
    # GO sequence library: /home/contrera/codigo/cvs/barleyGO/_grassesAthaGOdb.faa
    #
    # running BLAST jobs...
    #
    # output file with GO terms associated to input sequences: test.fna.go (3 sequences)
    #
    # output file with IBSC genes matched in input sequences: test.fna.ISBC.list (2 sequences)


Note that output *test.fna.ISBC.list* can be further used to calculate GO enrichment, with a command such as:

$ perl barleyGO.pl -f test.fna.ISBC.list


