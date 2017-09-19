#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include "myGOparser.h"

// http://bioinfoperl.blogspot.com.es/2011/08/decodificando-la-gene-ontology-cc.html

// gcc -o _get_go_annotation get_go_annotation.cpp myGOparser.cpp -lstdc++ -static

using namespace std;

#define DEFFLATFILE "gene_ontology.1_2.obo";

int main(int argc, char *argv[])
{
	string flatfilename = DEFFLATFILE;
	map <string, GOnode> parsedGO;
	string input_term, GOfilename, line, annot; 
	int option_char,n_of_records = 0;

	// parse arguments
	string usage = "# usage: \n\n";
	usage += " -h this message\n";
	usage += " -t GO term, such as -t GO:0046983\n";
	usage += " -f file with GO terms, 1 per line\n\n";
	
	if(argc == 1){	fprintf(stderr,"%s",usage.c_str());	exit(-1);	}	
	
	while ( (option_char=getopt(argc,argv, "t:f:")) != -1)
		switch (option_char)
		{  
         case 't': input_term = optarg; break;
			case 'f': GOfilename = optarg; break;
			
         default: fprintf(stderr,"%s",usage.c_str());	exit(-1);
      }
	
	printf("# parsing GO flat file ... ");
	n_of_records = parse_GO_flat_file(flatfilename, parsedGO);
	printf("done (%d records)\n\n",n_of_records);
	
	if(GOfilename.empty() == false)
	{
		std::ifstream GOstream(GOfilename.c_str());
		if(GOstream.is_open())
		{
			while(getline(GOstream, line))
			{
    			std::istringstream iss(line);
				if(!(iss >> input_term)) break; 
				annot = get_full_annotation_for_term(input_term,parsedGO);
				printf(">term: %s\n%s\n\n",input_term.c_str(),annot.c_str());
			}
		}
		else
		{
			fprintf(stderr,"# need a valid input file with GO terms, exit\n");
			exit(-2);
		}	
	}
	else
	{
		annot = get_full_annotation_for_term(input_term,parsedGO);
		printf(">term: %s\n%s\n\n",input_term.c_str(),annot.c_str());
	}	
	
	exit(0);
}
