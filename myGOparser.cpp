// Bruno Contreras-Moreira, 2011 EEAD/CSIC
#include <string.h>
#include <map>
#include "myGOparser.h"

using namespace std;

int parse_GO_flat_file(string &filename, map <string, GOnode> &GO)
{
	// 1) parse flat GO file
	FILE * gofile = fopen(filename.c_str(),"r");
	if(gofile == NULL)
	{
		printf("# parse_GO_flat_file : cannot read %s, exit ...\n",
			filename.c_str());
		return 0;
	}
	
	int n_of_records = 0;
	char line[LINE], field[FIELD];
	string term,name,alt_id,parent;
	map <string, string> synonym;
	while( fgets(line, LINE, gofile) != NULL )
	{
    	/*[Term]
		id: GO:0006355
		name: regulation of transcription, DNA-dependent
		namespace: biological_process
		alt_id: GO:0032583
		is_a: GO:0010468 ! regulation of gene expression */
		
		if(strncmp(line,"id:",3) == 0)
		{ 
			sscanf(line,"id: %s", (char *) &field);
			term = field; 
			n_of_records++;
		}
		else if(strncmp(line,"name:",5) == 0) 
		{
			strncpy(field,line+6,FIELD-1);
			field[strcspn (field,"\n")] = '\0'; // chomp
			name = field; 
			GO[term].name = name;
		}
		else if(strncmp(line,"alt_id:",7) == 0)
		{
			sscanf(line,"alt_id: %s",(char *) &field);
			alt_id = field;
			synonym[alt_id] = term; 
		}
		else if(strncmp(line,"is_a:",4) == 0)
		{
			sscanf(line,"is_a: %s",(char *) &field);
			parent = field;	
			GO[term].parents.push_back(parent);
		}
	}
	fclose(gofile);
	
	// 2) link synonims
	map <string, string>::iterator syn;
	for(syn = synonym.begin(); syn != synonym.end(); ++syn ) 
		GO[syn->first] = GO[syn->second];
	
	return n_of_records;
}

string get_full_annotation_for_term(string &term, map <string, GOnode> &GO)
{
	string annot, pterm;
	vector <string>::iterator pit;
	
	if(GO.find(term) == GO.end())
	{
		annot = "[cannot find GO term ";
		annot += term;
		annot += "]";
	}
	else
	{
		if(GO[term].parents.size() == 0)
		{
			annot += GO[term].name;
			annot += "(";
			annot += term;
			annot += ") | ";
		}
		else
		{
			for(pit=GO[term].parents.begin();pit != GO[term].parents.end();++pit)
			{ 	
				annot += GO[term].name;
				annot += "(";
				annot += term;
				annot += "), ";
				annot += get_full_annotation_for_term(*pit,GO);
			
			}
		}	
	}
	
	return annot;
}
