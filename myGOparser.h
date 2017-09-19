/* Bruno Contreras-Moreira, 2011 EEAD/CSIC */  

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>

using namespace std;

#define LINE 400
#define FIELD 200

struct GOnode 
{
	string name;
	vector <string> parents;	
};

int parse_GO_flat_file(string &filename, map <string, GOnode> &GO);

string get_full_annotation_for_term(string &term, map <string, GOnode> &GO);
