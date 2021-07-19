#include <iostream>
#include <fstream>
#include <string> 
#include <algorithm>
#include <ctype.h>


#include <R.h>
#include <Rdefines.h>


using namespace std;

///char* nucs= "ATCGN";

extern "C" {
	
	unsigned char GetCode(char prev_n, char next_n) {
		 unsigned char code=10;
		 prev_n= toupper(prev_n);
		 next_n= toupper(next_n);
		 
			switch (prev_n) 
				{
				case 'A':
					switch (next_n) {
					case 'A':
						code = 0;
						break;
					case 'T':
						code = 2;
						break;
					case 'C':
						code = 5;
						break;
					case 'G':
						code= 4;
						break;
					default:
						code = 10;
						break;
					}
					break;
				case 'T':
					switch (next_n) {
					case 'A':
						code = 1;
						break;
					case 'T':
						code = 0;
						break;
					case 'C':
						code = 6;
						break;
					case 'G':
						code= 3;
						break;
					default:
						code = 10;
						break;
					}
					break;
				case 'C':
					switch (next_n) {
					case 'A':
						code = 3;
						break;
					case 'T':
						code = 4;
						break;
					case 'C':
						code = 9;
						break;
					case 'G':
						code= 7;
						break;
					default:
						code = 10;
						break;
					}
					break;
				case 'G':
					switch (next_n) {
					case 'A':
						code = 6;
						break;
					case 'T':
						code = 5;
						break;
					case 'C':
						code = 8;
						break;
					case 'G':
						code= 9;
						break;
					default:
						code = 10;
						break;
					}
					break;
				default:
					code = 10;
					break;
				}
		
		return (unsigned char) (code+1);			
	}

	SEXP readDinuc(SEXP filename, SEXP maxsize) {
		
		int nprotect=0;
		
		PROTECT(filename = AS_CHARACTER(filename));nprotect++;
		PROTECT(maxsize = AS_INTEGER(maxsize));nprotect++;
		
		
		const char* seqfile;
		int nsize= INTEGER_VALUE(maxsize);
		seqfile = STRING_VALUE(filename);		
				
		SEXP result;
		SEXP dinuc;
		SEXP seqsize;
		SEXP resnames;
	
		PROTECT(dinuc = NEW_INTEGER(nsize)); nprotect++;
		PROTECT(seqsize = NEW_INTEGER(1)); nprotect++;
		PROTECT(result = NEW_LIST(2)); nprotect++;
		
		string tmpstr= seqfile;		
		int ext=tmpstr.find_last_of(".fa");
		if (ext >= tmpstr.size()) {
			cout << "Error reading:" << tmpstr << std::endl;
			UNPROTECT(nprotect);
			return result;  /* file not found */
		}		
		Rprintf("%s\n", seqfile);
		ifstream infile(seqfile);
		
		if (!infile) {
			cout << "Cannot open :" << seqfile << std::endl;
			UNPROTECT(nprotect);
			return result;  /* file not found */	
		}
		
				
		std::string line;
		
		char prev_char='N';	
		bool firstline= true;
		int  linenum = 0;
		unsigned char cc;
		
		cc=(unsigned char) 11;

		int pos = 0;
		INTEGER(dinuc)[pos++]= cc;
	#ifdef __DEBUG__	
		cout << (unsigned int) cc << " ";
	#endif 
		//outfile << cc;
		
		
		while (std::getline(infile, line)) {
			if (firstline) {
				cout << "Reading " << seqfile << "["<< line<<"]" << endl;
				firstline=false;
			} else  {
				for (int ii=0; ii<line.size(); ii++) {
					if (ii==0) {
						cc=GetCode(prev_char, line[ii]);
					} else {
						cc=GetCode(line[ii-1],line[ii]);					
					}
					prev_char = line[line.size()-1];  //The last character
	#ifdef __DEBUG__	
					cout << (unsigned int) cc << " ";
	#endif
					INTEGER(dinuc)[pos++]= cc;
				}
	#ifdef __DEBUG__
				cout << "[" << prev_char << "]" << endl;
	#endif
			}
			linenum++;
		}
		
		
		infile.close();
		
		INTEGER(seqsize)[0] = pos;
				
		SET_ELEMENT(result,0, dinuc);
		SET_ELEMENT(result,1, seqsize);		

					
		PROTECT(resnames= NEW_CHARACTER(2)); nprotect++;
		SET_STRING_ELT(resnames,0,mkChar("dinuc"));
		SET_STRING_ELT(resnames,1,mkChar("length"));
		SET_NAMES(result, resnames);
			
		UNPROTECT(nprotect);
		return result;
	}		

}  // extern "C"
