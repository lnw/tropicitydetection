#include <vector>
#include <iostream>
#include <fstream>
#include <list>
#include <iterator>
#include <regex>

using namespace std;


int main() {
smatch match;
regex alkaa("<DataArray");
regex loppuu("</DataArray");

fstream vti ("jvec-minimal.vti");
string vtiline;
while (getline (vti, vtiline)) {
	if (regex_search(vtiline, match, alkaa)) {
		cout<<"********************ALKAAA\n";
		cout<<vtiline<<"\n";
	}
	else if (regex_search(vtiline, match, loppuu)) {
		cout<<"********************LOPPUU\n";
		cout<<vtiline<<"\n";
	}
	else {
		//cout << "NOMATCH:    " << vtiline;
		cout<<vtiline<<"\n";
	}
}

cout<<"ENDING \n";
vti.close();
return 0;
}
