//
//  main.cpp
//  pdb_dump_seqs
//
//  Created by Jonathan deWerd on 5/2/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <utility>
#include <map>
#include <algorithm>
#include "PDB.h"
using namespace std;

int main(int argc, const char * argv[]) {
	const char* infile_name = NULL;
	for (int i=1; i<argc; i++) {
		if (!infile_name) infile_name=argv[i];
	}
	if (!infile_name) {
		cout<<"Usage: "<<argv[0]<<" PDB_file_to_dump_seqs_from.pdb\n";
		exit(0);
	}
	PDB::PDBFile pdb(infile_name);
	map<pair<string,char>, map<int,string>> SegChain_to_seq;
	for (PDB::Record* r : pdb.records) {
		PDB::AtomRecord* ar = dynamic_cast<PDB::AtomRecord*>(r);
		if (ar==nullptr) continue;
		SegChain_to_seq[make_pair(string(ar->segment),ar->chain)][ar->resid] = ar->resName;
	}
	for (auto it : SegChain_to_seq) {
		pair<string,char> SegChain = it.first;
		map<int,string>& seq = it.second;
		int idxmax = 0;
		for (auto it : seq) { idxmax = max(it.first,idxmax); }
		char* linear_seq = (char*)malloc(idxmax+1);
		linear_seq[idxmax] = 0;
		for (int i=0; i<idxmax; i++) {
			auto it = seq.find(i);
			if (it==seq.end())
				linear_seq[i] = '-';
			else
				linear_seq[i] = it->second[0];
		}
		cout<<">Seg"<<SegChain.first<<"Chain"<<SegChain.second<<endl;
		cout<<linear_seq<<endl;
		free(linear_seq);
	}
    return 0;
}

