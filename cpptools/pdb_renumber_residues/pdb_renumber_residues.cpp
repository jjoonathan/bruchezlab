//
//  main.cpp
//  pdb_renumber_residues
//
//  Created by Jonathan deWerd on 5/3/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#include <iostream>
#include <utility>
#include <fstream>
#include <algorithm>
#include "PDB.h"
using namespace std;

int main(int argc, const char * argv[]) {
	const char* infile_name = NULL;
	const char* outfile_name = NULL;
	for (int i=1; i<argc; i++) {
		if (!infile_name) infile_name=argv[i];
		else if (!outfile_name) outfile_name=argv[i];
	}
	if (!infile_name) {
		cout<<"Usage: "<<argv[0]<<" PDB_file_to_renumber.pdb out.pdb\n";
		exit(0);
	}
	PDB::PDBFile pdb(infile_name);
	map<pair<string,char>, map<int,PDB::AtomRecord*>> SegChain_resid_to_atom;
	for (PDB::Record* r : pdb.records) {
		PDB::AtomRecord* ar = dynamic_cast<PDB::AtomRecord*>(r);
		if (ar==nullptr) continue;
		SegChain_resid_to_atom[make_pair(string(ar->segment),ar->chain)][ar->resid] = ar;
	}
	for (auto it : SegChain_resid_to_atom) {
		pair<string,char> SegChain = it.first;
		map<int,PDB::AtomRecord*>& seq = it.second;
		int idxmax = 0;
		for (auto it : seq) { idxmax = max(it.first,idxmax); }
		int newidx = 1;
		for (int i=1; i<=idxmax; i++) {
			auto it = seq.find(i);
			if (it!=seq.end()) {
				PDB::AtomRecord* ar = it->second;
				ar->resid = newidx++;
			}
		}
	}
	ofstream out(outfile_name);
	out<<pdb;
    return 0;
}

