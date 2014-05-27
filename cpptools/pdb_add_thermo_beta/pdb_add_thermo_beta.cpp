//
//  main.cpp
//  pdb_add_thermo_beta
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
		const int wlen = 13;
		int idxmax = 0;
		for (auto it : seq) { idxmax = max(it.first,idxmax); }
		for (int i=1; i<=idxmax; i++) {
			double gt = 0;
			PDB::AtomRecord* bar = seq[i];
			if (!bar) continue;
			for (int j=0; j<wlen; j++) {
				PDB::AtomRecord* ar = seq[i+j];
				if (ar) {
					const char* rn = ar->resName;
					if (!strcmp(rn,"a")) gt += 2;
					if (!strcmp(rn,"t")) gt += 2;
					if (!strcmp(rn,"A")) gt += 2;
					if (!strcmp(rn,"T")) gt += 2;
					if (!strcmp(rn,"u")) gt += 2;
					if (!strcmp(rn,"U")) gt += 2;
					if (!strcmp(rn,"g")) gt += 4;
					if (!strcmp(rn,"G")) gt += 4;
					if (!strcmp(rn,"c")) gt += 4;
					if (!strcmp(rn,"C")) gt += 4;
				}
			}
			bar->tempFactor = 123;
		}
	}
	ofstream out(outfile_name);
	out<<pdb;
    return 0;
}

