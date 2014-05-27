//
//  main.cpp
//  pdb_dump_seqs
//
//  Created by Jonathan deWerd on 5/2/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//
//  Takes a field list and PDB file as arguments

#include <iostream>
#include <cstdlib>
#include <utility>
#include <map>
#include <algorithm>
#include <iomanip>
#include <climits>
#include <cstring>
#include <cassert>
#include <fstream>
#include "PDB.h"
using namespace std;

int main(int argc, char * argv[]) {
	// Parse Arguments
	char* infile_name = NULL;
	bool useseg = 0;
	vector<string> fields = {"resid","resName","beta"};
	for (int i=1; i<argc; i++) {
		if (!strcmp(argv[i],"-useseg")) {
			useseg = 1;
			continue;
		}
		if (!strcmp(argv[i],"-f")) {
			assert(i!=argc-1 /* -f option needs an arg list */);
			fields.clear();
			char* tokstate;
			for (char* word = strtok_r(argv[i+1], ", ", &tokstate);
				 word;
				 word=strtok_r(nullptr, ", ", &tokstate)) {
				fields.push_back(word);
			}
		}
		if (!infile_name) infile_name=argv[i];
	}
	if (!infile_name) {
		cout<<"Usage: "<<argv[0]<<" [-useseg] [-f field1,field2,...] PDB_file_to_dump_seqs_from.pdb\n";
		cout<<"    Outputs each chain as a file <root>.<chain>.fasTSV\n";
		cout<<"    -useseg: split output files into <seg>-<chain> vs just <chain>\n";
		cout<<"    -f fields (default resid,resName,beta):\n";
		size_t nfields = PDB::AtomRecord::fields_and_descs.size()/2;
		for (int i=0; i<nfields; i++) {
			const string& field_name = PDB::AtomRecord::fields_and_descs[2*i];
			const string& field_desc = PDB::AtomRecord::fields_and_descs[2*i+1];
			cout<<"      "<<setw(8)<<field_name<<"  "<<field_desc<<endl;
		}
		exit(0);
	}
	
	// Read the PDB records into a (seg,chain)->(resid,AtomRecord*) map
	PDB::PDBFile pdb(infile_name);
	map<pair<string,char>, map<int,PDB::AtomRecord*>> SegChain_to_seq;
	for (PDB::Record* r : pdb.records) {
		PDB::AtomRecord* ar = dynamic_cast<PDB::AtomRecord*>(r);
		if (ar==nullptr) continue;
		map<int,PDB::AtomRecord*>& seq = SegChain_to_seq[make_pair(string(ar->segment),ar->chain)];
		seq[ar->resid] = ar;
	}
	
	// Loop through entries in the (seg,chain) map
	for (auto it : SegChain_to_seq) {
		// Open the file to which we are sending the data of this (seg,chain)
		const string& seg = it.first.first;
		char chain = it.first.second;
		char fname[2048];
		strcpy(fname, infile_name);
		size_t fname_len = strlen(fname);
		char* ext = fname+fname_len-4;
		if (strcasecmp(ext, ".pdb")) {
			throw invalid_argument("Input PDB file ought to have a .pdb extension");
		}
		if (useseg) {
			sprintf(ext,".%s-%c.fasTSV",seg.c_str(),chain);
		} else {
			sprintf(ext,".%c.fasTSV",chain);
		}
		ofstream of(fname);
		size_t num_fields = fields.size();
		for (int i=0; i<num_fields; i++) {
			of << fields[i] << ((i==num_fields-1)?"\n":"\t");
		}
		
		// Loop through the residues in this (seg,chain), output the records
		// to the file we just opened
		map<int,PDB::AtomRecord*>& seq = it.second;
		int resid_max=0, resid_min=INT_MAX;
		for (auto it : seq) {
			resid_max = max(it.first,resid_max);
			resid_min = min(it.first,resid_min);
		}
		for (int i=resid_min; i<=resid_max; i++) {
			auto it = seq.find(i);
			if (it==seq.end()) {
				// Residue id i has no corresponding AtomRecord.
				// Current behavior: ignore
			} else {
				for (int i=0; i<num_fields; i++) {
					of << it->second->field(fields[i]) << ((i==num_fields-1)?"\n":"\t");
				}
			}
		}
	}
    return 0;
}

