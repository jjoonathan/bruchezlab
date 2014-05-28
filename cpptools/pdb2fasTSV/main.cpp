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
#define PATH_LEN 2048
using namespace std;

int main(int argc, char * argv[]) {
	// Parse Arguments
	const char* infile_name = nullptr;
	char out_base[PATH_LEN]; out_base[0]=0;
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
		if (!strcmp(argv[i],"-d")) {
			assert(i!=argc-1 /* -d option needs a directory */ );
			strcpy(out_base, argv[i+1]);
			size_t l = strlen(out_base);
			if (out_base[l]!='/') {
				out_base[l] = '/';
				out_base[l+1] = '\0';
			}
		}
		if (!infile_name) infile_name=argv[i];
	}
	if (!infile_name) {
		cout<<"Usage: "<<argv[0]<<" [-useseg] [-f field1,field2,...] PDB_file_to_dump_seqs_from.pdb\n";
		cout<<"    Outputs each chain as a file <root>.<chain>.fasTSV\n";
		cout<<"    -d: output directory. Default: same as input PDB.\n";
		cout<<"    -useseg: split output files into <seg>-<chain>. Default: just <chain>.\n";
		cout<<"    -f fields. Default: resid,resName,beta\n";
		size_t nfields = PDB::AtomRecord::fields_and_descs.size()/2;
		for (int i=0; i<nfields; i++) {
			const string& field_name = PDB::AtomRecord::fields_and_descs[2*i];
			const string& field_desc = PDB::AtomRecord::fields_and_descs[2*i+1];
			cout<<"      "<<setw(8)<<field_name<<"  "<<field_desc<<endl;
		}
		exit(0);
	}
	if (out_base[0]=='\0') {
		// No explicit output dir from -d => we infer it from the .pdb
		strncpy(out_base,infile_name,PATH_LEN);
		char* slash = strrchr(out_base,'/');
		if (slash) {
			slash[1] = '\0';
		} else {
			strcpy(out_base,"./"); //pdb didn't have a dirname => use cwd
		}
	}
	// At this point, out_base is "/path/to/out/dir/"
	// We now append the input file name so that for input.pdb
	// out_base is "/path/to/out/dir/input."
	char* end_of_outbase = out_base+strlen(out_base);
	char infile_basename[PATH_LEN];
	char* infile_slash = strrchr(infile_name,'/');
	if (infile_slash) {
		strcpy(infile_basename,infile_slash+1);
	} else {
		strcpy(infile_basename, infile_name);
	}
	char* infile_dot = strrchr(infile_basename,'.');
	if (infile_dot) infile_dot[0] = '\0';
	strcpy(end_of_outbase, infile_basename);
	// At this point, out_base should look like "/path/to/out/dir/input"
	
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
		char fname[PATH_LEN];
		strncpy(fname, out_base, PATH_LEN);
		char* ext = fname+strlen(fname);
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
	exit(0);
    return 0;
}

