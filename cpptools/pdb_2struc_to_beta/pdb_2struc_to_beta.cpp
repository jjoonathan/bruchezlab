//
//  main.cpp
//  2s2beta
//
//  Created by Jonathan deWerd on 4/28/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <set>
#include <string>
#include "PDB.h"
using namespace std;

set<string> noncanonical;
bool is_canonical(char* rnaview_call) {
	if (!strcmp(rnaview_call,"stacked")) return true;
	if (!strcmp(rnaview_call,"+/+")) return true;
	if (!strcmp(rnaview_call,"-/-")) return true;
	if (!strcmp(rnaview_call,"W/W")) return true;
	noncanonical.insert(string(rnaview_call));
	return false;
}

float bond_type_for_call(char* rnaview_call) {
	if (is_canonical(rnaview_call)) return 2.0;
	return 1.0;
}

int main(int argc, const char * argv[]) {
	ifstream smap("/tmp/cerR.txt");
	smap.ignore(100,'\n');
	smap.ignore(100,'\n');
	char tmp[100];
	char bond_call[100];
	string internalid_respair;
	map<pair<char,int>,float> bond_type; // pair(chain,resid) --> bond_type
	while (1) {
		if (!smap.good()) break;
		smap>>internalid_respair;
		if (internalid_respair=="END_base-pair") break;
		if (internalid_respair=="Summary") break;
		if (internalid_respair[internalid_respair.length()-1]!=',')
			break;
		char lChain; smap>>ws; smap.getline(tmp,100,':'); lChain = tmp[0];
		int lResid; smap>>ws>>tmp; lResid = atoi(tmp);
		smap>>bond_call;
		int rResid; smap>>tmp; rResid = atoi(tmp);
		char rChain; smap>>ws; smap.getline(tmp,100,':'); rChain = tmp[0];
		smap.ignore(100,'\n');
		bond_type[make_pair(lChain, lResid)] = bond_type_for_call(bond_call);
		bond_type[make_pair(rChain, rResid)] = bond_type_for_call(bond_call);
		cout<< lChain <<";"<< lResid << endl;
		cout<< rChain <<";"<< rResid << endl;
	}
	cout<<"----------------------\n";
	PDB::PDBFile pdb("/tmp/cerR.pdb");
	for (PDB::Record* r : pdb.records) {
		PDB::AtomRecord* ar = dynamic_cast<PDB::AtomRecord*>(r);
		if (ar==nullptr) continue;
		float beta = 0;
		auto bondcall_it = bond_type.find(make_pair(ar->chain,ar->resid));
		cout<< ar->chain <<":"<< ar->resid << endl;
		if (bondcall_it != bond_type.end()) {
			beta = bondcall_it->second;
		}
		ar->tempFactor = beta;
	}
	ofstream pdbo("/tmp/cerR_wc.pdb");
	pdbo<<pdb;
    return 0;
}

