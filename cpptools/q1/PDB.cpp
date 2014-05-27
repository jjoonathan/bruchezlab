//
//  PDB.cpp
//  q1
//
//  Created by Jonathan deWerd on 4/24/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#include "PDB.h"
#include <stdexcept>
#include <map>
#include <fstream>
#include <string>
using namespace std;

PDB::PDBFile::PDBFile(string fname) {
	ifstream f(fname);
	if (!f.is_open()) {
		throw logic_error("Couldn't open PDB file for reading.");
	}
	this->init(f);
}

PDB::PDBFile::PDBFile(istream& str) {
	this->init(str);
}

PDB::PDBFile::PDBFile(istream&& str) {
	this->init(str);
}

void PDB::PDBFile::init(istream& stream) {
	string l;
	while(true) {
		getline(stream,l);
		l.resize(80, ' ');
		if (!stream.good()) break;
		PDB::Record* rec = Record::fromText(l);
		if (rec) records.push_back(rec);
	}
}

bool PDB::PDBFile::check_lens(bool fix_lens) {
	bool good = true;
	for (PDB::Record* r : records) {
		good = good && r->check_lens(fix_lens);
	}
	return good;
}

ostream& operator<<(ostream& o, PDB::PDBFile& f) {
	for (auto r : f.records) {
		o<<r<<"\n";
	}
	return o;
}