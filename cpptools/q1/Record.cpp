//
//  Record.cpp
//  q1
//
//  Created by Jonathan deWerd on 4/24/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#include "Record.h"
#include "AtomRecord.h"
#include "MiscRecord.h"
#include <map>
#include <cstdlib>
#include <sstream>
#include <cctype>
using namespace std;

void PDB::copy_field(char* dst, string src, int firstchar, int lastchar) {
	const char* cstr = src.c_str();
	strncpy(dst,cstr+firstchar-1-6, lastchar-firstchar+1);
	dst[lastchar-firstchar+1] = 0;
}

PDB::Record::Record(string& name)
: mName(name) {}

string PDB::Record::text() {
	ostringstream oss;
	this->writeTo(oss);
	return oss.str();
}

bool PDB::Record::check_lens(bool fix_lens) {
	bool good = true;
	if (mName.length() != 6) {
		good = false;
		if (fix_lens) mName.resize(6,' ');
	}
	return good;
}





// Initialize the map from PDB 6-letter type identifiers to
// init routines
static map<string,PDB::Record*(*)(string& name, string& content)> type2init;

static bool type2init_inited = false;
static void init_type2init() {
	if (type2init_inited) return;
	type2init["ATOM  "] = PDB::AtomRecord::newRecord;
	type2init_inited = true;
}

PDB::Record* PDB::Record::fromText(string& text) {
	init_type2init();
	if (text.size()<6) {
		return nullptr;
	}
	string name = text.substr(0,6);
	string content = text.substr(6,text.length()-6+1);
	auto type_factory_it = type2init.find(name);
	if (type_factory_it == type2init.end()) {
		return new PDB::MiscRecord(name, content);
	} else {
		return (*type_factory_it->second)(name, content);
	}
	return nullptr;
}







ostream& operator<<(ostream& o, PDB::Record* r) {
	r->writeTo(o);
	return o;
}

ostream& operator<<(ostream& o, PDB::Record& r) {
	r.writeTo(o);
	return o;
}
