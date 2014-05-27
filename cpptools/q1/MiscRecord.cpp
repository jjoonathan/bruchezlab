//
//  Record.cpp
//  q1
//
//  Created by Jonathan deWerd on 4/24/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#include "MiscRecord.h"
#include <map>
#include <cstdlib>
#include <sstream>
using namespace std;

PDB::MiscRecord::MiscRecord(string& name, string& content)
: PDB::Record(name), mContent(content) {}


PDB::Record* PDB::MiscRecord::newRecord(string& name, string& content) {
	return new PDB::MiscRecord(name,content);
}

ostream& PDB::MiscRecord::writeTo(ostream& out) {
	out<<mName<<mContent;
	return out;
}

bool PDB::MiscRecord::check_lens(bool fix_lens) {
	bool good = true;
	if (mContent.size()!=80-6) {
		good = false;
		if (fix_lens) mContent.resize(80-6,' ');
	}
	return good && Record::check_lens(fix_lens);
}
