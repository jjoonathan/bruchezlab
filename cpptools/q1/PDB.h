//
//  PDB.h
//  q1
//
//  Created by Jonathan deWerd on 4/24/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#ifndef __q1__PDB__
#define __q1__PDB__

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include "Record.h"
#include "MiscRecord.h"
#include "AtomRecord.h"

namespace PDB {
	using namespace std;
	
	struct PDBFile {
		vector<Record*> records;
		PDBFile(string fname);
		PDBFile(istream& str);
		PDBFile(istream&& str);
		bool check_lens(bool fix_lens);
	private:
		void init(istream& str);
	};
	
}

std::ostream& operator<<(std::ostream&, PDB::PDBFile&);

#endif /* defined(__q1__PDB__) */
