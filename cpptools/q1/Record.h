//
//  Record.h
//  q1
//
//  Created by Jonathan deWerd on 4/24/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#ifndef __q1__Record__
#define __q1__Record__

#include <iostream>
#include <string>

namespace PDB {
	using namespace std;

	struct Record {
		string mName;
		virtual ostream& writeTo(ostream&)=0;
		Record(string& name);
		string text();
		static Record* fromText(string& text);
		bool check_lens(bool fix_lens=true); // Makes sure all values fit in PDB-defined field widths
	};
	
	//Copy 1-based-indexed characters firstchar THROUGH lastchar from
	//src (which is assumed to have had 6 characters chopped off the front)
	//to dst and append a terminating \0 to dst.
	//   Example: extract the "Residue Name" field, defined in PDB 3.3 as:
	//     COLUMNS    DATA TYPE    FIELD       DEFINITION
	//      18-20       string      resName     Residue Name
	//   like so:
	//      char temp[4];
	//      copy_field(temp, stuff_after_record_name, 18, 20);
	void copy_field(char* dst, string src, int firstchar, int lastchar);
}

std::ostream& operator<<(std::ostream&, PDB::Record*);
std::ostream& operator<<(std::ostream&, PDB::Record&);

#endif /* defined(__q1__Record__) */
