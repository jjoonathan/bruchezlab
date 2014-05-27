//
//  Record.h
//  q1
//
//  Created by Jonathan deWerd on 4/24/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#ifndef __q1__MiscRecord__
#define __q1__MiscRecord__

#include <iostream>
#include <string>
#include "Record.h"

namespace PDB {
	using namespace std;
	
	struct MiscRecord : Record {
		string mContent;
		ostream& writeTo(ostream&);
		MiscRecord(string& name, string& content);
		static Record* newRecord(string& name, string& content);
		bool check_lens(bool fix_lens=true); // Makes sure all values fit in PDB-defined field widths
	};

}

#endif /* defined(__q1__MiscRecord__) */
