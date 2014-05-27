//
//  MultiAtomRecord.h
//  q1
//
//  Created by Jonathan deWerd on 5/2/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#ifndef __q1__MultiAtomRecord__
#define __q1__MultiAtomRecord__

#include <iostream>
#include "Record.h"

namespace PDB {
	using namespace std;
	
	struct MultiAtomRecord : Record {
		string mContent;
		ostream& writeTo(ostream&);
		MiscRecord(string& name, string& content);
		static Record* newRecord(string& name, string& content);
		bool check_lens(bool fix_lens=true); // Makes sure all values fit in PDB-defined field widths
	};
	
}

#endif /* defined(__q1__MultiAtomRecord__) */
