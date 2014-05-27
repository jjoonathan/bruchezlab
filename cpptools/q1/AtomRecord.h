//
//  Record.h
//  q1
//
//  Created by Jonathan deWerd on 4/24/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#ifndef __q1__AtomRecord__
#define __q1__AtomRecord__

#include <iostream>
#include <string>
#include <vector>
#include "Record.h"

namespace PDB {
	using namespace std;
	
	struct AtomRecord : Record {
		//mName=="ATOM  "
		int serial; //5 digit Atom serial number
		char atomName[5]; // Atom name
		char altLoc;  // Alternate location indicator
		char resName[4]; // Residue name
		char chain;    // Chain identifier
		int resid;      // Residue sequence number
		char iCode;      // Code for insertion of residues
		float x,y,z;     // Location in angstroms
		float occupancy; // Occupancy
		float tempFactor; // Beta factor
		char segment[5]; // Segment identifier
		char element[3];    // Element symbol, right-justified
		char charge[3];  // Charge at the atom
		string text();
		AtomRecord(string& name, string& content);
		static Record* newRecord(string& name, string& content);
		ostream& writeTo(ostream&);
		bool check_lens(bool fix_lens=true); // Makes sure all values fit in PDB-defined field widths
		
		// Interface to dynamically access properties
		const static vector<string> fields_and_descs;
		vector<string> fields(); // Lists fields that can dynamically be fetched
		string field(string fieldName); // Fetches field dynamically
	};
	
}

#endif /* defined(__q1__AtomRecord__) */
