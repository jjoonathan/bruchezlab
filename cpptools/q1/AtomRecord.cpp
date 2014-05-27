//
//  Record.cpp
//  q1
//
//  Created by Jonathan deWerd on 4/24/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#include "AtomRecord.h"
#include <map>
#include <cstdlib>
#include <sstream>
#include <cassert>
#include <cmath>
#include <cctype>
#include <exception>
using namespace std;

const vector<string> PDB::AtomRecord::fields_and_descs = {
    "serial", "5 Digit Atom Serial Number",
    "atomName", "Name of the atom type (e.g. CA for alpha carbon)",
    "resName", "Residue name (e.g. K for lysine)",
    "resid", "Residue sequence number (e.g. 1,2,3 for MGD in N-MGDAG...-C)",
    "x", "X location in angstroms.",
    "y", "Y location in angstroms.",
    "z", "Z location in angstroms.",
    "beta", "Temperature factor (or other data)",
    "seg", "Segment Name (stripped of whitespace)",
    "chain", "Chain Name",
    "ele", "Element",
    "chg", "Charge"
};

vector<string> PDB::AtomRecord::fields() {
    size_t n = fields_and_descs.size()/2;
    vector<string> ret(n);
    for (int i=0; i<n; i++) {
        ret[i] = fields_and_descs[i*2];
    }
    return ret;
}

string PDB::AtomRecord::field(string fieldName) {
    char buf[100];
    if (fieldName=="serial") {
        sprintf(buf,"%i",serial);
        return string(buf);
    }
    if (fieldName=="atomName") return string(atomName);
    if (fieldName=="resName") return string(resName);
    if (fieldName=="resid") {sprintf(buf,"%i",resid); return string(buf);}
    if (fieldName=="x") {sprintf(buf,"%f",x); return string(buf);}
    if (fieldName=="y") {sprintf(buf,"%f",y); return string(buf);}
    if (fieldName=="z") {sprintf(buf,"%f",z); return string(buf);}
    if (fieldName=="beta") {sprintf(buf,"%f",tempFactor); return string(buf);}
    if (fieldName=="seg") return string(segment);
    if (fieldName=="chain") {buf[0]=chain; buf[1]='\0'; return string(buf);}
    if (fieldName=="ele") return string(element);
    if (fieldName=="chg") return string(charge);
    throw invalid_argument(fieldName);
}

// Remove leading and trailing spaces from C-string s
void trim_str(char* s) {
    char* dst = s;
    while (isspace(*s)) s++;
    while (*s) *(dst++) = *(s++);
    *dst = ' ';
    while (isspace(*dst)) *(dst--)='\0';
}

PDB::AtomRecord::AtomRecord(string& name, string& content) : Record(name) {
	if (content.length()<(80-6)) {
        //throw invalid_argument("line too short");
        content.resize(80-6,' ');
    }
	char temp[10];
	copy_field(temp,content,7,11); serial=atoi(temp); //5 digit Atom serial number
	copy_field(atomName,content,13,16); // Atom name
	altLoc = content[17-6-1];  // Alternate location indicator
	copy_field(resName,content,18,20); // Residue name
	chain = content[22-6-1];    // Chain identifier
	copy_field(temp,content,23,26); resid=atoi(temp);  // Residue sequence number
	iCode = content[27-6-1];      // Code for insertion of residues
	copy_field(temp,content,31,38); x=atof(temp); //x (Angstroms)
	copy_field(temp,content,39,46); y=atof(temp); //x (Angstroms)
	copy_field(temp,content,47,54); z=atof(temp); //x (Angstroms)
	copy_field(temp,content,55,60); occupancy=atof(temp); // Occupancy
	copy_field(temp,content,61,66); tempFactor=atof(temp); // Beta factor
	copy_field(segment,content,73,76); // Segment name
	copy_field(element,content,77,78);    // Element symbol, right-justified
	copy_field(charge,content,79,80);  // Charge at the atom
    trim_str(segment);
    trim_str(element);
    trim_str(charge);
    trim_str(resName);
    trim_str(atomName);
    if (serial==99999)
        cerr<<"5 digit atom serial# wrap isn't implemented";
    if (resid==9999)
        cerr<<"4 digit resid wrap isn't implemented";
}

PDB::Record* PDB::AtomRecord::newRecord(string& name, string& content) {
	return new PDB::AtomRecord(name,content);
}


// Returns false if field was too long or was fixed
static bool check_or_fix(bool fix, int& field, int maxlen) {
	if (field<rint(pow(10,maxlen))) return true; // All good; field has length <= len
	if (fix) field %= int(rint(pow(10,maxlen)));
	return false;
}
// Returns false if field was too long or was fixed
static bool check_or_fix(bool fix, char* field, int max_len) {
	assert(max_len<15);
	size_t len = strnlen(field, 15);
	if (len<=max_len) return true; // All good; field has correct length
    strncpy(field,"                ",max_len);
    field[max_len] = '\0';
	return false;
}
bool PDB::AtomRecord::check_lens(bool fix) {
	bool good = true;
	if (!check_or_fix(fix, serial, 5))
		good = false;
	if (!check_or_fix(fix, atomName, 4))
		good = false;
	//char                       altLoc
	if (!check_or_fix(fix, resName, 3))
		good = false;
	//char						 chain
	if (!check_or_fix(fix, resid, 4))
		good = false;
	//char					     iCode
	//float                      x,y,z,occupancy,tempFactor
	if (!check_or_fix(fix, segment, 4))
		good = false;
	if (!check_or_fix(fix, element, 2))
		good = false;
	if (!check_or_fix(fix, charge, 2))
		good = false;
	return good;
}


ostream& PDB::AtomRecord::writeTo(ostream& out) {
    char tmp[100];
	if (!check_lens(true)) {
		cerr<<"Field exceeded bounds";
	}
	out<<right<<fixed; // Right align, fixed notation,
	assert(mName.size()==6); out<<mName;
    out.width(5); out<<serial<<" ";
	out.width(4); out<<right<<atomName;
	out.width(1); out<<altLoc;
	out.width(3); out<<left<<resName<<" ";
	out.width(1); out<<chain;
	out.width(4); out<<right<<resid;
	out.width(1); out<<iCode<<"   ";
	sprintf(tmp,"%8.3f",x); out<<tmp;
	sprintf(tmp,"%8.3f",y); out<<tmp;
	sprintf(tmp,"%8.3f",z); out<<tmp;
	sprintf(tmp,"%6.2f",occupancy); out<<tmp;
	sprintf(tmp,"%6.2f",tempFactor); out<<tmp;
    out<<"      ";
	out.width(4); out<<left<<segment;
	out.width(2); out<<right<<element;
	out.width(2); out<<charge;
	return out;
}
