//
//  main.cpp
//  q1
//
//  Created by Jonathan deWerd on 4/23/14.
//  Copyright (c) 2014 a.b.c. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "PDB.h"

void usage_die(const char* argv[]) {
	std::cout<<"Usage: "<<argv[0]<<" -join|-split in.pdb out.pdb dat.txt\n";
	std::cout<<"    dat.txt  file into which information lost in join is dumped\n";
	exit(0);
}

int main(int argc, const char * argv[]) {
	using namespace std;
	if (argc<2) usage_die(argv);
	bool join = !strcmp(argv[1],"-join");
	bool split = !strcmp(argv[1],"-split");
	if (argc!=5 || !(join^split)) usage_die(argv);
	PDB::PDBFile pdb(argv[2]);
	
	if (join) {
		ofstream dat(argv[4]);
		dat<<"newRID\tresName\tatmNm\tsegment\tchain\tserial\toldRID\n";
		int oldResid=0, newResid=0;
		for (PDB::Record* r : pdb.records) {
			PDB::AtomRecord* ar = dynamic_cast<PDB::AtomRecord*>(r);
			if (ar==nullptr) continue;
			if (ar->resid != oldResid) {
				oldResid = ar->resid;
				newResid++;
			}
			dat                  << newResid
			               <<"\t"<< ar->resName
			               <<"\t"<< ar->atomName
			               <<"\t"<< ar->segment
			               <<"\t"<< ar->chain
			               <<"\t"<< ar->serial
			               <<"\t"<< ar->resid
			               <<"\n";
			strcpy(ar->segment,"");
			ar->chain = '1';
			ar->resid = newResid;
		}
		cerr<<"Lengths are good: "<<pdb.check_lens(false)<<"\n";
	}
	
	ofstream of(argv[3]);
	of<<pdb;
    return 0;
}

