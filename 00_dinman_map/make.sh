#!/bin/bash
set -x
mkdir -p out || exit

# Dump LSU sequences to out/3u5h.{5,7,8}.fasTSV
pdb2fasTSV in/3u5h.pdb -d out || exit

# Dump SSU sequences to out/3u5f.
pdb2fasTSV in/3u5f.pdb -d out || exit

# Dump the Dinman data for Empty80S to out/Dinman{25S,18S,5.8S}.fasTSV
cat>/tmp/dump_dinman.py <<HERE_DOC
import pandas as pd
import numpy as np
names = ['25S rRNA', '18S rRNA', '5.8S rRNA']
ofnames = ['out/Dinman%s.fasTSV'%(n.split()[0]) for n in names]
tbls = [pd.read_excel('in/ALL_COMPLEXES_reactivity_values.xlsx',s) for s in names]
tbls = [tbl[['seqnum','seq','Empty 80S']].dropna() for tbl in tbls]
for tbl,ofname in zip(tbls,ofnames):
	tbl['Empty 80S'] = np.log10(tbl['Empty 80S'])
	tbl['Empty 80S'] = np.convolve(tbl['Empty 80S'], np.ones(4), 'same')
	tbl.columns = ['resid','resName','Empty 80S']
	tbl.to_csv(ofname,sep='\t')
HERE_DOC
python /tmp/dump_dinman.py || exit

# Merge the PDB sequences with the Dinman data
joinTSV out/Dinman18S.fasTSV  out/3u5f.6.fasTSV -o out/18S.fasTSV || exit
joinTSV out/Dinman25S.fasTSV  out/3u5h.5.fasTSV -o out/25S.fasTSV || exit
joinTSV out/Dinman5.8S.fasTSV out/3u5h.8.fasTSV -o out/5.8S.fasTSV || exit

# Create new PDB files with the altered beta factors
fasTSV2pymol -n 3u5h.5-resName -i 3u5h.5-resid -b "Dinman25S-Empty 80S" out/25S.fasTSV in/3u5h.pdb -o out/3u5h_.pdb || exit
fasTSV2pymol -n 3u5h.8-resName -i 3u5h.8-resid -b "Dinman5.8S-Empty 80S" out/5.8S.fasTSV out/3u5h_.pdb -o out/3u5h__.pdb || exit
fasTSV2pymol -n 3u5f.6-resName -i 3u5f.6-resid -b "Dinman18S-Empty 80S" out/18S.fasTSV in/3u5f.pdb -o out/3u5f_.pdb || exit

# Combine the 3 RNA strands into one pymol session
cat>/tmp/coalesce_pymol.py <<HERE_DOC
import pymol, os, sys, re
from pymol import cmd, stored
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()
cmd.load('out/3u5h__.pdb')
cmd.load('out/3u5f_.pdb')
cmd.remove('hetatm')
cmd.save('out/annotated_rsome.pse')
HERE_DOC
python /tmp/coalesce_pymol.py || exit
