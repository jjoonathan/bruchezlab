#!/usr/bin/env python
from glob import glob
import pymol, os, sys, re
import pandas as pd
import numpy as np
from math import *
from pymol import cmd, stored

# Parse args
fasTSV_fnames, pdb_fname, objname = [], None, None
bcol, qcol, rncol, ridcol = "beta", None, "resName", "resid"
for arg in sys.argv[1:]:
    if arg=="-o":
        objname = 1.23
    elif objname==1.23:
        objname = arg
    elif arg=="-b":
        bcol = 1.23
    elif bcol==1.23:
        bcol=arg
    elif arg=="-q":
        qcol=1.23
    elif qcol==1.23:
        qcol=arg
    elif arg=="-n":
        rncol=1.23
    elif rncol==1.23:
        rncol=arg
    elif arg=="-i":
        ridcol=1.23
    elif ridcol==1.23:
        ridcol=arg
    elif arg[-7:] == ".fasTSV":
        fasTSV_fnames.append(arg)
    elif arg[-4:] in (".pdb", ".pse"):
        pdb_fname = arg
    else:
        raise ValueError("Can't figure out what to do with arg "+arg)
fasTSVbad = not fasTSV_fnames or any([not os.path.exists(p) for p in fasTSV_fnames])
pdbbad = not pdb_fname or not os.path.exists(pdb_fname)
if 1.23 in (bcol,qcol,rncol,ridcol) or (not bcol and not qcol):
    raise ValueError("Incomplete -b or -q argument")
if fasTSVbad or pdbbad:
    print (bcol,qcol,rncol,ridcol)
    print "fasTSVbad: "+str(fasTSVbad)
    print "pdbbad: "+str(pdbbad)
    print "Usage: %s data.fasTSV data2.fasTSV structure.{pdb,pse}" % sys.argv[0]
    print "   Action: copies beta values to the structure"
    print "   Output: structure_.pdb"
    print "   Option: -o objname (only applies changes to object o in a PSE file)"
    print "   Option: -b name of the column from which to set b-factors (default: 'beta')"
    print "   Option: -q name of the column from which to set q-factors (default: none)"
    print "   Option: -q name of the column with seq letters (default: 'resName')"
    print "   Option: -i name of the column with residue IDs (default: 'resid')"
    print "   Requires pymol be built+installed from source (module pymol)."
    exit(1)
file_chain = None
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()
cmd.load(pdb_fname)

for tsv_name in fasTSV_fnames:
    tsv = pd.read_table(tsv_name)
    stored.resid2b = {}
    stored.resid2q = {}
    stored.resid2n = {}
    stored.mismatch = 0
    stored.assigned = 0
    for num,row in tsv.iterrows():
        rid = row[ridcol]
        if bcol:
            try:
                val = log10(row[bcol])
            except:
                val = -1
            if not np.isfinite(val):
                val = -1
            stored.resid2b[rid] = max(val,-1)
        if qcol:
            val = row[qcol]
            if not np.isfinite(val):
                val = 0
            stored.resid2q[rid] = val
        stored.resid2n[rid] = row[rncol]
    m = re.match(r'(.+)\.(.)\.fasTSV',tsv_name)
    if m:
        file_chain = m.group(2)
    selectors = []
    if objname not in (None, 1.23):
        selectors.append(objname)
    if file_chain:
        selectors.append("chain " + file_chain)
    sel = ' and '.join(selectors)
    print "Selector: "+sel
    mm = "; stored.mismatch += int(stored.resid2n.get(int(resi),resn)!=resn)"
    asgnd = "; stored.assigned += 1"
    sample = "; stored.b=b; stored.resi=resi; stored.resn=resn"
    setb = "b=stored.resid2b.get(int(resi),b)"
    setq = "; q=stored.resid2q.get(int(resi),q)"
    cmd.alter(sel,setb+setq+mm+asgnd+sample)
    #mismatches = 0
    #assigned = 0
    #for obj_name in cmd.get_names("objects"):
    #    for atom in cmd.get_model(obj_name).atom:
    #        if file_chain and atom.chain!=file_chain:
    #            continue
    #        if objname and obj_name!=objname:
    #            continue
    #        resid = int(atom.resi)
    #        pdb_resn = atom.resn
    #        tsv_resn = stored.resid2n.get(resid,None)
    #        new_b = stored.resid2b.get(resid,None)
    #        if new_b != None:
    #            assigned += 1
    #            if pdb_resn!=tsv_resn:
    #                mismatches += 1
    #            atom.b = 999.0;#str(new_b)
    #    # would use load_model here to get this to work
    print "DONE, mismatches=%i, assigned=%i"%(stored.mismatch,stored.assigned)

cmd.save(os.path.splitext(pdb_fname)[0]+"_"+os.path.splitext(pdb_fname)[1])

