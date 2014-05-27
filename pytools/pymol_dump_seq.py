#!/usr/bin/env python
from glob import glob
import pymol, os, sys, re
from pymol import cmd, stored

if len(sys.argv)<2:
	print "Usage: %s in.pdb [out.fasta] [out.2fast]"
	print "   Requires clustalw2 to be in the PATH."
	print "   Requires pymol be built+installed from source (module pymol)."
infile, ofn, of2n = None, None, None
object_sel_regex = None
for arg in sys.argv[1:]:
	if object_sel_regex==1.23:
		object_sel_regex = arg
	elif arg=="-obj":
		object_sel_regex = 1.23
	elif infile==None:
		infile = arg
	elif ofn==None:
		ofn = arg
	elif of2n==None:
		of2n = arg

infile = argv[1]
ofn = argv[2] if (len(sys.argv)>=3) else (os.path.splitext(infile)[0] + ".fasta")
of2n = argv[3] if (len(sys.argv)>=4) else (os.path.splitext(infile)[0] + ".2fast")



pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

color_name_to_idx = dict(cmd.get_color_indices())
c2r_arr = """
deeppurple	5
purpleblue	-4
blue	-2
skyblue	-0.9
cyan	-0.7
green	-0.5
palegreen	-0.3
white	0
paleyellow	0.3
yellow	0.5
brightorange	0.7
orange	0.9
deepsalmon	2
tv_red	4
red	5
gray	NaN
""".split()
color_idx_to_reactivity = {}
for name, reactivity in zip(c2r_arr[0::2],c2r_arr[1::2]):
    idx = color_name_to_idx[name]
    color_idx_to_reactivity[idx] = float(reactivity)


of = open(ofn,'w')
of2 = open(of2n,'w')
cmd.load(infile)
for obj_name in cmd.get_names("objects"):
	if object_sel_regex not in (None, 1.23):
		if not re.match(object_sel_regex,obj_name):
			continue
	of.write(">%s\n"%(obj_name,))
	of2.write(">%s\n"%(obj_name,))
	atoms = cmd.get_model(obj_name).atom
	print dir(atoms)
	print dir(atoms[0])
	exit(0)
	colors = []  # List of indices
	stored.pse_extract_color_list = colors
	cmd.iterate(sel, "stored.pse_extract_color_list.append(color)")
	reactivities = map(lambda x: color_idx_to_reactivity.get(x,float("NaN")), colors)
	seq = ''.join(atom.resn for atom in atoms)
	of.write(seq+"\n")
	of_seqonly.write(seq+"\n")
	of.write(' '.join([str(r) for r in reactivities])+"\n")
cmd.delete("all")
of.close()
