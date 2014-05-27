#!/usr/bin/env python
from glob import glob
import pymol, os, sys
from pymol import cmd, stored

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


pymol_sessions = glob("in/*.pse")
#pymol_sessions = pymol_sessions[0:1]
for pse in pymol_sessions:
    session_name = os.path.splitext(os.path.basename(pse))[0]
    of = open("intermediate/%s.2fast"%(session_name,),'w')
    of_seqonly = open("intermediate/%s.fasta"%(session_name,),'w')
    print pse, session_name
    cmd.load(pse)
    obj_names = cmd.get_names("objects")
    rna_obj_names = []
    rna_obj_names.extend(filter(lambda name: 'LSU' in name, obj_names))
    rna_obj_names.extend(filter(lambda name: 'SSU' in name, obj_names))
    rna_obj_names.extend(filter(lambda name: 'RRNA' in name, obj_names))
    rna_obj_names = set(rna_obj_names)
    for obj_name in rna_obj_names:
        of.write(">%s\n"%(obj_name,))
        of_seqonly.write(">%s\n"%(obj_name,))
        sel = obj_name+" and n. C5"
        atoms = cmd.get_model(sel).atom
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
