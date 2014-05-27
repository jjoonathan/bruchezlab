#!/usr/bin/env python
from glob import glob
import os, sys, re
import os.path
from util import *


fn_fasTSV = {}  # Map root to list of numbered TSV names
fn_fasta = []
fn_unk = []
fasTSV_re = re.compile('(.+?)(\\.([^.]+))?\\.fasTSV')
fasta_re = re.compile('(.+?)\\.fasta')
for fname in sys.argv[1:]:
	m = fasTSV_re.match(fname)
	if m:
		arr = fn_fasTSV.get(m.group(1),[])
		arr.append(fname)
		fn_fasTSV[m.group(1)] = arr
		continue
	m = fasta_re.match(fname)
	if m:
		fn_fasta.append(fname)
		continue
	fn_unk.append(fname)


sys.stderr.write('Converting fasTSV to fasta:\n')
for k,v in fn_fasTSV.iteritems():
	for fn in v:
		sys.stderr.write('\t'+fn+'\n')
		fasTSV_2_fasta(v)
sys.stderr.write('Converting fasta to fasTSV:\n')
for fn in fn_fasta:
	sys.stderr.write('\t'+fn+'\n')
	fasta_2_fasTSV(fn)
sys.stderr.write('Ignoring because extension was not .fasta or .fasTSV:\n')
for fn in fn_unk:
	sys.stderr.write('\t'+fn+'\n')
