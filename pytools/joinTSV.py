#!/usr/bin/env python
import pandas as pd
import numpy as np
import os,sys,re,tempfile,itertools,subprocess
from util import *

tsvs = []
output_fname = None
for arg in sys.argv[1:]:
    if output_fname==1.23:
        output_fname = arg
    elif arg=='-o':
        output_fname = 1.23
    else:
        tsvs.append(arg)

if len(tsvs)<2 or output_fname in (None,1.23):
    print "Usage: %s file1.fasTSV file2.fasTSV -o joined.fasTSV"%sys.argv[0]
    print "    The column with header 'resName' should contain a"
    print "    sequence. The sequences of each fasTSV will be aligned"
    print "    so that each row of the output corresponds to a position"
    print "    in the alignment and the union of all columns in the"
    print "    input files are present in the output under the naming convention"
    print "    <file1>-<col1>, <file1>-<col2>, etc."
    exit(1)

# Take all .tsv files passed in arguments and eliminate common prefix/suffixes
# to determine a unique identifier for each.
# tsvs = sys.argv[1:]
# min_path_len = min(len(filename) for filename in tsvs)
# common_prefix_len = 0
# for i in range(0,min_path_len):
#     firstchar = tsvs[0][i]
#     all_chars_in_col_eq = all(filename[i]==firstchar for filename in tsvs)
#     if all_chars_in_col_eq:
#         common_prefix_len = += 1
#     else:
#         break
# common_suffix_len = 0
# for i in range(0,min_path_len):
#     firstchar = tsvs[0][-i]
#     all_chars_in_col_eq = all(filename[-i]==firstchar for filename in tsvs)
#     if all_chars_in_col_eq:
#         common_suffix_len = += 1
#     else:
#         break
# if common_prefix_len==min_path_len:
#     raise ValueError("All TSV paths are equal!")
# if common_suffix_len>0:
#     uniq_ids = [path[common_prefix_len:-common_suffix_len] for path in tsvs]
# else:
#     uniq_ids = [path[common_prefix_len:] for path in tsvs]
# if len(set(uniq_ids))!=len(uniq_ids):
#     raise ValueError("Could not generate unique identifiers from TSV paths: "+str(uniq_ids))

# Take all .tsv files passed in as arguments and use the part before the first '.' as the identifier
# (i.e. the prefix that will be prepended to its columns to make them unique)
uniq_ids = []
for fname in tsvs:
    fname = os.path.basename(fname)
    dotpos = fname.find('.')
    uniq_ids.append(fname[0:dotpos] if dotpos!=-1 else fname[0:])

# Read in the TSV files as pandas dataframes, generate sequences
dfs = [pd.read_table(tsvname,dtype={'resid':np.int32}) for tsvname in tsvs]
seqs = [''.join(df['resName']) for df in dfs]
hdrs = [str(n) for n in range(len(seqs))]

# Perform the alignment using clustalw2
alnIN_hndl, alnIN_name = tempfile.mkstemp(suffix=".fasta")
alnOUT_hndl, alnOUT_name = tempfile.mkstemp(suffix=".fasta")
alnIN_hndl = os.fdopen(alnIN_hndl,'w')
write_fasta(alnIN_hndl, hdrs, seqs)
alnIN_hndl.close()
os.close(alnOUT_hndl)
flags = ['-outorder=input','-output=fasta','-matrix=pam','-type=DNA']
subprocess.call(['clustalw2','-infile='+alnIN_name,'-outfile='+alnOUT_name]+flags)
alnHDR,alnSEQ = read_fasta(alnOUT_name)
os.unlink(alnIN_name)
os.unlink(alnOUT_name)

# Add an alignment index column to each data frame
for hdr,seq in zip(alnHDR,alnSEQ):
    unaln_seq = seqs[int(hdr)]
    df = dfs[int(hdr)]
    uniq_id = uniq_ids[int(hdr)]
    aln_col = pd.Series(np.ones(len(df.index),dtype=np.int32)*(-1),index=df.index)
    df.columns = [uniq_id+'-'+col for col in df.columns]
    i = 0
    for aln_i in range(len(seq)):
        if seq[aln_i]=='-':
            pass
        elif seq[aln_i]==unaln_seq[i]:
            aln_col[i] = aln_i
            i += 1
        else:
            raise ValueError("Mismatched Sequences!")
    df['aln'] = aln_col
    df.set_index('aln',inplace=True)

# Put the results into one big DF, write it
bigdf = dfs[0].join(dfs[1:],how='outer')
intcols = filter(lambda x: 'resid' in x, bigdf.columns)
betacols = filter(lambda x: 'beta' in x, bigdf.columns)
bigdf[intcols] = bigdf[intcols].fillna(-1).astype(np.int32)
bigdf[betacols] = bigdf[betacols].fillna(0)
bigdf.to_csv(output_fname,sep='\t')