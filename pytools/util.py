from glob import glob
import os, sys, re
import pandas as pd
import numpy as np
import os.path


def read_fasta(fname):
    """ Turns the file
            >a
            SEQA
            >b
            SEQB
        into [['a','b'],['SEQA','SEQB']] """
    in_file = open(fname, 'r')
    hdrs, seqs = [], []
    hdr = None
    seqlines = []
    while True:
        ln = in_file.readline()
        if (ln == '' or ln[0] == '>') and hdr:
            hdrs.append(hdr)
            seqs.append(''.join(seqlines))
            hdr, seqlines = None, []
        if ln == '':
            break
        ln = ln.strip()
        if ln[0] == '>':
            hdr = ln[1:]
        else:
            seqlines.append(ln)
    in_file.close()
    return hdrs, seqs


def write_fasta(fname, hdrs, seqs):
    """ Turns [['a','b'],['SEQA','SEQB']] into the file
            >a
            SEQA
            >b
            SEQB
        where the file location is determined by fname. """
    if type(fname)==str:
        of = open(fname, 'w')
    else:
        of = fname
    for h, s in zip(hdrs, seqs):
        of.write('>%s\n' % h)
        of.write('%s\n' % s)
    of.close()


def read_fasTSV(fname):
    """ Turns the file
           col1	col2
           1	a
           2	b
        into {'col1':['1','2'], 'col2':['a','b']}. """
    t = pd.read_table(fname, dtype={'num': np.int32, 'seq': np.char})
    return t


def fasta_2_fasTSV(fname):
    """ Reads the sequences in the fasta file fname and outputs each
    sequence to a separate file, rootname.seqname1.fasTSV, rootname.seqname2.fasTSV,
    etc. """
    root, ext = os.path.splitext(fname)
    hdrs, seqs = read_fasta(fname)
    for h, s, i in zip(hdrs, seqs, range(len(hdrs))):
        if '.' in h:
            sys.stderr.write('\tWARNING: "." in header "%s"-> ambiguity in fname\n' % h)
        of = open('%s.%s.fasTSV' % (root, h), 'w')
        of.write('num\tseq\n')
        for c, j in zip(s, range(len(s))):
            of.write('%i\t%s\n' % (j + 1, c))
        of.close()


def fasTSV_2_fasta(fasTSVs):
    """ Reads the files listed in the array fasTSVs and outputs
    the sequences to splitext(fasTSVs[0])[0]+'.fasta'. In other words, converts
    rootname.seq1name.fasTSV, ..., rootname.seqNname.fasTSV to rootname.fasta. """
    hdrs, seqs = [], []
    num_xtract_re = re.compile('(.+?)(\\.([^.]+))?\\.fasTSV')
    fn_parts = [num_xtract_re.match(fname).groups() for fname in fasTSVs]
    if None in fn_parts:
        raise ValueError('fasTSV file name must have the format <prefix>.<seq>.fasTSV!')
    if not all([fn[0] == fn_parts[0][0] for fn in fn_parts]):
        raise ValueError('All fasTSV files in a cohort must have the same prefix!')
    fasTSVparts = [(fn[2], fn[0], fn[1]) for fn in fn_parts]
    for seqname, root, dotseqname in fasTSVparts:
        ifile = open(root + dotseqname + '.fasTSV', 'r')
        hdr = ifile.readline().strip().split('\t')
        seq_col_idx = hdr.index('seq')
        entries = [l.strip().split('\t') for l in ifile.readlines()]
        if len(entries[-1]) != len(entries[-2]):
            entries = entries[:-1]
        entries = [e[seq_col_idx] for e in entries]
        hdrs.append(seqname)
        seqs.append(''.join(entries))
    write_fasta(fasTSVparts[0][1] + '.fasta', hdrs, seqs)
