
from Bio.Seq import Seq
import numpy as np
import itertools
from os.path import exists,dirname,join
from ete3 import Tree

def remove_3rd_foramplicons(record,return_seq=False,start=0):
    # start should be 0-coordinates
    # inplace remove the 3rd, number 2 is the 3rd codon
    ### amplicon sequence has three possible situations.
    # original: 0120120120
    # [1:] generates 120120120
    # [2:] generates 20120120
    seq = str(record.seq)
    if start %3 ==0:
        two_partitions = [seq[::3],seq[1::3]]
    elif start %3 ==1:
        two_partitions = [seq[2::3],seq[::3]]
    elif start %3 ==2:
        two_partitions = [seq[1::3],seq[2::3]]
    final_seq = ''.join([''.join(_) for _ in zip(*two_partitions)])
    if return_seq:
        return final_seq
    record.seq = Seq(final_seq)  # :-2 mean remove the lefted stop codon
    return record

def map_two_aln(r_short,r_long):
    """function for mapping the short sequence (usually amplicon sequence) to the long sequences (usually the full gene / alignments)
    * PLEASE make sure the short one can be fund in the long sequences

    It will return four index.
    1. index of the start at the long sequence (with gaps)
    2. index of the end at the long sequence (with gaps)
    3. index of the start at the long sequence (without gaps)
    4. index of the end at the long sequence (without gaps)
    """
    r_l = ''.join([_ for _ in str(r_long.seq).upper() if _!='-'])
    r_s = ''.join([_ for _ in str(r_short.seq).upper() if _!='-'])
    s = r_l.index(r_s)
    e = s+len(r_s)
    pos = [_ for _,bp in enumerate(r_long.seq) if bp !='-']
    new_s = [bp for _,bp in enumerate(r_long.seq) if bp !='-']
    pos_d = dict(zip(range(len(new_s)),pos))
    ns,ne = pos_d[s],pos_d[e]
    return ns,ne,s,e

def intervals_extract(iterable):
    """function for extracting the interval from a series of numbers
    for example:
    1,2,3,4,5,6,7, 10,11,12,13,14,15,19,25
    should output [(1,7),(10,15),(19),(25)]
    only ouput the continuous intervals
    """
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

def compare_array(larray,i1,i2):
    """
    Function for comparing two sequences found in a alignment file.
    larray: alignment in array format
    i1: index of the first sequence
    i2: index of the second sequence

    output:
    1. number of different bp
    2. total number of bp can be compare (after removing gaps)
    3. todo:
    4. todo:
    """
    array = larray[[i1,i2],:].T
    # remove all gaps
    parray = array[~(array=='-').all(1),:]
    compare = parray[np.isin(parray,list('actg')).all(1),:]
    s = np.where(np.isin(parray,list('actg')).all(1))[0]  # list of comparable bases
    interval_s = list(intervals_extract(list(s)))
    ug_inside = 0
    for _i1,_i2 in zip(interval_s,interval_s[1:]):
        num_gaps = _i2[0]-_i1[-1]
        ug_inside+=num_gaps
    total_num = compare.shape[0]
    unique_gap = parray.shape[0]-compare.shape[0]
    same_num = compare[compare[:,0]==compare[:,1],:].shape[0]
    diff_num = total_num-same_num
    return diff_num, total_num,unique_gap,ug_inside


import pandas as pd
from os.path import dirname,realpath,abspath
from Bio import SeqIO
import copy,os
from subprocess import check_call

def translate_amplicons(infasta,
                        dbfasta=None,
                        tblfile=None,
                        outdir=None):
    """
    first, blastx search infasta against dbfasta. No need to search all. Stop early.
    second, extract the frame and translate the infasta into protein sequences.

    """

    infasta = realpath(abspath(infasta))
    if tblfile is None and dbfasta is not None:
        ofile = f"{dirname(infasta)}/blastxSearch.tbl"
        cmd = f"blastx -query {infasta} -subject {dbfasta} -outfmt 6 |head -n100 > {ofile}"
        check_call(cmd,shell=1)
    else:
        ofile = tblfile
    if outdir is None:
        outdir = dirname(infasta)
    _df = pd.read_csv(ofile,sep='\t',header=None)
    _df.columns = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle qframe".split(' ')
    _df = _df.sort_values('evalue').groupby('qaccver').head(1)
    s2frame = dict(zip(_df['qaccver'],_df['qframe']))
    records = {_.id:_ for _ in SeqIO.parse(infasta,'fasta')}
    new_seqs = []
    for name,r in records.items():
        f = int(s2frame[name])
        seq = copy.deepcopy(r)
        seq.seq = seq[f-1:].translate().seq
        new_seqs.append(seq)
    n = infasta.split('/')[-1].rsplit('.',1)[0]
    with open(f"{outdir}/{n}_prot.faa",'w') as f1:
        SeqIO.write(new_seqs,f1,'fasta-2line')



def build_tree(aln,
               oprefix='',
               dtype='protein',
               excluded_gids = []):
    if not exists(dirname(oprefix)):
        os.makedirs(dirname(oprefix))
    if dtype=='protein':
        treepath = f"{oprefix}_prot.raw.newick"
        cmd = f"FastTreeMP -lg -gamma {aln} > {treepath} "
        os.system(cmd)
    elif dtype=='nucl':
        treepath = f"{oprefix}_nucl.raw.newick"
        cmd = f"FastTreeMP -gtr -gamma {aln} > {treepath} "
        os.system(cmd)

    gene_tre = Tree(treepath)
    gene_tre.resolve_polytomy()
    gene_tre.prune([_ 
                    for _ in gene_tre.get_leaf_names() 
                    if _ not in excluded_gids])
    for l in gene_tre.get_leaves():
        l.name = l.name.split('_')[0]
    c = 1
    for n in gene_tre.traverse():
        if not n.is_leaf():
            n.name = f"INode{c}"
            c+=1
    treepath = f"{oprefix}_nucl.renamed.newick"            
    gene_tre.write(outfile=treepath,format=3)
    