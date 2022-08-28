"""
Perform phylogenetic placement

Taken the niH analysis as example

"""
import copy
import os
import re
from collections import defaultdict
from glob import glob
from os.path import *
from subprocess import check_call
from ete3 import Tree
from Bio import SeqIO
from Bio.Data import IUPACData
from Bio.Seq import Seq
from tqdm import tqdm
import pandas as pd
import json
import numpy as np
tmp_d = IUPACData.ambiguous_dna_values

## following data can be found in HK server

os.chdir('/home-user/lling/5_Amplicon/PolF_Soil_seq_data/')

def renamed_tree(tre):
    count = 0
    for n in tre.traverse():
        if not n.name:
            if not str(int(n.support)) or str(int(n.support)) == '1':
                n.name = 'I%s' % (count)
            else:
                n.name = 'I%s_S%s' % (count, str(int(n.support)))
            count += 1
    return tre
## format the reference tree  (it could be generated with complete reference gene instead of trimmed sequences.)
ref_newick = "/home-user/lling/5_Amplicon/built_nifH_db/merged_nifH.newick"
tre = Tree(ref_newick)
tre.resolve_polytomy()
tre = renamed_tree(tre)
tre.write(outfile='./ref_nifH_resolved.newick',format=3)
ref_newick = '/home-user/lling/5_Amplicon/PolF_Soil_seq_data/ref_nifH_resolved.newick'

## get the branch lengths from the reference tree, the percentile95 will be taken as a threshold
n2plength = {n.name:n.dist for n in tre.traverse()}
all_dist = list(n2plength.values())
np.mean(all_dist),np.max(all_dist)
percentile95 = np.percentile(all_dist,95)

## trim the alignment
ref_aln = "/home-user/lling/5_Amplicon/built_nifH_db/merged_nifH.aln"
pre_sequences = list(SeqIO.parse(ref_aln,'fasta'))
primer_set = {"f": "TGCGAYCCSAARGCBGACTC",
              "r": "ATSGCCATCATYTCRCCGGA"}

def get_region(seq, f, r,include_primer=False):
    f_len, r_len = len(f), len(r)
    f = ''.join([l if len(tmp_d.get(l, l)) ==
                 1 else f"[{tmp_d.get(l,l)}]" for l in f])
    r = str(Seq(r).reverse_complement())
    f = ''.join([l if len(tmp_d.get(l, l)) ==
                 1 else f"[{tmp_d.get(l,l)}]" for l in f])
    r = ''.join([l if len(tmp_d.get(l, l)) ==
                 1 else f"[{tmp_d.get(l,l)}]" for l in r])
    if include_primer:
        matched_seq = re.search(f'({f}[ATCG]*{r})', seq)
    else:
        matched_seq = re.search(f'{f}([ATCG]*){r}', seq)
    if matched_seq is None:
        return None
    # return the sequence removing the primer
    return matched_seq.groups()[0]

def trim_aln_with_seq(ref_seq, regional_seq):
    # 0-coordinated, thus left closed and right closed
    aft2ori_idx = {}
    ori_i = aft_i = 0
    for s in ref_seq:
        if s.upper() not in 'ACGTN':
            ori_i += 1
        else:
            aft_i += 1
            ori_i += 1
        aft2ori_idx[aft_i] = ori_i
    ori_idx = list(range(0, len(ref_seq)))
    nogap_refseq = ''.join([_.upper()
                            for _ in ref_seq if _.upper() in 'ACGTN'])
    aft_idx = list(range(0, len(nogap_refseq)))

    idx = nogap_refseq.find(regional_seq.upper())
    if idx == -1:
        raise IOError
    start = aft2ori_idx[idx]
    end = aft2ori_idx[idx+len(regional_seq)]
    return start, end

matched_ = defaultdict(int)
for r in tqdm(pre_sequences):
    pos = [_ for _,bp in enumerate(r.seq) if bp !='-']
    new_s = [bp for _,bp in enumerate(r.seq) if bp !='-']
    pos_d = dict(zip(range(len(new_s)),pos))
    region = get_region(''.join(new_s),
                        primer_set['f'],
                        primer_set['r'],
                        include_primer=True)
    if not region:
        continue
    s, e = trim_aln_with_seq(''.join(new_s), region)
    ns = pos_d[s]
    ne = pos_d[e]
    matched_[(ns,ne)] +=1

### manual check whether it has only one start & end. Otherwise it might have some problems
assert len(matched_)==1
s,e = list(matched_.keys())[0]
print(s,e)

### align the reference sequences using the trimmed sequences
trimmed_seqs = []
for r in pre_sequences:
    sub_r = copy.deepcopy(r)[s:e]
    sub_r.seq = Seq(''.join([_ for _ in sub_r.seq if _!='-']))
    trimmed_seqs.append(sub_r)
with open('./ref_trimmed_nifH.fasta','w') as f1:
    SeqIO.write(trimmed_seqs,f1,'fasta-2line')
os.system(f"mafft --auto ./ref_trimmed_nifH.fasta > ./ref_trimmed_nifH.aln")
ref_phy = '/home-user/lling/5_Amplicon/PolF_Soil_seq_data/final_reference.phy'
records = list(SeqIO.parse('ref_trimmed_nifH.aln','fasta'))
n1,n2 = len(records),len(records[0].seq)
with open(ref_phy,'w') as f1:
    f1.write(f"{n1} {n2}\n")
    for r in records:
        f1.write(f"{r.id}{' '*10}{r.seq}\n")

## validations
# os.system(f"cat ./ref_trimmed_nifH.fasta ../PolF_Soil_seq_data/Brady_nifH.fasta > ./ALL_together.fasta; mafft --auto ./ALL_together.fasta > ./ALL_together.aln")
# os.system(f"FastTreeMP ./ALL_together.aln > ./ALL_together.newick")
# os.system(f"mkdir trees; iqtree -nt 15 -m MFP -redo -mset GTR,HKY85,K80 -s ./ALL_together.aln -pre ./trees/ALL_together -fast ")
########
os.chdir('/home-user/lling/5_Amplicon/PolF_Soil_seq_data/assignments')
infa = '/home-user/lling/5_Amplicon/PolF_Soil_seq_data/output_dir/dada2_output/rep.fa'
if exists('./papara_log.aln'):
    os.system(f"rm ./papara_*.aln")
cmd1 = f"papara -t {ref_newick} -s {ref_phy} -q {infa} -r -n aln -j 30"
check_call(cmd1,shell=1)
cmd2 = f"epa-ng --split {ref_phy} papara_alignment.aln --redo"
check_call(cmd2,shell=1)
cmd3 = f"epa-ng --ref-msa ./reference.fasta --tree {ref_newick} --query ./query.fasta -T 30 -w ./  --model 'GTR+F+R10' --redo"
check_call(cmd3,shell=1)
# ./reference.fasta and query.fasta are generated

from for_software.for_EPA.parse_jplace import parse_tree_with_edges
guppy_exe = "/home-user/thliao/download/pplacer-Linux-v1.1.alpha19/guppy "
jplace = './epa_result.jplace'
cmd = f"{guppy_exe} to_csv {jplace} > {dirname(jplace)}/jplace.csv"
check_call(cmd,shell=1)
cmd = f"guppy edpl {jplace} --csv -o {dirname(jplace)}/edpl.csv"
check_call(cmd,shell=1)
edpl_df = pd.read_csv(f"{dirname(jplace)}/edpl.csv",index_col=0,header=None)
df = pd.read_csv(f"{dirname(jplace)}/jplace.csv")
df = df.sort_values('like_weight_ratio',ascending=False).groupby('name').head(1)
df.index = df['name']
df = df.reindex(edpl_df.index)
df.loc[:,'EDPL'] =  edpl_df[1]
df = df.reindex(columns=['edge_num', 'like_weight_ratio','distal_length', 'pendant_length','EDPL'])
obj = json.load(open(jplace))
used_tree = obj["tree"]
tree,node_name2edge_num = parse_tree_with_edges(used_tree)
e2n = {v:k for k,v in node_name2edge_num.items()}
df.loc[:,'node'] = [e2n[_] for _ in df['edge_num']]
df.loc[:,'parental length'] =  [n2plength[_] for _ in df['node']]
final_df = df.loc[(df['pendant_length']<=percentile95) & (df['EDPL']<=0.1),:]
validated_seqs = list(final_df.index)

### classification method for assigning each ASV to different groups according to manual check.s
def get_num(tree, nodes):
    t = tree.get_common_ancestor(nodes)
    return [_.name for _ in t.traverse()]
def classification_criteria(tree):
    name2node = {_.name:_ for _ in tree.traverse()}
    all_names = list(name2node)

    FL1_group = []
    for target_n in ['I164_S0']:
        FL1_group+=list([_.name for _ in name2node[target_n].traverse()])
    _a = list([_.name for _ in name2node['I211_S0'].traverse()])
    FL1_group = set(FL1_group).difference(set(_a))

    pb_group = list(set(all_names).difference([_.name for _ in name2node['I132_S0'].traverse()]))

    FL2_group = []
    for target_n in ['I161_S0']:
        FL2_group+=list([_.name for _ in name2node[target_n].traverse()])

    Sym_group =  list([_.name for _ in name2node['I154_S0'].traverse()])
    Sym_group = set(Sym_group).difference(set(FL1_group))

    criteria = {"pb":pb_group,
                "FL1":FL1_group,
                "FL2":FL2_group,
                'Sym':Sym_group}
    return criteria

criteria = classification_criteria(tre)
seq2type = {}
for type_name, nodes in criteria.items():
    sub_df = final_df.loc[final_df['node'].isin(nodes)]
    names = list(sub_df.index)
    seq2type.update({n: type_name for n in names})

count_tab = f"/home-user/lling/5_Amplicon/PolF_Soil_seq_data/output_dir/dada2_output/profiling.csv"
count_df = pd.read_csv(count_tab, sep='\t', index_col=0)
count_df = count_df.T
count_df = count_df.loc[validated_seqs,:]
sample2total_reads = count_df.sum(0)
ratio_df = count_df/sample2total_reads * 100
sub_ratio_df = ratio_df.loc[ratio_df.index.isin(seq2type), :]
sub_ratio_df.loc[:, 'type'] = [seq2type[s] for s in sub_ratio_df.index]
type2ratio = sub_ratio_df.groupby('type').sum().T
type2ratio.columns = [f"{_} (%)" for _ in type2ratio.columns]
type2ratio.loc[:, 'reads inferred as nifH'] = count_df.sum(0).reindex(type2ratio.index)

sample2stats = pd.read_csv(dirname(count_tab)+'/profiling_stats.csv',sep='\t',index_col=0)
sample2stats = sample2stats.reindex(type2ratio.index)
sample2stats = sample2stats.reindex(columns=['input','filtered','denoised','merged','non-chimeric'])
sample2stats.columns=['input','After QC','dada2-denoised', 'dada2-merged','dada2-ASV']
type2ratio = type2ratio.join(sample2stats)
type2ratio.to_csv('./nifH_abundances.tsv',sep='\t',index=1)


