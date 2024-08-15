import os,sys
from os.path import abspath,dirname,exists
from subprocess import check_call
from glob import glob

from bin.ncbi_convertor import NCBI_convertor
import pandas as pd
from Bio import SeqIO

def run_cmd(cmd, dry_run=False, log_file=None, **kwargs):
    outstream = None
    if type(log_file) == str:
        if os.path.isfile(log_file):
            if os.path.getsize(log_file) != 0:
                outstream = open(log_file, 'a')
        if outstream is None:
            valid_path(log_file,check_ofile=True)
            outstream = open(log_file, 'w')
    elif log_file is None:
        outstream = sys.stdout
    else:
        outstream = log_file

    executable = "/usr/bin/zsh"
    if not os.path.exists(executable):
        executable = "/bin/bash"
    print(cmd, file=outstream)
    outstream.flush()
    if not dry_run:
        check_call(cmd,
                   shell=True,
                   executable=executable,
                   stdout=outstream,
                   stderr=outstream,
                   **kwargs)
        outstream.flush()

def get_validate_path(pth):
    if not pth.startswith('/'):
        pth = './' + pth
    pth = abspath(pth)
    return pth

def valid_path(in_pth,
               check_size=False,
               check_dir=False,
               check_glob=False,
               check_odir=False,
               check_ofile=False):
    if type(in_pth) == str:
        in_pths = [in_pth]
    else:
        in_pths = in_pth[::]
    for in_pth in in_pths:
        in_pth = os.path.abspath(os.path.realpath(in_pth))
        if in_pth is None:
            continue
        if check_glob:
            query_list = glob(in_pth)
            if not query_list:
                raise Exception('Error because of input file pattern %s' % in_pth)
        if check_dir:
            if not os.path.isdir(in_pth):
                raise Exception("Error because %s doesn't exist" % in_pth)
        if check_size:
            if os.path.getsize(in_pth) <= 0:
                raise Exception("Error because %s does not contain content." % in_pth)
        if check_odir:
            if not os.path.isdir(in_pth):
                os.makedirs(in_pth, exist_ok=True)
        if check_ofile:
            odir_file = os.path.dirname(in_pth)
            if not os.path.isdir(odir_file):
                os.makedirs(odir_file, exist_ok=True)
    return True

def get_dir_path(path,num=1):
    path = os.path.abspath(os.path.realpath(path))
    for _ in range(num):
        path = os.path.dirname(path)
    if path == "/":
        raise Exception("reach root path....")
    return path



def anno_repotu(rep_otu,db,name):
    ofile = f"{dirname(rep_otu)}/db_based_filter/{name}"
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    cmd = f"blastx -query {rep_otu} -db {db} -num_threads 10 -max_target_seqs 1000000 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle qframe' -evalue 1e-3 -max_hsps 1000000 -out {ofile} "
    if not exists(ofile):
        return cmd
    else:
        return "ls " + ofile
    
def assign_columns(df):
    if len(df.columns) == 16:
        return 'saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname staxids stitle qframe'.split(' ')
    if len(df.columns) == 14:
        return 'saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle qframe'.split(' ')

def get_positive(rep_otu,targets,outgroups):
    convertor = NCBI_convertor([],"nuccore")  
    postive_ofile = f"{dirname(rep_otu)}/db_based_filter/positive_656G.tbl"
    negative_ofile = f"{dirname(rep_otu)}/db_based_filter/negative_swissprot.tbl"
    _d1 = pd.read_csv(postive_ofile,sep='\t',index_col=0,header=None)
    _d1.columns = assign_columns(_d1)
    _d2 = pd.read_csv(negative_ofile,sep='\t',index_col=0,header=None)
    _d2.columns = assign_columns(_d2)
    df = pd.concat([_d1,_d2],axis=0)
    df = df.sort_values('bitscore')

    df.loc[:,'name'] = df.index
    df = df.loc[(df['bitscore']>100) ,:]
    sorted_df = df.sort_values('bitscore',ascending=False)

    subdf = sorted_df.loc[~sorted_df['staxids'].isna(),
                            ['saccver','staxids']]
    subject_id2taxid = dict(zip(subdf['saccver'],subdf['staxids']))
    subject_id2taxid = {k:str(v).split(';')[0]
                        for k,v in subject_id2taxid.items()}
    id2taxon = convertor.get_taxons_from_tid(tids=subject_id2taxid)
    tax_df = pd.DataFrame.from_dict(id2taxon, orient='index')
    allowed_genomes = set(list(id2taxon)+list(outgroups)+list(targets))
    for _ in outgroups:
        tax_df.loc[_,'genus'] = 'Sister/outgroup of Ruegeria'
    for _ in targets:
        tax_df.loc[_,'genus'] = 'Ruegeria'
    sorted_df = sorted_df.loc[[True
                                if _.split('_')[0] in allowed_genomes else False
                                for _ in sorted_df['saccver'] ],:]
    top10_df = sorted_df.groupby('name').head(1)
    asv2others = {}
    asv2ruegeria = set()
    for asv,row in top10_df.iterrows():
        if row['saccver'].split('_')[0] in targets:
            asv2ruegeria.add(asv)
        else:
            v = row['saccver']
            if v not in tax_df.index:
                v = v.split('_')[0]
            asv2others[asv] = tax_df.loc[v,'genus']
    rawseqs = [_ for _ in SeqIO.parse(rep_otu,'fasta') ]
    seqs = [_ for _ in rawseqs if _.id in asv2ruegeria]
    with open(f'{dirname(rep_otu)}/Positive_rep.fasta','w') as f1:
        SeqIO.write(seqs,f1,'fasta-2line')
    with open(f'{dirname(rep_otu)}/Negative_matched.tsv','w') as f1:
        for k,v in asv2others.items():
            f1.write(f"{k}\t{v}\n")

