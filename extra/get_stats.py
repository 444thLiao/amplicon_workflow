

from os.path import *
from glob import glob
import os
from Bio import SeqIO
from collections import defaultdict
def get_num(infile,format='fasta'):
    if format=='fastq':
        num = int(os.popen(f"zgrep -c '^+$' {infile}").read().strip().split(' ')[-1])
    elif format=='fasta':
        num = int(os.popen(f"zgrep -c '^>' {infile}").read().strip().split(' ')[-1])
    elif format=='merged-fasta':
        num = defaultdict(int)
        for r in SeqIO.parse(infile,'fasta'):
            num[r.id.split(';')[0].rsplit('_',1)[0]] += 1
    elif format=='merged-fastq':
        num = defaultdict(int)
        for r in SeqIO.parse(infile,'fastq'):
            num[r.id.split(';')[0].rsplit('_',1)[0]] += 1
    return num
import pandas as pd
def get_usearch_stats(idir):
    info_df = pd.DataFrame()
    clean_d = f"{idir}/OTU_pipelines/preprocessed/after_QC/"
    for f in glob(clean_d+'/*_R1.clean.fq.gz'):
        sample_name = f.split('/')[-1].split('_R1')[0]
        n = get_num(f,'fastq')
        info_df.loc[sample_name,'clean data/reads'] = n
    joined_d = f"{idir}/OTU_pipelines/preprocessed/joined_reads/"
    for f in glob(joined_d+'/*.fq.gz'):
        sample_name = f.split('/')[-1].split('.fq.gz')[0]
        n = get_num(f,'fastq')
        info_df.loc[sample_name,'joined/reads'] = n
    # num = get_num(f'{idir}/OTU_pipelines/preprocessed/merged.fastq','merged-fastq')
    # for s,v in num.items():
    #     info_df.loc[s,'merged/reads'] = v
    num = get_num(f"{idir}/USEARCH/filtered.fa",'merged-fasta')
    for s,v in num.items():
        info_df.loc[s,'filtered/reads'] = v
    otu = f'{idir}/USEARCH/zotu/zotutab_raw.txt'
    num = pd.read_csv(otu,sep='\t',index_col=0).sum(0).to_dict()
    for s,v in num.items():
        info_df.loc[s,'Final/reads'] = v
        
    return info_df

for idir in glob('/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/results/*'):
    if os.path.isdir(idir):
        info_df = get_usearch_stats(idir)
        sn = idir.split('/')[-1]
        info_df.to_csv(f'./usearch_stats_{sn}.tsv',sep='\t',index=1)