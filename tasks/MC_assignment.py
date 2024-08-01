"""
Novel part for MC(population) delienation.

"""

import luigi
from tasks.basic_tasks import base_luigi_task
from os.path import dirname,join
from config import soft_db_path
from toolkit import run_cmd,valid_path
import os
from os.path import *

from bin.ncbi_convertor import NCBI_convertor
import pandas as pd
from Bio import SeqIO


def anno_repotu(rep_otu,db,name):
    ofile = f"{dirname(rep_otu)}/db_based_filter/{name}"
    if not exists(dirname(ofile)):
        os.makedirs(dirname(ofile))
    cmd = f"blastx -query {rep_otu} -db {db} -num_threads 10 -max_target_seqs 1000000 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle qframe' -evalue 1e-3 -max_hsps 1000000 -out {ofile} "
    if not exists(ofile):
        return cmd
    else:
        return
    
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
    all_ASV = list(set(sorted_df.index))
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
    return 




class get_MC_assignment(base_luigi_task):
    analysis_type = luigi.Parameter(default="otu")
    def requires(self):
        self.get_config()
        kwargs = self.get_kwargs()
        required_tasks = {}
        antype = str(self.analysis_type).lower().split(',')
        if len(set(antype).intersection(set(["all", "otu", "deblur", "dada2",'usearch','qc','pdada2','pdeblur','pusearch'])))==0:
            raise Exception("analysis type must be one of the `all,otu,deblur,dada2,pdada2,pdeblur,pusearch`")
        
        if  "otu"  in antype or  "all" in antype:
            from tasks.for_otu import vsearch_otutable
            required_tasks["otu"] = vsearch_otutable(**kwargs)
        if "dada2" in antype or  "all" in antype:
            from tasks.for_dada2 import dada2_summarize
            required_tasks["dada2"] = dada2_summarize(**kwargs)
        if "deblur" in antype or  "all" in antype:
            from tasks.for_deblur import deblur_summarize
            required_tasks["deblur"] = deblur_summarize(**kwargs)
        if "usearch" in antype or  "all" in antype:
            from tasks.for_usearch import usearch_OTU_table,usearch_zOTU_table,usearch_denoise,usearch_OTU
            required_tasks["usearch_OTU"] = usearch_OTU_table(**kwargs)
            required_tasks["usearch_OTU_rep"] = usearch_OTU(**kwargs)
            required_tasks["usearch_zOTU"] = usearch_zOTU_table(**kwargs)
            required_tasks["usearch_zOTU_rep"] = usearch_denoise(**kwargs)
        if "pusearch" in antype or  "all" in antype:
            from tasks.for_Pusearch import Pusearch_zOTU_table
            required_tasks["run_Pusearch"] = Pusearch_zOTU_table(**kwargs)                 
        if "pdada2" in antype or  "all" in antype:
            from tasks.for_Pdada2 import format_Pdada2
            required_tasks["run_Pdada2"] = format_Pdada2(**kwargs)     
        if "pdeblur" in antype or  "all" in antype:
            from tasks.for_Pdeblur import custom_deblur_parse
            required_tasks["run_Pdeblur"] = custom_deblur_parse(**kwargs)               
        return required_tasks

    def output(self):
        ofiles = []
        for k, f in self.input().items():
            if type(f) == list:
                ofile = f[0].path
            else:
                ofile = f.path
            odir = dirname(ofile)
            
            ofiles.append(join(odir,
                               "Positive_rep.fasta"))
        ofiles = list(set(ofiles))
        valid_path(ofiles,check_ofile=1)
        return [luigi.LocalTarget(_)
                for _ in ofiles]


    def run(self):
        if self.dry_run:
            pass

        target_species = self.get_config_params('target_species')
        outgroup_species = self.get_config_params('outgroup_species')
        ref_genedb = self.get_config_params('ref_genedb')
        targets = open(target_species).read().strip().split('\n')
        outgroups = open(outgroup_species).read().strip().split('\n')
        for k, f in self.input().items():
            if k in ['usearch_OTU','usearch_zOTU']:continue
            if type(f) == list:
                ofile = f[0].path
                otu_table = f[1].path
            else:
                ofile = f.path
                otu_table = ''
            odir = dirname(ofile)
            c = anno_repotu(ofile,
                        db='/mnt/home-db/pub/protein_db/swissprot/swissprot',
                        name='negative_swissprot.tbl')
            if not exists(join(odir,'negative_swissprot.tbl')):

                os.system(c)
            c = anno_repotu(ofile,
                            db=ref_genedb,
                            name='positive_656G.tbl')
            if not exists(join(odir,'positive_656G.tbl')):
                os.system(c)
            get_positive(ofile,targets,outgroups)
            



