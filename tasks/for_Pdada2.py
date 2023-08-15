import sys
from os.path import join, dirname

sys.path.insert(0, dirname(dirname(__file__)))

from config import soft_db_path,luigi,run_cmd, valid_path,fileparser
from tasks.basic_tasks import base_luigi_task
from tasks.for_preprocess import QC_aft_screened,QC_trimmomatic
from os.path import *
vsearch = soft_db_path.vsearch_pth
rdp_gold = soft_db_path.rdp_gold_pth
map_pl = soft_db_path.map_pl_pth

class run_Pdada2(base_luigi_task):
    mission = 'Pdada2'
    def requires(self):
        df = fileparser(self.tab)
        R1_dict = df.R1
        R2_dict = df.R2
        kwargs = self.get_kwargs()
        
        tasks = {}
        for sid in R1_dict.keys():
            if self.screen:
                tasks[sid] = QC_aft_screened(sampleid=sid,
                                PE1=R1_dict[sid],
                                PE2=R2_dict[sid],
                                **kwargs)
            else:
                tasks[sid] = QC_trimmomatic(sampleid=sid,
                                PE1=R1_dict[sid],
                                PE2=R2_dict[sid],
                                **kwargs)            
        return tasks
    
    
    def output(self):
        ofiles = list(map(luigi.LocalTarget,
                          [join(self.odir,
                                "%s_output" % self.mission,
                                "rep.fasta"),
                           join(self.odir,
                                "%s_output" % self.mission,
                                "ASV_table.tsv"),
                           join(self.odir,
                                "%s_output" % self.mission,
                                "ASV_stats.txt")]
                          ))
        return ofiles

    def run(self):
        valid_path(self.output()[0].path, check_ofile=1)
        dada2_script = realpath(join(dirname(__file__),'..','static','dada2.R'))
        data_dir = dirname(list(self.input().values())[0][0].path)
        odir = dirname(self.output()[0].path)
        metadata = self.tab
        cmd = f"""Rscript {dada2_script} {metadata} {data_dir} {odir}"""
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)




class format_Pdada2(base_luigi_task):
    mission = 'Pdada2'
    def requires(self):
        kwargs = self.get_kwargs()
        return run_Pdada2(**kwargs)
    
    def output(self):
        ofiles = list(map(luigi.LocalTarget,
                          [join(self.odir,
                                "%s_output" % self.mission,
                                "rep.fasta"),
                           join(self.odir,
                                "%s_output" % self.mission,
                                "ASV_table.tsv"),
                           join(self.odir,
                                "%s_output" % self.mission,
                                "ASV_stats.txt"),
                           join(self.odir,
                                "%s_output" % self.mission,
                                "parsed.ok")]
                          ))
        return ofiles

    def run(self):
        ######## unify the format of the output
        from Bio import SeqIO
        import pandas as pd
        rep_fasta = self.output()[0].path
        tab = self.output()[1].path
        
        seqs = list(SeqIO.parse(rep_fasta,'fasta'))
        df = pd.read_csv(tab,sep='\t')
        for r in seqs:
            r.id = 'ASV' + '{n:06}'.format(n=int(r.id.replace('sq','').replace('OTU','').split(';')[0]))
        seq2name = {str(r.seq):r.id for r in seqs}
        df.index = [seq2name[_] for _ in df['Sequence']]
        _ = df.pop('sum')
        _ = df.pop('Sequence')
        with open(rep_fasta,'w') as f1:
            SeqIO.write(seqs,f1,'fasta-2line')
        df.to_csv(tab,sep='\t',index=1,index_label=df.index.name)  
        status = self.output()[3].path
        with open(status,'w') as f1:
            f1.write('ok')
        