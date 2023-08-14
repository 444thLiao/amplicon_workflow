"""
The qiime2 implemented deblur or the original deblur can not contain underscore in sample_name.
Thus, we use the original deblur to run the deblur analysis.



refer to https://github.com/biocore/deblur

"""
# odir = '/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/results/amplicon_ATP5B/20230811'

# input_fq = f'{odir}/OTU_pipelines/preprocessed/merged.fastq'

import sys
from os.path import join, dirname, basename,exists
sys.path.insert(0, dirname(dirname(__file__)))

from config import *


vsearch = soft_db_path.vsearch_pth
rdp_gold = soft_db_path.rdp_gold_pth
map_pl = soft_db_path.map_pl_pth
trimmomatic_jar = soft_db_path.trimmomatic_jar
trimmomatic_dir = dirname(trimmomatic_jar)
fastqc_path = soft_db_path.fastqc_path
multiac_path = soft_db_path.multiqc_path

from tasks.for_preprocess import joined_reads

class custom_deblur(base_luigi_task):
    mission = 'deblur'
    def requires(self):
        df = fileparser(self.tab)
        R1_dict = df.R1
        R2_dict = df.R2
        kwargs = self.get_kwargs()
        tasks = {}
        for sid in R1_dict.keys():
            tasks[sid] = joined_reads(sampleid=sid,
                                      PE1=R1_dict[sid],
                                      PE2=R2_dict[sid],
                                      screen=False,
                                      **kwargs)
        return tasks

    def output(self):
        odir = join(str(self.odir),
                    "%s_output" % self.mission)
        ofile = join(odir,
                     'all.biom')
        valid_path(ofile,check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        joined_dir = dirname(list(self.input().items())[0].path)
        output_dir = dirname(self.output().path)
        trunclen = self.get_config_params(('deblur_custom_args','trunclen'))
        deblur_ref = self.get_config_params(('deblur_custom_args','refdb'))
        jobs_to_start  = self.get_config_params(('deblur_custom_args','jobs_to_start'))
        cmd = f"deblur workflow --seqs-fp {joined_dir} --output-dir {output_dir} --trim-length {trunclen} --pos-ref-fp {deblur_ref} -w -a 0 -O {jobs_to_start} --log-file {output_dir}/deblur.log"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)

                
import pandas as pd
class custom_deblur_parse(base_luigi_task):
    mission = 'deblur'
    def requires(self):
        kwargs = self.get_kwargs()
        return custom_deblur(**kwargs)

    def output(self):
        ofiles = list(map(luigi.LocalTarget,
                          [join(self.odir,
                                "%s_output" % self.mission,
                                "rep.fasta"),
                           join(self.odir,
                                "%s_output" % self.mission,
                                "ASV_table.tsv"),]
                          ))
        return ofiles

    def run(self):
        infile = self.input().path
        ofile = self.ouput()[1].path
        cmd = f"biom convert -i {infile} -o {ofile} --to-tsv"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        rep_fa = join(dirname(infile),'all.seqs.fa')
        tmp = list(SeqIO.parse(rep_fa,'fasta'))
        seqmap = {}
        for idx,record in enumerate(tmp):
            record.id = record.name = record.description = ''
            record.id = 'ASV' + '{n:06}'.format(n=idx,)
            seqmap[str(record.seq)] = record.id
        c = open(ofile)
        next(c)
        names = next(c).strip('#').strip().split('\t')
        d = pd.read_csv(ofile,sep='\t',index_col=0,names=names,comment='#')
        d.index = [seqmap[_] for _ in d.index]
        d.to_csv(ofile,sep='\t',index=1)
        with open(self.ouput()[0].path,'w') as f1:
            SeqIO.write(tmp,f1,'fasta-2line')

        