import sys
from os.path import join, dirname

sys.path.insert(0, dirname(dirname(__file__)))

from config import soft_db_path,luigi,run_cmd, valid_path,fileparser
from tasks.basic_tasks import base_luigi_task, tabulate_seq
from static.q2_function import convert2otutab, convert2seq
from tasks.for_preprocess import multiqc
import pandas as pd
vsearch = soft_db_path.vsearch_pth
rdp_gold = soft_db_path.rdp_gold_pth
map_pl = soft_db_path.map_pl_pth


class run_dada2(base_luigi_task):
    mission = 'dada2'

    def requires(self):
        from tasks import import_data
        tasks = {}
        tasks['import'] = import_data(tab=self.tab,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           log_path=self.log_path)
        kwargs = self.get_kwargs()
        tasks["fastqc_after"] = multiqc(status='after',
                                         **kwargs)
        return tasks
    
    
    def output(self):
        ofiles = list(map(luigi.LocalTarget,
                          [join(self.odir,
                                "%s_output" % self.mission,
                                "profiling.qza"),
                           join(self.odir,
                                "%s_output" % self.mission,
                                "rep.qza"),
                           join(self.odir,
                                "%s_output" % self.mission,
                                "profiling_stats.qza")]
                          ))
        return ofiles

    def run(self):
        valid_path(self.output()[0].path, check_ofile=1)
        _infile = self.input()['fastqc_after'].path.replace('.html','_data/multiqc_fastqc.txt')
        _df = pd.read_csv(_infile,sep='\t',index_col=0)
        #_df.loc[:,'max seq len'] = [int(v.split('-')[1].strip()) for v in _df['avg_sequence_length'] ]
        r1_names = [_ for _ in _df.index if '_R1' in _]
        r2_names = [_ for _ in _df.index if '_R2' in _]
        
        r1_min_len = int(_df.loc[r1_names,'avg_sequence_length'].min())
        r2_min_len = int(_df.loc[r2_names,'avg_sequence_length'].min())
        
        extra_str = self.batch_get_params('dada2_args')
        extra_str += ' --p-trunc_len_f' + ' ' + str(r1_min_len)
        extra_str += ' --p-trunc_len_r' + ' ' + str(r2_min_len)
        # _d = config.dada2_args
        # _d['trunc_len_f'] = r1_min_len
        # _d['trunc_len_r'] = r2_min_len
        # extra_str = ''
        # for p, val in _d.items():
        #     p = p.replace('_', '-')
        #     if val is True:
        #         extra_str += ' --p-%s' % p
        #     elif val is not None and val is not False:
        #         extra_str += ' --p-%s %s ' % (p, val)
        
        cmd = """{qiime2_p} dada2 denoise-paired --i-demultiplexed-seqs {input_file} --o-representative-sequences {rep_seq} --o-table {profiling_tab} --o-denoising-stats {stats_file} --verbose""".format(
            qiime2_p=self.get_params('qiime2_p'),
            input_file=self.input()['import'].path,
            profiling_tab=self.output()[0].path,
            rep_seq=self.output()[1].path,
            stats_file=self.output()[2].path)

        cmd += extra_str
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


############################################################

class view_rep_seq(tabulate_seq):
    def requires(self):
        kwargs = dict(tab=self.tab,
                      odir=self.odir,
                      dry_run=self.dry_run,
                      log_path=self.log_path)
        return run_dada2(**kwargs)
    def output(self):
        return luigi.LocalTarget(self.input()[1].path.replace(".qza", ".qzv"))


class dada2_summarize(luigi.Task):
    tab = luigi.Parameter()
    odir = luigi.Parameter()
    dry_run = luigi.BoolParameter()
    log_path = luigi.Parameter(default=None)

    def requires(self):
        kwargs = dict(tab=self.tab,
                      odir=self.odir,
                      dry_run=self.dry_run,
                      log_path=self.log_path)
        required_tasks = {}
        required_tasks["dada2"] = run_dada2(**kwargs)
        required_tasks["view_rep_seq"] = view_rep_seq(**kwargs)
        return required_tasks

    def output(self):
        otutab, rep, stats = [_.path
                              for _ in self.input()["dada2"]]
        ofiles = [otutab.replace('.qza','.csv'),
                  rep.replace('.qza','.fa'),
                  stats.replace('.qza','.csv'),
                  ]
        return [luigi.LocalTarget(_) for _ in ofiles]

    def run(self):
        otutab, rep, stats = [_.path
                              for _ in self.input()["dada2"]]
        convert2otutab(otutab,
                       otutab.replace('.qza','.csv'))
        convert2seq(rep,
                    rep.replace('.qza','.fa'))
        convert2otutab(stats,
                       stats.replace('.qza','.csv'))
