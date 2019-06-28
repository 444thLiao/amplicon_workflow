import sys
from os.path import join, dirname

sys.path.insert(0, dirname(dirname(__file__)))
import luigi

from config import soft_db_path
from tasks.basic_tasks import base_luigi_task, tabulate_seq
from toolkit import run_cmd, valid_path
import config.default_params as config
from static.q2_function import convert2otutab, convert2seq

vsearch = soft_db_path.vsearch_pth
rdp_gold = soft_db_path.rdp_gold_pth
map_pl = soft_db_path.map_pl_pth


class run_dada2(base_luigi_task):
    mission = 'dada2'

    def requires(self):
        from tasks import import_data
        return import_data(tab=self.tab,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           log_path=self.log_path)

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
        extra_str = ''
        for p, val in config.dada2_args.items():
            p = p.replace('_', '-')
            if val is True:
                extra_str += ' --p-%s' % p
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p, val)

        cmd = """{qiime2_p} dada2 denoise-paired --i-demultiplexed-seqs {input_file} --o-representative-sequences {rep_seq} --o-table {profiling_tab} --o-denoising-stats {stats_file}""".format(
            qiime2_p=config.qiime2_p,
            input_file=self.input().path,
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
