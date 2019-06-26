import sys
from os.path import join, dirname

sys.path.insert(0, dirname(dirname(__file__)))
import luigi

from config import soft_db_path
from tasks.basic_tasks import base_luigi_task,tabulate_seq,visulize_seq
from toolkit import run_cmd, valid_path
import config.default_params as config
from luigi_workflow.static.q2_function import convert2otutab, convert2seq,convert2table


vsearch = soft_db_path.vsearch_pth
rdp_gold = soft_db_path.rdp_gold_pth
map_pl = soft_db_path.map_pl_pth


class joined_fastq(base_luigi_task):
    def requires(self):
        from tasks import import_data
        return import_data(tab=self.tab,
                           odir=self.odir,
                           dry_run=self.dry_run,
                           log_path=self.log_path, )

    def output(self):
        odir = join(str(self.odir),
                    "q2_pipelines")
        ofile = join(odir,
                     'joined_seq.qza')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        extra_str = ''
        for p, val in config.join_params.items():
            if val is True:
                extra_str += ' --p-%s' % p
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p, val)
        cmd = "{qiime2_p} vsearch join-pairs --i-demultiplexed-seqs {input_file} --o-joined-sequences {ofile}".format(
            qiime2_p=config.qiime2_p,
            input_file=self.input().path,
            ofile=self.output().path) + extra_str
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path,
                    dry_run=False)


class qualityFilter(base_luigi_task):
    def requires(self):
        return joined_fastq(tab=self.tab,
                            odir=self.odir,
                            dry_run=self.dry_run,
                            log_path=self.log_path)

    def output(self):
        return luigi.LocalTarget(self.input().path.replace("joined_seq",
                                                           "joined_qc_seq"))

    def run(self):
        extra_str = ''
        for p, val in config.qc_joined_params.items():
            p = p.replace('_', '-')
            if val is True:
                extra_str += ' --p-%s' % p
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p, val)
        cmd = "{qiime2_p} quality-filter q-score-joined --i-demux {input_qza} --o-filtered-sequences {output_seq} --o-filter-stats {output_stats}".format(
            qiime2_p=config.qiime2_p,
            input_qza=self.input().path,
            output_seq=self.output().path,
            output_stats=self.output().path.replace(
                '.qza', '-stats.qza'))
        # path of stats used in `share_tasks.summaried_tasks`
        cmd += extra_str
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False, )


class run_deblur(base_luigi_task):
    mission = "deblur"

    def requires(self):
        return qualityFilter(tab=self.tab,
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
        for p, val in config.deblur_args.items():
            p = p.replace('_', '-')
            if val is True:
                extra_str += ' --p-%s' % p
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p, val)

        cmd = """{qiime2_p} deblur denoise-16S --i-demultiplexed-seqs {input_file} --o-representative-sequences {rep_seq} --o-table {profiling_tab} --o-stats {stats_file}""".format(
            qiime2_p=config.qiime2_p,
            input_file=self.input().path,
            rep_seq=self.output()[1].path,
            profiling_tab=self.output()[0].path,
            stats_file=self.output()[2].path, )
        cmd += extra_str
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)




############################################################
class joined_summarize(visulize_seq):
    def requires(self):
        return joined_fastq(tab=self.tab,
                            odir=self.odir,
                            dry_run=self.dry_run,
                            log_path=self.log_path)

class view_rep_seq(tabulate_seq):
    def requires(self):
        kwargs = dict(tab=self.tab,
                      odir=self.odir,
                      dry_run=self.dry_run,
                      log_path=self.log_path)
        return run_deblur(**kwargs)


class deblur_summarize(luigi.Task):
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
        required_tasks["deblur"] = run_deblur(**kwargs)
        required_tasks["view_rep_seq"] = view_rep_seq(**kwargs)
        required_tasks["joined"] = joined_summarize(**kwargs)
        required_tasks["joined_qc"] = qualityFilter(**kwargs)
        return required_tasks

    def output(self):
        otutab, rep, stats = [_.path
                              for _ in self.input()["deblur"]]
        ofiles = [otutab.replace('.qza','.csv'),
                  rep.replace('.qza','.fa'),
                  stats.replace('.qza','.csv'),
                  ]
        return [luigi.LocalTarget(_) for _ in ofiles]

    def run(self):
        otutab, rep, stats = [_.path
                              for _ in self.input()["deblur"]]
        convert2otutab(otutab,
                       otutab.replace('.qza','.csv'))
        convert2seq(rep,
                    rep.replace('.qza','.fa'))
        convert2table(stats,
                       stats.replace('.qza','.csv'))