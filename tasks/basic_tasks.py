import luigi

from toolkit import run_cmd
from os.path import *
import os,sys

class base_luigi_task(luigi.Task):
    odir = luigi.Parameter()
    tab = luigi.Parameter(default=None)
    dry_run = luigi.BoolParameter(default=False)
    log_path = luigi.Parameter(default=None)
    screen = luigi.Parameter(default=False)
    config = luigi.Parameter(default=False)
    
    def get_log_path(self):
        base_log_path = self.log_path
        if base_log_path is not None:
            return base_log_path

    def get_kwargs(self):
        kwargs = dict(odir=self.odir,
                      tab=self.tab,
                      dry_run=self.dry_run,
                      log_path=self.log_path)
        return kwargs

    def get_config(self):
        if not self.config:
            self.config_p = f"{dirname(dirname(__file__))+'/config/default_params.py'}"
        else:
            self.config_p = self.config
        
        os.makedirs(join(self.odir,'tmp_import'))
        os.system(f"cat {self.config_p} > {join(self.odir,'tmp_import','__init__.py')}")
        sys.insert(0,join(self.odir,'tmp_import'))
        import tmp_import as config_params
        
        
    def get_params(self,arg):
        
        return self.config_params.get(arg,'')
    
    def batch_get_params(self,):
        extra_str = ''
        for p, val in self.config_params.join_params.items():
            if val is True:
                extra_str += ' --p-%s' % p
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p, val)        
class visulize_seq(base_luigi_task):
    """
    mainly for visualizing SingleFastqFormat qza
    """

    def output(self):
        return luigi.LocalTarget(self.input().path.replace(".qza", ".qzv"))

    def run(self):
        cmd = "qiime demux summarize --i-data {input_f} --o-visualization {output_f}".format(
            input_f=self.input().path,
            output_f=self.output().path)
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False, )


class tabulate_seq(base_luigi_task):
    """
    mainly for visualizing representative sequence
    """

    def output(self):
        return luigi.LocalTarget(self.input()[1].path.replace(".qza", ".qzv"))

    def run(self):
        cmd = "qiime feature-table tabulate-seqs --i-data {input_f} --o-visualization {output_f}".format(
            input_f=self.input()[1].path,
            output_f=self.output().path)
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False, )
