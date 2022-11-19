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
    config_params = luigi.Parameter(default=False)
    
    def get_log_path(self):
        base_log_path = self.log_path
        if base_log_path is not None:
            return base_log_path

    def get_kwargs(self):
        kwargs = dict(odir=self.odir,
                      tab=self.tab,
                      dry_run=self.dry_run,
                      log_path=self.log_path,
                      config=self.config,
                      config_params=self.config_params)
        return kwargs

    def get_config(self):
        self.config_p = f"{dirname(dirname(__file__))+'/config/default_params.py'}"
        sys.path.insert(0,dirname(dirname(__file__))+'/config')
        import config.default_params as config_params        
        if self.config:
            self.config_p = self.config
            os.makedirs(join(self.odir,'tmp_import'),exist_ok=True)
            os.system(f"cat {self.config_p} > {join(self.odir,'tmp_import','__init__.py')}")
            sys.path.insert(0,realpath(self.odir))
            
            import tmp_import as new_params
            for aparam in dir(new_params):
                if '__'  in aparam or aparam not in dir(config_params): continue
                if type(getattr(config_params,aparam))==dict:
                    getattr(config_params,aparam).update( getattr(new_params,aparam))
                else:
                    setattr(config_params,aparam,getattr(new_params,aparam))
            os.system(f"rm -r {join(self.odir,'tmp_import')}")
        self.config_params = config_params
        
        
    def get_config_params(self,arg):
        if type(arg) == str:
            if arg in dir(self.config_params):
                return getattr(self.config_params,arg)
        else:
            return getattr(self.config_params,arg[0])[arg[1]]
            
    def batch_get_config_params(self,name,new_key={}):
        params_dict = getattr(self.config_params,name)
        for k,v in new_key.items():
            params_dict[k] = v

        extra_str = ''
        for p, val in params_dict.items():
            if val is True:
                extra_str += ' --p-%s' % p.replace('_','-')
            elif val is not None and val is not False:
                extra_str += ' --p-%s %s ' % (p.replace('_','-'), val)        
        return extra_str
    
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
