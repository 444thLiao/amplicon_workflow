import luigi

from toolkit import run_cmd


class base_luigi_task(luigi.Task):
    odir = luigi.Parameter()
    tab = luigi.Parameter(default=None)
    dry_run = luigi.BoolParameter(default=False)
    log_path = luigi.Parameter(default=None)

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
