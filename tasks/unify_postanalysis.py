import luigi
from tasks import base_luigi_task
from os.path import dirname,join
from config import soft_db_path
from toolkit import run_cmd,valid_path

class get_tree(base_luigi_task):
    analysis_type = luigi.Parameter(default="otu")

    def requires(self):
        kwargs = dict(odir=self.odir,
                      tab=self.tab,
                      dry_run=self.dry_run,
                      log_path=self.log_path)

        antype = str(self.analysis_type).lower()
        if antype not in ["all", "otu", "deblur", "dada2"]:
            raise Exception("analysis type must be one of the `all,otu,deblur,dada2`")
        required_tasks = {}
        if antype == "otu" or antype == "all":
            from tasks.for_otu import vsearch_otutable
            required_tasks["otu"] = vsearch_otutable(**kwargs)
        if antype == "dada2" or antype == "all":
            from tasks.for_dada2 import dada2_summarize
            required_tasks["dada2"] = dada2_summarize(**kwargs)
        if antype == "deblur" or antype == "all":
            from tasks.for_deblur import deblur_summarize
            required_tasks["deblur"] = deblur_summarize(**kwargs)
        return required_tasks

    def output(self):
        ofiles = []
        for k, f in self.input().items():
            odir = dirname(f[1].path)
            ofiles.append(join(odir,"rep.tree"))
        valid_path(ofiles,check_ofile=1)
        return [luigi.LocalTarget(_) for _ in ofiles]


    def run(self):
        if self.dry_run:
            pass

        for k, f in self.input().items():
            odir = dirname(f[0].path)
            if k == "otu":
                cmd = "{usearch} -cluster_agg {rep_fa} -treeout {otree}".format(
                    usearch=soft_db_path.usearch_pth,
                    rep_fa=join(odir,"OTU_rep.fasta"),
                    otree=join(odir,"rep.tree"),
                )
                run_cmd(cmd,log_file=self.get_log_path(),dry_run=self.dry_run)
            else:
                cmd = "{usearch} -cluster_agg {rep_fa} -treeout {otree}".format(
                    usearch=soft_db_path.usearch_pth,
                    rep_fa=join(odir,f[1].path),
                    otree=join(odir,"rep.tree"),
                )
                run_cmd(cmd,log_file=self.get_log_path(),dry_run=self.dry_run)
