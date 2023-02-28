import luigi
from tasks.basic_tasks import base_luigi_task
from os.path import dirname,join
from config import soft_db_path
from toolkit import run_cmd,valid_path

class get_tree(base_luigi_task):
    analysis_type = luigi.Parameter(default="otu")

    def requires(self):
        self.get_config()
        kwargs = self.get_kwargs()
        required_tasks = {}
        
        
        antype = str(self.analysis_type).lower().split(',')
        if len(set(antype).intersection(set(["all", "otu", "deblur", "dada2",'usearch','qc','pdada2'])))==0:
            raise Exception("analysis type must be one of the `all,otu,deblur,dada2`")
        if 'qc' in antype:
            from tasks.for_preprocess import multiqc
            required_tasks["fastqc_before"] = multiqc(status='before',
                                                    **kwargs
                                                    )
            required_tasks["fastqc_after"] = multiqc(status='after',
                                                    **kwargs
                                                    )
            return required_tasks
        #####
        
        if  "otu"  in antype or  "all" in antype:
            from tasks.for_otu import vsearch_otutable
            required_tasks["otu"] = vsearch_otutable(**kwargs)
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
        if "pdada2" in antype or  "all" in antype:
            from tasks.for_Pdada2 import run_Pdada2
            required_tasks["run_Pdada2"] = run_Pdada2(**kwargs)      
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
                               "rep.tree"))
        ofiles = list(set(ofiles))
        valid_path(ofiles,check_ofile=1)
        return [luigi.LocalTarget(_)
                for _ in ofiles]


    def run(self):
        if self.dry_run:
            pass
        for k, f in self.input().items():
            # if not '_rep' in k:
            #     continue
            if type(f) == list:
                ofile = f[0].path
            else:
                ofile = f.path
            odir = dirname(ofile)
            cmd = "{usearch} -cluster_agg {rep_fa} -treeout {otree}".format(
                    usearch=soft_db_path.usearch_pth,
                    rep_fa=ofile,
                    otree=join(odir,"rep.tree"),
                )
            run_cmd(cmd,log_file=self.get_log_path(),dry_run=self.dry_run)
