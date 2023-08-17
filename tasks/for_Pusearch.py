from os.path import join,dirname
import sys
sys.path.insert(0,dirname(dirname(__file__)))
import luigi

from config import soft_db_path
from tasks.basic_tasks import base_luigi_task
from toolkit import run_cmd, valid_path

usearch = soft_db_path.usearch_pth


from tasks.for_otu import vsearch_derep,vsearch_filter


## denoise part
class Pusearch_denoise(base_luigi_task):
    def requires(self):
        kwargs = self.get_kwargs()
        return vsearch_derep(**kwargs)

    def output(self):
        odir = join(str(self.odir),"PUSEARCH")
        zotu_rep = join(odir, 'rep.fasta')
        zotu_mapping = join(odir, 'unoise3.txt')
        valid_path([zotu_rep,zotu_mapping], check_ofile=1)
        return [luigi.LocalTarget(zotu_rep),luigi.LocalTarget(zotu_mapping)]
    
    def run(self):
        derep_fa = self.input().path
        zotu_rep = self.output()[0].path
        zotu_mapping = self.output()[1].path
        valid_path([zotu_rep], check_ofile=1)
        cmd = f"{usearch} -unoise3 {derep_fa} -zotus {zotu_rep} -tabbedout {zotu_mapping}"
        

        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)
        else:
            run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())
            from Bio import SeqIO
            r = list(SeqIO.parse(zotu_rep,'fasta'))
            for _ in r:
                _.id = _.id.replace('Zotu','OTU')
                _.name= _.description = ''
            with open(zotu_rep.replace('.fasta','.relabelled.fasta'),'w') as f1:
                SeqIO.write(r,f1,'fasta-2line')


class usearch_zOTU_table(base_luigi_task):
    def requires(self):
        kwargs = self.get_kwargs()
        required_task = {}
        required_task["filtered"] = vsearch_filter(**kwargs)
        required_task["zOTU_rep"] = Pusearch_denoise(**kwargs)
        return required_task

    def output(self):
        odir = join(str(self.odir),"PUSEARCH")
        ofile = join(odir, 'ASV_table.tsv')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        filtered_fa = self.input()["filtered"].path
        zotu_rep_fa = self.input()["zOTU_rep"][0].path.replace('.fasta','.relabelled.fasta')
        zotu_otu_table = self.output().path
        valid_path(zotu_otu_table, check_ofile=1)
        cmd = f"{usearch} -otutab {filtered_fa} -otus {zotu_rep_fa} -otutabout {zotu_otu_table}"
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)
        else:
            run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())
