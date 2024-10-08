from os.path import join,dirname
import sys
sys.path.insert(0,dirname(dirname(__file__)))
import luigi

from config import soft_db_path
from tasks.basic_tasks import base_luigi_task
from toolkit import run_cmd, valid_path

usearch = soft_db_path.usearch_pth


class usearch_filter(base_luigi_task):
    def requires(self):
        from tasks.for_preprocess import merged_reads
        kwargs = self.get_kwargs()
        return merged_reads(**kwargs)

    def output(self):
        odir = join(str(self.odir),
                    "USEARCH",)
        ofile = join(odir,
                     'filtered.fa')
        valid_path(ofile,check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        merged_fq = self.input().path
        filtered_fa = self.output().path
        valid_path(filtered_fa, check_ofile=1)
        maxee = self.get_config_params('maxee')
        cmd = f"{usearch} -fastq_filter {merged_fq} -fastq_maxee {maxee} -fastaout {filtered_fa}"
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class usearch_derep(base_luigi_task):

    def requires(self):
        kwargs = self.get_kwargs()
        return usearch_filter(**kwargs)

    def output(self):
        odir = join(str(self.odir), "USEARCH",)
        ofile = join(odir, 'derep.fa')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        filtered_fa = self.input().path
        derep_fa = self.output().path
        cmd = f"{usearch} -fastx_uniques {filtered_fa} -fastaout {derep_fa} -sizeout"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class usearch_OTU(base_luigi_task):
    def requires(self):
        kwargs = self.get_kwargs()
        return usearch_derep(**kwargs)
    def output(self):
        odir = join(str(self.odir), "USEARCH","OTU")
        ofile = join(odir, 'otus_rep.fa')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        derep_fa = self.input().path
        OTU_rep_fa = self.output().path
        valid_path(OTU_rep_fa, check_ofile=1)
        cmd = f"{usearch} -cluster_otus {derep_fa} -otus {OTU_rep_fa} -relabel OTU"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class usearch_OTU_table(base_luigi_task):
    def requires(self):
        kwargs = self.get_kwargs()
        required_task = {}
        required_task["filtered"] = usearch_filter(**kwargs)
        required_task["OTU_rep"] = usearch_OTU(**kwargs)
        return required_task

    def output(self):
        odir = join(str(self.odir),"USEARCH","OTU")
        ofile = join(odir, 'otu_raw.tab')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        filtered_fa = self.input()["filtered"].path
        OTU_rep_fa = self.input()["OTU_rep"].path
        otu_raw_table = self.output().path

        cmd = f"{usearch} -otutab {filtered_fa} -otus {OTU_rep_fa} -otutabout {otu_raw_table}"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)

## denoise part
class usearch_denoise(base_luigi_task):
    get_positive = luigi.BoolParameter(default=False)
    def requires(self):
        kwargs = self.get_kwargs()
        return usearch_derep(**kwargs)

    def output(self):
        odir = join(str(self.odir),"USEARCH",'zotu')
        zotu_rep = join(odir, 'zotus.fa')
        zotu_mapping = join(odir, 'unoise3.txt')
    
        valid_path([zotu_rep,zotu_mapping], check_ofile=1)
        if not self.get_positive:
            return [luigi.LocalTarget(zotu_rep),luigi.LocalTarget(zotu_mapping)]
        else:
            positive_rep = join(odir, 'Positive_rep.fasta')
            return [luigi.LocalTarget(zotu_rep),luigi.LocalTarget(zotu_mapping),luigi.LocalTarget(positive_rep)]
    
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
            with open(zotu_rep.replace('.fa','.relabelled.fa'),'w') as f1:
                SeqIO.write(r,f1,'fasta-2line')
                
        if self.get_positive:
            self.get_positive_seqs(self.output()[0].path.replace('.fa','.relabelled.fa'))

            
            
class usearch_zOTU_table(base_luigi_task):
    
    def requires(self):
        kwargs = self.get_kwargs()
        required_task = {}
        required_task["filtered"] = usearch_filter(**kwargs)
        required_task["zOTU_rep"] = usearch_denoise(**kwargs)
        return required_task

    def output(self):
        odir = join(str(self.odir),"USEARCH",'zotu')
        ofile = join(odir, 'zotutab_raw.txt')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        filtered_fa = self.input()["filtered"].path
        zotu_rep_fa = self.input()["zOTU_rep"][0].path.replace('.fa','.relabelled.fa')
        zotu_otu_table = self.output().path
        valid_path(zotu_otu_table, check_ofile=1)
        cmd = f"{usearch} -otutab {filtered_fa} -otus {zotu_rep_fa} -otutabout {zotu_otu_table}"
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)
        else:
            run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())


if __name__ == '__main__':
    luigi.run()
