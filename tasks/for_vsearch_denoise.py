from os.path import join,dirname
import sys
sys.path.insert(0,dirname(dirname(__file__)))
import luigi,os
from config import soft_db_path

from tasks.basic_tasks import base_luigi_task
from toolkit import run_cmd, valid_path
from Bio import SeqIO
from tasks.for_otu import vsearch_derep,vsearch_filter


vsearch = soft_db_path.vsearch_pth
rdp_gold = soft_db_path.rdp_gold_pth
map_pl = soft_db_path.map_pl_pth

# num_threads = config.fq_screen_thread

# https://github.com/mooreryan/rhode_island_16s/blob/b2006eeb14a87c6ea09079776271536532e44d3c/scripts/make_asv_table/3_make_asv_table.Rmd
class vsearch_denoise(base_luigi_task):
    def requires(self):
        kwargs = self.get_kwargs()
        return vsearch_derep(**kwargs)

    def output(self):
        odir = join(str(self.odir),"VSEARCH_denosing")
        ofile = join(odir, "unoised.fa")
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        derep_fa = self.input()[0].path
        unoise_fa = self.output().path
        valid_path([unoise_fa], check_ofile=1)
        cmd = f"""
        {vsearch} \
  --threads {self.get_config_params('num_threads')} \
  --cluster_unoise {derep_fa} \
  --centroids {unoise_fa} \
  --relabel "ASV" \
  --fasta_width 0 \
  --sizeout
  """
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_uchime(base_luigi_task):
    def requires(self):
        kwargs = self.get_kwargs()
        return vsearch_denoise(**kwargs)

    def output(self):
        odir = join(str(self.odir),"VSEARCH_denosing")
        ofile = join(odir, "dechimer.fa")
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        input = self.input()[0].path
        output = self.output().path
        valid_path([output], check_ofile=1)
        cmd = f"""
        {vsearch} \
  --uchime3_denovo {input} \
  --nonchimeras {output} \
  --fasta_width 0 \
  --xsize
  """
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_denoise_otutab(base_luigi_task):
    def requires(self):
        kwargs = self.get_kwargs()
        required_task = {}
        required_task["filter"] = vsearch_filter(**kwargs)
        required_task["nochimer"] = vsearch_uchime(**kwargs)
        return required_task

    def output(self):
        odir = join(str(self.odir),"VSEARCH_denosing")
        ofile = join(odir, "asv.otu")
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        filtered_fa = self.input()['filter'].path
        dechimer_rep_fa = self.input()['nochimer'].path
        output = self.output().path
        valid_path([output], check_ofile=1)
        cmd = f"""
        {vsearch} \
  --usearch_global {filtered_fa} \
  --db {dechimer_rep_fa} \
  --id 0.97 \
  --otutabout {output}
  --threads 20
  """
  # use 0.97 to retrieve raw reads
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)
