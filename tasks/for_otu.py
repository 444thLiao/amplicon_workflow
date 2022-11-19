from os.path import join,dirname
import sys
sys.path.insert(0,dirname(dirname(__file__)))
import luigi,os
from config import soft_db_path
from tasks.basic_tasks import base_luigi_task
from toolkit import run_cmd, valid_path
from Bio import SeqIO

vsearch = soft_db_path.vsearch_pth
rdp_gold = soft_db_path.rdp_gold_pth
map_pl = soft_db_path.map_pl_pth


# trunclen = config.vesearch_args['trunclen']
# cluster_ratio = config.vesearch_args['cluster_ratio']


#
#
# class fq_screen(luigi.Task):
#     def requires(self):
#         return preprocess()
#
#     def output(self):
#         pass
#
#     def run(self):
#         cmd = "fastq_screen {input} --outdir {outdir} --nohits --aligner bowtie2".format(input=fqs, outdir=outdir)
#
#
# class fq_merge(luigi.Task):
#     def requires(self):
#         return fq_screen()
#
#     def output(self):
#         pass
#
#     def run(self):
#
#         def task_merge(indir, dirpath):
#             # dirpath could be 1. screen_dir or 2. joined_seq_path
#             indir = os.path.join(indir, dirpath, '*.fastq.gz')
#             odir = os.path.join(indir, merged_dir)
#             os.makedirs(odir, exist_ok=True)
#             opath = os.path.join(odir, merged_file)
#             # splitted by sample id; merge sequencing file joined after join_pairs.py
#             all_fqs = sorted(glob(os.path.join(indir, '*.gz')))
#             with open(opath, 'w') as f1:
#                 count = 0
#                 for fq in tqdm(all_fqs):
#                     stream = SeqIO.parse(gzip.open(fq, 'rt'), format='fastq')
#                     sid = os.path.basename(fq).replace('fastq.gz', '')
#                     for read in stream:
#                         read.id = read.description = read.name = ''
#                         read.id = '{sid}_{num};barcodelabel={sid}'.format(sid=sid,
#                                                                           num=str(count))
#                         count += 1
#                         SeqIO.write(read, f1, format='fastq')


# class input_for_vsearch_pp(luigi.ExternalTask):
#     odir = luigi.Parameter()
#     dry_run = luigi.BoolParameter()
#
#     def output(self):
#         odir = join(str(self.odir),
#                             'preprocessed')
#         ofile = join(odir, 'merged.fastq')
#         return luigi.LocalTarget(ofile)


class vsearch_filter(base_luigi_task):
    def requires(self):
        from tasks import merged_reads
        return merged_reads(tab=self.tab,
                            odir=self.odir,
                            dry_run=self.dry_run,
                            log_path=self.log_path)

    def output(self):
        odir = join(str(self.odir),
                    "OTU_pipelines",
                    "preprocessed",
                    'QC')
        ofile = join(odir,
                     'filtered.fasta')
        valid_path(ofile,check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        # cmd = "vsearch --fastx_filter {input} --fastq_maxee 1.0 --fastq_trunclen xx --fastaout {output}".format(input=input,output=output)
        merged_fa = self.input().path
        filtered_fa = self.output().path
        total_len = int(os.popen(f"grep -c '^+$' {merged_fa}").read().strip())
        length_dis = [len(_.seq) for _ in SeqIO.parse(merged_fa,'fastq')]
        
        trunclen = self.get_params('trunclen')
        valid_path(filtered_fa, check_ofile=1)
        cmd = f"{vsearch} --fastx_filter {merged_fa} --fastq_maxee 1.0 --fastq_trunclen {trunclen} --fastaout {filtered_fa}"
        run_cmd(cmd,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_derep(base_luigi_task):

    def requires(self):
        return vsearch_filter(tab=self.tab,odir=self.odir,
                              dry_run=self.dry_run,
                              log_path=self.log_path)

    def output(self):
        odir = join(str(self.odir), "OTU_pipelines",'derep')
        ofile = join(odir, 'derep.fa')

        o_uc = ofile.replace('.fa', '.uc')
        valid_path([ofile,o_uc], check_ofile=1)
        return [luigi.LocalTarget(ofile),
                luigi.LocalTarget(o_uc)]

    def run(self):
        filtered_fa = self.input().path
        derep_fa = self.output()[0].path
        derep_uc = self.output()[1].path
        cmd = f"{vsearch} --derep_fulllength {filtered_fa} --output {derep_fa} -sizeout --fasta_width 0 --uc {derep_uc}"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_rmSingle(base_luigi_task):

    def requires(self):
        return vsearch_derep(tab=self.tab,odir=self.odir, dry_run=self.dry_run,
                             log_path=self.log_path)

    def output(self):
        odir = join(str(self.odir),"OTU_pipelines", 'derep')
        ofile = join(odir, 'sorted_d2_fa.fa')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        derep_fa = self.input()[0].path
        sorted_d2_fa = self.output().path
        valid_path([sorted_d2_fa], check_ofile=1)
        cmd = f"{vsearch} --sortbysize  {derep_fa} --output {sorted_d2_fa} --minsize 1"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_pre_cluster(base_luigi_task):
    def requires(self):
        return vsearch_derep(tab=self.tab,odir=self.odir, dry_run=self.dry_run, log_path=self.log_path)

    def output(self):
        odir = join(str(self.odir),"OTU_pipelines", 'dechimera')
        ofile = join(odir, 'all.preclustered.fasta')

        o_uc = ofile.replace('.fasta', '.uc')
        valid_path([ofile, o_uc], check_ofile=1)
        return [luigi.LocalTarget(ofile),
                luigi.LocalTarget(o_uc)]

    def run(self):
        # cmd = """vsearch --cluster_size splited/derep.fa \
        # 	--id xxx --sizeout --fasta_width 0 \
        # 	--uc v_analysis/all.preclustered.uc \
        # 	--centroids v_analysis/all.preclustered.fasta"""
        derep_fa = self.input()[0].path
        precluster_uc = self.output()[1].path
        precluster_fa = self.output()[0].path
        cluster_ratio = self.get_params('cluster_ratio')
        cmd = f"{vsearch} --cluster_size {derep_fa} --id {cluster_ratio} --sizeout --fasta_width 0 --uc {precluster_uc} --centroids {precluster_fa}"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_uchime_deno(base_luigi_task):
    def requires(self):
        return vsearch_pre_cluster(tab=self.tab,odir=self.odir, dry_run=self.dry_run, log_path=self.log_path)

    def output(self):
        odir = join(str(self.odir),"OTU_pipelines", 'dechimera')
        ofile = join(odir, 'all.denovo.nonchimeras.fasta')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    # return denovo_nonchimer_fa
    def run(self):
        # cmd2 = """vsearch --uchime_deno v_analysis/all.preclustered.fasta \
        # 	--sizein --sizeout --fasta_width 0 \
        # 	--nonchimeras v_analysis/all.denovo.nonchimeras.fastaq"""
        precluster_fa = self.input()[0].path
        denovo_nonchimer_fa = self.output().path
        valid_path([denovo_nonchimer_fa], check_ofile=1)
        cmd = f"{vsearch} --uchime_deno {precluster_fa} --sizein --sizeout --fasta_width 0 --nonchimeras {denovo_nonchimer_fa}"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_uchime_ref(base_luigi_task):
    def requires(self):
        return vsearch_uchime_deno(tab=self.tab,odir=self.odir, dry_run=self.dry_run,
                                   log_path=self.log_path)

    def output(self):
        odir = join(str(self.odir),"OTU_pipelines", 'dechimera')
        ofile = join(odir, 'all.ref.nonchimeras.fasta')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        # cmd3 = """vsearch --uchime_ref v_analysis/all.denovo.nonchimeras.fastaq \
        # 	--db /home/liaoth/data2/project/16s_pipelines/microbiome_utils/vsearch_pipeliens/rdp_gold.fa --sizein --sizeout --fasta_width 0 \
        # 	--nonchimeras v_analysis/all.ref.nonchimeras.fasta"""
        denovo_nonchimer_fa = self.input().path
        ref_nonchimer_fa = self.output().path
        valid_path(ref_nonchimer_fa, check_ofile=1)
        cmd = f"{vsearch} --uchime_ref {denovo_nonchimer_fa} --db {rdp_gold} --sizein --sizeout --fasta_width 0 --nonchimeras {ref_nonchimer_fa}"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_map1(base_luigi_task):
    def requires(self):
        required_task = {}
        required_task["derep"] = vsearch_derep(tab=self.tab,odir=self.odir, dry_run=self.dry_run, log_path=self.log_path)
        required_task["precluster"] = vsearch_pre_cluster(tab=self.tab,odir=self.odir, dry_run=self.dry_run, log_path=self.log_path)
        required_task["uchime_ref"] = vsearch_uchime_ref(tab=self.tab,odir=self.odir, dry_run=self.dry_run, log_path=self.log_path)
        return required_task

    def output(self):
        odir = join(str(self.odir),"OTU_pipelines", 'dechimera')
        ofile = join(odir, 'all.nonchimeras.derep.fasta')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        # cmd4 = """perl /home/liaoth/data2/project/16s_pipelines/microbiome_utils/vsearch_pipeliens/map.pl splited/derep.fa v_analysis/all.preclustered.uc v_analysis/all.ref.nonchimeras.fasta > v_analysis/all.nonchimeras.derep.fasta"""
        derep_fa = self.input()["derep"][0].path
        precluster_uc = self.input()["precluster"][1].path
        ref_nonchimer_fa = self.input()["uchime_ref"].path
        nonchimera_derep_fa = self.output().path
        valid_path(nonchimera_derep_fa, check_ofile=1)
        cmd = f"perl {map_pl} {derep_fa} {precluster_uc} {ref_nonchimer_fa} > {nonchimera_derep_fa}"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_map2(base_luigi_task):

    def requires(self):
        required_task = {}
        required_task["map1"] = vsearch_map1(tab=self.tab,odir=self.odir, dry_run=self.dry_run, log_path=self.log_path)
        required_task["derep"] = vsearch_derep(tab=self.tab,odir=self.odir, dry_run=self.dry_run, log_path=self.log_path)
        return required_task

    def output(self):
        odir = join(str(self.odir),"OTU_pipelines", 'dechimera')
        ofile = join(odir, 'all.nonchimeras.fasta')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        # cmd5 = """perl /home/liaoth/data2/project/16s_pipelines/microbiome_utils/vsearch_pipeliens/map.pl splited/derep.fa splited/all.derep.uc v_analysis/all.nonchimeras.derep.fasta > v_analysis/all.nonchimeras.fasta"""
        derep_fa = self.input()["derep"][0].path
        derep_uc = self.input()["derep"][1].path
        nonchimera_derep_fa = self.input()["map1"].path
        nonchimera_fa = self.output().path
        valid_path(nonchimera_fa, check_ofile=1)
        cmd = f"perl {map_pl} {derep_fa} {derep_uc} {nonchimera_derep_fa} > {nonchimera_fa}"
        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in [self.output()]:
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_cluster(base_luigi_task):
    def requires(self):
        return vsearch_map2(tab=self.tab,odir=self.odir, dry_run=self.dry_run, log_path=self.log_path)

    def output(self):
        odir = join(str(self.odir),"OTU_pipelines", 'otu_output')
        ofile = join(odir, "OTU_rep.fasta")
        o_uc = ofile.replace('.fasta', '.uc')
        valid_path([ofile,o_uc], check_ofile=1)
        return [luigi.LocalTarget(ofile),
                luigi.LocalTarget(o_uc)]

    def run(self):
        # cmd = """vsearch --cluster_size v_analysis/all.nonchimeras.fasta --id 0.97 \
        # 	--sizein --sizeout --fasta_width 0 \
        # 	--uc v_analysis/all.clustered.uc \
        # 	--relabel OTU --centroids v_analysis/all.otus.fasta"""

        nonchimera_fa = self.input().path
        rep_fa = self.output()[0].path
        cluster_uc = self.output()[1].path
        valid_path([rep_fa, cluster_uc], check_ofile=1)
        cluster_ratio = self.get_params('cluster_ratio')
        cmd = f"{vsearch} --cluster_size {nonchimera_fa} --id {cluster_ratio} --sizein --sizeout --fasta_width 0 --uc {cluster_uc} --relabel OTU --centroids {rep_fa}"

        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class vsearch_otutable(base_luigi_task):
    def requires(self):
        required_task = {}
        required_task["filter"] = vsearch_filter(tab=self.tab,
                                                 odir=self.odir,
                                                 dry_run=self.dry_run,
                                                 log_path=self.log_path)
        required_task["cluster"] = vsearch_cluster(tab=self.tab,
                                                   odir=self.odir,
                                                   dry_run=self.dry_run,
                                                   log_path=self.log_path)
        return required_task

    def output(self):
        odir = join(str(self.odir),"OTU_pipelines", 'otu_output')
        ofile = join(odir, "raw_OTU.tab")

        mapfile = join(odir, "raw_OTU.uc")
        return [luigi.LocalTarget(ofile),
                luigi.LocalTarget(mapfile)]

    def run(self):
        # cmd2 = """vsearch --usearch_global splited/filtered_uparsed.fa --db v_analysis/all.otus.fasta --strand plus --id 0.97 --uc v_analysis/map.txt --otutabout v_analysis/otu_raw.tab    """
        filtered_fa = self.input()['filter'].path
        rep_fa = self.input()['cluster'][0].path

        map_output = self.output()[1].path
        raw_otutab = self.output()[0].path
        valid_path([map_output, raw_otutab], check_ofile=1)
        cluster_ratio = self.get_params('cluster_ratio')
        cmd = f"{vsearch} --usearch_global {filtered_fa} --db {rep_fa} --strand plus --id {cluster_ratio} --uc {map_output} --otutabout {raw_otutab}"

        run_cmd(cmd, dry_run=self.dry_run, log_file=self.get_log_path())

        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)



if __name__ == '__main__':
    luigi.run()
