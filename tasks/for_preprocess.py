from subprocess import check_call
import sys
from os.path import join, dirname, basename,exists

sys.path.insert(0, dirname(dirname(__file__)))
import luigi
from config import soft_db_path
from config import default_params
from tasks.basic_tasks import base_luigi_task
from toolkit import run_cmd, valid_path
from input_parser import fileparser
from Bio import SeqIO
import gzip

vsearch = soft_db_path.vsearch_pth
rdp_gold = soft_db_path.rdp_gold_pth
map_pl = soft_db_path.map_pl_pth
trimmomatic_jar = soft_db_path.trimmomatic_jar
trimmomatic_dir = dirname(trimmomatic_jar)


class QC_trimmomatic(base_luigi_task):
    sampleid = luigi.Parameter()
    odir = luigi.Parameter()
    PE1 = luigi.Parameter()
    PE2 = luigi.Parameter(default=None)
    prefix = luigi.Parameter(default="OTU")

    # def requires(self):
    #     from luigi_workflow.main import workflow
    #     return workflow()
    #     pass
    #
    def output(self):
        odir = join(str(self.odir),
                    "%s_pipelines" % self.prefix,
                    "preprocessed",
                    "after_QC", )

        ofile_name1 = join(str(odir),
                           "{}_R1.clean.fq.gz".format(str(self.sampleid)))
        ofile_name2 = join(str(odir),
                           "{}_R2.clean.fq.gz".format(str(self.sampleid)))
        valid_path(ofile_name1, check_ofile=1)
        
        if not exists(self.PE2):
            return [luigi.LocalTarget(ofile_name1)]
        else:
            return [luigi.LocalTarget(ofile_name1),
                luigi.LocalTarget(ofile_name2)]

    def run(self):
        valid_path(self.output()[0].path, check_ofile=1)
        # auto make output dir
        sample_name = str(self.sampleid)
        # if len(self.input()) == 2:
        #     input1 = self.input()[0].path
        #     input2 = self.input()[1].path
        # else:
        input1 = self.PE1
        input2 = self.PE2

        if not exists(self.PE2):
            cmdline = f"java -jar {trimmomatic_jar} SE -threads {default_params.trimmomatic_thread} {input1} {self.output()[0].path} ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
            
        else:
            cmdline = "java -jar {trimmomatic_jar} PE -threads {thread} {input1} {input2} {ofile1} {outdir}/{PE1_id}.unpaired.fq.gz {ofile2} {outdir}/{PE2_id}.unpaired.fq.gz ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50".format(
            trimmomatic_jar=trimmomatic_jar,
            trimmomatic_dir=trimmomatic_dir,
            input1=input1,
            input2=input2,
            PE1_id=sample_name + "_R1",
            PE2_id=sample_name + "_R2",
            ofile1=self.output()[0].path,
            ofile2=self.output()[1].path,
            outdir=dirname(self.output()[0].path),
            thread=default_params.trimmomatic_thread)
        
        run_cmd(cmdline,
                dry_run=self.dry_run,
                log_file=self.get_log_path())
        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class screen_false(base_luigi_task):
    sampleid = luigi.Parameter()
    odir = luigi.Parameter()
    PE1 = luigi.Parameter()
    PE2 = luigi.Parameter(default=None)

    def requires(self):
        return QC_trimmomatic(sampleid=self.sampleid,
                              PE1=self.PE1,
                              PE2=self.PE2,
                              odir=self.odir,
                              log_path=self.log_path,
                              dry_run=self.dry_run)

    def output(self):
        odir = join(str(self.odir),
                    "OTU_pipelines",
                    "preprocessed",
                    "screened", )

        ofile_name1 = join(str(odir),
                           "{}_R1.fq.gz".format(self.sampleid))
        ofile_name2 = join(str(odir),
                           "{}_R2.fq.gz".format(self.sampleid))
        valid_path(ofile_name1, check_ofile=1)
        if exists(self.PE2):
        
            return [luigi.LocalTarget(ofile_name1),
                luigi.LocalTarget(ofile_name2)]
        else:
            return [luigi.LocalTarget(ofile_name1)]

    def run(self):
        fq_screen = soft_db_path.fq_screen
        if exists(self.PE2):
            clean_pe1 = self.input()[0].path
            clean_pe2 = self.input()[1].path
            infiles = ' '.join([clean_pe1, clean_pe2])
        else:
            clean_pe1 = self.input()[0].path
            clean_pe2 = ''
            infiles = clean_pe1
            
        outdir = join(dirname(dirname(self.output()[0].path)),
                      "_screened_cache")
        valid_path(outdir, check_odir=1)
        # anyone could be ok
        threads = default_params.fq_screen_thread
        cmdline = f"{fq_screen} {infiles} --outdir {outdir} --nohits --aligner bowtie2 --threads {threads}"
        run_cmd(cmdline,
                log_file=self.get_log_path(),
                dry_run=self.dry_run,
                )
        ############################################################
        # renamed
        name1 = basename(clean_pe1).replace('.fq.gz', '')
        name2 = basename(clean_pe2).replace('.fq.gz', '')
        filtered_r1 = join(outdir,
                           "{}.tagged_filter.fastq.gz".format(name1))
        filtered_r2 = join(outdir,
                           "{}.tagged_filter.fastq.gz".format(name2))
        
        for filtered_f, output_target in zip([filtered_r1, 
                                              filtered_r2],
                                             self.output()):
            opath = output_target.path.replace('.gz', '')
            with open(opath, 'w') as new_file:
                stream = SeqIO.parse(gzip.open(filtered_f, 'rt'), format='fastq')
                for read in stream:
                    _cache = str(read.description)
                    read.id = read.description = read.name = ''
                    _cache = _cache.rpartition('#FQST')[0]
                    read.id = _cache
                    SeqIO.write(read, new_file, format='fastq')
            run_cmd("gzip -f %s" % opath,
                    dry_run=self.dry_run,
                    log_file=self.get_log_path())

        if self.dry_run:
            for _o in self.output():
                run_cmd("touch %s" % _o.path, dry_run=False)


class QC_aft_screened(QC_trimmomatic):

    def requires(self):
        return screen_false(sampleid=self.sampleid,
                            PE1=self.PE1,
                            PE2=self.PE2,
                            odir=self.odir,
                            log_path=self.log_path,
                            dry_run=self.dry_run)

    def output(self):
        odir = join(str(self.odir),
                    "OTU_pipelines" ,
                    "preprocessed",
                    "screened_QC", )

        ofile_name1 = join(odir,
                           "{}_R1.clean.fq.gz".format(str(self.sampleid)))
        ofile_name2 = join(odir,
                           "{}_R2.clean.fq.gz".format(str(self.sampleid)))
        valid_path(ofile_name1, check_ofile=1)
        if exists(self.PE2):
            return [luigi.LocalTarget(ofile_name1),
                luigi.LocalTarget(ofile_name2)]
        else:
            return [luigi.LocalTarget(ofile_name1)]


#
# class screen_summary(base_luigi_task):
#     def requires(self):
#         df = fileparser(self.tab)
#         R1_dict = df.R1
#         R2_dict = df.R2
#         tasks = {}
#         for sid in R1_dict.keys():
#             tasks[sid] = screen_false(sampleid=sid,
#                                       PE1=R1_dict[sid],
#                                       PE2=R2_dict[sid],
#                                       odir=self.odir,
#                                       log_path=self.log_path,
#                                       dry_run=self.dry_run)
#         return tasks
#
#     def output(self):
#         new_data_input = join(str(self.odir),
#                               "preprocessed",
#                               "screened",
#                               "new_input.tsv")
#         valid_path(new_data_input,check_ofile=1)
#         return luigi.LocalTarget(new_data_input)
#
#     def run(self):
#         if self.dry_run:
#             run_cmd("touch %s" % self.output().path, dry_run=False)
#         total_df = pd.DataFrame(columns=["R1", "R2"])
#         for sid in self.input():
#             filtered_r1 = self.input()[sid][0].path
#             filtered_r2 = self.input()[sid][1].path
#             opaths = []
#             for f in [filtered_r1, filtered_r2]:
#                 opath = f.replace(".tagged_filter.fastq.gz",
#                                   ".fq")
#                 with open(opath, 'w') as new_file:
#                     stream = SeqIO.parse(gzip.open(f, 'rt'), format='fastq')
#                     for read in stream:
#                         _cache = read.description
#                         read.id = read.description = read.name = ''
#                         _cache = _cache.rpartition('#FQST')[0]
#                         read.description = _cache
#                         SeqIO.write(read, new_file, format='fastq')
#                 run_cmd("gzip -f %s" % opath,
#                         dry_run=self.dry_run,
#                         log_file=self.get_log_path())
#                 opaths.append(opath + '.gz')
#             total_df = total_df.append(pd.DataFrame(columns=total_df.columns,
#                                                     index=[sid],
#                                                     data=[opaths]))
#         total_df.index.name = "Sample ID"
#         total_df.to_csv(self.output().path, index=1, sep='\t')
#


class joined_reads(base_luigi_task):
    sampleid = luigi.Parameter()
    odir = luigi.Parameter()
    PE1 = luigi.Parameter()
    PE2 = luigi.Parameter(default=None)
    
    def requires(self):
        if self.screen:
            return QC_aft_screened(sampleid=self.sampleid,
                               PE1=self.PE1,
                               PE2=self.PE2,
                               odir=self.odir,
                               log_path=self.log_path,
                               dry_run=self.dry_run)
        else:
            return QC_trimmomatic(sampleid=self.sampleid,
                               PE1=self.PE1,
                               PE2=self.PE2,
                               odir=self.odir,
                               log_path=self.log_path,
                               dry_run=self.dry_run)
    def output(self):
        odir = str(self.odir)

        ofile_name = join(str(odir),
                          "OTU_pipelines",
                          "preprocessed",
                          "joined_reads",
                          "{}.fq.gz".format(str(self.sampleid)))
        valid_path(ofile_name, check_ofile=1)
        return luigi.LocalTarget(ofile_name)

    def run(self):
        from static.join_pairs import _join_pairs_w_command_output
        clean_pe1 = self.input()[0].path
        if exists(self.PE2):
            clean_pe2 = self.input()[1].path
            if not self.dry_run:
                _join_pairs_w_command_output(fwd_fp=clean_pe1,
                                            rev_fp=clean_pe2,
                                            fastq_out=self.output().path.replace('.gz', ''),
                                            log_file=self.log_path)
            else:
                run_cmd("touch %s" % self.output().path, dry_run=False),
        else:
            cmd = f"ln -s `realpath {clean_pe1}` {self.output().path}"
            run_cmd(cmd,
                    dry_run=self.dry_run,
                    log_file=self.get_log_path())

class merged_reads(base_luigi_task):

    def requires(self):
        df = fileparser(self.tab)
        R1_dict = df.R1
        R2_dict = df.R2
        tasks = {}
        for sid in R1_dict.keys():
            tasks[sid] = joined_reads(sampleid=sid,
                                      PE1=R1_dict[sid],
                                      PE2=R2_dict[sid],
                                      odir=self.odir,
                                      screen=False,
                                      log_path=self.log_path,
                                      dry_run=self.dry_run)
        return tasks

    def output(self):
        odir = join(str(self.odir),
                    "OTU_pipelines",
                    'preprocessed')
        ofile = join(odir, 'merged.fastq')
        valid_path(ofile, check_ofile=1)
        return luigi.LocalTarget(ofile)

    def run(self):
        if self.dry_run:
            run_cmd("touch %s" % self.output().path, dry_run=False)

        with open(self.output().path, 'w') as f1:
            count = 0
            for sid in self.input():
                joined_fq = self.input()[sid].path
                try:
                    stream = SeqIO.parse(gzip.open(joined_fq, 'rt'), format='fastq')
                    for read in stream:
                        read.id = read.description = read.name = ''
                        read.id = '{sid}_{num};barcodelabel={sid}'.format(sid=sid,
                                                                        num=str(count))
                        count += 1
                        SeqIO.write(read, f1, format='fastq')
                except :
                    print(f" {joined_fq} is mal-formatted, deleting. Please re-run the whole pipeline")
                    run_cmd(f"rm {joined_fq}",dry_run=self.dry_run,log_file=self.get_log_path())


############################################################
# For qiime2 preprocessing.

class import_data(base_luigi_task):

    def requires(self):
        df = fileparser(self.tab)
        R1_dict = df.R1
        R2_dict = df.R2
        tasks = {}
        for sid in R1_dict.keys():
            tasks[sid] = QC_aft_screened(sampleid=sid,
                                         PE1=R1_dict[sid],
                                         PE2=R2_dict[sid],
                                         odir=self.odir,
                                         log_path=self.log_path,
                                         dry_run=self.dry_run,
                                         )
        return tasks

    def output(self):
        odir = join(str(self.odir),
                    "q2_pipelines")
        ofile = join(odir,
                     'raw_data.qza')
        return luigi.LocalTarget(ofile)

    def run(self):
        import config.default_params as config
        collect_params = dict(ids=[],
                              r1_files=[],
                              r2_files=[])
        for sid in self.input():
            r1_path = self.input()[sid][0].path
            r2_path = self.input()[sid][1].path
            collect_params["ids"].append(sid)
            collect_params["r1_files"].append(r1_path)
            collect_params["r2_files"].append(r2_path)

        from static.utils import write_manifest
        opath = self.output().path.replace('.qza','')
        write_manifest(opath=opath,
                       **collect_params)
        cmd = "{qiime2_p} tools import \
          --type 'SampleData[PairedEndSequencesWithQuality]' \
          --input-path {manifest} \
          --output-path {prefix}.qza \
          --input-format PairedEndFastqManifestPhred33".format(
            qiime2_p=config.qiime2_p,
            manifest=opath,
            prefix=opath, )
        run_cmd(cmd,dry_run=self.dry_run,log_file=self.get_log_path())

