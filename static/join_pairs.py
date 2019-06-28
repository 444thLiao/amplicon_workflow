import sys
from os.path import dirname, join

sys.path.insert(0, dirname(dirname(__file__)))
from static.utils import assign_work_pool
import pandas as pd

from config.default_params import *
from config import soft_db_path
from toolkit import run_cmd


def cal(args):
    func, args = args
    return func(*args)


def join_pairs(manifest,
               n_thread,
               truncqual: int = join_params['truncqual'],
               minlen: int = join_params['minlen'],
               maxns: int = join_params['maxns'],
               allowmergestagger: bool = join_params['allowmergestagger'],
               minovlen: int = join_params['minovlen'],
               maxdiffs: int = join_params['maxdiffs'],
               minmergelen: int = join_params['minmergelen'],
               maxmergelen: int = join_params['maxmergelen'],
               maxee: float = join_params['maxee'],
               qmin: int = join_params['qmin'],
               qminout: int = join_params['qminout'],
               qmax: int = join_params['qmax'],
               qmaxout: int = join_params['qmaxout'],
               qascii: int = join_params['qascii'],
               ):
    dir_path = dirname(manifest)
    os.makedirs(join(dir_path, joined_seq_path), exist_ok=True)
    manifest = pd.read_csv(manifest,
                           header=0, comment='#')

    id_to_fps = manifest.pivot(index='sample-id',
                               columns='direction',
                               values='absolute-filepath')
    all_args = []

    for i, (sample_id, (fwd_fp, rev_fp)) in enumerate(id_to_fps.iterrows()):
        fastq_out = join(dir_path, joined_seq_path, sample_id + '.fastq')
        all_args.append((_join_pairs_w_command_output,
                         (fwd_fp, rev_fp, fastq_out, truncqual, minlen, maxns, allowmergestagger,
                          minovlen, maxdiffs, minmergelen, maxmergelen, maxee, qmin, qminout,
                          qmax, qmaxout, qascii)
                         ))

    assign_work_pool(cal, all_args, num_thread=n_thread)


def _join_pairs_w_command_output(fwd_fp,
                                 rev_fp,
                                 fastq_out,
                                 truncqual: int = join_params['truncqual'],
                                 minlen: int = join_params['minlen'],
                                 maxns: int = join_params['maxns'],
                                 allowmergestagger: bool = join_params['allowmergestagger'],
                                 minovlen: int = join_params['minovlen'],
                                 maxdiffs: int = join_params['maxdiffs'],
                                 minmergelen: int = join_params['minmergelen'],
                                 maxmergelen: int = join_params['maxmergelen'],
                                 maxee: float = join_params['maxee'],
                                 qmin: int = join_params['qmin'],
                                 qminout: int = join_params['qminout'],
                                 qmax: int = join_params['qmax'],
                                 qmaxout: int = join_params['qmaxout'],
                                 qascii: int = join_params['qascii'],
                                 log_file=None):
    # this function exists only to simplify unit testing
    cmd = [soft_db_path.vsearch_pth,
           '--fastq_mergepairs', fwd_fp,
           '--reverse', rev_fp,
           '--fastqout', fastq_out,
           '--fastq_ascii', str(qascii),
           '--fastq_minlen', str(minlen),
           '--fastq_minovlen', str(minovlen),
           '--fastq_maxdiffs', str(maxdiffs),
           '--fastq_qmin', str(qmin),
           '--fastq_qminout', str(qminout),
           '--fastq_qmax', str(qmax),
           '--fastq_qmaxout', str(qmaxout),
           ]
    if truncqual is not None:
        cmd += ['--fastq_truncqual', str(truncqual)]
    if maxns is not None:
        cmd += ['--fastq_maxns', str(maxns)]
    if minmergelen is not None:
        cmd += ['--fastq_minmergelen', str(minmergelen)]
    if maxmergelen is not None:
        cmd += ['--fastq_maxmergelen', str(maxmergelen)]
    if maxee is not None:
        cmd += ['--fastq_maxee', str(maxee)]
    if allowmergestagger:
        cmd.append('--fastq_allowmergestagger')
    run_cmd(' '.join(cmd), dry_run=False, log_file=log_file)
    run_cmd(' '.join(['gzip', '-f', fastq_out]), dry_run=False, log_file=log_file)

# if __name__ == '__main__':
#     from utils import parse_param
#     import argparse
#
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-p", "--parameter",
#                         help="input file contains all parameters, template is place at the %s. This is a python script actually, so just use python code to write your code." % join(
#                             dirname(abspath(__file__)), 'param.template'), default=join(dirname(__file__), 'param.template'))
#     parser.add_argument("-nt", "--num_thread",
#                         help="thread")
#
#     args = parser.parse_args()
#     parameter_f = args.parameter
#     nt = args.num_thread
#     parse_param(parameter_f, g=globals())
#
#     r1_files = sorted(get_files(indir, r1_format))
#     r2_files = sorted(get_files(indir, r2_format))
#     ids = [re.findall(idpattern, basename(_))[0].strip(string.punctuation) for _ in r1_files]
#
#     if demux_on:
#         print("Demux is on, first start the demultiplexing process.")
#         demux_dict['seqfile1'] = r1_files
#         demux_dict['seqfile2'] = r2_files
#         if not_overwrite_demux and isfile(demux_stats):
#             print('demuxed files exist, pass it')
#         else:
#             r1_files, r2_files, ids, stats = demux_main(**demux_dict)
#             stats_df = pd.DataFrame.from_dict(stats, orient='index')
#             stats_df.loc[:, 'remaining reads/%'] = stats_df.loc[:, "remaining reads"] / stats_df.iloc[:, 1:4].sum(1)
#             stats_df.to_csv(demux_stats, index=True)
#         print("demultiplexing has been complete, indir will change to", demux_dir_samples)
#     if not isfile(opath):
#         write_manifest(opath=opath,
#                        ids=ids,
#                        r1_files=r1_files,
#                        r2_files=r2_files, )
#
#     os.makedirs(odir, exist_ok=True)
#     ## 跑命令
#     join_pairs(opath, int(nt), **join_params)
#     print("预处理完成,完成原始序列评估 与 joined, 去污染,去chimera,fix orientation")

# python3 /home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/preprocessing/join_pairs.py -p '/home/liaoth/data2/16s/shanghai_152/param.template' -nt 5
