import re
import string

from config.default_params import *
from pipelines.demux import main as demux_main
from pipelines.pipelines import selective_p, tax_assign_qiime2, g_tree
from utils import *


def preprocess(indir):
    r1_files = sorted(get_files(indir, r1_format))
    r2_files = sorted(get_files(indir, r2_format))
    ids = [re.findall(idpattern, os.path.basename(_))[0].strip(string.punctuation) for _ in r1_files]

    if demux_on:
        print("Demux is on, first start the demultiplexing process.")
        demux_dict['seqfile1'] = r1_files
        demux_dict['seqfile2'] = r2_files
        if not_overwrite_demux and os.path.isfile(demux_stats):
            print('demuxed files exist, pass it')
        else:
            r1_files, r2_files, ids, stats = demux_main(**demux_dict)
            stats_df = pd.DataFrame.from_dict(stats, orient='index')
            stats_df.loc[:, 'remaining reads/%'] = stats_df.loc[:, "remaining reads"] / stats_df.iloc[:, 1:4].sum(1)
            stats_df.to_csv(demux_stats, index=True)
        print("demultiplexing has been complete, indir will change to", demux_dir_samples)

    if not os.path.isfile(opath):
        write_manifest(opath=opath,
                       ids=ids,
                       r1_files=r1_files,
                       r2_files=r2_files, )
    print("manifest has been output to", opath)
    # 准备序列的输入
    raw_seq = import_data_with_manifest(opath)
    # 将序列信息导入qiime的环境,可另存为qza

    raw_seq_eval_vis = seq_eval(raw_seq,
                                n=n)
    raw_seq_eval_vis.save(os.path.join(odir,
                                       raw_seq_vis_path))

    join_params.update(qc_joined_params)
    join_params['raw_seq'] = raw_seq
    join_params['n'] = n
    joined_seq, joined_seq_eval_vis, \
    joined_qc_seq, joined_qc_eval_vis, joined_qc_stats = join_seqs(raw_seq,
                                                                   **join_params
                                                                   )
    os.makedirs(os.path.join(odir, 'preprocess'), exist_ok=True)
    raw_seq.save(os.path.join(odir, 'preprocess',
                              raw_seq_path))
    joined_qc_seq.save(os.path.join(odir, 'preprocess',
                                    joined_seq_path))
    joined_qc_eval_vis.save(os.path.join(odir, 'preprocess',
                                         joined_qc_seq_vis_path))
    joined_seq.save(os.path.join(odir, 'preprocess',
                                 joined_qc_seq_vis_path))
    joined_seq_eval_vis.save(os.path.join(odir, 'preprocess',
                                          joined_seq_vis_path))

    joined_qc_stats.view(pd.DataFrame).to_csv(os.path.join(odir, 'preprocess',
                                                           joined_qc_stats_tab_path), index=True)

    return raw_seq, joined_qc_seq


def run_pipelines(p, pipelines_args):
    pre_ = 'sOTU'
    tab, rep, tab_vis, seq_vis, stats_df = selective_p(p, pipelines_args)

    p_tab = tab.view(pd.DataFrame)
    name_dict = dict(zip(p_tab.columns, [pre_ + str(_) for _ in range(p_tab.shape[1])]))
    p_tab.columns = [name_dict[_] for _ in p_tab.columns]

    seq_f = glob(str(rep._archiver.data_dir) + '/*')
    if len(seq_f) != 1:
        raise Exception
    seq_f = seq_f[0]

    os.makedirs(os.path.join(odir, p), exist_ok=True)
    # output part
    tab.save(os.path.join(odir, p,
                          profiled_tab_path.format(prefix=p)))
    p_tab.to_csv(os.path.join(odir, p,
                              profiled_tab_path.format(prefix=p)),
                 sep='\t' if profiled_tab_path.endswith('.tab') else ',',
                 index=True)
    mv_seq(seq_f,
           opath=os.path.join(odir, p,
                              representative_sequence_path.format(prefix=p)),
           name_dict=name_dict)
    rep.save(os.path.join(odir, p,
                          representative_sequence_path.format(prefix=p)))
    tab_vis.save(os.path.join(odir, p,
                              profiled_tab_vis_path.format(prefix=p)))
    seq_vis.save(os.path.join(odir, p,
                              representative_sequence_vis_path.format(prefix=p)))
    stats_df.to_csv(os.path.join(odir, p,
                                 process_stats_path.format(prefix=p)))

    return tab, p_tab, rep


def after_otu(args):
    # assigning taxonomy, perform alpha,beta diversity
    tax_tab = tax_assign_qiime2(**args)
    rooted_tree = g_tree(**args)


if __name__ == '__main__':
    from utils import parse_param
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("pipelines", help="Which kinds of pipelines you want to perform. \
                             [deblur|dada2|none]",
                        type=str, choices=['deblur', 'dada2', 'none'])
    parser.add_argument("-p", "--parameter",
                        help="input file contains all parameters, template is place at the %s. This is a python script actually, so just use python code to write your code." % os.path.join(
                            os.path.dirname(os.path.abspath(__file__)), 'param.template'), default=os.path.join(os.path.dirname(__file__), 'param.template'))

    args = parser.parse_args()
    p = args.pipelines
    parameter_f = args.parameter

    os.makedirs(odir, exist_ok=True)
    parse_param(parameter_f, g=globals())
    ## 跑命令
    raw_seq, joined_qc_seq = preprocess(indir)
    print("预处理完成,完成原始序列评估 与 joined, 去污染,去chimera,fix orientation")
    if p == 'none':
        exit("Completed...")
    pipelines_args['deblur_input'] = joined_qc_seq
    pipelines_args['dada2_input'] = raw_seq
    tab, p_tab, rep = run_pipelines(p, pipelines_args)
    # after_otu_args['rep'] = rep
    # after_otu(after_otu_args)

    # python main.py deblur -p '/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/param.template2'
