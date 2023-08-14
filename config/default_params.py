# 输入
import os
from .soft_db_path import qiime_env

qiime2_p = "source activate %s; qiime " % qiime_env

############################################################
# 输出文件名集合
raw_seq_path = 'raw_data'
joined_seq_path = 'joined_seq'

raw_seq_vis_path = 'unjoined_seq_eval'
joined_seq_vis_path = 'joined_seq_eval'
joined_qc_seq_vis_path = 'joined_qc_seq_eval'
joined_qc_stats_tab_path = 'joined_qc_stats.csv'

profiled_tab_path = 'profiling.tab'
profiled_tab_vis_path = 'profiling'
representative_sequence_path = 'rep.fa'
representative_sequence_vis_path = 'rep'
process_stats_path = 'profiling_stats.csv'

root_tree_path = 'rep_rooted_tree.tab'
tax_tab = 'rep_sintax.tab'

############################################################
# demux部分
demux_on = False
not_overwrite_demux = True
id_col = 'SampleID'
fb_col = 'Forward_Barcode'
rb_col = 'Reverse_Barcode'
fp_col = 'Forward_Primer'
rp_col = 'Reverse_Primer'

############################################################
trimmomatic_thread = 5
fq_screen_thread = 5
# ################
# # demux部分
# demux_on = False
# not_overwrite_demux = True
# demux_dict = dict(
#     metadata=os.path.join(_p_dir, 'test', 'metadata.tab'),
#     id_col='SampleID',
#     fb_col='Forward_Barcode',
#     rb_col='Reverse_Barcode',
#     fp_col='Forward_Primer',
#     rp_col='Reverse_Primer',
#     attempt_read_orientation=True,
#     output_dir_pre=demux_dir_pre,
#     output_dir_samples=demux_dir_samples,
#     num_thread=0  # all threads
# )

# 序列评估 可视化部分参数
n = 10000
# join 部分参数
join_params = dict(
    truncqual=None,
    minlen=1,
    maxns=None,
    allowmergestagger=True,
    minovlen=10,
    maxdiffs=10,
    minmergelen=None,
    maxmergelen=None,
    maxee=None,
    qmin=0,
    qminout=0,
    qmax=41,
    qmaxout=41,
    #qascii=33
)
# join 序列评估
qc_joined_params = dict(
    min_quality=4,  # default
    quality_window=3,  # default
    min_length_fraction=0.75,  # default
    max_ambiguous=0,  # default
)

vesearch_args = dict(trunclen=240,
                     cluster_ratio=0.97
                     )

deblur_custom_args = dict(trunclen=240,
                          refdb='',
                          jobs_to_start=7,
                     )


deblur_args = dict(
    # deblur
    trim_length=220,
    sample_stats=True,
    mean_error=0.005,  # default
    indel_prob=0.01,  # default
    indel_max=3,  # default
    min_reads=10,  # default
    min_size=2,  # default
    jobs_to_start=7,
    hashed_feature_ids=True,  # default
)

# pipeliens args
dada2_args = dict(
    # dada2
    trunc_len_f=230,
    trunc_len_r=230,
    n_threads=0,  # all threads
    trunc_q=2,  # default
    n_reads_learn=1000000,  # default
    max_ee_f=2.0,  # default
    max_ee_r=2.0,  # default
)

# usearch part
maxee = 1.0

# Not implement yet
## after profiling
# after_otu_args = dict(
#     # classify
#     classifier_pth='/home/liaoth/data2/gg-13-8-99-nb-classifier.qza',
#     n_jobs=-2,
#     confidence=-1,
#     read_orientation=None,
#     # mafft tree
#     n_threads=0,
#     mask_max_gap_frequency=1.0,
#     mask_min_conservation=0.4,
# )
