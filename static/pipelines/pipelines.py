import qiime2
from qiime2.plugins.feature_classifier.methods import classify_sklearn
from qiime2.plugins.feature_table.visualizers import summarize as ft_vis
from qiime2.plugins.feature_table.visualizers import tabulate_seqs as ft_tabseq
from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree

from pipelines.dada2_pipelines import dada2_pipelines
from pipelines.deblur_pipelines import deblur_pipelines
from utils import *


def selective_p(p, args):
    if p not in ['dada2',
                 'deblur',
                 'otu']:

        return
    elif p == 'dada2':
        print("Perform dada2 pipelines.....")
        tab, rep, stats = dada2_pipelines(**args)
        dada2_tab_vis = ft_vis(tab)[0]
        dada2_seq_vis = ft_tabseq(rep)[0]
        dada2_stats_df = stats.view(qiime2.Metadata).to_dataframe()

        for c in list(dada2_stats_df.columns[1:]):
            now_loc = dada2_stats_df.columns.get_loc(c)
            new_name = '{} (%)'.format(c)
            dada2_stats_df.insert(now_loc + 1, new_name, dada2_stats_df.loc[:, c] / dada2_stats_df.iloc[:, 0] * 100)
        return tab, rep, dada2_tab_vis, dada2_seq_vis, dada2_stats_df

    elif p == 'deblur':
        print("Perform deblur pipelines.....")
        tab, rep, stats = deblur_pipelines(**args)
        deblur_tab_vis = ft_vis(tab)[0]
        deblur_seq_vis = ft_tabseq(rep)[0]
        deblur_stats_df = stats.view(pd.DataFrame)

        for c in [_ for _ in deblur_stats_df.columns[1:] if _.startswith('reads-')]:
            now_loc = deblur_stats_df.columns.get_loc(c)
            new_name = '{} (%)'.format(c)
            deblur_stats_df.insert(now_loc + 1, new_name, deblur_stats_df.loc[:, c] / deblur_stats_df.iloc[:, 0] * 100)
        return tab, rep, deblur_tab_vis, deblur_seq_vis, deblur_stats_df

    else:
        pass


def tax_assign_qiime2(rep=None,
                      classifier_pth=None,
                      n_jobs=-2,
                      confidence=-1,
                      read_orientation=None,  # default
                      **kwargs
                      ):
    classifier = Artifact.load(classifier_pth)
    tax_tab = classify_sklearn(rep=rep,
                               classifier=classifier,
                               n_jobs=n_jobs,
                               confidence=confidence,
                               read_orientation=read_orientation,  # default
                               )[0]

    return tax_tab


def g_tree(rep=None,
           n_threads=0,
           mask_max_gap_frequency=1.0,
           mask_min_conservation=0.4,
           ):
    aligned, masked_aligned, tree, r_tree = align_to_tree_mafft_fasttree(sequences=rep,
                                                                         n_threads=n_threads,  # all
                                                                         mask_max_gap_frequency=mask_max_gap_frequency,
                                                                         mask_min_conservation=mask_min_conservation, )

    return r_tree
