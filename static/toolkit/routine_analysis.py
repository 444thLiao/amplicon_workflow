import os
import sys

sys.path.insert(0,
                os.path.dirname(os.path.dirname(__file__)))
from toolkit.summarize_tax_report import summarize_tax
import argparse

USEARCH = '/home/liaoth/tools/usearch'
RDP_DB = '~/data2/rdp_16s_v16.fa'
draw_PD = "/home/liaoth/project/16s_pipelines/microbiome_utils/Visualization/draw_PD.py"
draw_stack_dis = "/home/liaoth/data2/project/16s_pipelines/microbiome_utils/Visualization/draw_stack_bar_plus.py"


def regular_analysis(OTU_table, rep_fa, outputdir, draw_pd=False):
    """
    OTU_TABLE is normalize and filtered.
    :param OTU_table:
    :param rep_fa:
    :param outputdir:
    :return:
    """
    if not os.path.isdir(outputdir):
        os.system('mkdir -p %s;' % outputdir)

    os.system('mkdir -p %s;' % (outputdir + '/beta_diversity'))
    os.system('mkdir -p %s;' % (outputdir + '/alpha_diversity'))
    os.system('mkdir -p %s;' % (outputdir + '/rarefaction'))
    os.system('mkdir -p %s;' % (outputdir + '/taxonomy_report'))
    print('%s -cluster_agg %s -treeout %s' % (USEARCH,
                                              rep_fa,
                                              os.path.join(outputdir, 'otus.tree')))
    os.system('%s -cluster_agg %s -treeout %s' % (USEARCH,
                                                  rep_fa,
                                                  os.path.join(outputdir, 'otus.tree')))
    os.system("%s -alpha_div %s -output '%s'" % (USEARCH,
                                                 OTU_table,
                                                 os.path.join(outputdir + '/alpha_diversity', 'alpha.txt')))
    os.system("%s -beta_div %s -tree %s -filename_prefix '%s'" % (USEARCH,
                                                                  OTU_table,
                                                                  os.path.join(outputdir, 'otus.tree'),
                                                                  outputdir + '/beta_diversity/'))
    os.system("%s -sintax %s -db %s -strand both -tabbedout %s -sintax_cutoff 0.8" % (USEARCH,
                                                                                      rep_fa,
                                                                                      RDP_DB,
                                                                                      os.path.join(outputdir, 'sintax.txt')))
    os.system("%s -alpha_div_rare %s -output %s/rare.txt" % (USEARCH,
                                                             OTU_table,
                                                             outputdir + '/rarefaction'))

    tax_summary_out_file = os.path.join(outputdir + '/taxonomy_report')

    summarize_tax(os.path.join(outputdir, 'sintax.txt'),
                  OTU_table,
                  tax_summary_out_file,
                  level='gpf')

    os.system('python3 %s %s' % (draw_stack_dis, os.path.join(outputdir + '/taxonomy_report')))

    if draw_pd:
        os.system("python3 %s -t %s -i %s -o %s -M shannon,observed_otus,faith_pd" % (draw_PD,
                                                                                      os.path.join(outputdir, 'otus.tree'),
                                                                                      OTU_table, outputdir + '/rarefaction'))

    os.system('cp %s %s' % (OTU_table, outputdir + '/'))
    os.system('cp %s %s' % (rep_fa, outputdir + '/otu_rep_seq.fasta'))


# regular_analysis('/home/liaoth/data2/16s/shandong/16s_pipelines/v_analysis_dechimera/duplicate_sample/redo_way/otu_norm_filtered_40k_s2.tab',
#                  '/home/liaoth/data2/16s/shandong/16s_pipelines/v_analysis_dechimera/duplicate_sample/redo_way/otus_rep.fa',
#                  '/home/liaoth/data2/16s/shandong/16s_pipelines/v_analysis_dechimera/duplicate_sample/redo_way/routine_analysis/')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='otutab', type=str, required=True,
                        help="OTU table path. Make sure you have filtered.")
    parser.add_argument('-fa', dest='fasta', type=str, required=True,
                        help="OTU represent sequence fasta file path.")
    parser.add_argument('-o', '--output', dest='odir', type=str, required=True,
                        help="output Dir path")
    parser.add_argument('-rc', dest='rarefaction', action='store_true',
                        help="if you want to draw a rarefaction curve(which will take long time.But it will have a progress bar to display its.)")
    args = parser.parse_args()

    otu_tab = os.path.abspath(args.otutab)
    fasta = os.path.abspath(args.fasta)
    odir = os.path.abspath(args.odir)
    draw_rc = args.rarefaction

    regular_analysis(otu_tab, fasta, odir, draw_pd=draw_rc)
