import os
import sys

sys.path.insert(0,
                os.path.dirname(os.path.dirname(__file__)))

import argparse
from toolkit.summarize_tax_report import summarize_tax
import pandas as pd


USEARCH = os.popen(f"which usearch").read().strip()
RDP_DB = '/home-user/thliao/db/16S/RDP/rdp_16s_v18.fa'
draw_stack_dis = f"{os.path.dirname(__file__)}/draw_stack_bar_plus.py"
draw_rarefraction  = f"{os.path.dirname(__file__)}/rarefraction.py "

draw_PD = "/home/liaoth/project/16s_pipelines/microbiome_utils/Visualization/draw_PD.py"

def regular_analysis(OTU_table, rep_fa, outputdir, draw_pd=False,simple=True):
    """
    OTU_TABLE is normalize and filtered.
    :param OTU_table:
    :param rep_fa:
    :param outputdir:
    :return:
    """
    if not os.path.isdir(outputdir):
        os.system('mkdir -p %s;' % outputdir)
    first_row = open(OTU_table).read().strip().split('\n')[0]
    if not first_row.startswith('#OTU'):
        c = open(OTU_table).read().split('\n')
        OTU_table = OTU_table.rsplit('.')[0]+'_renamed.csv'
        fr = '\t'.join(['#OTU'] + c[0].split('\t')[1:])
        contents = [fr] + c[1:]
        with open(OTU_table,'w') as f1:
            f1.write('\n'.join(contents))
    df = pd.read_csv(OTU_table,sep='\t',index_col=0)
    if df.shape[0]>df.shape[1]:
        df.to_csv(OTU_table,sep='\t',index=1,index_label=df.index.name)
    else:
        pass
    
    os.system('mkdir -p %s;' % (outputdir + '/beta_diversity'))
    os.system('mkdir -p %s;' % (outputdir + '/alpha_diversity'))
    os.system('mkdir -p %s;' % (outputdir + '/rarefaction'))
    os.system('mkdir -p %s;' % (outputdir + '/taxonomy_report'))
    
    output_tree = os.path.join(outputdir, 'otus.tree')

    sintax = os.path.join(outputdir, 'sintax.txt')
    if not simple:
        os.system('%s -cluster_agg %s -treeout %s' % (USEARCH,
                                                  rep_fa,
                                                  output_tree))
    os.system("%s -alpha_div %s -output '%s'" % (USEARCH,
                                                 OTU_table,
                                                 os.path.join(outputdir + '/alpha_diversity', 'alpha.txt')))
    os.system("%s -beta_div %s -tree %s -filename_prefix '%s'" % (USEARCH,
                                                                  OTU_table,
                                                                  output_tree,
                                                                  outputdir + '/beta_diversity/'))
    if not simple:
        os.system("%s -sintax %s -db %s -strand both -tabbedout %s -sintax_cutoff 0.8" % (USEARCH,
                                                                                      rep_fa,
                                                                                      RDP_DB,
                                                                                      sintax))
    os.system("%s -alpha_div_rare %s -output %s/rare.txt" % (USEARCH,
                                                             OTU_table,
                                                             outputdir + '/rarefaction'))
    if not simple:
        tax_summary_out_file = os.path.join(outputdir + '/taxonomy_report')

        summarize_tax(sintax,
                        OTU_table,
                        tax_summary_out_file,
                        level='gpf')

        os.system('python3 %s %s' %
                    (draw_stack_dis, os.path.join(outputdir + '/taxonomy_report')))
    os.system('python3 %s %s' %
                (draw_rarefraction, os.path.join(outputdir + '/rarefaction/rare.txt')))
    
    if draw_pd:
        os.system("python3 %s -t %s -i %s -o %s -M shannon,observed_otus,faith_pd" % (draw_PD,
                                                                                    os.path.join(
                                                                                        outputdir, 'otus.tree'),
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
                        help="if you want to draw a rarefaction curve(which will take long time. But it will have a progress bar to display its.)")
    parser.add_argument('-s', dest='simple', action='store_true',
                        help="no annotation. only run alpha and rarefaction")    
    args = parser.parse_args()

    otu_tab = os.path.abspath(args.otutab)
    fasta = os.path.abspath(args.fasta)
    odir = os.path.abspath(args.odir)
    draw_rc = args.rarefaction
    simple = args.simple
    regular_analysis(otu_tab, fasta, odir, draw_pd=draw_rc,simple=simple)

    # python ~/software/16s_workflow/static/toolkit/routine_analysis.py -i ./profiling.csv -fa rep.fa -o ./routine
