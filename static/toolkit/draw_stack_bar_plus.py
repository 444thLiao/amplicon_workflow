"""
Stacked bar plot for summarizing the taxonomic classification of 16S analysis

"""

from sys import settrace
import plotly
import random
import os
from plotly import graph_objs as go
import seaborn as sns
import pandas as pd
from glob import glob


def read_data(fp):
    try:
        first = pd.read_csv(fp)
        if len(first.columns) == 1 or len(first.columns) == 0:
            first = pd.read_table(fp, index_col=0)
        if len(first.columns) == 0:
            print('Unresolved delimiter.')
    except:
        first = pd.read_excel(fp)

    return first


def parse_data(tax, tax_otu_sub, colors):
    bucket = []
    if tax != 'family':
        judge = False
    else:
        judge = True
    
    tax_otu_sub = tax_otu_sub.reindex(columns=tax_otu_sub.sum(0).sort_values().index)
    tax_otu_sub = tax_otu_sub.sort_index()
    
    for col_name, col in tax_otu_sub.iteritems():
        # loop it and raw each OTU(taxed)
        c = list(tax_otu_sub.columns).index(col_name)
        if c >=len(colors):
            c = c-c//len(colors)*len(colors)
        else:
            c = colors[c]
        try:
            bucket.append(go.Bar(
            x=col.index,
            y=col/tax_otu_sub.sum(1),  # Normalization with each sample.
            visible=judge,
            name=col_name,
            # name will display at legend.
            
            marker=dict(color=c,
                        line=dict(width=1, color='#FFFFFF')
                        # use white line to separate the bar, it can use it to make it more clear.
                        )
        ))
        except:
            import pdb;pdb.set_trace()
    return bucket

def generate_html(result_dfs, ofile):
    tax2dis = {}
    for tax_name, (tax_otutab_norm, tax_otutab_sum) in zip('pfg', result_dfs):
        tax_otu = tax_otutab_sum
        tax_otu.loc[:, 'sum_all'] = tax_otu.sum(1)
        # sum all samples reads in each OTU(taxed).
        tax_otu.sort_values('sum_all', ascending=False, inplace=True)
        # Use total reads to sort in order to make big OTU been drawn above.
        tax_otu_sub = tax_otu.iloc[:, :-1]
        # deleted the sum_all columns.

        # processing lot of colors
        data = []
        colors = sns.husl_palette(len(tax_otu.columns))
        random.shuffle(colors)

        colors = colors.as_hex()

        tax2dis[tax_name] = parse_data(tax_name, tax_otu_sub, colors)

    family_shows = [True] * len(tax2dis['f']) + [False] * \
        len(tax2dis['g']) + [False] * len(tax2dis['g'])
    genus_shows = [False] * len(tax2dis['f']) + [True] * \
        len(tax2dis['g']) + [False] * len(tax2dis['g'])
    phylum_shows = [False] * len(tax2dis['f']) + [False] * \
        len(tax2dis['g']) + [True] * len(tax2dis['g'])

    updatemenus = list([
        dict(type="buttons",
             active=1,
             buttons=list([
                 dict(label='phylum',
                      method='restyle',
                      args=['visible', phylum_shows]),
                 dict(label='family',
                      method='restyle',
                      args=['visible', family_shows]),
                 dict(label='genus',
                      method='restyle',
                      args=['visible', genus_shows]),
                 # dict(label = 'OTU',
                 #  method = 'restyle',
                 #  args = ['visible', OTU_shows])
             ]),
             )
    ])

    layout = go.Layout(
        title='Taxonomic distribution Plus',
        barmode='stack',
        # height=2000,
        updatemenus=updatemenus

    )
    fig = go.Figure(data=tax2dis['f']+tax2dis['g']+tax2dis['p'], layout=layout)
    fig.write_html(ofile, include_plotlyjs='cdn')


def generate_html(tax_otutabs):
    for tax_otutab in tax_otutabs:
        tax_names = os.path.basename(tax_otutab).split('_')[-1].split('.')[0]

        tax_otu = read_data(tax_otutab)
        #tax_otu = tax_otu.loc[tax_otu.sum(0)!=0,:]
        # tax_otu_sub = tax_otu_sub.loc[:, sorted(list(tax_otu_sub.columns))]
        #new_cols = [i for i in sorted(list(tax_otu.columns)) if '-N' in i] + [i for i in sorted(list(tax_otu.columns)) if '-T' in i]
        #tax_otu = tax_otu.loc[:,new_cols]
        tax_otu.loc[:, 'sum_all'] = tax_otu.sum(1)
        # sum all samples reads in each OTU(taxed).
        tax_otu.sort_values('sum_all', inplace=True)
        # Use total reads to sort in order to make big OTU been drawn above.
        tax_otu_sub = tax_otu.iloc[:, :-1]
        # deleted the sum_all columns.
        #tax_otu_sub = tax_otu_sub.loc[:, sort_way(tax_otu_sub.columns)]
        # sort the columns according to the header.
        tax_otus_index = tax_otu_sub.index

        # processing lot of colors
        colors = sns.husl_palette(len(tax_otus_index))
        random.shuffle(colors)
        # shuffle the colors to make it divergence
        colors = colors.as_hex()

        if tax_names == 'family':
            
            family_dis = parse_data(tax_names, tax_otu_sub, colors)
        elif tax_names == 'genus':
            genus_dis = parse_data(tax_names, tax_otu_sub, colors)
        elif tax_names == 'phylum':
            phylum_dis = parse_data(tax_names,  tax_otu_sub, colors)
        elif tax_names == 'OTU':
            OTU_dis = parse_data(tax_names,  tax_otu_sub, colors)
        else:
            print(tax_names)
            # exit()
    family_shows = [True] * len(family_dis) + [False] * \
        len(genus_dis) + [False] * len(phylum_dis)
    genus_shows = [False] * len(family_dis) + [True] * \
        len(genus_dis) + [False] * len(phylum_dis)
    phylum_shows = [False] * len(family_dis) + [False] * \
        len(genus_dis) + [True] * len(phylum_dis)
    OTU_shows = [False] * len(family_dis) + [False] * \
        len(genus_dis) + [False] * len(phylum_dis)
    updatemenus = list([
        dict(type="buttons",
             active=1,
             buttons=list([
                 dict(label='phylum',
                      method='restyle',
                      args=['visible', phylum_shows]),
                 dict(label='family',
                      method='restyle',
                      args=['visible', family_shows]),
                 dict(label='genus',
                      method='restyle',
                      args=['visible', genus_shows]),
                 # dict(label = 'OTU',
                 #  method = 'restyle',
                 #  args = ['visible', OTU_shows])
             ]),
             )
    ])

    layout = go.Layout(
        title='Taxonomic distribution Plus',
        barmode='stack',
        # height=2000,
        updatemenus=updatemenus

    )

    fig = go.Figure(data=family_dis+genus_dis+phylum_dis, layout=layout)
    fig.write_html(input_dir+'/taxonomic_distribution.html')
    # plotly.offline.plot(fig, filename=)


if __name__ == '__main__':
    import sys
    # Config
    metadata = None
    input_dir = sys.argv[1]

    tax_otutabs = [glob(f"{input_dir}/*_{t}.txt")[0]
                   for t in ['family', 'genus', 'phylum']]
    tax_levels = [os.path.basename(i).split(
        '_')[-1].split('.')[0] for i in tax_otutabs]
    generate_html(tax_otutabs)