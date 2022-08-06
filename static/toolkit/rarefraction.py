import pandas as pd
import plotly.graph_objects as go
from os.path import *
import sys

if __name__=='__main__':
    infile = sys.argv[1]
    #infile = 'analysis/rarefaction/rare.txt'
    if '/' not in infile:
        infile = './'+infile
    _df = pd.read_csv(infile,sep='\t',index_col=0,low_memory=False)
    fig = go.Figure()
    for s,col in _df.iteritems():
        fig.add_scatter(x=_df.index,y=list(col),mode='markers+lines',
                        name=s)
    fig.write_html(dirname(infile)+'/rare.html')
    fig.write_image(dirname(infile)+'/rare.png')