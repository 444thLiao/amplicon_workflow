from qiime2 import Artifact
import qiime2
from os.path import join,abspath,exists
from glob import glob
from Bio import SeqIO
import pandas as pd
def convert2otutab(infile,ofile):
    artifact = Artifact.load(infile)
    dada2_df = artifact.view(qiime2.Metadata).to_dataframe()
    dada2_df.index = dada2_df.index.astype(str)
    ofile_path = abspath(ofile)
    dada2_df.to_csv(ofile_path,sep='\t',index=1)
    return ofile_path

def convert2table(infile,ofile):
    artifact = Artifact.load(infile)
    stats = artifact.view(pd.DataFrame)
    stats.index = stats.index.astype(str)
    ofile_path = abspath(ofile)
    stats.to_csv(ofile_path,sep='\t',index=1)
    return ofile_path


def convert2seq(infile,ofile):
    artifact = Artifact.load(infile)
    path = artifact._archiver.data_dir

    path = glob(join(path,'*'))
    print(path)
    if path:
        path = path[0]
    assert exists(path)
    ofile_path = abspath(ofile)
    print(path)
    with open(ofile_path,'w') as f1:
        instream = SeqIO.parse(path,format='fasta')
        for record in instream:
            SeqIO.write(record,f1,format="fasta-2line")
    return ofile_path