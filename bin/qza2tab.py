#!/home-user/thliao/anaconda3/envs/qiime2-2020.8/bin/python

import sys
from os.path import join, dirname

sys.path.insert(0, dirname(dirname(__file__)))
from static.q2_function import convert2otutab, convert2seq

  
if __name__ == "__main__":
    infile = sys.argv[1]
    
    convert2otutab(infile,
                    infile.replace('.qza','.csv'))
        