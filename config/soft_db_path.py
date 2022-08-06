from os.path import *
from os import popen
import os,sys
    
def get_project_root():
    from os.path import dirname
    project_root = dirname(dirname(__file__))
    return project_root
def env_exe(name):

    bin_dir = dirname(sys.executable)
    f = join(bin_dir,name)
    if exists(f):
        return f
    f = popen(f'which {name} 2> /dev/null').read().strip('\n')
    return f

project_root = get_project_root()
############################################################
# exe path
############################################################

project_root_path = dirname(dirname(__file__))
qiime_env = "qiime2-2020.8"
vsearch_pth = env_exe('vsearch') if env_exe('vsearch') else  '/home-user/thliao/software/vsearch-2.15.1-linux-x86_64/bin/vsearch'
rdp_gold_pth = "/home-user/thliao/db/16S/RDP/rdp_16s_v18.fa"
map_pl_pth = str(join(project_root_path, 'static', "map.pl"))

trimmomatic_jar = '/home-user/thliao/software/Trimmomatic-0.39/trimmomatic-0.39.jar'
fq_screen = "/home-user/thliao/software/FastQ-Screen-0.14.1/fastq_screen"
usearch_pth = "/home-user/software/usearch/usearch"

fastqc_path = env_exe('fastqc') if env_exe('fastqc') else "/home-user/thliao/anaconda3/envs/wgs/opt/fastqc-0.11.8/fastqc"
multiqc_path = env_exe('multiqc') if env_exe('multiqc') else "/home-user/thliao/anaconda3/bin/multiqc" 
