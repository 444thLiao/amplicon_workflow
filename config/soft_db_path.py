from os.path import dirname, join

project_root_path = dirname(dirname(__file__))
qiime_env = "qiime2-2019.4"
vsearch_pth = '/home/liaoth/tools/vsearch-2.6.2/bin/vsearch'
rdp_gold_pth = "/home/liaoth/data2/project/16s_pipelines/microbiome_utils/vsearch_pipeliens/rdp_gold.fa"
map_pl_pth = str(join(project_root_path, 'static', "map.pl"))

trimmomatic_jar = '/home/liaoth/tools/Trimmomatic-0.36/trimmomatic-0.36.jar'
fq_screen = "/home/liaoth/tools/fastq_screen_v0.11.1/fastq_screen"
usearch_pth = "/home/liaoth/tools/usearch"
