



odir = '/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/results/amplicon_ATP5B/ppo'

cmd = f"qiime deblur denoise-16S --i-demultiplexed-seqs {odir}/q2_pipelines/joined_qc_seq.qza --o-representative-sequences {odir}/deblur_output/rep.qza --o-table {odir}/deblur_output/profiling.qza --o-stats {odir}/deblur_output/profiling_stats.qza --p-trim-length 220  --p-sample-stats --p-mean-error 0.005  --p-indel-prob 0.01  --p-indel-max 3  --p-min-reads 10  --p-min-size 2  --p-jobs-to-start 7  --p-hashed-feature-ids"

os.system("unzip joined_qc_seq.zip -d joined/")
f = glob('joined/*/data/MANIFEST')[0]
nrows = []
for r in open(f).readlines():
    rows = r.strip('\n').split(',')
    if len(rows)==3 and rows[0]!='sample-id':
        c = [rows[0].replace('_','-')]
        c+= [rows[1],rows[2]]
        nrows.append(','.join(c))
    else:
        nrows.append(','.join(rows))
with open(f,'w') as f1:
    f1.write('\n'.join(nrows))

os.system("cd joined/ && zip -f ../joined_qc_seq.zip {f.replace('joined/','')}")
os.system("mv joined_qc_seq.zip joined_qc_seq.qza")

cmd = f"qiime deblur denoise-16S --i-demultiplexed-seqs joined_qc_seq.qza --o-representative-sequences ../deblur_output/rep.qza --o-table ../deblur_output/profiling.qza --o-stats ../deblur_output/profiling_stats.qza --p-trim-length 500  --p-sample-stats --p-mean-error 0.005  --p-indel-prob 0.01  --p-indel-max 3  --p-min-reads 10  --p-min-size 2  --p-jobs-to-start 7  --p-hashed-feature-ids"

# c = {"ACPI":"#c3deae",
#      #"Lhcr":"#fdd4fe",
#  'Lhcr1':"#fdd4fe",
#  'Lhcr2':"#fdd4fe",
#  'Lhcr3':"#fdd4fe",
#  'Lhcr4':"#fdd4fe",
#  'Lhcr5':"#fdd4fe",
#  'Lhcr6':"#fdd4fe",
#  'Lhcr7':"#fdd4fe",     
#      "FCPI":"#dbe8f4",
#      "LHCI":"#fae3dd",
#      'RED':"#fdd4fe"
#      }
# from api_tools.itol_func import to_color_range
# text = to_color_range({t:t.split('_')[0].split('-')[0] for t in tips},
#                       c
#                       )