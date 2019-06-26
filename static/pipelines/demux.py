import gzip
import re

from Bio import SeqIO
from tqdm import tqdm
from glob import glob
from config.default_params import *
from utils import data_parser, assign_work_pool
import itertools
from subprocess import check_call
def p_sto(col, stodge, is_id=False):
    iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
             'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
             'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}
    for v in col:
        if not is_id:
            v = re.compile(''.join([iupac[_.upper()] for _ in v]))
        stodge.append(v)


def rev_p_sto(col, stodge, is_id=False):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'R': 'Y', 'Y': 'R',
                  'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                  'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}

    iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
             'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
             'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}
    for v in col:
        if not is_id:
            v = re.compile(''.join([iupac[complement[_.upper()]] for _ in v[::-1]]))
        stodge.append(v)


def parse_metadata(metadata,
                   id_col='',
                   fb_col='',
                   rb_col='',
                   fp_col='',
                   rp_col='',
                   ):
    ids = []
    fb = []
    rb = []
    fp = []
    rp = []
    col2sto = dict(zip([id_col, fb_col, rb_col, fp_col, rp_col],
                       [ids, fb, rb, fp, rp]))
    metadf = data_parser(metadata, ft='csv')

    for col, sto in col2sto.items():
        if col and col in metadf.columns:
            column = metadf.loc[:, col]
            if col == id_col:
                p_sto(column, sto, is_id=True)
            else:
                p_sto(column, sto)
    return ids, fb, rb, fp, rp


def process_pair(read1_data,
                 read2_data,
                 forward_primers=None,
                 reverse_primers=None,
                 length_bc=None,
                 attempt_read_orientation=False,
                 ):
    # init
    is_reversed = False
    is_close_circle = False
    bc1_end = None
    bc2_end = None
    read1 = read2 = None
    fp_end = rp_end = None
    # stats = {"reversed_read":0,
    #          "single primer":0,
    #          "no primer":0,}
    stats = [0, 0, 0]
    r1_seq = str(read1_data.seq)
    r2_seq = str(read2_data.seq)

    ## force re orientation
    for curr_primer in forward_primers:
        matched_str = curr_primer.search(r1_seq)
        if matched_str is not None:
            read1 = read1_data
            read2 = read2_data
            bc1_end = fp_start = matched_str.start()
            fp_end = matched_str.end()
            break
        matched_str = curr_primer.search(r2_seq)
        if matched_str is not None:
            read1 = read2_data
            read2 = read1_data
            is_reversed = True
            bc1_end = fp_start = matched_str.start()
            fp_end = matched_str.end()
            break
    # Check reverse primers if forward primers not found
    for curr_primer in reverse_primers:
        if is_reversed:
            matched_str = curr_primer.search(r1_seq)  # 左闭右开
            if matched_str is not None:
                read1 = read2_data
                read2 = read1_data
                bc2_end = rp_start = matched_str.start()
                rp_end = matched_str.end()

                is_close_circle = True
                stats[0] = 1
                break
        elif bc1_end is not None and not is_reversed:
            matched_str = curr_primer.search(r2_seq)
            if matched_str is not None:
                read1 = read1_data
                read2 = read2_data
                bc2_end = rp_start = matched_str.start()
                rp_end = matched_str.end()
                is_close_circle = True
                break
    # if read1_data.id == 'ERR1750995.619':
    #     import pdb;pdb.set_trace()
    if is_close_circle:
        # 双向都找到了对应的barcode, 或者闭合且单向的也找到了.
        if length_bc is None:
            bc1_start = 0
            bc2_start = 0
        else:
            bc1_start = bc1_end - length_bc
            bc2_start = bc2_end - length_bc
        bc1 = read1[bc1_start:bc1_end]
        remaining_read1 = read1[fp_end:]  # cut primer and barcode before it
        bc2 = read2[bc2_start:bc2_end]
        remaining_read2 = read2[rp_end:]
        total_bc = bc1 + bc2

    else:
        # import pdb;pdb.set_trace()
        # print("Wrong seq.",
        #       "id1:%s,id2:%s" % (read1_data.id, read2_data.id),
        #       "bc1_end: %s" % str(bc1_end),
        #       "bc2_end: %s" % str(bc2_end),
        #       "read1_data[:40]: %s" % str(read1_data[:40].seq),
        #       "read2_data[:40]: %s" % str(read2_data[:40].seq),
        #       )
        if (bc1_end is None) and (bc2_end is None):
            stats[2] = 1
        else:
            stats[1] = 1
        total_bc = read1_data[0:length_bc] + read2_data[0:length_bc]
        remaining_read1 = read1_data[length_bc:]
        remaining_read2 = read2_data[length_bc:]

    total_bc = str(total_bc.seq)
    if is_close_circle and attempt_read_orientation:
        # 如果闭环,且要fix方向,则正常输出
        output_read1 = remaining_read1
        output_read2 = remaining_read2
    elif not attempt_read_orientation and is_close_circle:
        # 如果闭环,且不想fix方向,那么,按照原来的方向
        # 如果is_reversed,说明read1/2互换了,而read1肯定到read1_output2file,所以此时,应该是
        #       read1_output2file 输出到output_fastq2,才能保证按照原来方向
        output_read1 = remaining_read2 if is_reversed else remaining_read1
        output_read2 = remaining_read1 if is_reversed else remaining_read2
    else:
        # 如果不闭环,那么说明这个切(primer/barcode)都切不了..算了,告辞..不要这条read了
        output_read1 = None
        output_read2 = None

    return output_read1, output_read2, total_bc, stats


def process_single():
    pass


def demux_files(seqfile1,
                seqfile2,
                fp,
                rp,
                ids,
                bc,
                length_bc,
                output_dir,
                not_overwrite_demux=False,
                attempt_read_orientation=False
                ):
    """
    对单个的序列文件,进行分库,并且输出到output_dir,保持名称不变(不进行压缩)
    因为对于该脚本而言,这只是个中间文件.
    分库的结果只会标记到每个read的head中,而不会生成sample为单位的序列文件
    :param str seqfile1:
    :param str seqfile2:
    :param list fp:
    :param list rp:
    :param list ids:
    :param list bc:
    :param int length_bc:
    :param str output_dir:
    :param bool not_overwrite_demux:
    :param bool attempt_read_orientation:
    :return:
    """
    bc2id = dict(zip(bc, ids))
    if len(set(bc)) != len(set(ids)):
        raise Exception("same barcode corresponding multiple ids.")
    name = str(os.path.basename(seqfile1)).partition('_1')[0]
    dir_name = os.path.join(output_dir,name)
    os.makedirs(dir_name,exist_ok=True)
    stats = {"reversed reads": 0,
             "single primer": 0,
             "no primer": 0,
             "missing barcode": 0,
             "remaining reads": 0}
    if len(set(fp)) != 1:
        print('WARNING!!!! Different primers, it may occur errors.')
    fp = list(set(fp))
    rp = list(set(rp))
    if seqfile2:
        if '.gz' in seqfile1:
            seqfile1_stream = gzip.open(seqfile1, 'rt')
            seqfile2_stream = gzip.open(seqfile2, 'rt')
        else:
            seqfile1_stream = open(seqfile1, 'r')
            seqfile2_stream = open(seqfile2, 'r')

        for read1, read2 in tqdm(zip(SeqIO.parse(seqfile1_stream, format='fastq'),
                                     SeqIO.parse(seqfile2_stream, format='fastq'))):
            remaining_read1, remaining_read2, total_bc, _s = process_pair(read1_data=read1,
                                                                          read2_data=read2,
                                                                          forward_primers=fp,
                                                                          reverse_primers=rp,
                                                                          length_bc=length_bc,
                                                                          attempt_read_orientation=attempt_read_orientation)
            stats["reversed reads"] += _s[0]
            stats["single primer"] += _s[1]
            stats["no primer"] += _s[2]
            if remaining_read1 is not None:
                tmp1 = remaining_read1.description
                tmp2 = remaining_read2.description
                remaining_read1.id = remaining_read1.name = remaining_read1.description = ''
                remaining_read2.id = remaining_read2.name = remaining_read2.description = ''

                sid = bc2id.get(total_bc, None)
                if sid is not None:
                    stats["remaining reads"] += 1
                    remaining_read1.id = sid + ' ' + tmp1
                    remaining_read2.id = sid + ' ' + tmp2
                    id_f1 = os.path.join(dir_name,sid+'_1.fastq')
                    id_f2 = os.path.join(dir_name, sid + '_2.fastq')
                    if not os.path.isfile(id_f1):
                        stream1 = open(id_f1, 'w')
                        stream2 = open(id_f2, 'w')
                    else:
                        stream1 = open(id_f1, 'a')
                        stream2 = open(id_f2, 'a')
                    SeqIO.write(remaining_read1, stream1, format='fastq')
                    SeqIO.write(remaining_read2, stream2, format='fastq')
                else:
                    stats["missing barcode"] += 1
            # else:
            #     print(_s, read1.id, total_bc, bc2id.get(total_bc, None))
    else:
        # 不想写....
        process_single()

    return name, stats


def _split(f1, f2, output_dir_pre, output_dir_samples):
    nf1 = os.path.join(output_dir_pre, os.path.basename(f1).strip('.gz'))
    nf2 = os.path.join(output_dir_pre, os.path.basename(f2).strip('.gz'))

    seqfile1_stream = open(nf1, 'r')
    seqfile2_stream = open(nf2, 'r')
    for read1, read2 in tqdm(zip(SeqIO.parse(seqfile1_stream, format='fastq'),
                                 SeqIO.parse(seqfile2_stream, format='fastq'))):
        sid = read1.id
        id_f1 = os.path.join(output_dir_samples, '%s_1.fastq' % sid)
        id_f2 = os.path.join(output_dir_samples, '%s_2.fastq' % sid)
        if not os.path.isfile(id_f1):
            stream1 = open(id_f1, 'w')
            stream2 = open(id_f2, 'w')
        else:
            stream1 = open(id_f1, 'a')
            stream2 = open(id_f2, 'a')
        SeqIO.write(read1, stream1, format='fastq')
        SeqIO.write(read2, stream2, format='fastq')


def cal(args):
    func, args = args
    return func(*args)

#
# def split_into_files(seqfile1, seqfile2, output_dir_pre, output_dir_samples, num_thread):
#     """
#     将多个demux_files输出的结果,根据每个read前面的id,生成到以sample为单位的序列文件中
#     :param list seqfile1:
#     :param list seqfile2:
#     :param str output_dir_pre:
#     :param str output_dir_samples:
#     :return:
#     """
#     all_args = [(_split, (f1, f2, output_dir_pre, output_dir_samples)) for f1, f2 in zip(seqfile1, seqfile2)]
#
#     assign_work_pool(cal, all_args, num_thread=num_thread)
#
#     f1_files = sorted([_ for _ in glob(os.path.join(output_dir_samples,
#                                                     '*_1.fastq')) if _ not in seqfile1])
#     f2_files = sorted([_ for _ in glob(os.path.join(output_dir_samples,
#                                                     '*_2.fastq')) if _ not in seqfile2])
#     ids = [str(os.path.basename(_)).split('_1')[0] for _ in f1_files]
#     return f1_files, f2_files, ids
def merge_files(cmd_template,sid,num):
    check_call(cmd_template.format(sid=sid,
                                   num=num))
    return 'complete'

def main(metadata,
         id_col='',
         fb_col='',
         rb_col='',
         fp_col='',
         rp_col='',
         seqfile1=(),
         seqfile2=(),
         output_dir_pre=None,
         output_dir_samples=None,
         not_overwrite_demux=False,
         attempt_read_orientation=False,
         num_thread=0
         ):
    ids, fb, rb, fp, rp = parse_metadata(metadata,
                                         id_col=id_col,
                                         fb_col=fb_col,
                                         rb_col=rb_col,
                                         rp_col=rp_col,
                                         fp_col=fp_col)
    file_stats = {}
    os.makedirs(output_dir_pre, exist_ok=True)
    os.makedirs(output_dir_samples, exist_ok=True)

    if type(seqfile1) == str:
        seqfile1 = [seqfile1]
    else:
        seqfile1 = sorted(seqfile1)
    if type(seqfile2) == str:
        seqfile2 = [seqfile2]
    else:
        seqfile2 = sorted(seqfile2)

    if len(set([len(_.pattern) for _ in fb])) == 1:
        length_bc = len(fb[0].pattern)
        if rb:
            bc = [p1.pattern + p2.pattern for p1, p2 in zip(fb, rb)]
        else:
            bc = [p1.pattern for p1 in fb]
        all_args = [(demux_files, (seq1,
                                   seq2,
                                   fp,
                                   rp,
                                   ids,
                                   bc,
                                   length_bc,
                                   output_dir_pre,
                                   not_overwrite_demux,
                                   attempt_read_orientation)) for seq1, seq2 in zip(seqfile1, seqfile2)]

        results = assign_work_pool(cal, all_args, num_thread=num_thread)
        for name, stats in results:
            file_stats[name] = stats

        print("Start mergeing reads from multiple sampels fastq......")
        all_sample_files = glob(os.path.join(output_dir_pre,'*','*_1.fastq'))
        unique_samples_name = set([os.path.basename(_).replace('_1.fastq','') for _ in all_sample_files])
        # concat start.
        cmd_template = "cat %s > %s" % (os.path.join(output_dir_pre, '*', '{sid}_{num}.fastq'),
                                        os.path.join(output_dir_samples, "{sid}_{num}.fastq"))
        all_args = [(merge_files, (cmd_template,
                                   sid,
                                   num,)) for sid, num in itertools.product(unique_samples_name, ['1','2'])]
        results = assign_work_pool(cal, all_args, num_thread=num_thread)
    else:
        # 多种barcode长度,忽略吧...
        ## 还不知道该咋办...不应该的说
        pass

    return f1_files, f2_files, ids, file_stats


if __name__ == '__main__':
    # path = os.path.abspath(__file__)


    path = '/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/pp/demux.py'

    metadata = os.path.join(os.path.dirname(path), '..', 'test', 'metadata.tab')
    id_col = 'SampleID'
    fb_col = 'Forward_Barcode'
    rb_col = 'Reverse_Barcode'
    fp_col = 'Forward_Primer'
    rp_col = 'Reverse_Primer'
    seqfile1, seqfile2 = sorted(glob('/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq_data2/test_seq*_1.fastq.gz')), \
                         sorted(glob('/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq_data2/test_seq*_2.fastq.gz'))
    # the order of glob output is random, be careful !!!!!!!!!!!!!!!!!!!!!!!!1
    f1_files, f2_files, ids, stats = main(metadata=metadata,
                                          id_col=id_col,
                                          rb_col=rb_col,
                                          fb_col=fb_col,
                                          rp_col=rp_col,
                                          fp_col=fp_col,
                                          seqfile1=seqfile1,
                                          seqfile2=seqfile2,
                                          output_dir_pre='/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq2_demux/',
                                          output_dir_samples='/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq2_samples/',
                                          attempt_read_orientation=True)
    for f, s in sorted(stats.items(), key=lambda x: x[0]):
        print('{:>12}  {:>12}  {:>12} {:>12}'.format(*s.keys()))
        print('{:>12}  {:>12}  {:>12} {:>12}'.format(*s.values()))

    ## performance analysis
    ## ~5000 read per second per file

#    import cProfile
#
#    cProfile.run("""main(metadata=metadata,
# id_col=id_col,
# rb_col=rb_col,
# fb_col=fb_col,
# rp_col=rp_col,
# fp_col=fp_col,
# seqfile1=seqfile1,
# seqfile2=seqfile2,
# output_dir='/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq2_demux/',
# attempt_read_orientation=True)
#    """, sort=2)
#
#    cProfile.run("""for i1,i2 in tqdm(zip(SeqIO.parse(gzip.open(seqfile1, 'rt'), format='fastq'),
#                                     SeqIO.parse(gzip.open(seqfile2, 'rt'), format='fastq'))): pass
#
#    """, sort=2)
