#################################################################################
####  利用多进程+协程，进行高效的demutiplex data
####  20190609
####  by tianhua liao
#################################################################################
import sys
from os.path import dirname, join, basename

sys.path.insert(0, dirname(dirname(dirname(__file__))))
import gzip
import re
import os
import config.default_params as config
from Bio import SeqIO
from tqdm import tqdm
from glob import glob
from static.utils import data_parser
import multiprocessing as mp
from threading import Semaphore
import click
import pandas as pd

root_path = os.path.abspath(__file__)
r1_suffix_format = "_R1.fastq"
r2_suffix_format = "_R2.fastq"


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
    # parse a metadata and return barcode & sample info
    ids = []
    fb = []
    rb = []
    fp = []
    rp = []
    col2sto = dict(zip([id_col, fb_col, rb_col, fp_col, rp_col],
                       [ids, fb, rb, fp, rp]))
    metadf = data_parser(os.path.abspath(metadata), ft='csv')

    for col, sto in col2sto.items():
        if col and col in metadf.columns:
            column = metadf.loc[:, col]
            if col == id_col:
                p_sto(column, sto, is_id=True)
            else:
                p_sto(column, sto)
    return ids, fb, rb, fp, rp


def get_total_bc(bc1, bc2, bc_type):
    if bc_type == "concat":
        total_bc = bc1 + bc2
    elif bc_type == "same":
        assert bc1 == bc2
        total_bc = bc1
    return total_bc


def exec_fun(args):
    func, kwargs = args
    return func(**kwargs)


def process_pair(read1_data,
                 read2_data,
                 forward_primers,
                 reverse_primers,
                 length_bc=None,
                 attempt_read_orientation=False,
                 barcode_type="concat",
                 ):
    """
    默认双端都有barcode，默认拼起来使用（也可以选其一）
    Given pair record from fastq and assign sample name to its
    返回
    去除了primer,barcode的R1,R2的reads(已经校正了方向)，barcode，统计信息
    统计信息只有1/0，分别是 [方向相反，单个primer，无primer]
    :param read1_data:
    :param read2_data:
    :param forward_primers:
    :param reverse_primers:
    :param length_bc:
    :param attempt_read_orientation:
    :param barcode_type: concat / same
    :return:
    """
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
    stats = [0, 0, 0, 0]
    r1_seq = str(read1_data.seq)
    r2_seq = str(read2_data.seq)

    ## force re orientation
    for curr_primer in forward_primers:
        # iter forward primer first
        matched_str = curr_primer.search(r1_seq)
        if matched_str is not None:
            # search the primer, if we find this, it mean r1 is r1. orientation is right.
            read1 = read1_data
            read2 = read2_data
            bc1_end = fp_start = matched_str.start()
            # start of forward primer is the end of barcode.
            fp_end = matched_str.end()
            break
            # if we find it, just break it.
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

    if is_close_circle:
        # pair的reads上都找到了对应的barcode
        if length_bc is None:
            # 如果barcode 长度没给定。
            bc1_start = 0
            bc2_start = 0
        else:
            bc1_start = bc1_end - length_bc
            bc2_start = bc2_end - length_bc
        bc1 = str(read1[bc1_start:bc1_end].seq)
        remaining_read1 = read1[fp_end:]
        # cut primer and barcode before it
        bc2 = str(read2[bc2_start:bc2_end].seq)
        remaining_read2 = read2[rp_end:]


    else:
        if (bc1_end is None) and (bc2_end is None):
            # both primer are missing
            stats[2] = 1
        else:
            # one of them is missing.
            stats[1] = 1
        bc1 = str(read1_data[0:length_bc].seq)
        bc2 = str(read2_data[0:length_bc].seq)
        # force to cut a barcode (mostly it is wrong)
        remaining_read1 = read1_data[length_bc:]
        remaining_read2 = read2_data[length_bc:]

    total_bc = get_total_bc(bc1, bc2, barcode_type)
    if is_close_circle:
        if attempt_read_orientation:
            # 如果要fix方向,则正常输出
            output_read1 = remaining_read1
            output_read2 = remaining_read2
        else:
            # 如果不想fix方向,那么,按照原来的方向
            # 如果is_reversed,说明read1/2互换了,而read1肯定到read1_output2file,所以此时,应该是
            #       read1_output2file 输出到output_fastq2,才能保证按照原来方向
            output_read1 = remaining_read2 if is_reversed else remaining_read1
            output_read2 = remaining_read1 if is_reversed else remaining_read2
    else:
        # 如果不闭环,那么说明这个切(primer/barcode)都切不了..算了,告辞..不要这条read了
        output_read1 = None
        output_read2 = None

    return output_read1, output_read2, total_bc, stats


def process_single(read1_data,
                 forward_primers,
                 length_bc=None,):
    # init
    bc1_end = None
    read1 = None
    fp_end = None
    # stats = {"reversed_read":0,
    #          "single primer":0,
    #          "no primer":0,}
    stats = [0, 0, 0, 0]
    r1_seq = str(read1_data.seq)

    ## force re orientation
    for curr_primer in forward_primers:
        # iter forward primer first
        matched_str = curr_primer.search(r1_seq)
        if matched_str is not None:
            # search the primer, if we find this, it mean r1 is r1. orientation is right.
            read1 = read1_data
            bc1_end = fp_start = matched_str.start()
            # start of forward primer is the end of barcode.
            fp_end = matched_str.end()
            break
            # if we find it, just break it.
        if matched_str is not None:
            continue
        
    if length_bc is None:
        # 如果barcode 长度没给定。
        bc1_start = 0
    else:
        bc1_start = bc1_end - length_bc
    if bc1_start is not None and bc1_end is not None:
        bc1 = str(read1[bc1_start:bc1_end].seq)
        remaining_read1 = read1[fp_end:]
        # cut primer and barcode before it

        output_read1 = remaining_read1

        return output_read1, bc1, stats
    return None,None,stats

def demux_files(seqfile1,
                fp,
                ids,
                bc,
                length_bc,
                output_dir,
                seqfile2=None,
                rp=None,
                fileid='',
                attempt_read_orientation=False,
                barcode_type='concat',
                num_thread=0
                ):
    """
    对单个的序列文件,进行分库,并且输出到output_dir,保持名称不变(不进行压缩)
    因为对于该脚本而言,这只是个中间文件.
    分库的结果只会标记到每个read的head中,而不会生成sample为单位的序列文件
    Performance(complete process):
        For test/seq_data2:
            with asyncio : 3700-4000it/s  (78.59s user 10.72s system 173% cpu 51.609 total)
            without asyncio: 3700-4000it/s(78.26s user 11.53s system 174% cpu 51.596 total)
        For real data:
            without asyncio: 2000984it [11:15, 2828.67it/s]

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
    if num_thread == 0 or num_thread == -1:
        num_thread = mp.cpu_count()

    if len(set(bc)) != len(set(ids)):
        raise Exception("same barcode corresponding multiple ids.")
    bc2id = dict(zip(bc,
                     ids))

    name = fileid if fileid else "Demux"

    stats = {"mix-orientation reads": 0,
             "single primer": 0,
             "no primer": 0,
             "unknown barcode": 0,
             "output reads": 0}

    if len(set(fp)) != 1:
        print('WARNING!!!! Different primers, it may occur errors.')
    fp = list(set(fp))
    if rp is not None:
        rp = list(set(rp))

    if seqfile2:
        # PE sequencing files input
        if '.gz' in seqfile1:
            seqfile1_stream = gzip.open(seqfile1, 'rt')
            seqfile2_stream = gzip.open(seqfile2, 'rt')
        else:
            seqfile1_stream = open(seqfile1, 'r')
            seqfile2_stream = open(seqfile2, 'r')

        # init a params but within a generator.
        params = ((process_pair,
                   dict(read1_data=read1,
                        read2_data=read2,
                        forward_primers=fp,
                        reverse_primers=rp,
                        length_bc=length_bc,
                        attempt_read_orientation=attempt_read_orientation,
                        barcode_type=barcode_type))
                  for read1, read2 in zip(SeqIO.parse(seqfile1_stream, format='fastq'),
                                          SeqIO.parse(seqfile2_stream, format='fastq')))
        # use semaphore to control memory because of a bug of multiprocessing
        # if you don't fix it, it will take more and more memory instead of 0 increased memory of used-dump pattern.
        # see https://stackoverflow.com/questions/40922526/memory-usage-steadily-growing-for-multiprocessing-pool-imap-unordered
        semaphore = Semaphore(100000)

        def loader_f(semaphore):
            for p in params:
                semaphore.acquire()
                yield p

        # # async usage code.
        # async def write2file(read1, read2, file1, file2):
        #     # 异步处理 write2file
        #     if not os.path.isfile(file1):
        #         stream1 = open(file1, 'w')
        #         stream2 = open(file2, 'w')
        #     else:
        #         stream1 = open(file1, 'a')
        #         stream2 = open(file2, 'a')
        #     SeqIO.write(read1, stream1, format='fastq')
        #     SeqIO.write(read2, stream2, format='fastq')
        #     stream1.close()
        #     stream2.close()
        # # use generator instead of list to avoid memory explosion
        # async def async_work_main(imap_it,semaphore):
        #     nonlocal stats
        #     for remaining_read1, remaining_read2, total_bc, _s in tqdm(imap_it):
        #         stats["mix-orientation reads"] += _s[0]
        #         stats["single primer"] += _s[1]
        #         stats["no primer"] += _s[2]
        #         if remaining_read1 is not None:
        #             tmp1 = remaining_read1.description
        #             tmp2 = remaining_read2.description
        #             remaining_read1.id = remaining_read1.name = remaining_read1.description = ''
        #             remaining_read2.id = remaining_read2.name = remaining_read2.description = ''
        #
        #             sid = bc2id.get(total_bc, None)
        #             if sid is not None:
        #                 # try to get corresponding sample id
        #                 stats["output reads"] += 1
        #                 remaining_read1.id = sid + ' ' + tmp1
        #                 remaining_read2.id = sid + ' ' + tmp2
        #
        #                 file_f1 = os.path.join(output_dir, sid + r1_suffix_format)
        #                 file_f2 = os.path.join(output_dir, sid + r2_suffix_format)
        #                 # 调用异步调用
        #                 await write2file(remaining_read1,
        #                                  remaining_read2,
        #                                  file_f1,
        #                                  file_f2)
        #             else:
        #                 stats["unknown barcode"] += 1
        #                 continue
        #         semaphore.release()
        # # End.async usage code.
        # loop = asyncio.new_event_loop()
        # print('with asyncio')
        # 准备一个协程的事件管理loop
        with mp.Pool(processes=num_thread) as thread_pool:
            # 准备一个线程池
            loader = loader_f(semaphore)
            imap_it = thread_pool.imap(exec_fun, loader)

            for remaining_read1, remaining_read2, total_bc, _s in tqdm(imap_it):
                stats["mix-orientation reads"] += _s[0]
                stats["single primer"] += _s[1]
                stats["no primer"] += _s[2]
                if remaining_read1 is not None:
                    tmp1 = remaining_read1.description
                    tmp2 = remaining_read2.description
                    remaining_read1.id = remaining_read1.name = remaining_read1.description = ''
                    remaining_read2.id = remaining_read2.name = remaining_read2.description = ''

                    sid = bc2id.get(total_bc, None)
                    if sid is not None:
                        # try to get corresponding sample id
                        stats["output reads"] += 1
                        remaining_read1.id = sid + ' ' + tmp1
                        remaining_read2.id = sid + ' ' + tmp2

                        file_f1 = os.path.join(output_dir, sid + r1_suffix_format)
                        file_f2 = os.path.join(output_dir, sid + r2_suffix_format)
                        if not os.path.isfile(file_f1):
                            stream1 = open(file_f1, 'w')
                            stream2 = open(file_f2, 'w')
                        else:
                            stream1 = open(file_f1, 'a')
                            stream2 = open(file_f2, 'a')
                        SeqIO.write(remaining_read1, stream1, format='fastq')
                        SeqIO.write(remaining_read2, stream2, format='fastq')
                        stream1.close()
                        stream2.close()

                        # asyncio.run(async_work(remaining_read1,
                        #                                    remaining_read2,
                        #                                    file_f1,
                        #                                    file_f2))
                    else:
                        stats["unknown barcode"] += 1
                semaphore.release()
        #     loop.run_until_complete(async_work_main(imap_it,semaphore))
        # loop.close()
    else:
        # todo: 处理single end的测序数据（这么少了。。。干脆不写了吧。。
        # SE sequencing files input
        if '.gz' in seqfile1:
            seqfile1_stream = gzip.open(seqfile1, 'rt')
        else:
            seqfile1_stream = seqfile1
        # init a params but within a generator.
        params = ((process_single,
                   dict(read1_data=read1,
                        forward_primers=fp,
                        length_bc=length_bc))
                  for read1 in SeqIO.parse(seqfile1_stream, format='fastq'))  
        semaphore = Semaphore(100000)

        def loader_f(semaphore):
            for p in params:
                semaphore.acquire()
                yield p
                
        with mp.Pool(processes=num_thread) as thread_pool:
            # 准备一个线程池
            loader = loader_f(semaphore)
            imap_it = thread_pool.imap(exec_fun, loader)

            for remaining_read1, total_bc, _s in tqdm(imap_it):
                stats["mix-orientation reads"] += _s[0]
                stats["single primer"] += _s[1]
                stats["no primer"] += _s[2]
                if remaining_read1 is not None:
                    tmp1 = remaining_read1.description
                    remaining_read1.id = remaining_read1.name = remaining_read1.description = ''

                    sid = bc2id.get(total_bc, None)
                    if sid is not None:
                        # try to get corresponding sample id
                        stats["output reads"] += 1
                        remaining_read1.id = sid + ' ' + tmp1

                        file_f1 = os.path.join(output_dir, sid +'.fastq')
                        if not os.path.isfile(file_f1):
                            stream1 = open(file_f1, 'w')
                        else:
                            stream1 = open(file_f1, 'a')
                        SeqIO.write(remaining_read1, stream1, format='fastq')
                        stream1.close()

                        # asyncio.run(async_work(remaining_read1,
                        #                                    remaining_read2,
                        #                                    file_f1,
                        #                                    file_f2))
                    else:
                        stats["unknown barcode"] += 1
                semaphore.release()                
    return name, stats


def get_paths(files_input):
    if type(files_input) == str:
        return list(sorted(glob(files_input)))
    else:
        return files_input


def generate_df():
    pass

@click.command()
@click.option("-m", '--metadata', "metadata", help='path of metadata')
@click.option("-o", '--output-dir', "output_dir", help='The directory you want to output to')
@click.option("-r1", '--r1-files', "r1_files", type=str,
              help='forward reads you want to demuplexed. If you pass wildcard, it should use quote to quote it.')
@click.option("-r2", '--r2-files', "r2_files", type=str,
              help='Revesed reads you want to demuplexed. If you pass wildcard, it should use quote to quote it.')
@click.option("-f", '--overwrite-demux', "overwrite_demux", is_flag=True, help='Overwrite the output or not')
@click.option("-r", '--fix-orientation/--no-fix', "attempt_read_orientation",
              default=True, help='fix orientation or not?  [default is True]')
@click.option("-p", '--num_thread', "num_thread", default=0, show_default=True,
              help='The number of thread you want to use. 0 mean all threads')
def main(metadata,
         id_col=config.id_col,
         fb_col=config.fb_col,
         rb_col=config.rb_col,
         fp_col=config.fp_col,
         rp_col=config.rp_col,
         r1_files=(),
         r2_files=(),
         output_dir=None,
         overwrite_demux=True,
         attempt_read_orientation=False,
         num_thread=0
         ):
    if overwrite_demux and os.path.isdir(output_dir):
        os.system("rm -fr %s" % output_dir)
    os.makedirs(output_dir, exist_ok=True)

    ids, fb, rb, fp, rp = parse_metadata(metadata,
                                         id_col=id_col,
                                         fb_col=fb_col,
                                         rb_col=rb_col,
                                         rp_col=rp_col,
                                         fp_col=fp_col)
    file_stats = {}

    seqfile1 = get_paths(r1_files)
    seqfile2 = get_paths(r2_files)
    if not seqfile1:
        raise Exception("No input sequencing files detect.")
    barcode_type = "concat"
    # todo: deal with this param. input? or get from metadata?
    if len(set([len(_.pattern) for _ in fb])) == 1:
        # 根据barcode的长度进行set，有且仅有1个
        length_bc = len(fb[0].pattern)
        if rb:
            # 如果存在reverse barcode
            bc = [p1.pattern + p2.pattern
                  for p1, p2 in zip(fb, rb)]
        else:
            bc = [p1.pattern for p1 in fb]
        for seq1, seq2 in zip(seqfile1,
                              seqfile2):
            # 大多情况下，只有少数的seqfile，所以我选择把多进程+协程放在demux_files里
            filebasename = str(os.path.basename(seq1)).split('.')[0]
            name, stats = demux_files(seqfile1=seq1,
                                      seqfile2=seq2,
                                      fp=fp,
                                      rp=rp,
                                      ids=ids,
                                      bc=bc,
                                      length_bc=length_bc,
                                      output_dir=output_dir,
                                      fileid=filebasename,
                                      attempt_read_orientation=attempt_read_orientation,
                                      barcode_type=barcode_type,
                                      num_thread=num_thread)
            file_stats[name] = stats
    else:
        # 多种barcode长度,忽略吧...
        ## 还不知道该咋办...不应该的说
        pass

    stats_df = pd.DataFrame.from_dict(file_stats, orient='index')
    stats_df.index.name = "ori file name"
    stats_df.to_csv(os.path.join(output_dir, "demux_stats.csv"), index=1)
    # generate a tab for 16s pipelines
    R1_files = glob(join(output_dir,
                         "*" + r1_suffix_format))
    R2_files = glob(join(output_dir,
                         "*" + r2_suffix_format))
    sid_r1 = set([basename(_).split(r1_suffix_format)[0]
                  for _ in R1_files])
    sid_r2 = set([basename(_).split(r1_suffix_format)[0]
                  for _ in R2_files])
    assert len(sid_r1) == len(sid_r2)
    from config import input_template_path
    columns = open(input_template_path).read().strip('\n').split('\t')
    metadata_df = pd.DataFrame(columns=columns, )
    for sid in sid_r1:
        sid = str(sid)
        metadata_df = metadata_df.append(pd.DataFrame(columns=metadata_df.columns, index=[0], data=[[sid,
                                                                                                     join(output_dir, sid + r1_suffix_format),
                                                                                                     join(output_dir, sid + r2_suffix_format)]]))
    metadata_df.to_csv(os.path.join(output_dir,
                                    "demux_data.tsv"),
                       index=0,
                       sep='\t')

    return stats_df


if __name__ == '__main__':
    main()

    # python3 pp/demux.py -m /home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/metadata.tab -o /home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq2_demux/ -r1 "/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq_data2/test_seq*_1.fastq.gz" -r2 "/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq_data2/test_seq*_2.fastq.gz" -f -r

    # from os.path import dirname
    #
    # # path = '/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/pp/demux.py'
    # metadata = os.path.join(dirname(dirname(root_path)), 'test', 'metadata.tab')
    # id_col = 'SampleID'
    # fb_col = 'Forward_Barcode'
    # rb_col = 'Reverse_Barcode'
    # fp_col = 'Forward_Primer'
    # rp_col = 'Reverse_Primer'
    # seqfile1, seqfile2 = sorted(glob(os.path.join(dirname(dirname(root_path)),
    #                                               'test/seq_data2/test_seq*_1.fastq.gz')
    #                                  )), \
    #                      sorted(glob('/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq_data2/test_seq*_2.fastq.gz'))
    # # the order of glob output is random, be careful !!!!!!!!!!!!!!!!!!!!!!!!1
    # stats = main(metadata=metadata,
    #              id_col=id_col,
    #              rb_col=rb_col,
    #              fb_col=fb_col,
    #              rp_col=rp_col,
    #              fp_col=fp_col,
    #              seqfile1=seqfile1,
    #              seqfile2=seqfile2,
    #              output_dir='/home/liaoth/data2/16s/qiime2_learn/gpz_16s_pipelines/test/seq2_demux/',
    #              attempt_read_orientation=True,
    #              overwrite_demux=True)
    # for f, s in sorted(stats.items(), key=lambda x: x[0]):
    #     print('{:>12}  {:>12}  {:>12} {:>12}'.format(*s.keys()))
    #     print('{:>12}  {:>12}  {:>12} {:>12}'.format(*s.values()))

    ## performance analysis
    ## ~15000 read per second per file

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
