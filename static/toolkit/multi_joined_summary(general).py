import argparse
import os
import re
from collections import defaultdict
from glob import glob

import pandas as pd
from tqdm import tqdm

"""
190429
written by: TianHua Liao

This script is tend to used after multiplt_joined.
There have four purpose in this script.
    1.remove unjoined fastq file incase mixed in after pipelines.
    2.rename the complex sample name into simpify one.
    3.remove the tags from fastq_screen tools in case affect after pipelines.
    4.calculate the reads maintain_ratio after joined process.

"""

default_fastq_suffix = '*.gz'
sample_name_pattern = "(.*).[Rf]"


def count_fq_reads(fq):
    num_reads = int(os.popen("zgrep -c '^+$' %s" % fq).read())
    return num_reads


def count_fa_reads(fq):
    num_reads = int(os.popen("zgrep -c '^>' %s" % fq).read())
    return num_reads


def joined_summary(joined_result_dir, raw_dir, pattern, output_csv):
    """
    :param joined_result_dir: the dir which is multiple_join_paired_ends.py outputed without last '/'.
    :param output_csv: the summary file path you want to storge.
    :return: None
    """
    joined_dict = {}
    for joined_seq in tqdm(glob(os.path.join(joined_result_dir, default_fastq_suffix))):
        sample_name = re.findall(pattern, os.path.basename(joined_seq))
        num_reads = count_fq_reads(joined_seq)
        joined_dict[sample_name] = num_reads
    raw_dict = defaultdict(list)
    for raw_seq in tqdm(glob(os.path.join(raw_dir, default_fastq_suffix))):
        sample_name = re.findall(pattern, os.path.basename(raw_seq))
        num_reads = count_fq_reads(raw_seq)
        raw_dict[sample_name].append(num_reads)
    raw_dict = {k: sum(v) / len(v) for k, v in raw_dict.items()}
    data = pd.DataFrame(dict(raw=raw_dict,
                             joined=joined_dict,
                             ))
    data.loc[:, "joined ratio(%)"] = (data['raw'] / data['joined']) * 100
    data.to_csv(output_csv, index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='jrd', type=str, required=True,
                        help="joined seq dir")
    parser.add_argument('-i2', '--input2', dest='rawd', type=str, required=True,
                        help="raw seq dir")
    parser.add_argument('-o', '--output', dest='oc', type=str, required=True,
                        help="csv output file path")
    parser.add_argument('-pattern', dest='pattern', type=str, required=False, default='(.*)',
                        help="name pattern will be stodge in summary and renamed with. [Default] is the ori name.")
    args = parser.parse_args()
    joined_dir = os.path.abspath(args.jrd)
    raw_dir = os.path.abspath(args.rawd)
    pattern = args.pattern
    if not pattern:
        pattern = sample_name_pattern
    output_csv = os.path.abspath(os.path.abspath(args.oc))

    if not os.path.isfile(output_csv):
        joined_summary(joined_dir, raw_dir, pattern=sample_name_pattern, output_csv=output_csv)
        print('Output file already exist, please make sure and delete it.')
    else:
        print("output file is existed")
