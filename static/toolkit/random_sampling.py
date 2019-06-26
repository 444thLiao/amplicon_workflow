import gzip
import os
from glob import glob

import numpy as np
from Bio import SeqIO
from tqdm import tqdm


def random_sampling(r1, r2, odir, n=5000):
    if '*' in r1:
        r1_files = sorted(glob(r1))
        r2_files = sorted(glob(r2))
    else:
        r1_files = [r1]
        r2_files = [r2]
    count = 1
    for f1, f2 in zip(r1_files, r2_files):
        print(f1)
        rids = np.random.choice(np.arange(1000000), size=n, replace=False)
        stodge = []
        if '.gz' in f1:
            # Bio.SeqIO.parse is faster than skbio.io.read
            seq1 = SeqIO.parse(gzip.open(f1, 'rt'), format='fastq')
            seq2 = SeqIO.parse(gzip.open(f2, 'rt'), format='fastq')
        else:
            seq1 = SeqIO.parse(f1, format='fastq')
            seq2 = SeqIO.parse(f2, format='fastq')
        for idx, reads in tqdm(enumerate(zip(seq1, seq2))):
            read1, read2 = reads
            if idx in rids:
                stodge.append((read1, read2))
            if len(stodge) >= n:
                break

        with open(os.path.join(odir,
                               'test_seq%s_1.fastq' % count), 'w') as f1:
            SeqIO.write((x[0] for x in stodge), f1, 'fastq')
        with open(os.path.join(odir,
                               'test_seq%s_2.fastq' % count), 'w') as f1:
            SeqIO.write((x[1] for x in stodge), f1, 'fastq')
        count += 1


if __name__ == '__main__':
    r1 = '/home/data_public/Guangdong/*_1.fastq.gz'
    r2 = '/home/data_public/Guangdong/*_2.fastq.gz'

    odir = '/home/liaoth/temp_/fortestdata/'
    random_sampling(r1, r2, odir, 10000)
