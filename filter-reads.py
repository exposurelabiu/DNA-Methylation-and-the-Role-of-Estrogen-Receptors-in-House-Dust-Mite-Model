'''Summary for the script'''
import sys
import argparse

import gzip
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser()
parser.add_argument('sample', help='Sample name')
args = parser.parse_args()

if __name__ == '__main__':

    good_insert = { 'CGATGG' : None, 'TGGCGA' : None, 'CGACGG' : None, 'CGGCGA' : None }

    records_f = FastqGeneralIterator(gzip.open(f'{args.sample}-clean-1.fq.gz', 'rt'))
    records_r = FastqGeneralIterator(gzip.open(f'{args.sample}-clean-2.fq.gz', 'rt'))

    f_handle = gzip.open(f'{args.sample}-filtered-1.fq.gz', 'wt')
    r_handle = gzip.open(f'{args.sample}-filtered-2.fq.gz', 'wt')

#    records_f = SeqIO.parse(gzip.open(f'{args.sample}-clean-1.fq.gz', 'rt'), 'fastq')
#    records_r = SeqIO.parse(gzip.open(f'{args.sample}-clean-2.fq.gz', 'rt'), 'fastq')
    
    for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in zip(records_f,records_r):

        f_key = f_seq[0:3]

        key = f_key + r_seq[0:3]

        if key not in good_insert:
            continue

        # trim trim trim
        # reading from synthesized strand - trim 5' end of insert
        if f_key == 'CGA':
            f_seq = f_seq[2:]
            f_q = f_q[2:]
            r_seq = r_seq[:-2]
            r_q = r_q[:-2]
        # reading from original strand - trim 3' end of insert
        else:
            f_seq = f_seq[:-2]
            f_q = f_q[:-2]
            r_seq = r_seq[2:]
            r_q = r_q[2:]

        f_handle.write(f'@{f_id}\n{f_seq}\n+\n{f_q}\n')
        r_handle.write(f'@{r_id}\n{r_seq}\n+\n{r_q}\n')
