'''Read a compressed CGmap and wig file and produce new versions with only CG sites'''
import sys
import argparse

from collections import defaultdict
import gzip
import csv

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser()
parser.add_argument('--cg-map', required=True, help='CGmap File')
parser.add_argument('--wig', required=True, help='WIG file')
parser.add_argument('--sample', required=True, help='Sample name')
args = parser.parse_args()

if __name__ == '__main__':

    keep = defaultdict(dict)

    # read in CG file, write out CG file
    with gzip.open(args.cg_map, 'rt', newline='') as tsvfile, gzip.open(f'{args.sample}-cg-methylation.CGmap.gz', 'wt', newline='') as outfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        for row in reader:
            if row[3] == 'CG':
                keep[row[0]][row[2]] = None
                writer.writerow(row)

    # read in WIG, filter
    with gzip.open(args.wig, 'rt', newline='') as file, gzip.open(f'{args.sample}-cg-methylation.wig.gz', 'wt', newline='') as outfile:
        reader = csv.reader(file, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        chr = None

        for row in reader:
            if row[0] == 'track type=wiggle_0':
                writer.writerow(row)

            elif row[0][0] == 'v':
                writer.writerow(row)
                chr = row[0].replace('variableStep chrom=','')

            elif row[0] in keep[chr]:
                writer.writerow(row)

