'''Parse CGmap files for all samples and combined them into a super file'''
import sys
import argparse

import csv
import gzip

from collections import defaultdict

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser()
parser.add_argument('samples', nargs='+', help='Samples to process')
args = parser.parse_args()


if __name__ == '__main__':

    site_map = defaultdict(dict)
    sample_map = defaultdict(lambda: defaultdict(dict))

    chromosomes = {}
    
    for sample in args.samples:
        with gzip.open(f'{sample}-cg-methylation.CGmap.gz', 'rt', newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            for row in reader:
                c = row[0]
                i = int(row[2])

                if c not in chromosomes:
                    chromosomes[c] = None

                # if we haven't seen this site before, remember it
                if i not in site_map[c]:
                    site_map[c][i] = [row[0], row[1], row[2]]

                # store the count and methylation for this sample
                sample_map[sample][c][i] = [row[6], row[7]]

    print('# ' + ' '.join(args.samples))

    for c in chromosomes:        
        for i in sorted(site_map[c].keys()):
            answer = site_map[c][i]

            for s in args.samples:
                answer += sample_map[s][c][i] if i in sample_map[s][c] else ['0','0']
            print('\t'.join(answer))

