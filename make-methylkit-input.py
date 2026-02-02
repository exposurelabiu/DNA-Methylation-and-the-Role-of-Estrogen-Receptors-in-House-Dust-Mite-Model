'''Takes a CG map file for a sample and a list of CG sites to use and
create a Methylkit-compatible file for those sites.'''
import sys
import argparse

from collections import defaultdict
import gzip
import csv

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser()
parser.add_argument('--cg-map', required=True, help='CGmap File')
parser.add_argument('--master', required=True, help='Master CG file')
args = parser.parse_args()

def read_master_list():
    '''Read in the list of CG sites to use'''
    
    sites = defaultdict(dict)

    with open(args.master, newline='') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            sites[row[0]][row[1]] = None
    return sites

if __name__ == '__main__':

    sites = read_master_list()

    print('chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT')

    # read in CG file, write out CG file
    with gzip.open(args.cg_map, 'rt', newline='') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            if row[2] in sites[row[0]]:
                
                strand = 'F' if row[1] == 'C' else 'R'
                c = 100 * (int(row[6]) / int(row[7]))
                t = 100 - c

                print(f'{row[0]}.{row[2]}\t{row[0]}\t{row[2]}\t{strand}\t{row[7]}\t{c:.2f}\t{t:.2f}')

