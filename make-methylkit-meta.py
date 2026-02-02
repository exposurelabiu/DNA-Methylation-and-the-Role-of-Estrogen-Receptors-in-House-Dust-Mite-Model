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
parser.add_argument('--field', required=True, help='Field to split on')
parser.add_argument('--control', required=True, help='Control condition')
parser.add_argument('--experiment', required=True, help='Experimental condition')
parser.add_argument('--filter-field', nargs='*', default=[], 
                    help='Field to filter on')
parser.add_argument('--filter-value', nargs='*', default=[],
                    help='Only include samples with this value')
parser.add_argument('--meta', required=True, help='Metafile to read')
args = parser.parse_args()

def read_meta():
    '''Read the meta file'''

    samples = []
    with open(args.meta, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            samples.append(row)
    return samples

def passes_filter(s, filter_count):    
    for i in range(filter_count):
        if s[args.filter_field[i]] != args.filter_value[i]:
            return False
    return True

if __name__ == '__main__':
    
    samples = read_meta()
    
    print('sample\tsex\tgenotype\ttreatment\texperiment\tfile')

    filter_count = len(args.filter_field)


    # parse samples and print entries
    for s in samples:

        if not passes_filter(s, filter_count):
            continue

        # are we filtering
#        if args.filter_field:
#            if s[args.filter_field] != args.filter_value:
#                continue
        
        exp = None
        if s[args.field] == args.control:
            exp = 0
        elif s[args.field] == args.experiment:
            exp = 1
        else:
            continue

        print(f's{s['sample']}\t{s['sex']}\t{s['genotype']}\t{s['treatment']}\t{exp}\t{s['sample']}-mk-input.tsv')
