'''Summary for the script'''
import sys
import argparse

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--bam', required=True, help='BAM file to filter')
args = parser.parse_args()

keep = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
         '13', '14', '15', '16', '17', '18', '19', 'X', 'Y' ]
KEEP = { k:None for k in keep }

if __name__ == '__main__':

    with pysam.AlignmentFile(args.bam, 'rb') as samfile:
        with pysam.AlignmentFile('-', 'w', template=samfile) as outfile:
            for read in samfile.fetch():
                if read.reference_name in KEEP:
                    outfile.write(read)
