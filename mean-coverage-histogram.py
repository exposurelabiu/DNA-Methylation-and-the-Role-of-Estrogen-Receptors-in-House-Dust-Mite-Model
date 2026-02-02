'''Parse CGmap files and generate histograms of site coverage depth'''
import sys
import argparse

import csv

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

parser = argparse.ArgumentParser()
parser.add_argument('--cg-file', required=True, help='File with CG data')
parser.add_argument('--output', required=True, help='Output PDF file')
args = parser.parse_args()

if __name__ == '__main__':

    with PdfPages(args.output) as pdf:

        means = []
        cvs = []

        all_count = 0
        count = 0

        with open(args.cg_file, newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')

            for row in reader:
                if row[0][0] == '#':
                    continue

                all_count += 1

                depth = []

                for i in range(4, len(row), 2):
                    depth.append(min(100,int(row[i])))

                mean = np.average(depth)
                std = np.std(depth)
                cv = std / mean if mean > 0 else 0
                    
                means.append(mean)
                
                if mean >= 20 and mean <= 90:
                    print(f'{row[0]}\t{row[2]}')
                    count += 1
                    cvs.append(cv)

        plt.figure()
        plt.clf()

        plt.hist(means, bins=50)
        plt.xlabel('CG Site Coverage Depth Mean')
        plt.ylabel('GG Site Count')
        plt.title(f'CG Coverage Mean for all samples (n={all_count})')
        
        pdf.savefig(bbox_inches="tight")
        plt.close()

        plt.figure()
        plt.clf()

        plt.hist(means, bins=50)
        plt.yscale('log')
        plt.xlabel('CG Site Coverage Depth Mean')
        plt.ylabel('Log GG Site Count')
        plt.title(f'CG Coverage Mean for all samples (n={all_count})')
        
        pdf.savefig(bbox_inches="tight")
        plt.close()
        
        plt.figure()
        plt.clf()
        
        plt.hist(cvs, bins=40)
        plt.xlabel('CG Site Coverage Coefficient of Variation')
        plt.ylabel('GG Site Count')
        plt.title(f'CG Coverage CV for 20 <= mean <= 90 (n={count})')
        
        pdf.savefig(bbox_inches="tight")
        plt.close()

        plt.figure()
        plt.clf()
        
        plt.hist(cvs, bins=40)
        plt.yscale('log')
        plt.xlabel('CG Site Coverage Coefficient of Variation')
        plt.ylabel('Log GG Site Count')
        plt.title(f'CG Coverage CV for 20 <= mean <= 90 (n={count})')
        
        pdf.savefig(bbox_inches="tight")
        plt.close()
