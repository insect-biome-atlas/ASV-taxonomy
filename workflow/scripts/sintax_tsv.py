#!/usr/bin/env python3
##### Program description #######
#
# Title: Subsetting Insect Biome ASVs 
#
# Author(s): Lokeshwaran Manoharan
#
#
#
# Description: The ASVs that are supposed to be of particular length (418 ± nx3 and between 403 and 418), possibly does not contain any stop codon in the right reading frame. 
#  
# List of subroutines: 
#
#
#
# Overall procedure: using functions and dictionaries in python
#
# Usage: create_q2_md_lanes.py <sample_info_file> <location_of_trimmed_files> <output_ss_file>
##################################

import re
import sys
import copy
import argparse

usage = '''This program takes the output of SINTAX and creates readable taxonomy file'''

parser = argparse.ArgumentParser(description=usage)


parser.add_argument(
    '-i', '--infile',
    dest='infile',
    metavar='INFILE',
    type=argparse.FileType('r'),
    help='SINTAX output',
    required=True
    )

parser.add_argument(
    '-o', '--outfile',
    dest='outfile',
    metavar='OUTFILE',
    type=argparse.FileType('w'),
    help='Final Taxonomy in tsv format',
    default=sys.stdout
    )

args = parser.parse_args()

def complete_list(some_list, target_len):
    return some_list[:target_len] + ['unclassified.' + some_list[-1]]*(target_len - len(some_list))

p0 = re.compile('>')
p1 = re.compile(' ')
p2 = re.compile('\,')
p3 = re.compile('\t')
p4 = re.compile('\;')

print('ASV', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'BOLD_bin', sep = '\t', file = args.outfile)

for line in args.infile:
    line = line.rstrip('\n')
    tmp_list = re.split(p3, line)
    if tmp_list[3] != '':
        annot = re.sub('.\:', '', tmp_list[3])
        tmp_list1 = re.split(p2, annot)
        if len(tmp_list1) != 8:
            tmp_list1 = complete_list(tmp_list1, 8)
            print(tmp_list[0], '\t'.join(tmp_list1), sep = '\t', file = args.outfile)
        else:
            print(tmp_list[0], '\t'.join(tmp_list1), sep = '\t', file = args.outfile)
    else:
        print(tmp_list[0], '\t'.join(['unclassified'] * 8), sep = '\t', file = args.outfile)

args.outfile.close()
args.infile.close()