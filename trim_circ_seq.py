#!/usr/bin/env python3
#
# Moritz Blumer | 2025-05-10
#
# This script trims trailing sequence from a linearized circular DNA assembly 
# by identifying and removing the region that matches the start of the sequence.

# import packages
import sys
import os
import re
import gzip

# CLI & help
global input_path, seed_length, output_path

if len(sys.argv) < 4:
    print(
        '\n   python trim_circ_seq.py <input_path> <seed_length> <output_path>\n\n\
        <input_path>          str  path to input FASTA\n\
        <seed_length>         int  number of bp at start of sequence to query\n\
        <output_path>         str  path to output FASTA ("-" for STDOUT)\n\
        ',
    file=sys.stderr,
    )
    sys.exit()

# fetch arguments    
_, input_path, seed_length, output_path = sys.argv

# convert optional paths to NONE bool if "NONE"
try:
    seed_length = int(seed_length)
except:
    print('\n[ERROR] <seed_length> must be an integer.', file=sys.stderr)
    sys.exit()

# check if input file was specified and exists
if not os.path.isfile(input_path):
    print('\n[ERROR] Input file does not exist.', file=sys.stderr)
    sys.exit()

# read input FASTA
read_func = gzip.open if input_path.endswith('.gz') else open
fasta = read_func(input_path, 'r')

# read header and sequence
header = fasta.readline().strip()
seq = fasta.read().replace('\n', '')

# infer search string
search_string = seq[0:seed_length]
match_lst = [m.start() for m in re.finditer(search_string, seq)]

if len(match_lst) == 1:
    print(f'[ERROR] no other matches found, consider reducing <seed_length>.')
    sys.exit(flush=True, file=sys.stderr)

elif len(match_lst) > 2:
    print(f'[ERROR] more than two matches found, consider increasing <seed_length>')
    sys.exit(flush=True, file=sys.stderr)
else:
    seq = seq[0:match_lst[1]]

# print output
if output_path == '-':
    print(f'{header}\n{seq}', file=sys.stdout)
else:
    read_func = gzip.open if output_path.endswith('.gz') else open
    with read_func(output_path, 'w') as out_fasta:
        out_fasta.write(f'{header}\n{seq}')
