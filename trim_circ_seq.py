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
global input_path, output_path, seed_length

if len(sys.argv) < 4:
    print(
        '\n   python trim_circ_seq.py <input_path> <seed_length> <output_path>\n\n\
        <input_path>          str  path to input FASTA\n\
        <output_path>         str  path to output FASTA ("-" for STDOUT)\n\
        <seed_length>         int  number of bp at start of sequence to query\n\
        ',
    file=sys.stderr,
    )
    sys.exit()

# fetch arguments    
_, input_path, output_path, seed_length = sys.argv

# convert optional paths to NONE bool if "NONE"
try:
    seed_length = int(seed_length)
except:
    print(
        '\n[ERROR] <seed_length> must be an integer.', 
        flush=True,
        file=sys.stderr,
    )
    sys.exit()

# check if input file was specified and exists
if not os.path.isfile(input_path):
    print(
        '\n[ERROR] {input_path} does not exist.',
        flush=True,
        file=sys.stderr,
    )
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
    print(
        f'\n[ERROR] {input_path}: no other matches found, consider reducing <seed_length>.',
        flush=True,
        file=sys.stderr,
    )
    sys.exit()

elif len(match_lst) > 2:
    print(
        f'\n[ERROR] {input_path}: more than two matches found, consider increasing <seed_length>',
        flush=True,
        file=sys.stderr,
    )
    sys.exit()
else:
    seq = seq[0:match_lst[1]]

# print output
if output_path == '-':
    print(f'{header}\n{seq}\n', file=sys.stdout)
else:
    read_func = gzip.open if output_path.endswith('.gz') else open
    with read_func(output_path, 'w') as out_fasta:
        out_fasta.write(f'{header}\n{seq}\n')
