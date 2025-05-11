#!/usr/bin/env python3
#
# Moritz Blumer | 2025-05-10
#
# This script trims trailing sequence from a linearized circular DNA assembly 
# by identifying and removing the region that matches the start of the sequence.
# Allows one mismatch.

# import packages
import sys
import os
import re
import gzip

# set k range
k_max = 50
k_min = 5

# CLI & help
global input_path, output_path

if len(sys.argv) < 3:
    print(
        '\n   python trim_circ_seq.py <input_path> <seed_length> <output_path>\n\n\
        <input_path>          str  path to input FASTA\n\
        <output_path>         str  path to output FASTA ("-" for STDOUT)\n\
        ',
    file=sys.stderr,
    )
    sys.exit()

# fetch arguments    
_, input_path, output_path = sys.argv

# check if input file was specified and exists
if not os.path.isfile(input_path):
    print(
        f'\n[ERROR] {input_path} does not exist.\n',
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

# initiate terminate
terminate = None

# try all k values from k_max to k_min
for k in range(k_min, k_max+1)[::-1]:
    
    # infer search string
    search_string = seq[0:k]
    match_lst = [m.start() for m in re.finditer(search_string, seq)]

    # if unique second match found, terminate
    if len(match_lst) == 2:
        seq = seq[0:match_lst[1]]
        terminate = True
        mismatches = 0
        break
    
    # if more than two unique matches are found this suggests a mismatch
    # (= difference between leading and trailing assembly). Try shifting the
    # search string to k + 1 and if finding a unique second match, trim off
    # the trailing sequence, offsetting again to the entire string incl. the
    # single mismatch
    if len(match_lst) > 2:
        skip = k+1
        for kk in range(k_min, k_max+1)[::-1]:
            # infer search string
            search_string = seq[skip:kk+skip]
            match_lst = [m.start() for m in re.finditer(search_string, seq)]

            # if unique second match found, terminate
            if len(match_lst) == 2:
                seq = seq[0:match_lst[1]-skip]
                terminate = True
                mismatches = 1
                break
        break

# if no unique second match could be found for the specified k range and 
# allowing for one mismatch, exit with error message
if not terminate:
    print(
        f'\n[ERROR] {input_path}: No trailing duplication found.\n',
        flush=True,
        file=sys.stderr,
    )
    sys.exit()

# print info
if mismatches == 0:
    print(
    f'[DONE] {input_path}: k={k}.',
    flush=True,
    file=sys.stderr,
    )
else:
    print(
    f'[DONE] {input_path}: k={kk+skip}; {mismatches} mismatch.',
    flush=True,
    file=sys.stderr,
    )

# print output
if output_path == '-':
    print(f'{header}\n{seq}\n', file=sys.stdout)
else:
    read_func = gzip.open if output_path.endswith('.gz') else open
    with read_func(output_path, 'w') as out_fasta:
        out_fasta.write(f'{header}\n{seq}\n')
