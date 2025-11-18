#!/usr/bin/env python3

# imports
import sys
from collections import Counter

# parameters
MAX_MISSINGNESS = 0.05

# define mappings
base_map_dct = {
    '.': 0,
    'A': 1,
    'a': 1,
    'C': 2,
    'c': 2,
    'G': 3,
    'g': 3,
    'T': 4,
    't': 4,
}

gt_map = {
    './.': [0, 0],
    './0': [0, 1],
    './1': [0, 2],
    './2': [0, 3],
    './3': [0, 4],
    '0/.': [1, 0],
    '0/0': [1, 1],
    '0/1': [1, 2],
    '0/2': [1, 3],
    '0/3': [1, 4],
    '1/.': [2, 0],
    '1/0': [2, 1],
    '1/1': [2, 2],
    '1/2': [2, 3],
    '1/3': [2, 4],
    '2/.': [3, 0],
    '2/0': [3, 1],
    '2/1': [3, 2],
    '2/2': [3, 3],
    '2/3': [3, 4],
    '3/.': [4, 0],
    '3/0': [4, 1],
    '3/1': [4, 2],
    '3/2': [4, 3],
    '3/3': [4, 4],
    '.|.': [0, 0],
    '.|0': [0, 1],
    '.|1': [0, 2],
    '.|2': [0, 3],
    '.|3': [0, 4],
    '0|.': [1, 0],
    '0|0': [1, 1],
    '0|1': [1, 2],
    '0|2': [1, 3],
    '0|3': [1, 4],
    '1|.': [2, 0],
    '1|0': [2, 1],
    '1|1': [2, 2],
    '1|2': [2, 3],
    '1|3': [2, 4],
    '2|.': [3, 0],
    '2|0': [3, 1],
    '2|1': [3, 2],
    '2|2': [3, 3],
    '2|3': [3, 4],
    '3|.': [4, 0],
    '3|0': [4, 1],
    '3|1': [4, 2],
    '3|2': [4, 3],
    '3|3': [4, 4],
}

if sys.argv[1] == '-h':
    print(
        'CHROM\tPOS\t.\tA\tC\tG\tT',
        file=sys.stdout,
    )

# parse lines from STDIN
for line in sys.stdin:
    if line.startswith('#'):
        continue

    # initiate allele counts
    #                   .  A  C  G  T
    allele_count_lst = [0, 0, 0, 0, 0]

    # parse VCF data rows
    chrom, pos, _, ref, alt, _, _, _, _, gt_fields = \
        line.rstrip('\n').split('\t', 9)

    # parse and encode alleles
    allele_lst = [ref] + alt.split(',')
    allele_lst = [ref if x == '.' else x for x in allele_lst]
    allele_lst = [0] + [base_map_dct[x] for x in allele_lst]

    # extract genotypes
    gt_lst = [g.partition(':')[0] for g in gt_fields.split('\t')]

    # count gts
    gt_count_dct = Counter(gt_lst)

    # convert to allele counts
    for gt, count in gt_count_dct.items():
        a_lst = gt_map[gt]
        a_1, a_2 = a_lst[0], a_lst[1]
        allele_count_lst[allele_lst[a_1]] += count
        allele_count_lst[allele_lst[a_2]] += count

    # write to STDOUT
    line = "\t".join([chrom, str(pos)] + [str(x) for x in allele_count_lst])
    print(
        line,
        file=sys.stdout,
    )

