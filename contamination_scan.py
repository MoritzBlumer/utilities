#!/usr/bin/env python
'''
Compute allelic depths from SAM/BAM/CRAM using samtools mpileup. Compute
fraction of total sites with reads supporting non-major alleles.
Requires samtools version >=1.10.
'''



## FILE INFO
__author__ = 'Moritz Blumer, 2026'
__email__ = 'lmb215@cam.ac.uk'



## CONFIG

# general
INFO_SIZE = 1000000

# samtools mpileup settings
MPILEUP_C = 0
MPILEUP_d = 250000
MPILEUP_q = 20
MPILEUP_Q = 13
MPILEUP_ff = 'UNMAP,SECONDARY,QCFAIL,DUP'


## SETUP

# packages
import gzip
import sys
import shutil
import argparse
import subprocess



## CLI

def cli():

    '''
    Parse command line arguments.
    '''

    parser = argparse.ArgumentParser(
        description='Calculate allelic depths and fraction of surplus alleles'
                    ' from SAM/BAM/CRAM.')

    # add arguments
    tab = '\u200B \u200B \u200B \u200B '
    parser.add_argument(
        'CRAM_PATH',
        type=str,
        help='Input CRAM file (single sample).',
    )
    parser.add_argument(
        'REF_PATH',
        type=str,
        help='Corresponding reference FASTA (must have .FAI index).',
    )
    parser.add_argument(
        'SAMPLE_ID',
        type=str,
        help='Identifier for the input sample.',
    )
    parser.add_argument(
        '-r', '--region',
        dest='region',
        required=False,
        metavar='\b',
        default=None,
        help=f'{tab}Genomic region in format chr:start-stop [default: None]',
    )
    parser.add_argument(
        '-p', '--ploidy',
        dest='ploidy',
        type=int,
        required=False,
        metavar='\b',
        default=2,
        help=f'{tab}Classify reads supporting any minor allel more than ploidy'
             ' as surplus [default: 2]',
    )
    parser.add_argument(
        '-d', '--depth_threshold',
        dest='depth_threshold',
        type=float,
        required=False,
        metavar='\b',
        default=8,
        help='Minimum number of total reads to use a site for genome-wide'
             ' statistic [default: 8]',
    )
    parser.add_argument(
        '-m', '--min_count',
        dest='min_count',
        type=float,
        required=False,
        metavar='\b',
        default=1,
        help=f'{tab}Minimum number of reads supporting surplus allel(s) to '
             ' count a site as having surplus allele in genome-wide statistic'
             ' (threshold is applied to all surplus alleles jointly, not per'
             ' allele if there are multiple) [default: 1]',
    )
    parser.add_argument(
        '-w', '--write_per_site_tsv',
        dest='per_site_tsv',
        required=False,
        metavar='\b',
        default=None,
        help='Specify path for per-site output TSV containing chrom, pos,'
             ' allelic depths, total depth, surplus reads and surpus'
             ' proportion',
    )

    # parse
    return parser.parse_args()



## FUNCTIONS

def parse_cram(
        SAMTOOLS,
        MPILEUP_C,
        MPILEUP_d,
        MPILEUP_q,
        MPILEUP_Q,
        MPILEUP_ff,
        CRAM_PATH,
        REF_PATH,
        region,
        per_site_tsv,
        ploidy,
        depth_threshold,
        min_count,
        ):

    '''
    Main function to parse SAM records and make all calculations.
    '''

    # construct samtools command
    if region:
        shell_command = f'{SAMTOOLS} mpileup' \
            f' {CRAM_PATH}' \
            f' -f {REF_PATH}' \
            f' -r {region}' \
            f' -C {MPILEUP_C}' \
            f' -d {MPILEUP_d}' \
            f' -q {MPILEUP_q}' \
            f' -Q {MPILEUP_Q}' \
            f' --ff {MPILEUP_ff}' \
            f' --no-output-ins --no-output-ins' \
            f' --no-output-del --no-output-del' \
            f' --no-output-ends'
    else:
        shell_command = f'{SAMTOOLS} mpileup' \
            f' {CRAM_PATH}' \
            f' -f {REF_PATH}' \
            f' -C {MPILEUP_C}' \
            f' -d {MPILEUP_d}' \
            f' -q {MPILEUP_q}' \
            f' -Q {MPILEUP_Q}' \
            f' --ff {MPILEUP_ff}' \
            f' --no-output-ins --no-output-ins' \
            f' --no-output-del --no-output-del' \
            f' --no-output-ends'

    # execute samtools mpileup
    process = subprocess.Popen(
        shell_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
        text=True,
    )

    # initiate counts
    pass_sites = 0
    surp_sites = 0

    # open output file if specified
    if per_site_tsv:
        output_file = gzip.open(per_site_tsv, 'wt')
        output_file.write(
            'chrom\tpos\tA\tT\tC\tG\tdepth\tsurplus_reads\tsurplus_fraction\n'
        )
    else:
        output_file = None

    # iterrate STDIN lines
    for i, line in enumerate(process.stdout):

        # parse line
        line = line.strip().split()
        chrom = line[0]
        pos = line[1]
        ref = line[2].lower()
        pile = line[4]

        # count ref
        r = pile.count(',') + pile.count('.')

        # count alleles
        c_dct = {
            'a': pile.lower().count('a'),
            't': pile.lower().count('t'),
            'c': pile.lower().count('c'),
            'g': pile.lower().count('g'),
        }

        # set ref
        c_dct[ref] = r

        # sort
        counts = sorted(
            [
                c_dct['a'],
                c_dct['t'],
                c_dct['c'],
                c_dct['g']
            ]
        , reverse=True)

        # calculate total depth
        depth = sum(counts)

        # count reads supporting surplus alleles
        surplus = sum(counts[ploidy:])

        # calculate surpluss
        surplus_frac = round(surplus / depth, 2) if depth > 0 else 'NA'

        # update counters
        if depth > depth_threshold:
            pass_sites +=1
            if surplus > min_count:
                surp_sites +=1

        # write AD data
        if per_site_tsv:
            output_file.write(
                f"{chrom}\t{pos}"
                f"\t{c_dct['a']}\t{c_dct['t']}\t{c_dct['c']}\t{c_dct['g']}"
                f"\t{depth}\t{surplus}\t{surplus_frac}\n",
            )

        # print info
        if i % INFO_SIZE == 0:

            print(
                f'[INFO] Processed {i} sites',
                file=sys.stderr,
            )

    # close output file
    if output_file:
        output_file.close()

    return [pass_sites, surp_sites]


## MAIN

def main():

    '''
    Main.
    '''

    # parse commend line arguments
    args = cli()

    # parse variable names
    CRAM_PATH = args.CRAM_PATH
    REF_PATH = args.REF_PATH
    SAMPLE_ID =  args.SAMPLE_ID
    region = args.region
    ploidy = args.ploidy
    depth_threshold = args.depth_threshold
    min_count = args.min_count
    per_site_tsv = args.per_site_tsv

    # fetch samtools from execution $PATH
    SAMTOOLS = shutil.which('samtools')

    # check samtools version
    samtools_v = subprocess.getoutput(
        f'{SAMTOOLS} --version'
    ).split(' ')[1].split('-')[0]
    if float(samtools_v.split('.')[0]) <= 1 \
        and float(samtools_v.split('.')[1]) < 10:
        print(
            f'[ERROR] samtools must be at least v1.10 but found v{samtools_v}',
            file=sys.stdout,
            flush=True,
        )
        sys.exit(1)

    # run
    stats = parse_cram(
        SAMTOOLS,
        MPILEUP_C,
        MPILEUP_d,
        MPILEUP_q,
        MPILEUP_Q,
        MPILEUP_ff,
        CRAM_PATH,
        REF_PATH,
        region,
        per_site_tsv,
        ploidy,
        depth_threshold,
        min_count,
        )

    # print output to STDOUT
    print(
        f'{SAMPLE_ID}\t{stats[0]}\t{stats[1]}\t{round(stats[1]/stats[0], 8)}',
        file=sys.stdout,
        flush=True
    )



## EXECUTE

if __name__ == "__main__":
    main()
