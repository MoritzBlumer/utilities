#!/usr/bin/env python
#
# Moritz Blumer | 2025-03-28
#
# Subset (optionally gzipped) BEAGLE file to a subset of samples and/or sites.
# BEAGLE is assumed to have marker, allele1, allele2 fields followed by three 
# GL/PL fields per sample.
# Sample order remains the same as in the input BEAGLE.


## FILE INFO
__author__ = 'Moritz Blumer, 2023'
__email__ = 'lmb215@cam.ac.uk'


## SETUP

#  import packages
import argparse
import gzip
import sys


## CLI

def cli():

    '''
    Parse command line arguments.
    '''

    global input_beagle_path, output_beagle_path, sample_file_path, subset, \
        drop_header

    parser = argparse.ArgumentParser(description="Subset BEAGLE file.")

    # add arguments
    parser.add_argument('input_beagle_path', type=str, help='Path to input \
        BEAGLE or BEAGLE.gz')
    parser.add_argument('output_beagle_path', type=str, help='Path to output \
        BEAGLE or BEAGLE.gz')
    parser.add_argument('-s', '--samples', type=str, help='Path to file with \
        one sample ID per line', required=False)
    parser.add_argument('-l', '--lines', type=int, help='Subset to every nth \
        line (incl. first line), specify e.g. 10', required=False)
    parser.add_argument('-d', '--drop_header', action='store_true', 
        help='do not write header line', required=False)

    # parse
    args = parser.parse_args()

    # reassign variable names
    input_beagle_path, output_beagle_path, sample_file_path, subset, \
        drop_header = args.input_beagle_path, args.output_beagle_path, \
        args.samples, args.lines, args.drop_header


## FUNCTIONS

def read_samples(sample_file_path):

    # read samples
    keep_sample_lst = []
    with open(sample_file_path, 'r') as sample_file:
        for sample in sample_file.readlines():
            keep_sample_lst.append(sample.strip())

    # ensure there are no duplicates
    if len(set(keep_sample_lst)) != len(keep_sample_lst):
        print(f'[ERROR] {sample_file_path} contains duplicate sample IDs.', 
            file=sys.stderr)
        sys.exit()

    return keep_sample_lst


def parse_beagle(input_beagle_path, output_beagle_path, sample_file_path, 
    subset, drop_header):
    
    '''
    parse BEAGLE file
    '''

    # parse samples if specified
    if sample_file_path:
        keep_sample_lst = read_samples(sample_file_path)

    # determine read function
    read_func = gzip.open if input_beagle_path.endswith('.gz') else open
    write_func = gzip.open if output_beagle_path.endswith('.gz') else open

    # iterrate file
    with read_func(input_beagle_path, 'rt') as input_beagle:
        with write_func(output_beagle_path, 'wt') as output_beagle:

            # subset sample mode:
            if sample_file_path:

                # parse header
                header = input_beagle.readline().strip().split('\t')

                # ensure all specified samples are in the header
                if not all(item in set(header[3:]) for item in keep_sample_lst):
                    print(f'[ERROR] {sample_file_path} contains sample IDs'
                        f' that are not in {input_beagle_path}', 
                        file=sys.stderr)
                    sys.exit()

                # get indices to keep
                keep_idx_lst = [0, 1, 2] + \
                    [i for i, x in enumerate(header) if x in keep_sample_lst]
                
                # write new header if not skipping
                if not drop_header:
                    header = '\t'.join([header[i] for i in keep_idx_lst])
                    output_beagle.write(header + '\n')

                # normal mode (write all sites)
                if not subset:
                    for line in input_beagle:
                        line = line.strip().split('\t')
                        line = [line[i] for i in keep_idx_lst]
                        output_beagle.write('\t'.join(line) + '\n')

                # subset mode
                else:
                    n = subset - 1
                    for line in input_beagle:
                        n += 1
                        if n == subset:
                            line = line.strip().split('\t')
                            line = [line[i] for i in keep_idx_lst]
                            output_beagle.write('\t'.join(line) + '\n')
                            n = 0

            # only lines mode
            else:
                
                if subset:

                    # parse header
                    header = input_beagle.readline().strip()

                    # write header if not skipping
                    if not drop_header:
                        output_beagle.write(header + '\n')

                    # write every nth line                
                    n = subset - 1
                    for line in input_beagle:
                        n += 1
                        if n == subset:
                            output_beagle.write(line)
                            n = 0
                
                else:
                    print(f'[ERROR] No action specified. Exiting.', 
                        file=sys.stderr)
                    sys.exit()



## MAIN

def main():

    # parse input arguments
    cli()

    # parse BEAGLE file
    parse_beagle(input_beagle_path, output_beagle_path, sample_file_path, \
        subset, drop_header)

if __name__ == '__main__':
    main()
