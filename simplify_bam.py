#!/usr/bin/env python
#
# Moritz Blumer | 2023-05-15
#
# Simplify CIGAR string and remove tags in SAM/BAM file to reduce size and IGV load time (using 
# pysam functions) and write read length to 9th field



## File info
__author__ = 'Moritz Blumer, 2023'
__email__ = 'lmb215@cam.ac.uk'


## Dependencies
import pysam
import sys, os


## Functions

def parse_arguments():
    '''
    Parse command line arguments & print help message if # of arguments is incorrect
    '''

    global input_path, output_path

    # print help message if incorrect number of arguments was specified
    if len(sys.argv) < 3:
        print(
            '\n   python simplify_bam.py <input_path> <output_path>\n\n\
            <input path>           str    path to input SAM/BAM\n\
            <output path>          str    path to output BAM\n',
        file=sys.stderr,
        )
        sys.exit()

    # fetch arguments    
    _, input_path, output_path = sys.argv

    # print info
    print('\n[INFO] Input file: ' + input_path, file=sys.stderr)
    print('\n[INFO] Output file: ' + output_path, file=sys.stderr)


def init_output():
    '''
    Open input file & output file and append new PG line to header
    '''
    # open input file
    in_bam = pysam.AlignmentFile(input_path, "rb")

    # extract input header
    header = in_bam.header.as_dict()

    # add new PG line
    header['PG'].append(
                        {'ID': 'bam_simplify.py',
                         'PN': 'bam_simplify.py',
                         'PP': 'bam_simplify.py',
                        },
                        )

    # open output file and add modified header
    out_bam = pysam.AlignmentFile(output_path, "wb", header=header)

    return in_bam, out_bam


def process_records(in_bam, out_bam):

    '''
    compensate alignment length difference between reference and query by adding# an insertion or
    deletion to the matching sequence (unless alignment length is equal)
    '''

    # iterate input BAM:
    for alignment in in_bam.fetch():

        # only modify alignments ≥ 0
        if not alignment.query_length == 0:

            # calculate difference between reference/query alignment length
            diff = alignment.reference_length - alignment.query_alignment_length

            # handle different scenarios to obtain compensation factor
            if diff == 0:
                cigar_match = [(0, alignment.reference_length)]
                compensation = [(0, 0)]
            elif diff > 0:
                compensation = [(2, abs(diff))]
                cigar_match = [(0, alignment.reference_length - abs(diff))]
            else:
                compensation = [(1, abs(diff))]
                cigar_match = [(0, alignment.reference_length)]
            
            # fetch hard-/ soft clipping
            start_soft_clipped = []
            end_soft_clipped = []
            start_hard_clipped = []
            end_hard_clipped = []
            
            # no hard clipping
            if alignment.cigartuples[0][0] == 4:
                start_soft_clipped = [alignment.cigartuples[0]]
            if alignment.cigartuples[-1][0] == 4:
                end_soft_clipped = [alignment.cigartuples[-1]]

            # hard clipping
            if alignment.cigartuples[0][0] == 5:
                start_hard_clipped = [alignment.cigartuples[0]]
                if alignment.cigartuples[1][0] == 4:
                    start_soft_clipped = [alignment.cigartuples[1]]
            if alignment.cigartuples[-1][0] == 5:
                end_hard_clipped = [alignment.cigartuples[-1]]
                if alignment.cigartuples[-2][0] == 4:
                    end_soft_clipped = [alignment.cigartuples[-2]]
            
            # compile new simplified CIGAR string
            new_cigartuples = start_hard_clipped + start_soft_clipped + cigar_match \
                + compensation + end_soft_clipped + end_hard_clipped
            alignment.cigartuples = new_cigartuples

            # ADDITIONALLY write read length into 9th field:
            start_hard_clipped = [(0, 0)] if start_hard_clipped == [] else start_hard_clipped
            end_hard_clipped = [(0, 0)] if end_hard_clipped == [] else end_hard_clipped
            alignment.template_length = start_hard_clipped[0][1] + alignment.query_length + end_hard_clipped[0][1]
            
            # remove query sequence
            alignment.query_sequence = None
            
            # remove tags
            for tag in alignment.get_tags():
                tag = tag[0]
                alignment.set_tag(tag,
                                  value=None,
                                  )
                 
            # alignment.set_tag('ql', value=alignment.query_alignment_length)

            # write output line
            out_bam.write(
                alignment
            )
            
    # close input file
    in_bam.close()

    #  close output file and index
    out_bam.close()
    pysam.index(output_path)

    # print exit message and exit
    print(
        '\n[INFO] Done\n',
        file=sys.stderr,
    )


## Main

def main():

    # parse input arguments
    parse_arguments()

    # initiate output file
    in_bam, out_bam = init_output()

    # process input records and write output
    process_records(in_bam, out_bam)


if __name__ == '__main__':
    main()
