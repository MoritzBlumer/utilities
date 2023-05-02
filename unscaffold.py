#!/usr/bin/env python
#
# Moritz Blumer | 2023-05-02
#
# Split FASTA sequences at gaps of specified size


## File info
__author__ = 'Moritz Blumer, 2023'
__email__ = 'lmb215@cam.ac.uk'


## Dependencies
import sys
import gzip


## Main

def parse_arguments():
    '''
    Parse command line arguments & print help message if # of arguments is incorrect
    '''

    global input_path, output_path, contig_prefix, gap_size, out_line_len

    # fetch arguments    
    _, input_path, output_path, contig_prefix, gap_size, out_line_len = sys.argv

    # print help message if incorrect number of arguments was specified
    if len(sys.argv) != 6:
        print(
            '   python unscaffold.py <input path> <output path> <contig prefix> <gap size>\
                                <output line length>\n\n\
            <input path>             str    path to uncompressed or gzipped input FASTA\n\
            <output path>            str    path to output FASTA (*.gz for gzipped output)\n\
            <contig prefix>          str    prefix to use for output contigs (e.g. "contig_")\n\
            <gap size>               int    gap size to split at (this should be the gap\n\
                                            size used by the scaffolding software)\n\
            <output line length>     int    line length for output FASTA, e.g. 100',
        file=sys.stderr,
        )

    # change str to int where appropriate
    gap_size, out_line_len = int(gap_size), int(out_line_len)


def collect_scaffolds(input_path):
    '''
    Read scaffolds from input file into list
    '''

    scaffold = ''
    scaffold_lst = []

    read_func = gzip.open if input_path.endswith('.gz') else open

    with read_func(input_path, 'rt') as fasta:

        next(fasta) # skip first header

        for line in fasta:

            if line.startswith('>'):

                scaffold_lst.append(scaffold)
                
                scaffold = ''

            else:
                
                scaffold += line.strip()

    scaffold_lst.append(scaffold)

    return scaffold_lst


def unscaffold(scaffold_lst, gap_size):
    '''
    Split scaffolds into contigs at stretches of Ns of the specified gap_size used by the 
    scaffolding software
    '''
    
    contig_lst = []

    for scaffold in scaffold_lst:
        
        contig = ''
        gap = False

        for i in scaffold:

            if gap:
                
                if i == 'N':
                    gap += 1

                else:
                                
                    if gap != gap_size:
                        contig += ''.join(['N'] * gap) + i

                    else:   
                        if len(contig) > 0: contig_lst.append(contig) 
                        contig = i
                    
                    gap = False

            elif i == 'N':
                gap = 1

            else:
                contig += i
        
        if gap:
            if gap != gap_size:
                contig += ''.join(['N'] * gap) + i

        if len(contig_lst) > 0:
            contig_lst.append(contig)
    
    return contig_lst


def sort_and_save(contig_lst, output_path, contig_prefix, out_line_len):
    '''
    Sort contigs by size, and write them to file using specified line length
    '''

    contig_lst.sort(key=len, reverse=True)

    write_func = gzip.open if output_path.endswith('.gz') else open
    with write_func(output_path, 'wt') as fasta:
        
        for idx, contig in enumerate(contig_lst):
            fasta.write('>' + str(contig_prefix) + str (idx) + '\n')
            for line in range(0, len(contig), out_line_len):
                    fasta.write(contig[line:line+out_line_len] + '\n')
                    

## Main

def main():
    
    # parse command line arguments
    parse_arguments()

    # read in sequences
    scaffold_lst = collect_scaffolds(input_path)

    # split sequences into contigs
    contig_lst = unscaffold(scaffold_lst, gap_size)

    # sort by size and write output
    sort_and_save(contig_lst, output_path, contig_prefix, out_line_len)

if __name__ == '__main__':
    main()
