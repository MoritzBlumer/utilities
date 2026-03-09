#!/usr/bin/env python3

"""
Type stratified haplotype pattern in a single winpca window.
"""


## IMPORT PACKAGES

import pandas as pd
import sys



## MAIN

def main():

    '''
    Main.
    '''

    # parse arguments
    if len(sys.argv) == 1:
        sys.exit()
    else:
        _, pc_tsv_path, pos, u_threshold, l_threshold = sys.argv

    # convert pos and thresholds to numbers
    try:
        pos, u_threshold, l_threshold = \
            float(pos), float(u_threshold), float(l_threshold)
    except:
        print(
            '[ERROR] pos, l_threshold & u_threshold must all be numbers but the'
            f' following values were providedwere: {pos}, {l_threshold} &'
            f'{u_threshold}',
            flush=True,
            file=sys.stderr,
        )
        sys.exit()

    # read data
    pc_df = pd.read_csv(
        pc_tsv_path,
        sep='\t',
    )

    # subset to focal window
    pc_df = pc_df.loc[pc_df['pos'] == pos]

    # check if anything was found
    if len(pc_df) != 1:
        print(
            f'[ERROR] No window found at pos {pos:.0f}.',
            flush=True,
            file=sys.stderr,
        )
        sys.exit()

    # transpose and reformat names
    pc_df = pc_df.set_index('pos').transpose()
    pc_df.index.name = 'sample_id'
    pc_df.columns.name = None
    pc_df.columns = ['PC']
    pc_df['typing'] = None

    # find top samples
    pc_df.loc[
        pc_df['PC'] > u_threshold
        , 'typing'] = 2

    # middle samples
    pc_df.loc[
        (pc_df['PC'] <= u_threshold) & (pc_df['PC'] >= l_threshold)
        , 'typing'] = 1

    # lower samples
    pc_df.loc[
        pc_df['PC'] < l_threshold
        , 'typing'] = 0

    # print output
    for i, row in pc_df.iterrows():
        print(
            f'{i}\t{row["PC"]}\t{row["typing"]}',
            file=sys.stdout,
        )



## EXECUTE
if __name__ == "__main__":
    try:
        main()
    except:
        print(
        '\nHELP MESSAGE'
        '\n\n'
        '   type_winpca_locus.py <pc_tsv> <window_pos> <upper_threshold> '
           '<lower_threshold>\n'
        '\n'
        '      <pc_tsv>            str   Path to WinPCA PC file'
                                        ' (.pc_*.tsv.gz)\n'
        '      <window_pos>        int   Target window position\n'
        '      <upper_threshold>  float  Upper PC threshold\n'
        '      <lower_threshold>  float  Lower PC threshold\n'
        '\n'
        '   STDOUT is <samples_id>, <pc_value>, <typing> as TSV with:\n'
        '      --> "2" if pc_value > upper_threshold\n'
        '      --> "1" if upper_threshold >= pc_value >= lower_threshold\n'
        '      --> "0" if pc_value < lower_threshold\n',
        file=sys.stderr,
        )
        sys.exit()
