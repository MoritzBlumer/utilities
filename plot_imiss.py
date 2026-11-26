#!/usr/bin/env python
'''
Simple script to plot missingness per chromosome from vcftools --missing-indv
'''



## FILE INFO
__author__ = 'Moritz Blumer, 2026'
__email__ = 'lmb215@cam.ac.uk'



## SETUP

# packages
import argparse
import glob
import os
import sys
import pandas as pd
import plotly.express as px



## CLI

def cli():

    '''
    Parse command line arguments.
    '''

    parser = argparse.ArgumentParser(
        description='Plot .imiss files')

    # add arguments
    tab = '\u200B \u200B \u200B \u200B '
    parser.add_argument(
        'prefix',
        type=str,
        help='Full shared prefix (incl. any . or / etc) of all imiss files',
    )
    parser.add_argument(
        'chroms',
        type=str,
        help='Comma-separated list of chromosomes, used to construct file'
             ' names and to determine plotting order',
    )
    parser.add_argument(
        'output_path',
        type=str,
        help='Path for output file, suffix determines file format'
             ' (.pdf/.png/.html)',
    )

    # parse
    return parser.parse_args()



## FUNCTIONS

def find_suffix(prefix, chrom_lst):

    '''
    Find file unique shared suffix
    '''

    start = f'{prefix}{chrom_lst[0]}'
    candidate_lst = [
        file.removeprefix(start) \
            for file in glob.glob(f'{glob.escape(start)}*.imiss')
    ]
    shared_lst = [
        candidate for candidate in candidate_lst \
            if not candidate[0].isalnum() \
            if all(os.path.isfile(f'{prefix}{chrom}{candidate}') \
                for chrom in chrom_lst)
    ]

    if len(shared_lst) != 1:
        print(
            f'[ERROR] No shared suffix found: {shared_lst}',
            file=sys.stderr,
        )
        sys.exit(1)

    return shared_lst[0]


def read_data(prefix, suffix, chrom_lst):

    '''
    Read data from multiple files
    '''

    # compile genome df
    imiss_df_lst = []
    for chrom in chrom_lst:
        file = f'{prefix}{chrom}{suffix}'
        imiss_df = pd.read_csv(
            file,
            sep='\t',
        )
        imiss_df['chrom'] = chrom
        imiss_df_lst.append(imiss_df)
        print(
            f'[INFO] Read file {file}',
            file=sys.stderr,
        )
    imiss_df = pd.concat(imiss_df_lst, ignore_index=True)

    return imiss_df


def plot(miss_df, chrom_lst):
    
    '''
    Plot
    '''

    # plot
    fig = px.strip(
        miss_df,
        x='chrom',
        y='F_MISS',
        hover_data=['INDV', 'N_MISS', 'N_DATA'],
        category_orders={'chrom': chrom_lst},
        template='simple_white',
    )
    fig.update_traces(marker=dict(size=3, opacity=0.4), jitter=0.4)
    fig.update_layout(width=1100, height=450, showlegend=False)
    
    return fig


def save_fig(fig, output_path):

    '''
    Write output figure
    '''

    if output_path.endswith('.html'):
        fig.write_html(output_path)
    else:
        fig.write_image(output_path)



## MAIN

def main():

    '''
    Main.
    '''

    # parse command line arguments
    args = cli()

    # parse variable names
    prefix = args.prefix
    output_path =  args.output_path
    chroms = args.chroms

    # fetch chromosome list
    chrom_lst = chroms.split(',')

    # find shared suffix
    suffix = find_suffix(prefix, chrom_lst)

    print(
        f'\n[INFO] Found suffix {suffix}\n',
        file=sys.stderr,
    )

    # read data
    imiss_df = read_data(
        prefix,
        suffix,
        chrom_lst,
    )

    # plot
    imiss_fig = plot(
        imiss_df,
        chrom_lst,
    )

    # save
    save_fig(imiss_fig, output_path)

    print(
        f'\n[INFO] Done\n',
        file=sys.stderr,
    )


## EXECUTE

if __name__ == "__main__":
    main()
