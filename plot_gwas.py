#!/usr/bin/env python
#
# Moritz Blumer | 2026-03-12
#
# Visualize GEMMA formatted GWAS output



## FILE INFO
__author__ = 'Moritz Blumer, 2026'
__email__ = 'lmb215@cam.ac.uk'



## CONFIG
COLOR_A = '#7d7d7d'
COLOR_B = '#363636'



## SETUP

# packages
import sys
import argparse
import pandas as pd
import numpy as np
import plotly.graph_objects as go



## CLI

def cli():

    '''
    Parse command line arguments.
    '''

    global gwas_path, faidx_path, output_prefix, p_val_column, plot_fmt, \
        alpha, fig_height, fig_width

    parser = argparse.ArgumentParser(description="Visualize pairwise genome \
        alignment(s) from PAF.")

    # add arguments
    parser.add_argument(
        'gwas_path',
        type=str,
        help='GEMMA-formatted input file.',
    )
    parser.add_argument(
        'faidx_path',
        type=str,
        help='Chromosome sizes tsv (e.g. FAIDX index).',
    )
    parser.add_argument(
        'output_prefix',
        type=str,
        help='Output prefix',
    )
    parser.add_argument(
        '-p', '--p_val_column',
        dest='p_val_column',
        required=False,
        metavar='\b',
        default='p_wald',
        help='P value column to use from GEMMA output [default: p_wald]',
    )
    parser.add_argument(
        '-f', '--format',
        dest='plot_fmt',
        required=False,
        metavar='\b',
        default='HTML',
        help='Output plot file format ("HTML", "PDF", "SVG" or "PNG"),'
             ' may also be a comma-separated list [default: "HTML"]',
    )
    parser.add_argument(
        '-a', '--alpha',
        dest='alpha',
        required=False,
        metavar='\b',
        default='0.05',
        help='Global significance level for Bonferroni correction [default: 0.05]',
    )
    parser.add_argument(
        '-y', '--fig_height',
        dest='fig_height',
        required=False,
        metavar='\b',
        default=500,
        help='Figure height [default: 500]',
    )
    parser.add_argument(
        '-x', '--fig_width',
        dest='fig_width',
        required=False,
        metavar='\b',
        default=1500,
        help='Figure width [default: 1500]',
    )

    # parse
    args = parser.parse_args()

    # reassign variable names
    gwas_path, faidx_path, output_prefix, p_val_column, plot_fmt, alpha, \
        fig_height, fig_width = \
        args.gwas_path, args.faidx_path, args.output_prefix, \
            args.p_val_column, args.plot_fmt, \
            float(args.alpha), int(args.fig_height), int(args.fig_width)



## MAIN

def main():

    '''
    Main.
    '''

    # parse commend line arguments
    cli()

    # read FAIDX
    chrom_sizes_dct = pd.read_csv(
        faidx_path,
        sep='\t',
        usecols=[0,1],
        index_col=0,
        names=['chrom', 'size'],
    ).to_dict(orient='index')

    # read GWAS
    gwas_df = pd.read_csv(
        gwas_path,
        sep='\t',
        usecols=[0,2,12,13,14],
    )

    # subset for selected p_val_column
    gwas_df = gwas_df[['chr', 'ps', p_val_column]]

    # rename columns
    gwas_df.columns = ['chrom', 'pos', p_val_column]

    # add chr prefix back in
    gwas_df['chrom'] = 'chr' + gwas_df['chrom'].astype(str)

    # calculate offsets
    chrom_lst = list(set(gwas_df['chrom']))
    offset = 0
    for chrom in chrom_sizes_dct.keys():
        if chrom in chrom_lst:
            chrom_sizes_dct[chrom]['offset'] = offset
            chrom_sizes_dct[chrom]['offset_mid'] = \
                offset + chrom_sizes_dct[chrom]['size']/2
            offset += chrom_sizes_dct[chrom]['size']
            chrom_sizes_dct[chrom]['offset_end'] = offset

    # add genome pos to gwas_df
    gwas_df['offset'] = [
        chrom_sizes_dct[chrom]['offset'] for chrom in gwas_df['chrom']
    ]
    gwas_df['genome_pos'] = gwas_df['offset'] + gwas_df['pos']
    gwas_df = gwas_df[['chrom', 'pos', 'genome_pos', p_val_column]]

    # calculate -log10(p_val_column)
    gwas_df['-log10_p'] = -np.log10(gwas_df[p_val_column])

    # plot
    fig = go.Figure()

    # collect chromosome labels
    chrom_label_lst = []

    # sort data by chromosome
    gwas_df['chrom'] = pd.Categorical(gwas_df['chrom'],
        categories=chrom_sizes_dct.keys(),
        ordered=True,
    )
    gwas_df = gwas_df.sort_values(['chrom'])

    # calculate significance threshold
    significance_threshold = -np.log10(alpha / len(gwas_df))

    # plot chromosome-wise
    for i, (chrom, chrom_df) in enumerate(
        gwas_df.groupby('chrom', observed=True)
        ):

        # add name
        chrom_label_lst.append(chrom)

        # alternating colors
        color = COLOR_A if (i+1) % 2 == 0 else COLOR_B

        # add points
        fig.add_trace(
            go.Scattergl(
                x=chrom_df['genome_pos'],
                y=chrom_df['-log10_p'],
                mode='markers',
                marker={
                    'size': 5,
                    'color': color,
                },
                customdata=chrom_df['chrom'].astype(str) + ":" + chrom_df['pos'].astype(str),
                hovertemplate="<b>%{customdata}</b><br>-log10(p): %{y:.2f}<extra></extra>"        )
        )

    # plot threshold
    fig.add_hline(
        y=significance_threshold,
        line_dash='dot',
        line_color='grey',
        line_width=1,
    )

    # formatting
    fig.update_layout(
        xaxis={
            'range': [0, chrom_sizes_dct[list(gwas_df['chrom'])[-2]]['offset_end']],
            'showgrid': False,
            'zeroline': False,
            'showline': True,
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
            'ticks': 'outside',
            'tickmode': 'array',
            'tickvals': [chrom_sizes_dct[x]['offset_mid'] for x in chrom_lst],
            'ticktext': chrom_lst,
        },
        yaxis={
            'showgrid': False,
            'zeroline': False,
            'showline': True,
            'linewidth': 1,
                'linecolor': 'black',
                'mirror': True,
                'ticks': 'outside'
            },
            height=fig_height,
            width=fig_width,
            title=None,
            xaxis_title='Genomic position',
            yaxis_title='-log10(p-value)',
            template='plotly_white',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            showlegend=False,
            font={
                'family': 'Arial',
                'size': 11,
                'color': 'black'
            },
        )

    # save figure
    for fmt in plot_fmt.split(','):
        if fmt in ['html', 'HTML']:
            fig.write_html(f'{output_prefix}.html')
        else:
            try:
                print(f'{output_prefix}.{fmt.lower()}')
                fig.write_image(f'{output_prefix}.{fmt.lower()}')
            except:
                print(
                    f'[INFO] {fmt.lower()} not supported – skipping',
                    file=sys.stderr,
                    flush=True,
                )
                sys.exit()



## EXECUTE

if __name__ == "__main__":
    main()
