#!/usr/bin/env python
'''
Visualize genome-wide association data
'''


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

    parser = argparse.ArgumentParser(description="Visualize pairwise genome \
        alignment(s) from PAF.")

    # add arguments
    parser.add_argument(
        'assoc_path',
        type=str,
        help='Input file with associations.',
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
        '-c', '--chrom_column',
        dest='chrom_column',
        required=False,
        metavar='\b',
        default='chrom',
        help='Column name that contains the chromosome [default: chrom]',
    )
    parser.add_argument(
        '-v', '--var_column',
        dest='var_column',
        required=False,
        metavar='\b',
        default='pos',
        help='Variant position [default: pos]',
    )
    parser.add_argument(
        '-p', '--p_val_column',
        dest='p_val_column',
        required=False,
        metavar='\b',
        default='p',
        help='P value column to use from GEMMA output [default: p]',
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
        '-l', '--log',
        dest='log',
        required=False,
        metavar='\b',
        default=True,
        help='"False" to plot raw values from p_val_column [default: True]',
    )
    parser.add_argument(
        '-a', '--alpha',
        dest='alpha',
        required=False,
        metavar='\b',
        default='0.01',
        help='Significance level for Bonferroni correction [default: 0.01]',
    )
    parser.add_argument(
        '-r', '--remove',
        dest='remove',
        required=False,
        metavar='\b',
        default=None,
        help='Remove n-th lowest quantile of associations to reduce plot size/'
             'complexity (applied after log unless disapled with --log False) '
             '[default: None]',
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
    return parser.parse_args()



## FUNCTIONS

def save_fig(fig, format_string, prefix):

    '''
    Save figure in specified formats.
    '''

    for fmt in format_string.split(','):
        if fmt in ['html', 'HTML']:
            fig.write_html(f'{prefix}.html')
        else:
            try:
                print(f'{prefix}.{fmt.lower()}')
                fig.write_image(f'{prefix}.{fmt.lower()}')
            except Exception as e:
                print(
                    f'[INFO] {fmt.lower()} not supported – skipping',
                    file=sys.stderr,
                    flush=True,
                )
                continue

## MAIN

def main():

    '''
    Main.
    '''

    # parse commend line arguments
    args = cli()

    # parse variable names
    assoc_path = args.assoc_path
    faidx_path = args.faidx_path
    output_prefix =  args.output_prefix
    chrom_column = args.chrom_column
    var_column = args.var_column
    p_val_column = args.p_val_column
    plot_fmt = args.plot_fmt
    log = args.log not in ['FALSE', 'False', False]
    alpha = float(args.alpha)
    remove = float(args.remove)
    fig_height = int(args.fig_height)
    fig_width = int(args.fig_width)

    # read FAIDX
    chrom_sizes_dct = pd.read_csv(
        faidx_path,
        sep='\t',
        usecols=[0,1],
        index_col=0,
        names=['chrom', 'size'],
    ).to_dict(orient='index')

    # read GWAS
    assoc_df = pd.read_csv(
        assoc_path,
        sep='\t',
    )

    # subset for selected p_val_column
    assoc_df = assoc_df[[chrom_column, var_column, p_val_column]]

    # calculate offsets
    chrom_lst = list(set(assoc_df[chrom_column]))
    offset = 0
    for chrom in chrom_sizes_dct.keys():
        if chrom in chrom_lst:
            chrom_sizes_dct[chrom]['offset'] = offset
            chrom_sizes_dct[chrom]['offset_mid'] = \
                offset + chrom_sizes_dct[chrom]['size']/2
            offset += chrom_sizes_dct[chrom]['size']
            chrom_sizes_dct[chrom]['offset_end'] = offset

    # add genome pos to assoc_df
    assoc_df['offset'] = [
        chrom_sizes_dct[chrom]['offset'] for chrom in assoc_df[chrom_column]
    ]
    assoc_df['genome_pos'] = assoc_df['offset'] + assoc_df[var_column]
    assoc_df = assoc_df[
        [chrom_column, var_column, 'genome_pos', p_val_column]
    ]

    # calculate -log10(p_val_column) & significance threshold
    if log:
        assoc_df['-log10_{p_val_column}'] = -np.log10(assoc_df[p_val_column])
        plot_col = '-log10_{p_val_column}'
        significance_threshold = -np.log10(alpha / len(assoc_df))
    else:
        plot_col = p_val_column

    # remove low associations to reduce plot size/complexity
    if remove is not None:
        r_threshold = assoc_df[plot_col].quantile(remove)
        assoc_df = assoc_df.loc[assoc_df[plot_col] >= r_threshold]
    else:
        r_threshold = None

    # plot
    fig = go.Figure()

    # collect chromosome labels
    chrom_label_lst = []

    # sort data by chromosome
    assoc_df[chrom_column] = pd.Categorical(assoc_df[chrom_column],
        categories=chrom_sizes_dct.keys(),
        ordered=True,
    )
    assoc_df = assoc_df.sort_values([chrom_column])

    # plot removal threshold
    if remove:
        fig.add_hline(
            y=r_threshold,
            line_dash='solid',
            line_color='grey',
            line_width=0.1,
        )
        fig.add_annotation(
                x=0,
                y=r_threshold,
                text=f"associations <{r_threshold} {plot_col} not shown",
                showarrow=False,
                xanchor='left',
                yanchor='top',
                yshift=-2,
                font={'size': 5, 'color': 'grey'},
                xref='x',
                yref='y',
            )

    # plot chromosome-wise
    for i, (chrom, chrom_df) in enumerate(
        assoc_df.groupby(chrom_column, observed=True)
        ):

        # add name
        chrom_label_lst.append(chrom)

        # alternating colors
        color = COLOR_A if (i+1) % 2 == 0 else COLOR_B

        # add points
        fig.add_trace(
            go.Scattergl(
                x=chrom_df['genome_pos'],
                y=chrom_df[plot_col],
                mode='markers',
                marker={
                    'size': 5,
                    'color': color,
                },
                customdata=chrom_df[chrom_column].astype(str) \
                    + ":" + chrom_df[var_column].astype(str),
                hovertemplate="<b>%{customdata}</b><br>" + plot_col \
                    + ": %{y:.2f}<extra></extra>",
            )
        )

    # plot significance threshold
    if log:
        fig.add_hline(
            y=significance_threshold,
            line_dash='dot',
            line_color='grey',
            line_width=1,
        )

    # formatting
    x_max = chrom_sizes_dct[list(assoc_df[chrom_column])[-2]]['offset_end']
    fig.update_layout(
        xaxis={
            'range': [0, x_max],
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
            yaxis_title=plot_col,
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
    save_fig(fig, plot_fmt, output_prefix)



## EXECUTE

if __name__ == "__main__":
    main()
