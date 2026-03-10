#!/usr/bin/env python
#
# Moritz Blumer | 2026-02-13
#
# Visualize pairwise genome alignment(s) from PAF



## FILE INFO
__author__ = 'Moritz Blumer, 2026'
__email__ = 'lmb215@cam.ac.uk'



## CONFIG

# mimimum alignment score
MIN_ALN_SCORE = 60

# ploting colors
SEQ_BORDER_COLOR = '#4f4f4f'
FWD_ALN_COLOR_1 = '#a3a3a3'
FWD_ALN_COLOR_2 = '#707070'
REV_ALN_COLOR = '#2c52db'
GAP_COLOR = '#f53665'

# plotting opacities
ALN_COLOR_OPACITY = 0.4
CHROM_COLOR_OPACITY = 0.7

# minimum display width of gaps relative to largest genome size
MIN_GAP_WIDTH_PCT = 0.01

# y axis boundaries for r, p, s sequences
P_Y_UPPER = 600
P_Y_LOWER = 400
R_Y_UPPER = 100
R_Y_LOWER = -100
S_Y_UPPER = -400
S_Y_LOWER = -600



## SETUP

#  import packages
import argparse
import sys
import re
import pandas as pd
import plotly.graph_objects as go



## CLI

def cli():

    '''
    Parse command line arguments.
    '''

    global p_paf_path, output_prefix, plot_fmt, seq_prefixes, s_paf_path, \
        r_fai_path, p_fai_path, s_fai_path, r_gaps_path, p_gaps_path, \
        s_gaps_path, assembly_names

    parser = argparse.ArgumentParser(description="Visualize pairwise genome \
        alignment(s) from PAF.")

    # add arguments
    parser.add_argument(
        'p_paf_path',
        type=str,
        help='Primary input PAF',
    )
    parser.add_argument(
        'output_prefix',
        type=str,
        help='Output prefix',
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
        '-p', '--seq_prefixes',
        dest='seq_prefixes',
        required=False,
        metavar='\b',
        default='chr,scf,SCAFFOLD_,atg_chr',
        help='Comma-separated list of prefixes to name-sort reference'
             ' sequences for plotting. Specifying "chr,scf" will sort as'
             ' chr[1,2,3..10,11..,chrM],scf[1,2,3..10,11...]. Default is'
             ' "chr,scf,SCAFFOLD_,atg_chr"',
    )
    parser.add_argument(
        '-s', '--s_paf_path',
        dest='s_paf_path',
        required=False,
        metavar='\b',
        default=None,
        help='Secondary input PAF',
    )
    parser.add_argument(
        '-a', '--r_fai_path',
        dest='r_fai_path',
        required=False,
        metavar='\b',
        default=None,
        help='Chromosome-sizes-file for reference genome (.fai)',
    )
    parser.add_argument(
        '-b', '--p_fai_path',
        dest='p_fai_path',
        required=False,
        metavar='\b',
        default=None,
        help='Chromosome-sizes-file for primary aligned assembly (.fai)',
    )
    parser.add_argument(
        '-c', '--s_fai_path',
        dest='s_fai_path',
        required=False,
        metavar='\b',
        default=None,
        help='Chromosome-sizes-file for secondary aligned assembly (.fai)',
    )
    parser.add_argument(
        '-x', '--r_gaps_path',
        dest='r_gaps_path',
        required=False,
        metavar='\b',
        default=None,
        help='Reference genome gap locations (chrom\tgap_start\tgap_end)',
    )
    parser.add_argument(
        '-y', '--p_gaps_path',
        dest='p_gaps_path',
        required=False,
        metavar='\b',
        default=None,
        help='Primary aligned assembly gap locations (same format)',
    )
    parser.add_argument(
        '-z', '--s_gaps_path',
        dest='s_gaps_path',
        required=False,
        metavar='\b',
        default=None,
        help='Secondary aligned assembly gap locations (same format)',
    )
    parser.add_argument(
        '-n', '--names',
        dest='assembly_names',
        required=False,
        metavar='\b',
        default=None,
        help='List of assembly names used to label output plot. Provide a'
             ' comma-separated list in this order:'
             ' reference_name,primary_name[,secondary_name',
    )

    
    # parse
    args = parser.parse_args()

    # reassign variable names
    p_paf_path, output_prefix, plot_fmt, seq_prefixes, s_paf_path, \
        r_fai_path, p_fai_path, s_fai_path, r_gaps_path, p_gaps_path, \
        s_gaps_path, assembly_names = \
    args.p_paf_path, args.output_prefix, args.plot_fmt, args.seq_prefixes, \
        args.s_paf_path, args.r_fai_path, args.p_fai_path, args.s_fai_path, \
        args.r_gaps_path, args.p_gaps_path, args.s_gaps_path, \
        args.assembly_names

# # def read_paf(input_paf_path):
# p_paf_path = '/Users/moritzblumer/Downloads/fLabFul222.P_vs_fAstCal68.P.paf.gz'
# s_paf_path = '/Users/moritzblumer/Downloads/fLabFul222.S_vs_fAstCal68.P.paf.gz'
# seq_prefixes = 'chr,scf,SCAFFOLD_,atg_chr'
# r_fai_path = '/Users/moritzblumer/Downloads/fAstCal68.P.fa.gz.fai'
# p_fai_path = '/Users/moritzblumer/Downloads/fLabFul222.P.fa.gz.fai'
# s_fai_path = '/Users/moritzblumer/Downloads/fLabFul222.S.fa.gz.fai'
# r_gaps_path = '/Users/moritzblumer/Downloads/fAstCal68.P.fa.gz.gaps.tsv'
# p_gaps_path = '/Users/moritzblumer/Downloads/fLabFul222.P.fa.gz.gaps.tsv'
# s_gaps_path = '/Users/moritzblumer/Downloads/fLabFul222.S.fa.gz.gaps.tsv'





## FUNCTIONS

def read_paf(paf_path):

    '''
    Read first 12 (core) columns from a PAF file.
    '''

    # read alignment data
    paf_df = pd.read_csv(
        paf_path,
        sep='\t',
        usecols=range(12),
        names=[
            'q_name',
            'q_len',
            'q_start',
            'q_end',
            'q_strand',
            'r_name',
            'r_len',
            'r_start',
            'r_end',
            'aln_res',
            'aln_len',
            'aln_qual',
        ],
        dtype={
            'q_name': str,
            'q_len': int,
            'q_start': int,
            'q_end': int,
            'q_strand': str,
            'r_name': str,
            'r_len': int,
            'r_start': int,
            'r_end': int,
            'aln_res': int,
            'aln_len': int,
            'aln_qual': int,
        }
    )

    return paf_df


def get_len_dct(input_df, r_q):
    '''
    Get dict {r_name: {'r_len': r_len}} or {q_name: {'q_len': q_len}} from 
    paf_df.
    '''
    len_df = input_df[
        [f'{r_q}_name', f'{r_q}_len']
    ].drop_duplicates().set_index(f'{r_q}_name')
    len_dct = len_df.to_dict(orient='index')
    
    return len_dct


def get_genome_size(seq_dct, q_r):
    '''
    Get the total size of a genome from a seq_dct.
    '''
    
    genome_size = 0
    for seq_name in seq_dct.keys():
        genome_size += seq_dct[seq_name][f'{q_r}_len']
    
    return genome_size


def natural_key(string):
    '''
    Split string into a list of integers and text chunks.
    '''

    key = [
        int(c) if c.isdigit() else c for c in re.split(r'(\d+)', string)
    ]

    return key


def v_sort(lst, prefix):
    '''
    Sort a prefixed list, numbers first, then sortings like 'M'.
    '''

    prefix_len = len(prefix)
    sort_lst = []

    for i in lst:
        if i.startswith(prefix):
            sort_lst.append(i[prefix_len:])

    out_lst = [
        f'{prefix}{x}' for x in sorted(sort_lst, key=natural_key)
        
        ]

    return out_lst


def fetch_q_rank_offset(aln_df, q_dct, r_name_order_lst):

    '''
    Infer query sequence plotting order based on the order of reference 
    sequences and the order of individual queries aligning to the same
    reference sequence.
    '''

    # iterrate over query sequences
    for q_name in q_dct.keys():

        # set default values
        best_r = 'NA'
        max_aln_len = 0
        best_r_aln_mid = 0

        # iterrate over sequences that have alignments with current query
        for r_name in r_name_order_lst:

            # check if there are alignments
            if q_name in set(aln_df.loc[aln_df['r_name'] == r_name]['q_name']):

                # get all pairwise alignments
                r_q_df = aln_df.loc[(aln_df['r_name'] == r_name) \
                        & (aln_df['q_name'] == q_name)]
                
                # get total lengths across all alignments for current ref
                aln_len_sum = sum(r_q_df['aln_len'])

                # if this is more than the current best reference
                if aln_len_sum > max_aln_len:

                    # set current ref as best_r and update max_aln_len
                    max_aln_len = aln_len_sum
                    best_r = r_name

                    # calculate midpoint
                    best_r_aln_mid = \
                        min(r_q_df['r_start']) \
                            + (max(r_q_df['r_end'] - min(r_q_df['r_start']))/2)

        # save best_r per query sequence after full pass of ref sequences into
        # q_dct
        q_dct[q_name]['best_r'] = best_r
        q_dct[q_name]['best_r_aln_mid'] = best_r_aln_mid

    # convert q_dct to df and sort by reference order
    q_rank_df = pd.DataFrame.from_dict(q_dct, orient='index')
    q_rank_df['best_r'] = pd.Categorical(
        q_rank_df['best_r'],
        categories=r_name_order_lst + ['NA'],
        ordered=True
    )
    q_rank_df = q_rank_df.sort_values(['best_r', 'best_r_aln_mid'])

    # infer plotting rank for each query sequence according to this order
    q_rank_df['rank'] = list(range(0, len(q_rank_df)))

    # infer offset based on rank and save both in q_dct
    q_offset = 0
    for q_name, row in q_rank_df.iterrows():
        q_dct[q_name]['rank'] = row['rank']
        q_dct[q_name]['offset'] = q_offset
        q_offset += row['q_len']
    
    return q_dct


def get_aln(aln_dct, aln_df, aln_type):
    
    '''
    Fetch alignment list and largest total alignment length per query.
    '''

    # iterrate through alignment file by query sequence
    for q_name in set(aln_df['q_name']):

        # get all alignments per query sequence
        q_name_df = aln_df.loc[aln_df['q_name'] == q_name]

        # get reference with longest alignment per query sequence
        best_r = \
            q_name_df.groupby('r_name')['aln_len'].sum().sort_values().index[-1]
        
        # subset to that query sequence
        q_r_df = q_name_df.loc[q_name_df['r_name'] == best_r].copy()

        # compute midpoint per alignment
        q_r_df['r_mid'] = \
            q_r_df['r_start'] + (q_r_df['r_end']- q_r_df['r_start'])/2
        
        # convert alignments to list
        aln_lst = []
        for _, row in q_r_df.iterrows():
            aln_lst.append(
                [
                    row['q_strand'],
                    row['q_start'],
                    row['q_end'],
                    row['r_start'],
                    row['r_end']
                ]
            )

        # add to aln_dct
        aln_dct[best_r][aln_type][q_name] = {
            'median_pos': q_r_df['r_mid'].median(),
            'aln_lst': aln_lst,
        }

    return aln_dct


def check_noaln_ref(r_name_lst, aln_dct, aln_type):
    '''
    Check if any reference sequences have no alignments and print info.
    '''
    r_name_no_aln_lst = [
    r_name for r_name in r_name_lst if aln_dct[r_name][aln_type] == {}
        ]
    _type = 'primary' if aln_type == 'p_aln' else 'secondary'
    if len(r_name_no_aln_lst) > 0:                                       
        print(
            f'[INFO] No major alignments found for reference sequences '
            f'{", ".join(r_name_no_aln_lst)} in {_type} alignments', 
            flush=True,
            file=sys.stderr,
        )

    return r_name_no_aln_lst


def read_gaps(gaps_path, seq_dct):

    '''
    Read gaps from file and add them to a sequence.
    '''

    gaps_df = pd.read_csv(
        gaps_path,
        sep='\t',
        names=['chrom', 'start', 'end'],
    )
    for seq_name in set(gaps_df['chrom']):
        gaps_lst = []
        seq_df = gaps_df.loc[gaps_df['chrom'] == seq_name]
        for _, row in seq_df.iterrows():
            start = row.iloc[1]
            end = row.iloc[2]
            size = abs(end - start)
            gaps_lst.append([start, end, size])
        seq_dct[seq_name]['gaps'] = gaps_lst
    
    return seq_dct


def plot_alignments(
    SEQ_BORDER_COLOR,
    FWD_ALN_COLOR_1,
    FWD_ALN_COLOR_2,
    REV_ALN_COLOR,
    GAP_COLOR,
    ALN_COLOR_OPACITY,
    CHROM_COLOR_OPACITY,
    MIN_GAP_WIDTH_PCT,
    P_Y_UPPER,
    P_Y_LOWER,
    R_Y_UPPER,
    R_Y_LOWER,
    S_Y_UPPER,
    S_Y_LOWER,
    aln_dct,
    s,
    r_dct,
    p_dct,
    s_dct,
    r_name_order_lst,
    p_offset,
    s_offset,
    max_genome_size,
    assembly_names,
    ):

    '''
    Plotting function.
    '''

    fig = go.Figure()


    ## ALIGNMENTS

    def plot_aln(
            fig,
            r_name,
            aln_type,
            q_dct,
            q_offset,
            r_offset,
            R_Y_UPPER,
            R_Y_LOWER,
            Q_Y_UPPER,
            Q_Y_LOWER,
            FWD_ALN_COLOR,
        ):

        '''
        Plot alignments.
        '''

        #for q_name in aln_dct[r_name][aln_type]['q_order_lst']:
        for q_name in aln_dct[r_name][aln_type].keys():

            # plot traces
            for aln in aln_dct[r_name][aln_type][q_name]['aln_lst']:
                
                # get coordinates
                q_strand = aln[0]
                q_start = aln[1] + q_offset + q_dct[q_name]['offset']
                q_end = aln[2] + q_offset + q_dct[q_name]['offset']
                r_start = aln[3] + r_offset
                r_end = aln[4] + r_offset

                # determine alignment color
                aln_color = FWD_ALN_COLOR

                # swap coordinates and color for inverted alignments
                if q_strand == '-':
                    aln_color = REV_ALN_COLOR
                    q_start, q_end = q_end, q_start

                # set display text
                text = f'{r_name}:{aln[3]}-{aln[4]} | {q_name}:{aln[1]}-{aln[2]}'
                
                # plot
                fig.add_trace(
                    go.Scatter(
                        x=[
                            q_start,
                            q_start,
                            q_end,
                            q_end,
                            r_end,
                            r_end,
                            r_start,
                            r_start,
                            q_start,
                        ], 
                        y=[
                            Q_Y_LOWER,
                            Q_Y_UPPER,
                            Q_Y_UPPER,
                            Q_Y_LOWER,
                            R_Y_UPPER,
                            R_Y_LOWER,
                            R_Y_LOWER,
                            R_Y_UPPER,
                            Q_Y_LOWER,
                        ],
                        fill='toself',
                        fillcolor=aln_color,
                        opacity=ALN_COLOR_OPACITY,
                        hoveron='fills',
                        mode='lines',
                        name=text,
                        legendgroup=f'{r_name}',
                        showlegend=False,
                        line=dict(width=0, color=aln_color),
                    )
                )

    # iterrate thrugh reference names in order
    r_offset = 0

    for i, r_name in enumerate(r_name_order_lst):

        # set even/uneven chrom color
        FWD_ALN_COLOR = FWD_ALN_COLOR_1 if i % 2 == 0 else FWD_ALN_COLOR_2

        # plot alignments
        plot_aln(
            fig,
            r_name,
            'p_aln',
            p_dct,
            p_offset,
            r_offset,
            R_Y_UPPER,
            R_Y_LOWER,
            P_Y_UPPER,
            P_Y_LOWER,
            FWD_ALN_COLOR,
        )
        if s:
            plot_aln(
                fig,
                r_name,
                's_aln',
                s_dct,
                s_offset,
                r_offset,
                R_Y_LOWER,
                R_Y_UPPER,
                S_Y_LOWER,
                S_Y_UPPER,
                FWD_ALN_COLOR,
            )

        # add reference offset to r_dct
        r_dct[r_name]['offset'] = r_offset

        # itteratively add sequence length to offset
        r_offset += aln_dct[r_name]['r_len']


    ## REFERENCE BOXES

    for i, r_name in enumerate(r_name_order_lst):

        # set even/uneven chrom color
        FWD_ALN_COLOR = FWD_ALN_COLOR_1 if i % 2 == 0 else FWD_ALN_COLOR_2

        # get coordinates
        r_offset = r_dct[r_name]['offset']
        r_start = r_offset
        r_end = r_dct[r_name]['r_len'] + r_offset

        # plot
        fig.add_trace(
            go.Scatter(
                x=[
                    r_start,
                    r_start,
                    r_end,
                    r_end,
                    r_start,
                ], 
                y=[
                    R_Y_LOWER,
                    R_Y_UPPER,
                    R_Y_UPPER,
                    R_Y_LOWER,
                    R_Y_LOWER,
                ],
                fill='toself',
                fillcolor=FWD_ALN_COLOR,
                opacity=CHROM_COLOR_OPACITY,
                mode='lines',
                name=f'{r_name}',
                legendgroup=f'{r_name}',
                showlegend=True,
                line=dict(width=1, color=SEQ_BORDER_COLOR),
            )
        )


    ## QUERY BOXES

    def plot_q_boxes(
            fig,
            q_dct,
            q_offset,
            r_name_order_lst,
            Q_Y_UPPER,
            Q_Y_LOWER
        ):
        
        '''
        Plot bounding boxes for query sequences
        '''
        for q_name in q_dct.keys():
            
            # set even/uneven chrom color
            r_name_order_ext_lst = r_name_order_lst + ['NA']
            r_index = r_name_order_ext_lst.index(q_dct[q_name]['best_r'])
            FWD_ALN_COLOR = FWD_ALN_COLOR_1 if r_index % 2 == 0 \
                else FWD_ALN_COLOR_2

            # get coordinates
            q_start = q_offset + q_dct[q_name]['offset']
            q_end = q_offset + q_dct[q_name]['offset'] + q_dct[q_name]['q_len']

            # plot
            fig.add_trace(
                go.Scatter(
                    x=[
                        q_start,
                        q_start,
                        q_end,
                        q_end,
                        q_start,
                    ], 
                    y=[
                        Q_Y_LOWER,
                        Q_Y_UPPER,
                        Q_Y_UPPER,
                        Q_Y_LOWER,
                        Q_Y_LOWER,
                    ],
                    fill='toself',
                    fillcolor=FWD_ALN_COLOR,
                    opacity=CHROM_COLOR_OPACITY,
                    mode='lines',
                    name=f'{q_name}',
                    legendgroup=f'{q_dct[q_name]["best_r"]}',
                    showlegend=False,
                    line=dict(width=1, color=SEQ_BORDER_COLOR),
                )
            )

    plot_q_boxes(fig, p_dct, p_offset, r_name_order_lst, P_Y_UPPER, P_Y_LOWER)
    if s:
        plot_q_boxes(fig, s_dct, s_offset, r_name_order_lst, S_Y_UPPER, S_Y_LOWER) 


    ## GAPS

    # add dummy trace for legend
    if any([r_gaps_path, p_gaps_path, s_gaps_path]):
        fig.add_trace(
            go.Scatter(
                x=[0, 0, 0, 0], 
                y=[0, 0, 0, 0], 
                fill='toself',
                fillcolor=GAP_COLOR,
                mode='lines',
                name=f'gaps',
                legendgroup=f'gaps',
                line=dict(width=0),
                showlegend=True,
            )
        )

    def plot_gaps(fig, seq_dct, offset, max_width, Y_UPPER, Y_LOWER, MIN_GAP_WIDTH_PCT, GAP_COLOR):

        min_width = max_width * MIN_GAP_WIDTH_PCT * 0.01

        for seq_name in seq_dct.keys():
            
            if 'gaps' in seq_dct[seq_name].keys():

                for gap in seq_dct[seq_name]['gaps']:

                    # get coordinates
                    g_start = gap[0] + offset + seq_dct[seq_name]['offset']
                    g_end = gap[1]  + offset + seq_dct[seq_name]['offset']
                    g_size = gap[2]

                    # adjust to minimum display width
                    if g_size < min_width:
                        g_midpoint = g_start + (g_end - g_start)/2
                        g_start = g_midpoint - min_width/2
                        g_end = g_midpoint + min_width/2
                    
                    # plot
                    fig.add_trace(
                        go.Scatter(
                            x=[
                                g_start,
                                g_start,
                                g_end,
                                g_end,
                                g_start,
                            ], 
                            y=[
                                Y_LOWER,
                                Y_UPPER,
                                Y_UPPER,
                                Y_LOWER,
                                Y_LOWER,
                            ],
                            fill='toself',
                            fillcolor=GAP_COLOR,
                            opacity=1,
                            mode='lines',
                            name=f'{seq_name}:{gap[0]}-{gap[1]}',
                            legendgroup=f'gaps',
                            showlegend=False,
                            line=dict(width=1, color=GAP_COLOR),
                        )
                    )

    plot_gaps(
        fig,
        r_dct,
        0,
        max_genome_size,
        R_Y_UPPER,
        R_Y_LOWER,
        MIN_GAP_WIDTH_PCT,
        GAP_COLOR,
    )

    plot_gaps(
        fig,
        p_dct,
        p_offset,
        max_genome_size,
        P_Y_UPPER,
        P_Y_LOWER,
        MIN_GAP_WIDTH_PCT,
        GAP_COLOR,
    )
    if s:
        plot_gaps(
            fig,
            s_dct,
            s_offset,
            max_genome_size,
            S_Y_UPPER,
            S_Y_LOWER,
            MIN_GAP_WIDTH_PCT,
            GAP_COLOR,
        )

    # ASSEMBLY NAMES
    if assembly_names:

        try:

            # get list
            assembly_names_lst = assembly_names.split(',')

            # reference
            fig.add_annotation(
                x=0-0.01*max_genome_size,
                y=R_Y_LOWER + (R_Y_UPPER-R_Y_LOWER)/2,
                text=assembly_names_lst[0],
                showarrow=False,
                xanchor='right',
                yanchor='middle',
                bgcolor=SEQ_BORDER_COLOR,
                bordercolor=SEQ_BORDER_COLOR,
                borderwidth=2,
                borderpad=4,
                font={'size': 20, 'color': 'white'},
            )

            # primary
            fig.add_annotation(
                x=p_offset-0.01*max_genome_size,
                y=P_Y_LOWER + (P_Y_UPPER-P_Y_LOWER)/2,
                text=assembly_names_lst[1],
                showarrow=False,
                xanchor='right',
                yanchor='middle',
                bgcolor=SEQ_BORDER_COLOR,
                bordercolor=SEQ_BORDER_COLOR,
                borderwidth=2,
                borderpad=4,
                font={'size': 20, 'color': 'white'},
            )

            # secondary
            if s:
                fig.add_annotation(
                    x=s_offset-0.01*max_genome_size,
                    y=S_Y_LOWER + (S_Y_UPPER-S_Y_LOWER)/2,
                    text=assembly_names_lst[2],
                    showarrow=False,
                    xanchor='right',
                    yanchor='middle',
                    bgcolor=SEQ_BORDER_COLOR,
                    bordercolor=SEQ_BORDER_COLOR,
                    borderwidth=2,
                    borderpad=4,
                    font={'size': 18, 'color': 'white'},
                )
        except():
            print(
                '[INFO] Failed to plot assembly names, please provide 2' \
                ' (or 3, if specifying a secondary PAF) assembly names as a'
                ' comma-separated list (see help message))'
            )


    ## FORMATTING

    # fig size
    fig_width = 1800
    fig_height = 500 if not s else 800

    # transparent background, hide axis labels, format legend text size and spacing
    fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        xaxis={'visible': False},
        yaxis={'visible': False},
        legend={'font': {'size': 15}, 'tracegroupgap': 0},
        width=fig_width, 
        height=fig_height,
        font_family='Arial',
    )

    return fig


def main():

    # parse command line arguments
    cli()

    # read P PAF
    p_paf_df = read_paf(p_paf_path)

    # set secondary alignment flag
    s = True if s_paf_path else None

    # read S PAF
    if s: s_paf_df = read_paf(s_paf_path)

    # extract reference names & sizes from P PAF
    r_dct = get_len_dct(p_paf_df, 'r')

    # [optional] extract reference names & sizes from S PAF and add anything 
    # missing to r_dct
    if s:
        if s: s_r_dct = get_len_dct(s_paf_df, 'r')
        for key in s_r_dct.keys():
            if key in r_dct.keys():
                if r_dct[key] != s_r_dct[key]:
                    print(
                        f'[ERROR] ref sizes for {key} do not match between the'
                        ' specified PAF files',
                        file=sys.stderr,
                        flush=True,
                    )
            else:
                r_dct[key] = s_r_dct[key]

    # [optional] extract names & sizes from ref FAIDX and add anything missing
    # to r_dct
    if r_fai_path:
        r_fai_df = pd.read_csv(
        r_fai_path,
        sep='\t',
        usecols=[0,1],
        names=['r_name', 'r_len'],
        )
        r_fai_dct = get_len_dct(r_fai_df, 'r')
        for key in r_fai_dct.keys():
            if key in r_dct.keys():
                if r_dct[key] != r_fai_dct[key]:
                    print(
                        f'[ERROR] ref sizes for {key} do not match between the'
                            ' specified PAF files and the chromosome file',
                        file=sys.stderr,
                        flush=True,
                    )
            else:
                r_dct[key] = r_fai_dct[key]

    # extract P names & sizes from P PAF
    p_dct = get_len_dct(p_paf_df, 'q')

    # [optional] extract names & sizes from P FAIDX and add anything missing to
    # p_dct
    if p_fai_path:
        p_fai_df = pd.read_csv(
        p_fai_path,
        sep='\t',
        usecols=[0,1],
        names=['q_name', 'q_len'],
        )
        p_fai_dct = get_len_dct(p_fai_df, 'q')
        for key in p_fai_dct.keys():
            if key in p_dct.keys():
                if p_dct[key] != p_fai_dct[key]:
                    print(
                        f'[ERROR] query sizes for {key} do not match between the'
                        f' specified PAF file ({p_paf_path}) and the chromosome'
                        f' file ({p_fai_path}).',
                        file=sys.stderr,
                        flush=True,
                    )
            else:
                p_dct[key] = p_fai_dct[key]

    # extract S names & sizes from S PAF
    s_dct = get_len_dct(s_paf_df, 'q') if s else None

    # [optional] extract names & sizes from S FAIDX and add anything missing to
    # s_dct
    if s and s_fai_path:
        s_fai_df = pd.read_csv(
        s_fai_path,
        sep='\t',
        usecols=[0,1],
        names=['q_name', 'q_len'],
        )
        s_fai_dct = get_len_dct(s_fai_df, 'q')
        for key in s_fai_dct.keys():
            if key in s_dct.keys():
                if s_dct[key] != s_fai_dct[key]:
                    print(
                        f'[ERROR] query sizes for {key} do not match between'
                        f' the specified PAF file ({s_paf_path}) and the'
                        f' chromosome file ({s_fai_path}).',
                        file=sys.stderr,
                        flush=True,
                    )
            else:
                s_dct[key] = s_fai_dct[key]

    # fetch genome sizes
    r_size = get_genome_size(r_dct, 'r')
    p_size = get_genome_size(p_dct, 'q')
    s_size = get_genome_size(s_dct, 'q') if s else 0

    # get list of reference sequence names
    r_name_lst = list(r_dct.keys())

    # sort reference sequence names for plotting
    # sort optionally according to seq_prefixes, otherwise alphanumerically
    if seq_prefixes:
        r_name_order_lst = []
        for prefix in seq_prefixes.split(','):
            r_name_order_lst += v_sort(r_name_lst, prefix)
        r_name_order_lst += sorted(
            [x for x in r_name_lst if x not in r_name_order_lst]
        )
    else:
        r_name_order_lst = sorted(r_name_lst)

    # ensure all sequences are still there
    missing_seq_lst = []
    for i in r_name_lst:
        if i not in r_name_order_lst:
            missing_seq_lst.append(i)
    if len(missing_seq_lst) > 0:
        print(
            f'[ERROR] Sequence(s) {", ".join(missing_seq_lst)} not captured by the'
            f' the specified prefixes {seq_prefixes}',
            flush=True,
            file=sys.stderr,
        ),
        sys.exit()
    if len(r_name_lst) != len(r_name_order_lst):
        print(
            f'ERROR Sth. is wrong with the v_sort() logic --> tell the author',
            flush=True,
            file=sys.stderr,
        )
        sys.exit()

    # filter low scoring alignments
    p_aln_df = p_paf_df.loc[p_paf_df['aln_qual'] >= MIN_ALN_SCORE]
    if s:
        s_aln_df = s_paf_df.loc[s_paf_df['aln_qual'] >= MIN_ALN_SCORE]

    # infer query sequence plotting order based on the order of reference 
    # sequences and the order of individual queries aligning to the same
    # reference sequence.
    p_dct = fetch_q_rank_offset(p_aln_df, p_dct, r_name_order_lst)
    if s:
        s_dct = fetch_q_rank_offset(s_aln_df, s_dct, r_name_order_lst)

    # start aln_dct (copy of r_dct)
    aln_dct = r_dct.copy()

    # add placeholder items for query alignments
    for r_name in aln_dct.keys():
        aln_dct[r_name]['p_aln'] = {}
    if s:
        for r_name in aln_dct.keys():
            aln_dct[r_name]['s_aln'] = {}

    # fetch query alignments from filtered alignment df
    aln_dct = get_aln(aln_dct, p_aln_df, 'p_aln')
    if s: aln_dct = get_aln(aln_dct, s_aln_df, 's_aln')

    # check if there are major alignments to all reference sequences, otherwise
    # print INFO message
    r_name_no_aln_dct = {}
    r_name_no_aln_dct['p_aln'] = check_noaln_ref(r_name_lst, aln_dct, 'p_aln')
    if s: r_name_no_aln_dct['s_aln'] = check_noaln_ref(r_name_lst, aln_dct, 's_aln')

    # determine plotting offsets
    p_offset = (r_size-p_size)/2
    s_offset = (r_size-s_size)/2 if s else 0

    # read in gaps if specified
    if r_gaps_path: r_dct = read_gaps(r_gaps_path, r_dct)
    if p_gaps_path: p_dct = read_gaps(p_gaps_path, p_dct)
    if s_gaps_path and s: s_dct = read_gaps(s_gaps_path, s_dct)

    # get largest genome size
    max_genome_size = max(r_size, p_size, s_size)

    # plot
    fig = plot_alignments(
        SEQ_BORDER_COLOR,
        FWD_ALN_COLOR_1,
        FWD_ALN_COLOR_2,
        REV_ALN_COLOR,
        GAP_COLOR,
        ALN_COLOR_OPACITY,
        CHROM_COLOR_OPACITY,
        MIN_GAP_WIDTH_PCT,
        P_Y_UPPER,
        P_Y_LOWER,
        R_Y_UPPER,
        R_Y_LOWER,
        S_Y_UPPER,
        S_Y_LOWER,
        aln_dct,
        s,
        r_dct,
        p_dct,
        s_dct,
        r_name_order_lst,
        p_offset,
        s_offset,
        max_genome_size,
        assembly_names,
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