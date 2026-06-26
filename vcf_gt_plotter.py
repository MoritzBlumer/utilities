#!/usr/bin/env python3

'''
Parse and plot SNP genotypes from VCF/BCF
'''


## FILE INFO
__author__ = 'Moritz Blumer, 2026'
__email__ = 'lmb215@cam.ac.uk'



## CONFIG

# settings
F_MISSING = 1                         # bcftools F_MISSING setting
MISSING_COLOR = '#e0e0e0'           # plot color for missing sites/GT calls
LOCATIONS_DEFAULT_COLOR = '#0000FF' # default locations color
LOCATIONS_DEFAULT_SPREAD = 100        # default locations vertical spread (%)
SNPS_MODE_N_TICKS = 6                 # number of X ticks in snps mode

# color codes
multi_color_dct = {
    0:  MISSING_COLOR, # .
    1:  '#E41A1C',   # AA
    2:  '#a38cab',   # AT/TA
    3:  '#d4c772',   # AC/CA
    4:  '#e3aada',   # AG/GA
    5:  '#377EB8',   # TT
    6:  '#80ada6',   # TC/CT
    7:  '#979cde',   # TG/GT
    8:  '#4DAF4A',   # CC
    9:  '#97c29e',   # CG/GC
    10: '#984EA3',   # GG
    11: '#ffffff',   # major GT
}
bw_color_dct = {
    0:  MISSING_COLOR, # .
    1:  '#000000',   # AA
    2:  '#808080',   # AT/TA
    3:  '#808080',   # AC/CA
    4:  '#808080',   # AG/GA
    5:  '#000000',   # TT
    6:  '#808080',   # TC/CT
    7:  '#808080',   # TG/GT
    8:  '#000000',   # CC
    9:  '#808080',   # CG/GC
    10: '#000000',   # GG
    11: '#ffffff',   # major GT
}
rw_color_dct = {
    0:  MISSING_COLOR, # .
    1:  '#ff0000',   # AA
    2:  '#ff8080',   # AT/TA
    3:  '#ff8080',   # AC/CA
    4:  '#ff8080',   # AG/GA
    5:  '#ff0000',   # TT
    6:  '#ff8080',   # TC/CT
    7:  '#ff8080',   # TG/GT
    8:  '#ff0000',   # CC
    9:  '#ff8080',   # CG/GC
    10: '#ff0000',   # GG
    11: '#ffffff',   # major GT
}



## SETUP

# packages
import sys
import shutil
import argparse
import subprocess
from collections import Counter
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# genotype endoding
gt_code_dct = {
    'AA': 1,
    'AT': 2,
    'TA': 2,
    'AC': 3,
    'CA': 3,
    'AG': 4,
    'GA': 4,
    'TT': 5,
    'TC': 6,
    'CT': 6,
    'TG': 7,
    'GT': 7,
    'CC': 8,
    'CG': 9,
    'GC': 9,
    'GG': 10,
    '..': 0,
    'A.': 0,
    '.A': 0,
    'T.': 0,
    '.T': 0,
    'C.': 0,
    '.C': 0,
    'G.': 0,
    '.G': 0,
}
code_gt_dct = {
    1:  'AA',
    2:  'AT',
    3:  'AC',
    4:  'AG',
    5:  'TT',
    6:  'TC',
    7:  'TG',
    8:  'CC',
    9:  'CG',
    10: 'GG',
    0:  '..',
}



## CLI

def cli():

    '''
    Parse command line arguments.
    '''

    parser = argparse.ArgumentParser(description="Parse and plot SNP \
        genotypes from VCF/BCF.")

    # add arguments
    tab = '\u200B \u200B \u200B \u200B '
    parser.add_argument(
        'vcf_path',
        type=str,
        help='input VCF/BCF with annotated PASS sites (plot_mode=<full> '
             ' requires VCF/BCF with invariant positions)',
    )
    parser.add_argument(
        'region',
        type=str,
        help='genomic region in format "chrom:start-end"')
    parser.add_argument(
        'ref_fasta_path',
        type=str,
        help='FASTA used as reference in VCF/BCF',
    )
    parser.add_argument(
        'sample_id_path',
        type=str,
        help='file with sample IDs to include; determines plotting order'
             ' (except if phenotypes are specified)',
    )
    parser.add_argument(
        'output_prefix',
        type=str,
        help='output prefix',
    )
    parser.add_argument(
        '-v', '--viz_mode',
        dest='viz_mode',
        type=str,
        required=False,
        metavar='\b',
        choices=['full', 'snps'],
        default='snps',
        help=f'{tab}<full> to display all sites accross region or <snps> to'
              ' display only variable SNPs [default: snps]',
    )
    parser.add_argument(
        '-c', '--color_mode',
        dest='color_mode',
        type=str,
        required=False,
        metavar='\b',
        choices=['multi', 'bw', 'rw'],
        default='multi',
        help=f'{tab}<multi> to use 10 distinct colors representing all'
              ' genotypes or <bw>/<rw> to use black-white/red-white to only'
              ' distinguish hom vs. het (polarized by major allele) [default:'
              '  multi]',
    )
    parser.add_argument(
        '-m', '--no_mask_major_gt',
        dest='no_mask_major_gt',
        required=False,
        default=False,
        action='store_true',
        help='set flag to disable default bahavior of masking the major GT '
             ' by setting its color to white (--> if set, all genotypes are '
             ' colored)'
    )
    parser.add_argument(
        '-p', '--phenotypes',
        dest='phenotypes',
        type=str,
        required=False,
        metavar='\b',
        default=False,
        help=f'{tab}TSV with phenotypes for all included samples which will'
              ' be used to sort samples on the x axis (overrides sample order'
              ' in <sample_id_path>; requires -n/--phenotype_name)',
    )
    parser.add_argument(
        '-n', '--phenotype_name',
        dest='phenotype_name',
        type=str,
        required=False,
        metavar='\b',
        default=False,
        help='column name in phenotypes file (requires -p/--phenotypes)',
    )
    parser.add_argument(
        '-a', '--annotate',
        dest='annotate',
        action='store_true',
        default=False,
        help='set to hover-annotate <phenotype_name>, but keep the order from '
             '<sample_id_path>',
    )
    parser.add_argument(
        '-r', '--rev_sample_order',
        dest='rev_sample_order',
        action='store_true',
        default=False,
        help='set to reverse sample plotting order (only effective with '
             '-p/-n)',
    )
    parser.add_argument(
        '-l', '--locations',
        dest='locations',
        type=str,
        required=False,
        metavar='\b',
        default=False,
        help=f'{tab}TSV with regions to highlight on top of the plot with'
              ' fields: (1) sequence ID, (2) start, (3) end, (4) [optional]'
              ' HEX code to specify color, (5) [optional] an integer [1-100]'
              ' specifying in per cent the vertical spread of the element on'
              ' the locations track',
    )
    parser.add_argument(
        '-f', '--format',
        dest='plot_fmt',
        required=False,
        metavar='\b',
        default='HTML',
        help=f'{tab}output plot file format (<HTML>, <PDF>, <SVG> or <PNG>),'
              ' may also be a comma-separated list [default: HTML]',
    )
    parser.add_argument(
        '-x', '--fig_width',
        dest='fig_width',
        required=False,
        metavar='\b',
        default=1400,
        help=f'{tab}figure width [default: 1400]',
    )
    parser.add_argument(
        '-y', '--fig_height',
        dest='fig_height',
        required=False,
        metavar='\b',
        default=600,
        help=f'{tab}figure height [default: 600]',
    )

    # parse
    return parser.parse_args()



##  FUNCTIONS

def read_phenotype_data(
        phenotypes,
        phenotype_name,
        sample_lst,
        annotate,
        rev_sample_order,
    ):

    '''
    read phenotypes and subset to sample list and target phenotype
    '''

    # read
    phenotypes_df = pd.read_csv(
        phenotypes,
        sep='\t',
    )

    # subset samples
    try:
        phenotypes_df = phenotypes_df.loc[
            phenotypes_df[phenotypes_df.columns[0]].isin(sample_lst)
        ]
    except:
        print(
            '[ERROR] Please provide phenotype data for all samples in sample' \
            ' list',
            file=sys.stderr,
        )
        sys.exit(1)

    # subset columns
    if phenotype_name not in phenotypes_df.columns:
        print(
            f'[ERROR] Specified phenotype column ({phenotype_name}) not'
            ' found',
            file=sys.stderr,
        )
        sys.exit(1)
    phenotypes_df = phenotypes_df[[phenotypes_df.columns[0], phenotype_name]]
    phenotypes_df.columns = ['sample_id', phenotype_name]

    # sort by phenotype
    if not annotate:
        phenotypes_df = phenotypes_df.sort_values(
            phenotype_name,
            ascending=True,
        )
    
    # reverse order if set
    if rev_sample_order:
        phenotypes_df = phenotypes_df.iloc[::-1]

    # update sample list
    sample_lst = phenotypes_df['sample_id'].to_list()

    return phenotypes_df, sample_lst


def read_locations_data(
        locations,
        chrom,
        start,
        end,
        locations_default_color,
        locations_default_spread,
    ):

    '''
    Read locations data and subset to focal region
    '''

    # read
    locations_df = pd.read_csv(
        locations,
        sep='\t',
        header=None,
    )

    # format
    if len(locations_df.columns) == 5:
        locations_df.columns = ['chrom', 'from', 'to', 'color', 'spread']
    if len(locations_df.columns) == 4:
        locations_df.columns = ['chrom', 'from', 'to', 'color']
        locations_df['spread'] = locations_default_spread
    elif len(locations_df.columns) == 3:
        locations_df.columns = ['chrom', 'from', 'to']
        locations_df['color'] = locations_default_color
        locations_df['spread'] = locations_default_spread
    elif len(locations_df.columns) < 3:
        print(
            '[ERROR] locations file (-l/--locations) must have at least'
            'columns',
            file=sys.stderr,
        )
        sys.exit(1)

    # drop locations outside focal region
    locations_df = locations_df.loc[
        (locations_df['chrom'] == chrom) & \
        (
            (locations_df['from'] >= start) | \
            (locations_df['to'] <= end)
        )
    ]

    # if no locations left, set locations to False and print message
    if locations_df.empty:
        locations = False
        print(
             '[INFO] no locations (-l/--locations) found overlapping focal'
            f' region {start}-{end}',
            file=sys.stderr,
        )

    return locations, locations_df


def run_bcftools(
        BCFTOOLS,
        vcf_path,
        chrom,
        start,
        end,
        sample_id_path,
        ref_fasta_path,
    ):

    '''
    Run bcftools subprocesss
    '''

    # construct bcftools command
    shell_command = f'{BCFTOOLS} view' \
                    f'  -S {sample_id_path}' \
                    f'  -f PASS' \
                    f'  -O u' \
                    f'  {vcf_path}' \
                    f'  {chrom}:{start}-{end}' \
                    f'| {BCFTOOLS} filter' \
                    f'  --SnpGap 5' \
                    f'  -O u' \
                    f'| {BCFTOOLS} view' \
                    f'  -V indels,mnps,ref,bnd,other' \
                    f'  -O u' \
                    f'| {BCFTOOLS} view' \
                    f"  -i 'F_MISSING<{F_MISSING}'" \
                    f'  -O u' \
                    f'| {BCFTOOLS} norm' \
                    f'  -a' \
                    f'  --fasta-ref {ref_fasta_path}' \
                    f'  -O u' \
                    f'| {BCFTOOLS} view' \
                    f'  -A' \
                    f'  -a' \
                    f'  -O v'

    # execute bcftools
    bcftools_process = subprocess.Popen(
        shell_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )

    return bcftools_process


def parse_vcf_data(
        bcftools_process,
        sample_lst,
        viz_mode,
    ):

    '''
    Parse bcftools output and return an array of integer-encoded genotypes
    '''

    # iterrate STDIN lines
    gt_idx = 0
    vcf_sample_idx_lst = []
    pos_lst = []
    gts_lst = []
    major_gt_lst = []

    # skip header and read variant sample order
    for line in bcftools_process.stdout:

        # extract sample order
        if line.startswith('#CHROM'):
            vcf_sample_lst = line.strip().split('\t')[9:]
            vcf_sample_idx_lst = [vcf_sample_lst.index(x) for x in sample_lst]
            break

    # parse variant lines
    prev_pos = None
    miss_pos_gts = [0 for x in sample_lst]
    for line in bcftools_process.stdout:

        # parse
        line = line.strip().split('\t')

        # skip sites with spanning indels
        if '*' in line[4]: continue

        # extract relevant information
        pos = int(line[1])
        ref = line[3]
        alt = line[4].split(',')
        fmt = line[8]
        sample_fields = line[9:]

        # add empty entries for missing POS
        if viz_mode == 'full' and prev_pos is not None:
            while pos > prev_pos+1:
                pos_lst.append(prev_pos+1)
                gts_lst.append(miss_pos_gts)
                major_gt_lst.append(-1)
                prev_pos +=1

        # extract GT field index
        gt_idx = fmt.split(':').index('GT')

        # build allele lookup
        alleles = [ref] + alt

        # parse per-sample genotypes
        gts = []
        for idx in vcf_sample_idx_lst:
            sample_field = sample_fields[idx]
            gt_str = sample_field.split(':')[gt_idx]
            a_1, a_2 = gt_str.replace('|', '/').split('/')
            if a_1 != '.': a_1 = alleles[int(a_1)]
            if a_2 != '.': a_2 = alleles[int(a_2)]
            gt = gt_code_dct[''.join([a_1, a_2])]
            gts.append(gt)

        # get major gt
        gts_nonmissing = [gt for gt in gts if gt != 0]
        gt_counts = Counter(gts_nonmissing)
        major_gt = gt_counts.most_common(1)[0][0]

        # append data
        # if SNPs, only add variable columns
        if viz_mode == 'snps':
            if len(set(gts_nonmissing)) > 1:
                pos_lst.append(pos)
                gts_lst.append(gts)
                major_gt_lst.append(major_gt)
        # if full add all column, and fill in missing data for missing columns
        if viz_mode == 'full':
            pos_lst.append(pos)
            gts_lst.append(gts)
            major_gt_lst.append(major_gt)

        # update prev_pos
        prev_pos = pos

    # convert GTs to arr and transpose such that POS is on x
    gt_arr = np.array(gts_lst, dtype=int).T

    # assess bcftools process
    return_code = bcftools_process.wait()
    if return_code != 0:
        print(
            f'[ERROR] bcftools failed with exit code {return_code}', 
            file=sys.stderr,
        )
        sys.exit(return_code)

    return gt_arr, major_gt_lst, pos_lst


def plot(
    start,
    end,
    gt_arr,
    major_gt_lst,
    pos_lst,
    sample_lst,
    phenotypes,
    phenotypes_df,
    phenotype_name,
    locations,
    locations_df,
    viz_mode,
    gt_color_dct,
    no_mask_major_gt,
    fig_width,
    fig_height,
    ):

    '''
    Genotype plotting function
    '''

    # for hover annotation, generate a complementary matrix of the decoded GTs
    gt_str_arr = np.vectorize(code_gt_dct.get)(gt_arr)

    # construct hover string array (optionally with phenotypes)
    sample_arr = np.array(sample_lst)[:, np.newaxis]
    pos_arr = np.array([f"{p:,}" for p in pos_lst], dtype=str)[np.newaxis, :]
    hover_matrix = np.char.add(sample_arr, '<br>')
    hover_matrix= np.char.add(hover_matrix, pos_arr)
    hover_matrix = np.char.add(hover_matrix, ": ")
    hover_matrix = np.char.add(hover_matrix, gt_str_arr)
    if phenotypes:
        phenotype_arr = \
            phenotypes_df[phenotype_name].to_numpy().astype(str)[:, np.newaxis]
        hover_matrix = np.char.add(hover_matrix, f'<br>{phenotype_name}: ')
        hover_matrix = np.char.add(hover_matrix, phenotype_arr)

    # convert major GTs to list
    major_gt_arr = np.array(major_gt_lst, dtype=int)

    # mask major GT (if not disabled)
    if not no_mask_major_gt:
        major_gt_mask = gt_arr == major_gt_arr[np.newaxis, :]
        gt_arr[major_gt_mask] = 11

    # construct color map
    # this a list of [value, color] lists which map the 12 gt_color_dct to the
    # 0-1 space
    colormap = []
    for i in range(12):
        colormap.append([i / 12, gt_color_dct[i]])
        colormap.append([(i + 1) / 12, gt_color_dct[i]])

    # initiate figure
    fig = go.Figure()

    # plot phenotypes on a magma scale (if provided)
    if phenotypes:

        # extract phenotypes
        phenotype_arr = phenotypes_df[phenotype_name].to_numpy()

        # create 1D array
        phenotype_matrix = phenotype_arr.reshape(-1, 1)

        # generate hover data
        pheno_hover = [
            [f'{sample_id}: {pheno_val}']
            for sample_id, pheno_val in zip(sample_lst, phenotype_arr)
        ]

        # plot
        fig.add_trace(
            go.Heatmap(
                z=phenotype_matrix,
                colorscale='magma',
                showscale=False,
                xgap=0, ygap=0,
                xaxis='x2',
                yaxis='y2',
                text=pheno_hover,
                hovertemplate='%{text}<extra></extra>',
            )
        )

    # plot genotype data
    fig.add_trace(
        go.Heatmap(
            z=gt_arr,
            x=list(range(gt_arr.shape[1])) if viz_mode == 'snps' else pos_lst,
            zmin=0,
            zmax=12,
            colorscale=colormap,
            showscale=False,
            xgap=0, ygap=0,
            xaxis='x', yaxis='y',
            customdata=pos_lst,
            text=hover_matrix,
            hovertemplate='%{text}<extra></extra>',
        )
    )

    # plot locations on top of plot (if provided)
    if locations:

        # iterrate locations
        for _, row in locations_df.iterrows():

            # extract info
            region_start = int(row['from'])
            region_end = int(row['to'])
            region_spread = float(row['spread'])

            # process spread
            height_fraction = region_spread / 100.0
            y_0_shape = 0.5 - (height_fraction / 2.0)
            y_1_shape = 0.5 + (height_fraction / 2.0)

            # viz_mode == 'full': use POS
            if viz_mode == 'full':
                x_0 = region_start
                x_1 = region_end

            # viz_mode == 'snps': map to idx
            else:

                # find left and right boundaries
                x_0 = 0
                for i, pos in enumerate(pos_lst):
                    if int(pos) >= region_start:
                        x_0 = i
                        break
                else:
                    x_0 = len(pos_lst) - 1
                x_1 = 0
                for i, pos in enumerate(pos_lst):
                    if int(pos) >= region_end:
                        x_1 = i
                        break
                else:
                    x_1 = len(pos_lst) - 1

            # pad for single POS or missing POS in snp mode
            if x_0 == x_1:
                x_0 -= 0.5
                x_1 += 0.5

            # plot
            fig.add_shape(
                type='rect',
                x0=x_0,
                y0=y_0_shape,
                x1=x_1,
                y1=y_1_shape,
                fillcolor=row['color'],
                line_width=1 if (x_0 == x_1 - 1 and viz_mode == 'full') else 0,
                line_color=row['color'],
                layer='above',
                xref='x',
                yref='y3',
                opacity=1,
            )

    # format
    fig.update_layout(
        template='simple_white',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        font_family='Arial',
        font_color='black',
        width=fig_width,
        height=fig_height,
        xaxis={
            'domain': [0.02, 1] if phenotypes else [0, 1],
            'range': [start, end] if viz_mode == 'full' else None,
            'autorange': viz_mode != 'full',
            'mirror': True,
            'showline': True,
            'linecolor': 'black',
            'ticks': 'outside',
            'title': 'Genomic position'
        },
        yaxis={
            'domain': [0, 0.94] if locations else [0, 1],
            'range': [-0.5, gt_arr.shape[0] - 0.5],
            'autorange': False,
            'mirror': True,
            'showticklabels': False,
            'ticks': '',
            'showline': True,
            'linecolor': 'black',
        },
    )

    if phenotypes:
        fig.update_layout(
            xaxis2={
                'domain': [0, 0.02],
                'showticklabels': False,
                'ticks': '',
                'showline': True,
                'linecolor': 'black',
                'mirror': True,
            },
            yaxis2={
                'matches': 'y',
                'domain': [0, 0.94] if locations else [0, 1],
                'range': [-0.5, gt_arr.shape[0] - 0.5],
                'anchor': 'x2',
                'autorange': False,
                'showticklabels': False,
                'ticks': '',
                'showline': True,
                'mirror': True,
                'linecolor': 'black',
                'title': f'Samples'
            },
        )

    # locations specific formatting
    if locations:
        fig.update_layout(
            yaxis3={
                'domain': [0.96, 1.0],
                'range': [0, 1],
                'showticklabels': False,
                'showline': False,
                'ticks': '',
                'zeroline': False
            }
        )

    # format x ticks
    if viz_mode == 'snps':

        # get list of tick indices
        tick_idx_lst = np.linspace(
            0,
            len(pos_lst) - 1,
            SNPS_MODE_N_TICKS,
            dtype=int,
        ).tolist()

        # get the corresponding display values
        tick_val_lst = [f'{pos_lst[i]:,}' for i in tick_idx_lst]

        # plot
        fig.update_layout(
            xaxis={
                'tickmode': 'array',
                'tickvals': tick_idx_lst,
                'ticktext': tick_val_lst,
                'showticklabels': True,
                'ticks': 'outside'
            }
        )

    # tick format
    fig.update_xaxes(
        tickformat='d' if viz_mode == 'full' else '',
        hoverformat=',d',
    )

    return fig


def save_fig(fig, format_string, prefix):

    '''
    Save figure in specified formats.
    '''

    for fmt in format_string.split(','):
        if fmt in ['html', 'HTML']:
            fig.write_html(f'{prefix}.html')
        else:
            try:
                fig.write_image(f'{prefix}.{fmt.lower()}')
            except:
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
    vcf_path = args.vcf_path
    region = args.region
    sample_id_path =  args.sample_id_path
    ref_fasta_path = args.ref_fasta_path
    output_prefix = args.output_prefix
    viz_mode = args.viz_mode
    color_mode = args.color_mode
    no_mask_major_gt = args.no_mask_major_gt
    phenotypes = args.phenotypes
    phenotype_name = args.phenotype_name
    annotate = args.annotate
    rev_sample_order = args.rev_sample_order
    locations = args.locations
    plot_fmt = args.plot_fmt
    fig_height = int(args.fig_height)
    fig_width = int(args.fig_width)

    # parse region
    chrom = region.split(':')[0]
    start = int(region.split(':')[1].split('-')[0])
    end = int(region.split(':')[1].split('-')[1])
    if start > end:
        start, end = end, start

    # handle unspecified arguments
    if not phenotypes:
        phenotypes_df = False
    if not locations:
        locations_df = False

    # read samples
    sample_lst = pd.read_csv(
        sample_id_path,
        sep='\t',
        header=None,
    )[0].to_list()

    # read phenotypes and subset to sample list and target phenotype
    # if specified
    if phenotypes:

        # print warning if no column name is specified
        if not phenotype_name:
            print(
                ' [ERROR] Please provide a phenotype column name in ' \
                f' {phenotypes} through -n/--phenotype_name',
                file=sys.stderr,
            )
            sys.exit(1)

        # read data
        phenotypes_df, sample_lst = read_phenotype_data(
            phenotypes,
            phenotype_name,
            sample_lst,
            annotate,
            rev_sample_order,
        )

    # read locations data if specified
    if locations:

        # read data
        locations, locations_df = read_locations_data(
            locations,
            chrom,
            start,
            end,
            LOCATIONS_DEFAULT_COLOR,
            LOCATIONS_DEFAULT_SPREAD,
            )

    # fetch samtools from execution $PATH
    try:
        BCFTOOLS = shutil.which('bcftools')
    except:
        print(
            '[ERROR] bcftools not found',
            file=sys.stderr,
            )
        sys.exit(1)

    # run bcftools
    bcftools_process = run_bcftools(
        BCFTOOLS,
        vcf_path,
        chrom,
        start,
        end,
        sample_id_path,
        ref_fasta_path,
    )

    # parse VCF data
    gt_arr, major_gt_lst, pos_lst = parse_vcf_data(
        bcftools_process,
        sample_lst,
        viz_mode,
    )

    # select color dict
    if color_mode == 'multi':
        gt_color_dct = multi_color_dct
    elif color_mode == 'bw':
        gt_color_dct = bw_color_dct
    elif color_mode == 'rw':
        gt_color_dct = rw_color_dct

    # plot genotypes
    fig = plot(
        start,
        end,
        gt_arr,
        major_gt_lst,
        pos_lst,
        sample_lst,
        phenotypes,
        phenotypes_df,
        phenotype_name,
        locations,
        locations_df,
        viz_mode,
        gt_color_dct,
        no_mask_major_gt,
        fig_width,
        fig_height,
    )

    # save figure
    save_fig(fig, plot_fmt, output_prefix)



## EXECUTE

if __name__ == "__main__":
    main()
