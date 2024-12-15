#!/usr/bin/env python3

"""
Perform/plot PCA.
"""

## IMPORT PACKAGES

import sys
import os
import argparse
import numpy as np
import pandas as pd
import gzip
import allel
import plotly.graph_objects as go

############################################
#  set default values in self.default_dct  #
############################################



## CLASSES

class CLI:
    '''
    Command line interface and argument parser.
    '''

    def __init__(self):

        # initiate argument parser
        self.parser = argparse.ArgumentParser(
            description='PCA v1.0',
            epilog='contact: lmb215@cam.ac.uk',
            formatter_class=argparse.RawDescriptionHelpFormatter
        )

        # add -v/--version
        self.parser.add_argument(
            '-v', '--version',
            action='version',
            version=f'PCA 1.0'
        )

        # initiate subparsers
        self.subparsers = self.parser.add_subparsers(
            dest='pca', required=True
        )

        # variables
        self.args_dct = None

        # allowed variant file suffixes
        self.variant_file_suffixes = [
            '.vcf', '.vcf.gz',
        ]

        self.plot_file_suffixes = [
            '.html', '.pdf', '.svg', '.png',
        ]

        self.default_dct = {
            'skip_monomorphic': True,
            'vcf_pass_filter': True,
            'gt_mean_impute': True,
            'plot_w': 1000,
            'plot_h': 800,
        }


    @staticmethod
    def shared_arguments(subparser):
        '''
        Arguments shared between all sub-commands.
        '''
        # define arguments that are shared across all subparsers
        subparser.add_argument(
            dest='prefix', metavar='<PREFIX>',
            help='Prefix for this run, i.e. prefix for all results generated in'
            ' this PCA analysis.')


    def pca(self):
        '''
        PCA on called genotypes(GT) using scikit-allel.
        '''

        # add subparser
        pca_parser = self.subparsers.add_parser(
            'pca', help='Perform windowed PCA on called genotypes (GT) or on'
            ' genotype likelihoods (GL/PL).'
        )

        # positional arguments
        pca_parser.add_argument(
            dest='variant_file_path', metavar='<VARIANT_FILE>', help='Path to'
            ' variant file (optionally gzipped VCF, TSV or BEAGLE; see'
            ' documentation for input file specifications).')
        self.shared_arguments(pca_parser)

        # optional arguments
        pca_parser.add_argument(
            '-m', '--min_maf', dest='min_maf', required=False, type=float,
            default=0.01, metavar='\b', help='Minor allele frequency'
            f' threshold [default: 0.01].')
        pca_parser.add_argument(
            '-n', '--n_pcs', dest='n_pcs', required=False, type=int,
            default=10, metavar='\b', help='Number of principal components'
            ' to compute [default: 10].')
        pca_parser.add_argument(
            '-r','--region', dest='region', metavar='\b', required=False,
            type=str, default=None, help='Genomic region in format '
            '"chrom:start-end".')
        pca_parser.add_argument(
            '-s', '--samples', dest='samples', required=False, metavar='\b',
            help='''Comma-separated list of samples to include or file with one'
            ' sample per line.''')


    def plot(self):
        '''
        Plot PC1, PC2, heterozygosity and per window stats for a specified input
        chromosome.
        '''

        # add subparser
        plot_parser = self.subparsers.add_parser(
            'plot', help='Plot PCA.'
        )

        # positional arguments
        self.shared_arguments(plot_parser)

        # optional arguments
        plot_parser.add_argument(
            '-p', '--plot_pcs', dest='plot_pcs', required=False,
            default='1,2', metavar='\b', help='Specify which PCs to plot'
            ' (comma-separated list) [default: 1,2].')
        plot_parser.add_argument(
            '-m', '--metadata', dest='metadata_path', required=False,
            metavar='\b', help='Path to metadata TSV where first column are'
            ' sample IDs. Additional columns will be used to annotate data in'
            ' HTML plot.')
        plot_parser.add_argument(
            '-g', '--groups', dest='color_by', required=False, metavar='\b',
            help='Metadata column for color-grouping. Requires -m/--metadata.')
        plot_parser.add_argument(
            '-c', '--colors', dest='hex_codes', required=False, metavar='\b',
            help='HEX codes (drop "#") for color groups. Check documentation'
            ' for formatting instructions. Requires -g/--groups.')
        plot_parser.add_argument(
            '-i', '--invert', dest='invert', required=False, metavar='\b',
            default=None, help='Specify whether to invert the values along the'
            ' x ("x") or y axis ("y") or both ("xy").')
        plot_parser.add_argument(
            '-f', '--format', dest='plot_fmt', required=False, metavar='\b',
            default='HTML', help='Output plot file format ("HTML",'
            ' "PDF", "SVG" or "PNG") [default: {HTML}].')
        plot_parser.add_argument(
            '--n', '--numeric', dest='numeric', required=False, 
            action='store_true', help='Set flag to apply a continuous color'
            ' scale (requires numerical data).')
        plot_parser.add_argument(
            '--r', '--reverse', dest='reverse', required=False, 
            action='store_true', help='Set flag to reverse the plotting'
            ' order of specified group (-g/--groups).')


    def parse_args(self):
        '''
        Parse and sanity check command-line arguments.
        '''

        # parse arguments
        args = self.parser.parse_args()

        # handle interdependent options
        if hasattr(args, 'color_by') and not hasattr(args, 'metadata_path'):
            self.parser.error(
                '-m/--metadata is required to infer -g/--groups.')
        if hasattr(args, 'hex_codes') and not hasattr(args, 'color_by'):
            self.parser.error(
                '-g/--groups is required if -c/--colors is set.')

        # check formatting
        if hasattr(args, 'variant_file_path'):
            if args.variant_file_path.lower().endswith('.gz'):
                suffix = '.' + '.'.join(args.variant_file_path.split('.')[-2:])
            else:
                suffix = '.' + '.'.join(args.variant_file_path.split('.')[-1:])
            if suffix.lower() not in self.variant_file_suffixes:
                self.parser.error(
                    f'{args.variant_file_path} has an unexpected file format'
                    + f' based on its suffix "{suffix}")')
        if hasattr(args, 'region') and not args.region is None:                                           
            if ':' not in args.region or '-' not in args.region:
                self.parser.error(
                    f'{args.region} is formatted incorrectly. Please make sure'
                    + ' to specify genomic coordinates correctly:'
                    + ' chrom:start-end (start and end must always be'
                    + ' specified).' )
        if hasattr(args, 'invert'):
            if args.invert not in ['x', 'y', 'xy', 'yx', None]:
                self.parser.error(
                    '-i/--invert must be one of "x", "y" or "xy"')

        # decompose/preprocess arguments
        if hasattr(args, 'plot_pcs'):
            args.plot_pcs = [f'PC{x}' for x in args.plot_pcs.split(',')]
        if hasattr(args, 'region'):
            if args.region is None:
                chrom=None
                start=None
                end=None
            else:
                chrom = args.region.split(':')[0]
                start = int(args.region.split(':')[1].split('-')[0])
                end = int(args.region.split(':')[1].split('-')[1])
        if hasattr(args, 'samples') and not args.samples is None:
            if ',' in args.samples:
                sample_lst = args.samples.split(',')
            else:
                if not os.path.exists(args.samples):
                    self.parser.error(
                        args.samples + ': file does not exist.')
                else:
                    sample_lst = []
                    with open(args.samples, 'r') as sample_file:
                        for line in sample_file:
                            sample_lst.append(line.strip().split('\t')[0])
        if hasattr(args, 'flip_windows') and not args.flip_windows is None:
            if ',' in args.flip_windows:
                flip_window_lst = args.flip_windows.split(',')
            else:
                if os.path.exists(args.flip_windows):
                    flip_window_lst = []
                    with open(args.flip_windows, 'r') as flip_windows_file:
                        for line in flip_windows_file:
                            flip_window_lst.append(line.strip().split('\t')[0])
                else:
                    flip_window_lst = [args.flip_windows]
        if hasattr(args, 'hex_codes'):
            if args.hex_codes:
                hex_code_dct = {}
                for i in args.hex_codes.split(','):
                    group, hex_code = i.split(':')[0], i.split(':')[1]
                    hex_code_dct[group] = '#' + hex_code
            else:
                hex_code_dct = None
        if hasattr(args, 'plot_fmt'):
            if ',' in args.plot_fmt:
                plot_fmt_lst = args.plot_fmt.split(',')
            else:
                plot_fmt_lst = [args.plot_fmt]
            plot_fmt_lst = [x.lower() for x in plot_fmt_lst]

        # further checks
        if hasattr(args, 'plot_fmt'):
            for fmt in plot_fmt_lst:
                if '.' + fmt not in self.plot_file_suffixes:
                    self.parser.error(
                        f'"{fmt}" is not supported as output format.')

        # check if files exist
        if hasattr(args, 'variant_file_path'):
            if not os.path.exists(args.variant_file_path):
                self.parser.error(
                    args.variant_file_path + ': file does not exist.')
        if hasattr(args, 'metadata_path') and args.metadata_path:
            if not os.path.exists(args.metadata_path):
                self.parser.error(
                    args.metadata_path + ': file does not exist.')

        # convert args to dict and add to class as instance variable
        self.args_dct = vars(args)

        # add derived/processed arguments
        if hasattr(args, 'region'):
            self.args_dct['chrom'] = chrom
            self.args_dct['start'] = start
            self.args_dct['end'] = end
        if hasattr(args, 'samples') and not args.samples is None:
            self.args_dct['sample_lst'] = sample_lst
        if hasattr(args, 'hex_codes'):
            self.args_dct['hex_code_dct'] = hex_code_dct
        if hasattr(args, 'plot_fmt'):
            self.args_dct['plot_fmt_lst'] = plot_fmt_lst
        if hasattr(args, 'flip_windows') and not args.flip_windows is None:
            self.args_dct['flip_window_lst'] = flip_window_lst
        if hasattr(args, 'plot_pcs'):
            self.args_dct['plot_pcs'] = args.plot_pcs
        if hasattr(args, 'n_pcs'):
            self.args_dct['n_pcs'] = args.n_pcs
        if hasattr(args, 'invert'):
            self.args_dct['invert'] = args.invert
        
        # add in settings from config
        self.args_dct['skip_monomorphic'] = self.default_dct['skip_monomorphic']
        self.args_dct['vcf_pass_filter'] = self.default_dct['vcf_pass_filter']
        self.args_dct['gt_mean_impute'] = self.default_dct['gt_mean_impute']
        self.args_dct['plot_w'] = self.default_dct['plot_w']
        self.args_dct['plot_h'] = self.default_dct['plot_h']

        # add sample_lst if not unset
        if not 'sample_lst' in self.args_dct:
            self.args_dct['sample_lst'] = None



class Log():

    def __init__(self):
        pass

    def newline(self):
        '''
        Print ERROR message to STDERR and exit.
        '''

        print('\n',
              file=sys.stderr,
              flush=True,
              )


    def info(self, message):
        '''
        Print INFO message to STDERR.
        '''

        print(f'[INFO] {message}.',
              file=sys.stderr,
              flush=True,
              )


    def error(self, message):
        '''
        Print ERROR message to STDERR and exit.
        '''

        print(f'[ERROR] {message}.',
              file=sys.stderr,
              flush=True,
              )
        sys.exit(1)

# instantiate logger
log = Log()

class PCA:

    '''
    Parse hard-called genotypes (GT) from VCF and apply a function.
    '''

    def __init__(self,
                 prefix,
                 variant_file_path,
                 file_fmt,
                 sample_lst,
                 n_pcs,
                 chrom, start, stop,
                 skip_monomorphic,
                 gt_mean_impute,
                 vcf_pass_filter,
                 min_maf,
                 ):

        # variant file/region
        self.prefix = prefix
        self.variant_file_path = variant_file_path
        self.file_fmt = file_fmt
        self.sample_lst = sample_lst
        self.n_pcs = n_pcs
        self.chrom = chrom
        self.start = start
        self.stop = stop

        # parameters
        self.skip_monomorphic = skip_monomorphic
        self.gt_mean_impute = gt_mean_impute
        self.vcf_pass_filter = vcf_pass_filter
        self.min_maf = min_maf

        # transient variables
        self.variants = None
        self.variant_arr = None
        self.pca_out = None

        # genotype encoding
        self.gt_code_dct = {
            '0/0':  0,
            '0|0':  0,
            '0/1':  1,
            '0|1':  1,
            '1/0':  1,
            '1|0':  1,
            '1/1':  2,
            '1|1':  2,
            './.': np.nan,
            '.|.': np.nan,
            '0/.': np.nan,
            '0|.': np.nan,
            './0': np.nan,
            '.|0': np.nan,
            '1/.': np.nan,
            '1|.': np.nan,
            './1': np.nan,
            '.|1': np.nan,
            '.':   np.nan,
        }


    def gt_mean_imputation(self):
        '''
        Mean impute missing GT calls.
        '''

        # calculate row/site means (ignoring NaN)
        row_mean_arr = np.nanmean(self.variant_arr, axis=1)

        # fetch indices of rows/sites WITH missing call(s)
        miss_mask_arr = np.isnan(self.variant_arr)

        # replace NaN with row/site means
        self.variant_arr[miss_mask_arr] = np.take(
            row_mean_arr, np.where(miss_mask_arr)[0]
        )


    def gt_drop_missing_sites(self):
        '''
        remove any site with at least one missing GT call.
        '''

        # fetch indices of rows/sites WITHOUT missing call(s)
        mask = ~np.isnan(self.variant_arr).any(axis=1)

        # drop those rows/sites
        self.variant_arr = self.variant_arr[mask]


    def gt_min_maf_filter(self):
        '''
        Drop SNPs with minor allele frequency below specified value.
        '''

        # count # called GTs per row (np.sum counts True), *2 because diploid
        allele_counts_arr = np.sum(~np.isnan(self.variant_arr), axis=1) * 2

        # count # alleles per row/site
        allele_sums_arr = np.nansum(self.variant_arr, axis=1)

        # calculate allel frequencies and multiple with -1 if AF > 0.5 (because
        # input data may not be polarized by major/minor allel)
        af_arr = allele_sums_arr / allele_counts_arr
        af_arr[af_arr > 0.5] = 1 - af_arr[af_arr > 0.5]

        # keep only sites where AF >= min_maf
        self.variant_arr = self.variant_arr[af_arr >= self.min_maf]


    def pca(self):
        '''
        Apply min_maf filter, drop rows/sites with missing call(s) OR mean
        impute, conduct PCA, but if (n_var < min_var_per_w) generate empty/dummy
        output instead.
        '''

        # min_maf filter
        if self.min_maf:
            self.gt_min_maf_filter()

        # count variants
        n_var = self.variant_arr.shape[0]

        # mean impute
        if self.gt_mean_impute:

            # mean impute
            self.gt_mean_imputation()

        # drop missing sites (& re-count)
        else:

            # drop missing sites
            self.gt_drop_missing_sites()

            # re-count count remaining variants
            n_var = self.variant_arr.shape[0]

        # if # variants passes specified threshold
        if n_var == 0:
            log.newline()
            log.error('No variants found.')

        # print info
        log.newline()
        log.info('Performing PCA ...')

        # perform PCA
        self.pca_out = allel.pca(
            self.variant_arr,
            n_components=self.n_pcs,
            copy=True,
            scaler='patterson',
            ploidy=2,
        )


    def variant_parser(self):
        '''
        Apply a target function to windows of variants (GT field) in an
        (optionally gzipped) VCF file.
        '''

        # open iput file
        read_func = gzip.open if self.variant_file_path.endswith('.gz') else open
        with read_func(self.variant_file_path, 'rt') as variant_file:

            # fetch sample ids from different header types
            for line in variant_file:
                if line.startswith('#CHROM'):
                    var_file_sample_lst = line.strip().split('\t')[9:]
                    break
        
        # use var_file_sample_lst (drop duplicates) if no samples specified
        if self.sample_lst is None:
            self.sample_lst = list(dict.fromkeys(var_file_sample_lst))

        # obtain index positions (returns first hit)
        sample_idx_lst = [
            var_file_sample_lst.index(x) for x in self.sample_lst
        ]

        # initiate variant list
        self.variants = []

        # traverse VCF
        with read_func(self.variant_file_path, 'rt') as variant_file:

            # if a region was specified
            if self.chrom:
                for line in variant_file:
                    if line.startswith('#'): continue
                    line = line.strip().split('\t')
                    q_chrom = line[0]
                    if q_chrom != self.chrom: continue
                    filter_field = line[6]
                    if filter_field != 'PASS' and self.vcf_pass_filter:
                        continue
                    pos = int(line[1])
                    if pos < self.start: continue
                    if pos > self.stop: break
                    gts = [line[9:][idx].split(':')[0] for idx in sample_idx_lst]
                    gts = [self.gt_code_dct[x] for x in gts]
                    if self.skip_monomorphic and len(set(gts)) == 1:
                        continue
                    self.variants.append([pos] + gts)

            # else use all variants in the file that pass filters
            else:
                for line in variant_file:
                    if line.startswith('#'): continue
                    line = line.strip().split('\t')
                    filter_field = line[6]
                    if filter_field != 'PASS' and self.vcf_pass_filter:
                        continue
                    pos = int(line[1])
                    gts = [line[9:][idx].split(':')[0] for idx in sample_idx_lst]
                    gts = [self.gt_code_dct[x] for x in gts]
                    if self.skip_monomorphic and len(set(gts)) == 1:
                        continue
                    self.variants.append([pos] + gts)
                    
            self.variants = [x for x in self.variants]
            self.variant_arr = np.array([x[1:] for x in self.variants],dtype=np.float32)
            self.pca()

        # print message
        log.newline()
        log.info('Writing output to files')
        log.newline()

        # write output
        pc_df = pd.DataFrame(self.pca_out[0])
        pc_df.columns = [f'PC{x}' for x in range(1, self.n_pcs+1)]
        pc_df['id'] = self.sample_lst
        pc_df = pc_df[['id'] + [f'PC{x}' for x in range(1, self.n_pcs+1)]]
        pc_df.to_csv(
            f'{self.prefix}.pc.tsv.gz',
            sep='\t',
            index=False, 
            compression='gzip'
        )
        pct_exp_df = pd.DataFrame(
            [
                [
                    round(self.pca_out[1].explained_variance_ratio_[x-1]*100, 2) \
                        for x in range(1, self.n_pcs+1)]
                ]
        )
        pct_exp_df.columns = pc_df.columns[1:]
        pct_exp_df.to_csv(
            f'{self.prefix}.ve.tsv.gz',
            sep='\t',
            index=False,
            compression='gzip',
        )



class Plot:
    '''
    Plot windowed PC and associated data.
    '''

    def __init__(self,
                 plot_pcs=None,
                 prefix=None,
                 metadata_path=None,
                 color_by=None,
                 hex_code_dct=None,
                 numeric=None,
                 reverse=None,
                 plot_w=None,
                 plot_h=None,
                 plot_fmt_lst=None,
                 invert=None,
                 ):

        # input variables
        self.plot_pcs = plot_pcs
        self.prefix = prefix
        self.metadata_path = metadata_path
        self.invert = invert
        self.color_by = color_by
        self.hex_code_dct = hex_code_dct
        self.numeric = numeric
        self.reverse = reverse
        self.plot_w = plot_w
        self.plot_h = plot_h
        self.plot_fmt_lst = plot_fmt_lst

        # instance variables
        self.pc_df = pd.DataFrame()
        self.ve_df  = pd.DataFrame()
        self.metadata_df = pd.DataFrame()
        self.group_id = None
        self.group_lst = None
        self.color_dct = None
        self.fig = None

        # define custom color scale
        colors = ['#000080', '#0000FF', '#800080', '#FF0000', '#FFFF00']
        thresholds = [0, 0.25, 0.5, 0.75, 1]

        # create the colorscale
        colorscale = [[threshold, color] for threshold, color in zip(thresholds, colors)]
                     
        # set color scheme
        self.color_scale = colorscale # or set inbuilt schemes like 'Plasma'

        # set allowed NA strings
        self.na_lst = [None, 'NA', 'na', 'NaN']


    def annotate(self):
        '''
        Annotate per-sample data with a metadata file if supplied, otherwise
        just reformat for plotting function.
        '''

        # fetch sample names and order (=data_df column names)
        sample_lst = list(self.pc_df['id'])

        # read metadata if provided, do sanity checks and use for annotation
        if self.metadata_path:

            # read metadata and print error message if there are non-unique IDs
            self.metadata_df = pd.read_csv(
                self.metadata_path, sep='\t', index_col=0, dtype=str
            )
            if len(self.metadata_df.index) != len(set(self.metadata_df.index)):
                log.newline()
                log.error('The provided metadata file contains non-unique'
                          ' sample IDs')
                log.newline()
            # subset and reorder metadata_df to match data_df individuals
            self.metadata_df = self.metadata_df.loc[
                self.metadata_df.index.intersection(sample_lst)
            ].reindex(sample_lst)

            # if individuals are missing in the metadata print error message
            if len(self.metadata_df) != len(sample_lst):
                log.newline()
                log.error('One or more sample IDs are missing in the'
                          ' provided metadata file')
                log.newline()

            # add metadata columns to data_df
            for column_name in self.metadata_df.columns:
                self.pc_df[column_name] = list(self.metadata_df[column_name])


    def set_colors(self):
        '''
        Parse per-sample plot color specifications and compile color_dct.
        '''

        # fetch group_id (=color_by) if specified, else default to 'id'
        self.group_id = self.color_by if not self.color_by is None else 'id'

        # get list of groups/ids
        self.group_lst = list(set(self.pc_df[self.group_id]))

        # set color_dct for continuous colors
        if self.numeric:
            import plotly.colors as p_colors
            # replace None (pandas) with np.nan (numpy), else through an error
            try:
                val_lst = np.array(
                    [np.nan if x in self.na_lst else x for x in self.group_lst],
                    dtype=np.float64
                )
            except:
                log.newline()
                log.error('Provided column (-g/--groups) contains non-numerical'
                          ' values')
                log.newline()
            # scale to 0-1
            norm_val_lst = (
                val_lst-np.nanmin(val_lst))/\
                (np.nanmax(val_lst) - np.nanmin(val_lst)
            )
            # compile color_dct with all values as keys, using again None
            # instead of np.nan
            self.color_dct = {}
            for val, norm_val in zip(self.group_lst, norm_val_lst):
                if val in self.na_lst:
                    self.color_dct[val] = 'lightgrey'
                else:
                    self.color_dct[val] = p_colors.sample_colorscale(
                        self.color_scale, [norm_val],
                    )[0]
            # adjust plotting order
            self.group_lst = \
                [
                    x for x in self.group_lst \
                        if self.na_lst
                ] \
                + sorted(
                [
                        x for x in self.group_lst \
                            if x not in self.na_lst
                ])

        # define colors based on plotly default colors or specified HEX codes;
        # print error messages if HEX codes are missing for specified groups
        elif self.hex_code_dct:
            self.color_dct = self.hex_code_dct
            if not all(x in self.color_dct.keys() for x in self.group_lst):
                log.newline()
                log.error('HEX codes missing for one or more groups')
                log.newline()
            else:
                # set color_dct keys as group_lst to set plotting order
                self.group_lst = list(self.color_dct.keys())

        # use defaultcolors
        else:
            import plotly.colors as p_colors
            def_col_lst = p_colors.DEFAULT_PLOTLY_COLORS
            self.color_dct = {
                self.group_lst[idx]: def_col_lst[idx % len(def_col_lst)] \
                    for idx in range(len(self.group_lst))
            }

        # reverse plotting order if specified
        if self.reverse:
            self.group_lst.reverse()

    @staticmethod
    def derive_ticks(values):
        data_range = max(values) - min(values)
        raw_step = data_range / 7
        magnitude = 10 ** int(np.floor(np.log10(raw_step)))
        step = max(5, round(raw_step / magnitude) * magnitude)
        if step % 10 != 0:
            step = int((step // 10 + 1) * 10)
        start_tick = int(step * np.floor(min(values) / step))
        end_tick = int(step * np.ceil(max(values) / step))
        ticks = list(range(start_tick, end_tick + step, step))
        while len(ticks) > 7:
            step *= 2
            start_tick = int(step * np.floor(min(values) / step))
            end_tick = int(step * np.ceil(max(values) / step))
            ticks = list(range(start_tick, end_tick + step, step))
        return ticks


    def savefig(self):
        '''
        Save figure in HTML and/or other (PDF, SVG, PNG) format(s).
        '''
        for fmt in self.plot_fmt_lst:
            title = f'{self.plot_pcs[0]}-{self.plot_pcs[0].replace("PC", "")}'
            if fmt == 'html':
                self.fig.write_html(f'{self.prefix}.{title}.{fmt}')
            else:
                self.fig.write_image(f'{self.prefix}.{title}.{fmt}')


    def plot(self):
        '''
        Plot per-sample values for one chromosome (e.g. PC 1) with small panel
        of per window values (e.g. PC 1 variance explained) on top.
        '''

        # LOAD & PREPARE DATA

        # load data
        self.pc_df = pd.read_csv(
            f'{self.prefix}.pc.tsv.gz', sep='\t'
        )[['id', self.plot_pcs[0], self.plot_pcs[1]]]
        self.ve_df = pd.read_csv(
            f'{self.prefix}.ve.tsv.gz', sep='\t'
        )

        # reflect one or both plot axes if specified
        if self.invert:
            if self.invert == 'x':
                 self.pc_df[self.plot_pcs[0]] =  self.pc_df[self.plot_pcs[0]]*-1
            if self.invert == 'y':
                 self.pc_df[self.plot_pcs[1]] =  self.pc_df[self.plot_pcs[1]]*-1
            if self.invert in ['xy', 'yx']:
                 self.pc_df[self.plot_pcs[0]] =  self.pc_df[self.plot_pcs[0]]*-1
                 self.pc_df[self.plot_pcs[1]] =  self.pc_df[self.plot_pcs[1]]*-1

        # annotate per-sample data if metadata were supplied
        self.annotate()

        # set per-sample plot colors
        self.set_colors()

        # figure setup
        self.fig = go.Figure()
    
        # compile hover data strings
        self.pc_df['hover_label'] = 'NA'

        # iterate through each individual
        for sample in set(self.pc_df['id']):

            # subset data
            sample_df = self.pc_df[self.pc_df['id'] == sample]

            # compile hover string
            hover_str = \
                ''.join(
                    [
                        f'<b>{c}</b>: {sample_df[c].iloc[0]}<br>' \
                            for c in sample_df.columns
                    ]
                )

            # add to pc_df
            self.pc_df.loc[self.pc_df['id'] == sample, 'hover_label'] = hover_str

        # derive scaling factor (using VE range and PC data ranges) and scale
        # the second specified PC
        ve_a = self.ve_df[self.plot_pcs[0]].iloc[0]
        ve_b = self.ve_df[self.plot_pcs[1]].iloc[0]
        ve_ratio = ve_b/ve_a
        pc_a_data = self.pc_df[f'{self.plot_pcs[0]}']
        pc_b_data = self.pc_df[f'{self.plot_pcs[1]}']
        pc_a_range = max(pc_a_data) - min(pc_a_data)
        pc_b_range = max(pc_b_data) - min(pc_b_data)
        range_ratio = pc_a_range / pc_b_range
        scale_factor = ve_ratio * range_ratio

        print(ve_ratio)
        print(range_ratio)
        print(scale_factor)
        self.pc_df[f'{self.plot_pcs[1]}_scaled'] = \
            self.pc_df[f'{self.plot_pcs[1]}'] * scale_factor

        # set show_legend to false if using a color scale
        show_legend = not self.numeric

        # plot each specified group (-g) separately (or 'id' if unspecified)
        for group in self.group_lst:

            # subset data to group
            group_df = self.pc_df[self.pc_df[self.group_id] == group].copy()

            # plot
            self.fig.add_trace(
                go.Scatter(
                    x=group_df[self.plot_pcs[0]],
                    y=group_df[f'{self.plot_pcs[1]}_scaled'],
                    hovertext=group_df['hover_label'],
                    hoverinfo='text',
                    hoverlabel=dict(font_size=8),
                    name=group,
                    legendgroup=group,
                    mode='markers',
                    marker=dict(color=self.color_dct[group]),
                    showlegend=show_legend,
                ),
            )
        # general layout
        self.fig.update_layout(
            template='simple_white',
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            font_family='Arial',
            font_color='black',
            width=self.plot_w,
            height=self.plot_h,
            legend=dict(font=dict(size=10)),
            )

        # format x axis
        self.fig.update_xaxes(
            scaleanchor="y",
            linewidth=1,
            side='bottom', 
            mirror=True,
            ticks='outside', tickfont=dict(size=10), tickformat=',.0f',
            title_font=dict(size=12),
            title=dict(
                text=f'<b>{self.plot_pcs[0].replace("PC", "PC ")}</b>'
                    f' ({ve_a} %)' , 
                standoff=1,
            ),
        )

        # format y axis
        self.fig.update_yaxes(
            linewidth=1,
            side='left', 
            mirror=True,
            ticks='outside', tickfont=dict(size=10),
            title_font=dict(size=12),
            title=dict(
                text=f'<b>{self.plot_pcs[1].replace("PC", "PC ")}</b>'
                     f' ({ve_b} %)' , 
                standoff=1,
            ),
        )

        # update axis ranges (± 5%) and y ticks/marks
        x_vals = self.pc_df[f'{self.plot_pcs[0]}']
        min_x = min(x_vals)-(0.05*(max(x_vals)-min(x_vals)))
        max_x = max(x_vals)+(0.05*(max(x_vals)-min(x_vals)))
        y_scaled_vals = self.pc_df[f'{self.plot_pcs[1]}_scaled']
        min_y = min(y_scaled_vals)-(0.05*(max(y_scaled_vals)-min(y_scaled_vals)))
        max_y = max(y_scaled_vals)+(0.05*(max(y_scaled_vals)-min(y_scaled_vals)))
        self.fig.update_xaxes(
            range=(min_x, max_x),
            constrain='domain'
        )

        # derive sensible ticks for pseudo-scaled axis
        y_tick_vals = self.derive_ticks(self.pc_df[f'{self.plot_pcs[1]}'])

        # apply to y axis
        self.fig.update_yaxes(
            range=(min_y, max_y),
            constrain='domain',
            tickmode='array',
            tickvals=[
                v * scale_factor for v in y_tick_vals
            ],
            ticktext=y_tick_vals,
        )

        # plot colorscale instead of per-sample legend for numeric metadata
        if self.numeric:
            val_lst = [
                float(x) for x in self.data_df[self.color_by] \
                    if not x in self.na_lst
            ]
            min_val = min(val_lst)
            max_val = max(val_lst)
            self.fig.update_layout(
                coloraxis=dict(
                    colorscale=self.color_scale,
                    cmin=min_val,
                    cmax=max_val,
                    colorbar=dict(
                        len=0.7,
                        thickness=10,
                        title=dict(
                            text=self.color_by,
                            font=dict(size=100),
                            side='right'
                        ),
                        tickvals=[min_val, max_val],
                        tickfont=dict(size=100),
                    ),
                ),
            )
            self.fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(
                    color=[0],
                    coloraxis='coloraxis'
                ),
                showlegend=False,
            ),
        )



## MAIN

def main():
    '''
    Execute
    '''

    # Parse CLI

    # instantiate, call subparsers & parse
    cli = CLI()
    cli.pca()
    cli.plot()
    cli.parse_args()
    args_dct = cli.args_dct

    # set mode
    mode = args_dct['pca']

    # mode: PCA

    if mode == 'pca':

        if args_dct['variant_file_path'].endswith('.gz'):
            file_fmt = args_dct['variant_file_path'].split('.')[-2].upper()
        else:
            file_fmt = args_dct['variant_file_path'].split('.')[-1].upper()

        pca = PCA(
            variant_file_path = args_dct['variant_file_path'],
            prefix=args_dct['prefix'],
            file_fmt=file_fmt,
            n_pcs=args_dct['n_pcs'],
            sample_lst=args_dct['sample_lst'],
            chrom=args_dct['chrom'],
            start=args_dct['start'],
            stop=args_dct['end'],
            skip_monomorphic=args_dct['skip_monomorphic'],
            gt_mean_impute=args_dct['gt_mean_impute'],
            vcf_pass_filter=args_dct['vcf_pass_filter'],
            min_maf=args_dct['min_maf'],
            )
        pca.variant_parser()

    # mode: plot

    if mode == 'plot':

        # print info
        log.newline()
        log.info('Creating PCA plot')

        # load module & instantiate
        #from modules.plot import Plot
        plot = Plot(plot_pcs=args_dct['plot_pcs'],
                    prefix=args_dct['prefix'],
                    metadata_path=args_dct['metadata_path'],
                    color_by=args_dct['color_by'],
                    hex_code_dct=args_dct['hex_code_dct'],
                    plot_w=args_dct['plot_w'],
                    plot_h=args_dct['plot_h'],
                    plot_fmt_lst=args_dct['plot_fmt_lst'],
                    numeric=args_dct['numeric'],
                    reverse=args_dct['reverse'],
                    invert=args_dct['invert']
        )
        plot.plot()
        plot.savefig()

    # print exit message
    log.info('Done')
    log.newline()


## EXECUTE
if __name__ == "__main__":
    main()


## JUPYTER NOTEBOOK USAGE
# import sys
# emulate command line
# command_line = 'INSERT COMMAND LINE'# -r chr1:1-10000000'
# command_line = 'pca.py plot TEST -m ../data/species_rep_metadata.tsv -g clade -c Mbuna:a120ed,AstCal:a2cd5b,Rhampho:8a4513,Diplo:ffa34e,Shallow:ff6247,Deep:4876ff,Utaka:006400'
# command_line = command_line.strip()
# sys.argv = command_line.split(' ')
