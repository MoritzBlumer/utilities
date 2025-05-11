#!/usr/bin/env python3

"""
Calculate allellic depth for mitochonrdial variants.
Since mitochondria are haploid, there should be only one allele per variant.
Exception: heteroplasmy (rare – https://en.wikipedia.org/wiki/Heteroplasmy)
"""

import sys, subprocess
import numpy as np
import pandas as pd
import io


## functions

def parse_command_line_arguments():
    '''
    parse command line arguments
    '''

    # declare all variables global
    global vcf_path, primary_ids_path, mito_name, mito_len, base_error_rate, \
        ad_output_path, stats_output_path

    # print help message if incorrect number of arguments was specified
    if len(sys.argv) != 8:
        print('\nUsage:', file=sys.stderr)
        print('\tpython extract_allelic_depths.py <vcf_path> <primary_ids_path> <mito_name> <mito_len> <base_error_rate> <ad_output_path> <stats_output_path>\n', file=sys.stderr)
        print('\t\t<vcf_path>\tstr\tpath to VCF/BCF that contains mitochondrion', file=sys.stderr)
        print('\t\t<primary_ids_path>\t\tstr\tpath to file with one primary_id per line to be included', file=sys.stderr)
        print('\t\t<mito_name>\t\tstr\tname of the mitochondrial scaffold', file=sys.stderr)
        print('\t\t<mito_len>\t\tstr\length of the mitochondrial sequence', file=sys.stderr)
        print('\t\t<base_error_rate>\tstr\tsequencing error rate, e.g. 0.001', file=sys.stderr)
        print('\t\t<ad_output_path>\tstr\tpath to main output TSV that will contain per position allelic depth data for all input samples', file=sys.stderr)
        print('\t\t<stats_output_path>\tstr\tpath to per-sample summary_stats TSV', file=sys.stderr)

    # fetch arguments
    _, vcf_path, primary_ids_path, mito_name, mito_len, base_error_rate, \
        ad_output_path, stats_output_path  = sys.argv
    mito_len = int(mito_len)
    base_error_rate = float(base_error_rate)


def zero_division(val_1, val_2):
    '''
    return 0 if dividing by 0
    '''

    return val_1/val_2 if val_2 != 0 else 0


def parse_row(row, query_primary_id):
    '''
    parse VCF rows and retrieve POS, all alleles (max 4 [ATCG], and ignoring 
    REF/ALT) and the corresponding allelic depths
    '''

    # extract 
    row = np.array(row)
    pos = row[0]
    alleles_arr = row[1:5]
    a_depths_arr = row[5:9]

    # sort ALT alleles and corresponding depths from largest to smalles depths
    alleles_arr = alleles_arr[np.flip(np.argsort(a_depths_arr))]
    a_depths_arr = a_depths_arr[np.flip(np.argsort(a_depths_arr))]

    # remove alleles (and corresponding allelic depth values) if: 1) depth = 1 
    # [only single read supports], 2) % supporting reads < than sequencing 
    # error rate
    # [would also remove indel alleles, but these are now ignored generally 
    # ignored the bcftools command]
    pass_idx_arr = \
        [
            idx for idx, (allele, depth) \
                in enumerate(zip(alleles_arr, a_depths_arr)) \
                    if (depth != 1) \
                        and (zero_division(depth, sum(a_depths_arr)) > base_error_rate) \
                            and len(allele) == 1]
    a_depths_arr = a_depths_arr[pass_idx_arr]
    alleles_arr = alleles_arr[pass_idx_arr]


    # [SANITY CHECK I] print error message if # alleles != # depths (should 
    # really not be possible...)
    if len(alleles_arr) != len(a_depths_arr):
        print('[ERROR] Number of alleles ≠ number of allelic depth values', 
              file=sys.stderr)
        sys.exit()


    # [SANITY CHECK II] make sure there are no other alleles than ATCG, else 
    # raise error and exit
    bases = ['A', 'T', 'C', 'G']
    for allele in alleles_arr:
        if allele not in bases:
            print('[ERROR] Unknown base found: ' + allele, file=sys.stderr)
            sys.exit()


    # return pd.Series of primary_id, POS, allele_1, allele_2, allele_3, 
    # allele_4, allelic_depth_1, allelic_depth_2, allelic_depth_3, 
    # allelic_depth_4, sum of all allelic depths, sum of all minor allelic 
    # depths, ratio of the latter two
    parsed_row = pd.Series(
        [query_primary_id, pos] \
         + list(alleles_arr) \
         + ['NaN' for x in list(range(0, 4-len(alleles_arr)))] \
         + list(a_depths_arr) \
         + [0 for x in list(range(0, 4-len(a_depths_arr)))] \
         + [sum(a_depths_arr)] \
         + [sum(a_depths_arr[1:])] \
         + [zero_division(sum(a_depths_arr[1:]), sum(a_depths_arr))])

    return parsed_row



def retrieve_variants(vcf_path, query_primary_id, mito_name):    
    '''
    Use bcftools view and bcftools query to extract allelic depth info for
    mitochondrial sequence
    '''

    # fetch allelic depths from VCF
    shell_command = "bcftools view -V indels -e 'INFO/INDEL!=0' -s " \
        + query_primary_id + ' ' + vcf_path + ' ' + mito_name \
            + " | bcftools query -f \
    '%POS\t%REF\t%ALT{0}\t%ALT{1}\t%ALT{2}[\t%AD{0}\t%AD{1}\t%AD{2}\t%AD{3}]\n'" 

    shell_command_out = subprocess.getoutput(shell_command)

    allelic_depths_df = pd.read_csv(
        io.StringIO(shell_command_out), 
        sep='\t', 
        header=None, 
        names=[
            'pos', 'a_1', 'a_2', 'a_3', 'a_4', 'ad_1', 'ad_2', 'ad_3', 'ad_4'
        ]
    )

    # replace '.' with 0 and adjust dtypes
    allelic_depths_df = \
        allelic_depths_df.replace('.', 0).astype(
            {
                'pos': int, 
                'a_1': str, 
                'a_2': str, 
                'a_3': str, 
                'a_4': str, 
                'ad_1': int, 
                'ad_2': int, 
                'ad_3': int, 
                'ad_4': int
            }
        )

    # make sure there are no duplicate rows (e.g. InDels being mistakingly read 
    # in and coinciding with a SNP)
    if len(allelic_depths_df.pos) != len(set(allelic_depths_df.pos)):
    
        print('[ERROR] Number of rows ≠ number of positions', file=sys.stderr)

        sys.exit()

    # re-format and then return as dataframe
    filtered_df = pd.DataFrame(
        allelic_depths_df.apply(parse_row, args=[query_primary_id], axis=1)
    )

    filtered_df.columns = [
        'primary_id', 
        'pos', 
        'a_1', 'a_2', 'a_3', 'a_4', 
        'ad_1', 'ad_2', 'ad_3', 'ad_4', 
        'cov_tot', 
        'cov_alt', 
        'cov_alt_pct'
    ]

    return filtered_df



def main():
    
    # fetch command line arguments
    parse_command_line_arguments()

    # read in query primary ids
    query_primary_ids_lst = [
        x.strip() for x in open(primary_ids_path, 'r').readlines()
    ]

    # initiate (a) list to hold the individual filtered per-sample dataframes 
    # that are later concatenated and (b) a summary_df that will be populated 
    # with per-sample information
    filtered_df_lst = []
    summary_stats_df = pd.DataFrame(
        columns=[
            'primary_id',
            'mito_n_variable_sites',
            'mito_cumulative_minor_bases',
            'mito_cumulative_total_bases',
            'mito_cumulative_minor_bases_pct'
        ]
    )

    # iterate input samples
    print(f'\n[INFO] Processing {len(query_primary_ids_lst)} samples\n', 
          file=sys.stderr)

    for idx, query_primary_id in enumerate(query_primary_ids_lst):

        # retrieve AD info for each poistion in the VCF for current sample
        sample_df = retrieve_variants(vcf_path, query_primary_id, mito_name)
        filtered_df_lst.append(sample_df)

        # obtain per-sample stats and append to summary_stats_df
        variable_sites = len(sample_df[sample_df['cov_alt'] > 0])
        cumulative_minor_bases = sum(sample_df['cov_alt'])
        cumulative_total_bases = sum(sample_df['cov_tot'])
        cumulative_minor_bases_pct = (
            cumulative_minor_bases/cumulative_total_bases
        )*100
        summary_stats_df = pd.concat(
            [
                summary_stats_df, 
                pd.DataFrame(
                    [
                        [query_primary_id, 
                         variable_sites, 
                         cumulative_minor_bases, 
                         cumulative_total_bases, 
                         cumulative_minor_bases_pct
                         ]
                    ],
                    columns = [
                        'primary_id', 
                        'mito_n_variable_sites', 
                        'mito_cumulative_minor_bases', 
                        'mito_cumulative_total_bases', 
                        'mito_cumulative_minor_bases_pct'
                    ]
                )
            ]
        )

        # print info
        print(f'[INFO] Processed {idx+1} samples', 
              file=sys.stderr, flush=True)

    # concatenate all sample_dfs and save as TSV
    filtered_df = pd.concat(filtered_df_lst)
    filtered_df.to_csv(ad_output_path, sep='\t', index=False)

    # write summary_stats to TSV
    summary_stats_df.to_csv(
        stats_output_path, 
        sep='\t', 
        index=False
    )

    # print info
    print('\n[INFO] Done', file=sys.stderr)


# EXECUTE
if __name__ == "__main__":
