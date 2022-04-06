# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/4/29


import os
import re
import sys
sys.path.append('/public1/home/stu_panhaoran/biosoft/pyhr/scripts')
import argparse
from path import Path
import pandas as pd
import richlog as rl
import cal_interval_dist as cid
from rich.traceback import install


__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20210926'
__version__ = 'v1.0'


def run_Ksdensity(Ks_file,limit,step,ks_name):
    ks_count_dict = cid.cal_Ks_valuecounts(Ks_file,limit)
    tree = cid.split_windows(limit,step)
    newtree = cid.cal_Ks_density(ks_count_dict,tree,limit)
    output = cid.convert_result(newtree,ks_name)
    return output


def read_Ksfile(ks_file_path):
    ks_file_list = [i for i in Path(ks_file_path).files() if i.endswith('KaKs.result')]
    return ks_file_list


def merge_ks_reads_mapping_df(ks_files,Upper_limit,step):
    file_dict = {}
    ks_file_name_list = []
    for i in ks_files:
        j = re.findall('([A-Za-z0-9\_]+)\.KaKs.result',i)[0]
        ks_file_name_list.append(j)
    rl.info_out('Read first Ks file:{}...'.format(ks_files[0]))
    first_file = run_Ksdensity(ks_files[0],Upper_limit,step,ks_file_name_list[0])
    for i in range(1,len(ks_files)):
        rl.info_out('Read Ks file:{}...'.format(ks_files[i]))
        file_dict[i] = run_Ksdensity(ks_files[i],Upper_limit,step,ks_file_name_list[i])
        if i < len(ks_files):
            first_file = pd.merge(first_file,file_dict[i],on='index',how='outer')
    choose_columns_number =  []
    for i in range(0,first_file.shape[1],2):
        choose_columns_number.append(i)
    output_ks_file = first_file.iloc[:,choose_columns_number]
    output_ks_file = output_ks_file.fillna(0)
    output_ks_file = output_ks_file.sort_values(by='index')
    return output_ks_file
    

def output_results(output_ks_file,output_name):
    output_ks_file.to_excel(output_name,header=True,index=False)


def main(args):
    install()
    Ks_file_path = args.path
    Upper_limit = args.limit
    step = args.step
    output_file = args.output
    
    ks_files = read_Ksfile(Ks_file_path)
    result = merge_ks_reads_mapping_df(ks_files,Upper_limit,step)
    output_results(result,output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Calculate the interval distribution of gene pairs with different synonymous substitution ratio''',
        usage="python {} -p path -1 Upper_limit -s step -o output.xlsx".format(sys.argv[0]),
        epilog= 'author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
    parser.add_argument('-p', '--path', required=True, help='Path to store ks file')
    parser.add_argument('-l', '--limit', required=True, type=float,help='Upper limit of ks reads mapping')
    parser.add_argument('-s', '--step',required=True,type=float,help="Input the step size")
    parser.add_argument('-o', '--output', required=True,help='The name of the output file(.xlsx)')
    args = parser.parse_args()
    main(args)
