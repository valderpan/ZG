#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/12/30


import sys
import pandas as pd

def parse_bamcov(covfile):
    Chr2cov = {}
    with open(covfile) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            line_list = line.split()
            if int(line_list[3]) != 0:
                if line_list[0] not in Chr2cov:
                    Chr2cov[line_list[0]] = 0
                    Chr2cov[line_list[0]] += int(line_list[2]) - int(line_list[1])
                else:
                    Chr2cov[line_list[0]] += int(line_list[2]) - int(line_list[1])
    return Chr2cov


def read_fai(file):
    df = pd.read_table(file,sep='\t',names=['seqid','len','uk1','uk2','uk3'])
    Chr2len = df.set_index('seqid').to_dict()['len']
    return Chr2len


def Output_res(covfile,Chr2cov,Chr2len):
    Gen_len = 0 ;Cov_len = 0
    with open('{}.stat'.format(covfile),'w') as w:
        w.write('{}\t{}\t{}\t{}\n'.format('seqid','Cov_len','Gen_cov','ratio'))
        for key in Chr2cov.keys():
            w.write('{}\t{}\t{}\t{}\n'.format(key,Chr2cov[key],Chr2len[key],round(Chr2cov[key]/Chr2len[key],5)*100))
            Gen_len += Chr2len[key]
            Cov_len += Chr2cov[key]
        w.write('{}\t{}\t{}\t{}'.format('Total',Cov_len,Gen_len,round(Cov_len/Gen_len,5)*100))


if __name__ == "__main__":
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==2:
        Cov_file = args[0]
        fai_file = args[1]
        Chr2cov = parse_bamcov(Cov_file)
        Chr2len = read_fai(fai_file)
        Output_res(Cov_file,Chr2cov,Chr2len)
    else:
        sys.exit(p.print_help())
