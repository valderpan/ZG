#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/8/12

'''
%prog <gff> <gff_type>

Convert gff format to bed format
Note: only the content after ID= is output
'''


import re
import sys
import pandas as pd


def read_gff(gff):
    header=['seqid','source','type','start','end','score','strand','phase','attributes']
    df = pd.read_table(gff,comment='#',names=header)
    prefix = re.findall('([0-9a-zA-Z\_\-\.]+)\.gff3',gff)[0]
    return df,prefix


def Convertgff2bed(gdf,Type):
    q_df = gdf[gdf['type'] == Type]
    ID = q_df['attributes'].str.split(';',expand=True).iloc[:,0].str.split('=',expand=True).iloc[:,1]
    bed = q_df.loc[:,['seqid','start','end']]
    bed = pd.merge(bed,ID,on=bed.index)
    bed.rename(columns={1:'ID'},inplace=True)
    return bed.iloc[:,1:]


def removectg(bed):
    cl_bed = bed[~bed['seqid'].str.contains('Contig')]
    return cl_bed


def output(bed,prefix,Type):
    bed.to_csv('{}.{}.bed'.format(prefix,Type),sep='\t',header=False,index=False)


if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 2:
        gff_file = args[0]
        Type = args[1]
        df,pre = read_gff(gff_file)
        bed = Convertgff2bed(df,Type)
        cl_bed = removectg(bed)
        output(cl_bed,pre,Type)
    else:
        sys.exit(p.print_help())