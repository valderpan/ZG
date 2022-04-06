#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/4/18


import sys
import re
import pandas as pd

def reads_GeneCounts(file):
    df = pd.read_table(file,sep='\t')
    return df


def countfamilies(GeneCountsdf):
    col_names = [col for col in GeneCountsdf][1:-1]
    species2OG = {}
    for colnames in col_names:
        species_name = re.findall('([0-9a-zA-Z\_\.]+)\.pep',colnames)[0]
        families_OG = GeneCountsdf[GeneCountsdf[colnames] !=0].iloc[:,0].tolist()
        with open(species_name+'.group.txt','w') as f:
            for OG in families_OG:
                f.write(OG+'\n')
        species2OG[species_name] = families_OG
    return species2OG


def Out_dataframe(speciesOG):
    speciesL = []
    for s in species2OG.keys():
        qdf = pd.DataFrame({s:species2OG[s]})
        speciesL.append(qdf)
    resdf = pd.concat(speciesL,axis=1)
    resdf.to_csv('Orthogroups.GeneCount2venn.tab',sep='\t',header=True,index=False)


if __name__ == '__main__':
    from optparse import OptionParser
    from rich.console import Console
    console = Console()
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args) == 1:
        count_file = args[0]
        with console.status("Working...", spinner="dots"):
            df = reads_GeneCounts(count_file)
            species2OG = countfamilies(df)
            Out_dataframe(species2OG)
    else:
        sys.exit(p.print_help())
