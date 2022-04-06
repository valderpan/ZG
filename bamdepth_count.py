#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/11/28


'''
%prog <depth.info> [depth.info.stat]

Statistics on the results of samtools depth

>>> samtools depth -aa bam > depth.info
>>> python %prog depth.info depth.info.stat

__author__ = 'Haoran Pan'
__mail__ = 'panpyhr@gmail.com'
__date__ = '20201128'
__version__ = 'v1.0'
'''


import sys


def depth_count(depthfile,outputfile):
    reads_num = 0
    coverage_depth = 0
    contig = ''
    result = open(outputfile,'w')
    with open(depthfile) as f:
        for line in f:
            line_list = line.split('\t')
            if reads_num == 0:
                coverage_depth += int(line_list[2].rstrip())
                reads_num += 1
                contig = line_list[0].rstrip()
            elif reads_num != 0 and line_list[0] == contig:
                coverage_depth += int(line_list[2].rstrip())
                reads_num += 1
            else:
                depth = float(coverage_depth) / float(reads_num)
                reads_num = 0
                coverage_depth = 0
                output = contig + '\t' + str(depth)
                result.write(output.rstrip()+'\n')
        last_depth = float(coverage_depth)/float(reads_num)
        result.write(contig+'\t'+str(last_depth)+'\n')
    result.close()


if __name__ == "__main__":
    from optparse import OptionParser
    from rich.traceback import install
    install()
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    if len(args)==2:
        de_file = args[0]
        out_file = args[1]
        depth_count(de_file,out_file)
    else:
        sys.exit(p.print_help())