#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/11/28

import sys
import gzip
import time
import Fontcolor
import argparse


__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20201207'
__version__ = 'v1.1'


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


def calcultate_sample_length(fq_gz_file):
    start = time.time()
    reads_num = 0
    line_index = 0
    reads_length = 0
    with gzip.open(fq_gz_file) as fq1:
        for line in fq1:
            line_index +=1
            if line_index %4 == 1:
                reads_num += 1
            elif line_index % 4 ==2:
                reads_length += len(line)
    end = time.time()
    Fontcolor.Tips_output('File {} reads num are {}'.format(str(fq_gz_file),reads_num))
    Fontcolor.Tips_output('File {} reads length are {}'.format(str(fq_gz_file),reads_length))
    Fontcolor.Tips_output('Time cost: {:.5f}'.format(end-start))
    return reads_length


def cal_genome_size(genome_fasta):
    fa_dict = {}
    with open(genome_fasta) as f:
        for line in f:
            if line.startswith('>'):
                fa_key = line.strip().split('>')[-1].split(' ')[0]
                fa_dict[fa_key] = ''
            else:
                fa_dict[fa_key] += line.strip().upper().replace('\n', '')
    genome_size = 0
    for key in fa_dict.keys():
        genome_size += len(fa_dict[key])
    Fontcolor.Tips_output('Genome size is {} bp'.format(genome_size))
    return genome_size


def run_step2(fq1file,fq2file,genomefile):
    fq1length = calcultate_sample_length(fq1file)
    fq2length = calcultate_sample_length(fq2file)
    genomesize = cal_genome_size(genomefile)
    std_depth = round((fq1length + fq2length) / genomesize)
    Fontcolor.Tips_output('Standard Sequencing Depth is {}X'.format(std_depth))
    return std_depth


def Fa2dict(fasta_file):
    fa_dict = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                fa_key = line.strip().split('>')[-1].split(' ')[0]
                fa_dict[fa_key] = ''
            else:
                fa_dict[fa_key] += line.strip().upper().replace('\n', '')
    return fa_dict


def ctg_phasing_base_depth(fasta_dict,depthcountfile,std_depth,Ploidy,result):
    depth_file = open(depthcountfile,'r')
    std_depth  = int(std_depth)
    ploidy = int(Ploidy)
    phasing_file=open('Phasing_record.txt','w')
    last_output = open(result,'w')
    for dep_line in depth_file:
            dep_line_list=dep_line.rstrip().split('\t')
            raw_depth=float(dep_line_list[1])/float(std_depth)
            depth=int(round(raw_depth))
            if depth > 1 :
                    phasing_file.write(dep_line.rstrip()+'\t'+str(raw_depth)+'\n')
                    if depth - 1 < ploidy:
                        for i in range(depth-1):
                                last_output.write('>'+dep_line_list[0]+'_'+str(i)+'\n')
                                last_output.write(fasta_dict[dep_line_list[0]]+'\n')
                    else:
                        for i in range(ploidy - 1):
                            last_output.write('>' + dep_line_list[0] + '_' + str(i) + '\n')
                            last_output.write(fasta_dict[dep_line_list[0]] + '\n')
    phasing_file.close()
    last_output.close()
    depth_file.close()


def main(args):
    depthfile = args.depth
    depthfileout = args.depthout
    Fontcolor.Log_output('-' * 10 + 'step 1 beginning !' + '-' * 10)
    depth_count(depthfile,depthfileout)
    Fontcolor.Log_output('-' * 10 + 'step 1 completed !' + '-' * 10)

    if not args.stdth:
        fq1_file = args.fastq1
        fq2_file = args.fastq2
        genome_file = args.genome
        Fontcolor.Log_output('-' * 10 + 'step 2 beginning !' + '-' * 10)
        std_depth = run_step2(fq1_file,fq2_file,genome_file)
        Fontcolor.Log_output('-' * 10 + 'step 2 completed !' + '-' * 10)
    else:
        Fontcolor.Log_output('-' * 10 + 'step 2 beginning !' + '-' * 10)
        std_depth = args.stdth
        Fontcolor.Tips_output('Input std depth is {}X'.format(std_depth))
        Fontcolor.Log_output('-' * 10 + 'step 2 completed !' + '-' * 10)

    Ploidy = args.ploidy
    output = args.output
    genome_file = args.genome
    Fontcolor.Log_output('-' * 10 + 'step 3 beginning !' + '-' * 10)
    fasta_dict = Fa2dict(genome_file)
    ctg_phasing_base_depth(fasta_dict,depthfileout,std_depth,Ploidy,output)
    Fontcolor.Log_output('-' * 10 + 'step 3 completed !' + '-' * 10)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Polyploid genome phasing:Use Second-generation sequencing reads mapping Third-generation sequencing genome''',
        usage="python {} -d depth -o depthout -f1 fq1 -f2 fq2 -g genome -p ploidy -O output\n\t\tor\n\t\tpython {} -d depth -o depthout -s std_depth -g genome -p ploidy -O output".format(sys.argv[0],sys.argv[0]),
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__, __mail__, __date__,
                                                                            __version__))
    parser.add_argument('-d', '--depth', required=True, help='Input samtools depth outputfile')
    parser.add_argument('-o', '--depthout', required=True, help='Output ctg depth result(.tab)')
    parser.add_argument('-f1', '--fastq1', required=None, help='Input sample1 file(fastq1.gz)')
    parser.add_argument('-f2', '--fastq2', required=None, help='Input sample2 file(fastq2.gz)')
    parser.add_argument('-s', '--stdth', required=None, help='Input the standard depth of known sequencing data, which is equal to the sequencing data/genome size(bp)')
    parser.add_argument('-g', '--genome', required=True, help='Input Thr Sequencing assmbly genome(.fasta)')
    parser.add_argument('-p', '--ploidy', required=True, help='Input Ploidy of species(eg:8)')
    parser.add_argument('-O', '--output', required=True, help='Output the fasta sequence of the contig that only needs to be phased(.fasta)')

    args = parser.parse_args()
    main(args)
