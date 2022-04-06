#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/4/19


import sys
from path import Path
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
sys.path.append('/share/home/stu_panhaoran/scripts')
import SetLog as log


__author__ = 'Haoran Pan'
__mail__ = 'haoran_pan@qq.com'
__date__ = '20210425'
__version__ = 'v1.2'



def get_msa(filepath):
    msa_files = [file for file in Path(filepath).files() if file.endswith('.msa')]
    return msa_files


def parse_msa(msa_files):
    gene2seq = {}
    for msa in msa_files:
        for seq in SeqIO.parse(msa,'fasta'):
            gene2seq[seq.id] = seq.seq
    return  gene2seq


def msa2phylip(gene2seq,species_list,file_type):
    species_length = []
    species2seq = {}
    for species in species_list:
        if not '/' in species:
            concat_seq = Seq('')
            for key in gene2seq.keys():
                if key.startswith(species):
                    concat_seq += gene2seq[key]
            species_length.append(len(concat_seq))
            species2seq[species] = concat_seq
        else:
            name1,name2 = species.split('/')
            concat_seq = Seq('')
            for key in gene2seq.keys():
                if key.startswith(name1):
                    concat_seq += gene2seq[key]
                elif key.startswith(name2):
                    concat_seq += gene2seq[key]
            species_length.append(len(concat_seq))
            species2seq[name1] = concat_seq

    length_uniq = list(set(species_length))
    if len(length_uniq) == 1:
        with open('SingleCopy.{}.phylip'.format(file_type),'w') as f:
            f.write('{}\t{}\n'.format(len(species_list),length_uniq[0]))
            for species in species_list:
                if '/' in species:
                    name1, name2 = species.split('/')
                    # print(name1,species2seq[name1])
                    f.write('{}    {}\n'.format(name1,species2seq[name1]))
                else:
                    # print(species,species2seq[species])
                    f.write('{}    {}\n'.format(species,species2seq[species]))
    else:
        species_length = [str(i) for i in species_length]
        print('\t'.join(species_length))
        log.error_out('Species super-gene length is inconsistent')


def main(args):
    msapath = args.path
    specieslist = args.species
    file_type = args.type
    msa_files = get_msa(msapath)
    gene2seq = parse_msa(msa_files)
    msa2phylip(gene2seq,specieslist,file_type)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawTextHelpFormatter,
        description='''Combine the results of multiple sequence alignments end to end to form a super-gene to generate a phylip file''',
        usage="python {} -p msa_path -s ID1 ID2 ID3..\n\nNote:ID is the ID of the gene in the genome of each species!!!\nGenerated file: SingleCopy.pep.phylip".format(sys.argv[0]),
        epilog= 'author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
    parser.add_argument('-p','--path',required=True,help='Input the path where the msa file is located')
    parser.add_argument('-t','--type',required=True,choices=['pep','cds'],help='Specify the type of input sequence')
    parser.add_argument('-s','--species',required=True,nargs='+',help='Input the gene ID prefix for each species genome')
    args = parser.parse_args()
    main(args)