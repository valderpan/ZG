#! /usr/bin/env Rscript
# @Author : Haoran Pan
# date: 2021/04/10

#=================================================
#LTR InsertTime calcultate
#=================================================

library(dplyr)
library(ggplot2)
library(viridis)
library(ggsci)
library(ggprism)
# option_list <- list(make_option(c("-i","--input"),type="character", default = NULL,help  = "Input ATAC insertSize file"),
#                     make_option(c("-o","--output"),type="character",default = NULL,help  = "Output the PDF file"),
#                     make_option(c("-n","--name"),type="character", default = NULL,help="Experiment name, which will be used to generate output plot title names"))

# args <- parse_args(OptionParser(option_list=option_list))
# setwd('E:\\潘浩然\\桌面\\脚本测试\\ltr')
# filelist <- dir(pattern = ".pass.list")
# allfile <- lapply(filelist, read.table)
# for (i in allfile) {
#   print(which(i))
#   # a <- i
#   # view(a)
#   break
# }

#df <- read.table('E:\\潘浩然\\桌面\\脚本测试\\ltr\\Sspon.genome.fasta.pass.list',sep='\t')
#df %>% mutate(MyA=V12/1000000) %>% ggplot(aes(MyA))+geom_density(adjust=1.5,fill="#69b3a2", color="#e9ecef", alpha=0.8)


setwd('E:\\潘浩然\\脚本测试\\20210417ltr')
df_ss <- read.table('Sspon.genome.fasta.pass.list',sep = '\t')
df_sb <- read.table('Sb.genome.fasta.pass.list',sep = '\t')
df_npx <- read.table('NpX.genome.fasta.pass.list',sep = '\t')
df_la <- read.table('LA.genome.fasta.pass.list',sep = '\t')

df_Ss <- df_ss %>% mutate(species='Ss',MyA=V12/1000000) %>% select(species,MyA)
df_Sb <- df_sb %>% mutate(species='Sb',MyA=V12/1000000) %>% select(species,MyA)
df_NpX <- df_npx %>% mutate(species='NpX',MyA=V12/1000000) %>% select(species,MyA)
df_La <- df_la %>% mutate(species='La',MyA=V12/1000000) %>% select(species,MyA)


df <- rbind(df_Ss,df_Sb,df_NpX,df_La)
#第一种配色
ggplot(df,aes(MyA,fill=species,color=species))+
  geom_density(adjust=1.5, alpha=0.7)+
  theme_bw()+scale_fill_viridis(discrete = T)+scale_color_viridis(discrete=TRUE)+
  xlab('Insertion time(MYA)')+ylab('Number of intact LTR-RTs')+
  theme_prism(border = TRUE, base_rect_size = 1) +
  coord_cartesian(clip = "off")+scale_x_continuous(breaks = seq(0,6,1))
#第二种配色
ggplot(df,aes(MyA,fill=species,color=species))+
  geom_density(adjust=1.5, alpha=0.7)+
  theme_bw()+scale_fill_viridis(discrete = T)+scale_fill_tron()+scale_color_tron()+
  xlab('Insertion time(MYA)')+ylab('Number of intact LTR-RTs')+
  theme_prism(border = TRUE, base_rect_size = 1) +
  coord_cartesian(clip = "off")+scale_x_continuous(breaks = seq(0,6,1))
