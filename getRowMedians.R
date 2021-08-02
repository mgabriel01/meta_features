#!/usr/bin/env Rscript


args <- commandArgs(TRUE)

suppressPackageStartupMessages(library(matrixStats))

suppressPackageStartupMessages(library(data.table))

my_file<-args[1]

my_start<-args[2]

my_table<-fread(my_file, header=F, stringsAsFactors=F,data.table=F,sep="\t")

my_table<-my_table[,my_start:ncol(my_table)]

my_table<-as.matrix(my_table)

#write.table(data.frame(rowMedians(my_table)),stdout(),col.names = F,row.names = F)

write.table(data.frame(rowMeans(my_table)),stdout(),col.names = F,row.names = F)
