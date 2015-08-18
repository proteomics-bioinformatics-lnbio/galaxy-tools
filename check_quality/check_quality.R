#!/usr/bin/env Rscript

# check_quality.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: April 2015
#
# Required packages to work: ("getopt")
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.

require(getopt, quietly=TRUE)
# Specification of how arguments will work
spec = matrix(c(
    'inputfilename', 'i', 1, "character",
    'outputfilename', 'o', 1, "character"
    ), byrow = TRUE, ncol=4);
opt <- getopt(spec);

#read table
table <- read.csv(file=opt$inputfilename, header=TRUE, fill=TRUE);
#remove blank rows
table_nozeros <- table[c(-1, which(rowSums(table[-1,-(1, 2)])>0)),-(1, 2)]
#TODO descobrir o que fazer com isso:
dim(table)
dim(table_nozeros)
#write out
write.csv(table_nozeros, file=opt$outputfilename, row.names=FALSE, col.names=FALSE)
