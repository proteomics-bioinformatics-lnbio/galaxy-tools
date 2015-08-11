#!/usr/bin/env Rscript

# marshall.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: August 2015
#
# Required packages to work: ("getopt")
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.
require('getopt', quietly=TRUE);
#define de options input that the code will have
opt = matrix(c(
    'inputfile_name', 'i', 1, 'character',
    'outputfile_name', 'o', 1, 'character'
),byrow=TRUE, ncol=4);

# parse de input
options = getopt(opt);

table <- read.delim(options$inputfile_name, header=TRUE, fill=TRUE);

regex.intensity <- "^Intensity[.]([^[:digit:]]+)[[:digit:]]+$";
regex.lfqintensity <- "^LFQ[.]intensity[.]([^[:digit:]]+)[[:digit:]]+$";
regex.spectral <- "^MS[.]MS[.]Count[.]([^[:digit:]]+)[[:digit:]]+$";
regex.proteinIDs <- "^Protein[.]IDs$"

column_names.intensity <- grep(regex.intensity, colnames(table), value=TRUE);
column_names.lfqintensity <- grep(regex.lfqintensity, colnames(table), value=TRUE);
column_names.spectral <- grep(regex.spectral, colnames(table), value=TRUE);
column_names.proteinIDs <- grep(regex.proteinIDs, colnames(table), value=TRUE);

categories.intensity <- gsub(regex.intensity, '\\1.intensity', column_names.intensity);
categories.lfqintensity <- gsub(regex.lfqintensity, '\\1.lfq.intensity.', column_names.lfqintensity);
categories.spectral <- gsub(regex.spectral, '\\1.speccount', column_names.spectral);

#function that inserts a row in the row_index
#insertRow <- function(existingdf, newrow, r_ind) {
#    existingdf[seq(r_ind+1, nrow(existingdf)+1),] <- existingdf[seq(r_ind, nrow(existingdf)),];
#    existingdf[r_ind,] <- newrow;
#    return(existingdf);
#}

#add a column
table <- cbind(table[,c(column_names.proteinIDs)], rep("", nrow(table)), table[,c(column_names.intensity, column_names.lfqintensity, column_names.spectral)])
colnames(table) <- c(column_names.proteinIDs,"Blank",column_names.intensity, column_names.lfqintensity, column_names.spectral)

table[,column_names.proteinIDs] <- sapply(table[,column_names.proteinIDs], as.character)

table <- rbind(table[0,], c("", "", categories.intensity, categories.lfqintensity, categories.spectral), table[seq(1, nrow(table)),])

#table <- data.frame(table[1,], rep("", ncol(table)), table[seq(2, ncol(table)),])

#table <- insertRow(table[,c(column_names.proteinIDs, column_names.intensity, column_names.lfqintensity, column_names.spectral)], c(NA, NA, categories.intensity, categories.lfqintensity, categories.spectral), 1);

output_handler <- file(options$outputfile_name, "w")
write.table(table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
