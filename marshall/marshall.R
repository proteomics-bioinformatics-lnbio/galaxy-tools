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
regex.id.uniprot.1 <- "^[OPQ][[:digit:]][[:upper:][:digit:]]{3}[[:digit:]]"
regex.id.uniprot.2 <- "^[A-NR-Z][[:digit:]]([[:upper:]][[:upper:][:digit:]]{2}[[:digit:]]){1,2}"
regex.id.ipi <- "^IPI[[:digit:]]+$"
regex.id.tair <- "^A[Tt][1-5MC][:word:][[:digit:]]+$"
regex.id.ensembl <- "^ENS[[:word:]]+$"
regex.id.refseq <- "^AC_|^N[CGTWSZMR]_|^X[MR]_|^[ANYXZ]P_"
regex.id.contaminant_reversed <- "^CON_|REV_"

column_names.intensity <- grep(regex.intensity, colnames(table), value=TRUE);
column_names.lfqintensity <- grep(regex.lfqintensity, colnames(table), value=TRUE);
column_names.spectral <- grep(regex.spectral, colnames(table), value=TRUE);
column_names.proteinIDs <- grep(regex.proteinIDs, colnames(table), value=TRUE);

categories.intensity <- gsub(regex.intensity, '\\1.intensity', column_names.intensity);
categories.lfqintensity <- gsub(regex.lfqintensity, '\\1.lfq.intensity.', column_names.lfqintensity);
categories.spectral <- gsub(regex.spectral, '\\1.speccount', column_names.spectral);

#add a blank column
table <- cbind(table[,c(column_names.proteinIDs)], rep("", nrow(table)), table[,c(column_names.intensity, column_names.lfqintensity, column_names.spectral)])
#rename the column names
colnames(table) <- c(column_names.proteinIDs,"Major.Protein.IDs",column_names.intensity, column_names.lfqintensity, column_names.spectral)
#set the protein IDs column as character
table[,column_names.proteinIDs] <- sapply(table[,column_names.proteinIDs], as.character)
#add a blank row
table <- rbind(table[0,], c("", "", categories.intensity, categories.lfqintensity, categories.spectral), table[seq(1, nrow(table)),])

output_handler <- file(options$outputfile_name, "w")
write.table(table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
