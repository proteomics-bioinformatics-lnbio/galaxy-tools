#!/usr/bin/env Rscript

# fold_change.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: October 2015
#
# Required packages to work: (getopt", "gtools", "utils")
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.
require('gtools', quietly=TRUE);
require('utils', quietly=TRUE);
require('getopt', quietly=TRUE);
opt <- matrix(c(
    'inputfile_name', 'i', 1, 'character',
    'type', 't', 1, 'character',
    'outputfile_name', 'o', 1, 'character'
),byrow=TRUE, ncol=4);
# xml wrapper argument reading
options <- getopt(opt);

#input table import
table <- read.delim(options$inputfile_name, header=TRUE, fill=TRUE);

#get the type to work with
if (options$type == "lfqlog2") {
  regexpr <- "([^[:digit:]]+).*[[:digit:]]+[.]lfq[.]intensity";
  code <- "LFQ";
} else if (options$type == "intensity") {
  regexpr <- "([^[:digit:]]+).*[[:digit:]]+[.]intensity";
  code <- "INT";
} else {
  regexpr <- "([^[:digit:]]+).*[[:digit:]]+[.]speccount";
  code <- "MS";
}
if (!(TRUE %in% grepl(regexpr, colnames(table)))) {
  sprintf("Error: No columns of type %s in input table", code);
}

# define the columns that will be taken in account for the foldchange
columns_names <- grep(regexpr, colnames(table), value=TRUE);

# here I extract the different experiment names in an array for easier
# manipulation, ordering them
experiment_names <- mixedsort(gsub("([^[:digit:]]+.*[[:digit:]]+)[.].*", "\\1",
                                    columns_names));

# extract from the experiment names all the different categories in the table
different_categories <- unique(gsub("([^[:digit:]]+).*", "\\1",
                                    experiment_names));

#routine to define which columns will be worked
i<-1;
columns <- list();
aux <- c();
for (cat in different_categories) {
  col <- columns_names[gsub(regexpr, "\\1", columns_names) == cat]
  aux <- c(aux, col);
  columns[[i]] <- col;
  i<-i+1;
}
# this is a filtered table to help with calculations
table_only_columns <- table[-1, aux]
#generate a combination of all different categories 2 by 2
aux <- combn(different_categories, 2)

for (j in 1:ncol(aux)) {
  table[paste0(code, ".fold.change.result.", aux[1, j], ".vs.", aux[2, j])] <- c("");
}
for (j in 1:length(different_categories)) {
  table[paste0(code, ".stddev.", different_categories[j])] <- c("");
}
# this loop computes the foldchange result for each row
# and adds it to a vector
i <- 1;
for (i in seq(1, nrow(table_only_columns))) {
  j <- 1;
  for (j in 1:ncol(aux)) {
    a <- as.numeric(table_only_columns[i, columns_names[gsub(regexpr, "\\1", columns_names) == aux[1, j]]]);
    b <- as.numeric(table_only_columns[i, columns_names[gsub(regexpr, '\\1', columns_names) == aux[2, j]]]);
    table[i+1, paste0(code, ".fold.change.result.", aux[1, j], ".vs.", aux[2, j])] <- foldchange(mean(a, na.rm=TRUE), mean(b, na.rm=TRUE));
  } #compute the standar deviation
  for (j in 1:length(different_categories)) {
    a <- as.numeric(table_only_columns[i, columns_names[gsub(regexpr, "\\1", columns_names) == different_categories[j]]]);
    table[i+1, paste0(code, ".stddev.", different_categories[j])] <- sd(a, na.rm=TRUE);
  }
}

# write out the table
output_handler <- file(options$outputfile_name, "w")
write.table(table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
