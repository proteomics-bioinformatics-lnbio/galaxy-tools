#!/usr/bin/env Rscript

# fold_change.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: April 2015
#
# Required packages to work: (getopt", "gtools", "utils")
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.
require('gtools', quietly=TRUE);
require('utils', quietly=TRUE);
source('../R_util/read-utils.R');
#define de options input that the code will have
args <- get_cmd_options(TRUE);

# define the columns that will be taken in account for the anova
columns_names <- grep(args$regexpr, colnames(args$table), value=TRUE);

# here I extract the different experiment names in an array for easier
# manipulation, ordering them
experiment_names <- mixedsort(gsub(".*[.]([^[:digit:]]+[[:digit:]]+).*", "\\1",
                                    columns_names));

# extract from the experiment names all the different categories in the table
different_categories <- unique(gsub("([^[:digit:]]+).*", "\\1",
                                    experiment_names));

i<-1;
columns <- list();
aux <- c();
for (cat in different_categories) {
  col <- columns_names[gsub(args$regexpr, "\\1", columns_names) == cat]
  aux <- c(aux, col);
  columns[[i]] <- col;
  i<-i+1;
}
# this is a filtered table to help with calculations
table_only_columns <- args$table[-1, aux]
aux <- combn(different_categories, 2)
# this loop computes the ttest result for each row
# and adds it to a vector
for (j in 1:ncol(aux)) {
  args$table[paste0(args$code, ".fold.change.result.", aux[1, j], ".vs.", aux[2, j])] <- NA;
}
for (j in 1:length(different_categories)) {
  args$table[paste0(args$code, ".stddev.", different_categories[j])] <- NA;
}
i <- 1;

for (i in seq(1, nrow(table_only_columns))) {
  j <- 1;
  for (j in 1:ncol(aux)) {
    a <- as.numeric(table_only_columns[i, columns_names[gsub(args$regexpr, "\\1", columns_names) == aux[1, j]]]);
    b <- as.numeric(table_only_columns[i, columns_names[gsub(args$regexpr, '\\1', columns_names) == aux[2, j]]]);
    args$table[i+1, paste0(args$code, ".fold.change.result.", aux[1, j], ".vs.", aux[2, j])] <- foldchange(mean(a), mean(b));
  }
  for (j in 1:length(different_categories)) {
    a <- as.numeric(table_only_columns[i, columns_names[gsub(args$regexpr, "\\1", columns_names) == different_categories[j]]]);
    args$table[i+1, paste0(args$code, ".stddev.", different_categories[j])] <- sd(a);
  }
}


# write out the table
output_handler <- file(args$options$outputfile_name, "w")
write.table(args$table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
