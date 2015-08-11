#!/usr/bin/env Rscript

# remove_zeros.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: May 2015
#
# Required packages to work: ("getopt")
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2015  Brazilian Center for Research in Energy and Materials
# All rights reserved.

source('../R_util/read-utils.R');

args <- get_cmd_options(TRUE);

colnames <- grep(args$regexpr, colnames(args$table), value=TRUE);

excluded_rows <- c();
for (row in seq(1, nrow(args$table[colnames]))) {
  if (grepl(args$table[row,colnames], 0) == TRUE) {
    excluded_rows <- c(excluded_rows, row);
  }
}

included_rows <- seq(1, nrow(args$table))[seq(1, nrow(args$table)) %in% excluded_rows == FALSE];

result <- args$table[included_rows,];

output_handler <- file(args$options$outputfile_name, "w")

write.table(result, file=output_handler, sep="\t", row.names=FALSE);

close(output_handler)
