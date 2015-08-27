#!/usr/bin/env Rscript

# wilcox.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: August 2015
#
# Required packages to work: ("getopt", "gtools")
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2015  Brazilian Center for Research in Energy and Materials
# All rights reserved.
require('gtools', quietly=TRUE);
require('getopt', quietly=TRUE);

# Specification of how arguments will work
spec = matrix(c(
    'inputfilename', 'i', 1, "character",
    'type', 't', 1, "character",
    'outputfilename', 'o', 1, "character"
    ), byrow = TRUE, ncol=4);
opt <- getopt(spec);

library('../r_utils/read_util.R')
library('../r_utils/write_util.R')

read_util <- read_util(opt);

i<-1;
columns <- list();
aux <- c();
for (cat in read_util$diff_cat) {
  col <- read_util$col_names[gsub(read_util$regex, "\\1", read_util$col_names) == cat]
  aux <- c(aux, col);
  columns[[i]] <- col;
  i<-i+1;
}
# this is a filtered read_util$table to help with calculations
table_only_columns <- read_util$table[-1, aux]

# this loop computes the ttest result for each row
# and adds it to a vector
i <- 2;
wtestresult <- c("");
wtestsignificant <- c("");
if (length(read_util$diff_cat) != 2) {
  print(sprintf("Can't calculate wilcox-test. Only two %s collumns allowed", read_util$code));
  q(1,save="no");
}

for (i in seq(2, nrow(table_only_columns)+1)) {
  # the wilcox-test arguments are the control values vector, the treatment values vector
  # and some extra arguments. mu=0 sets the hipothesis to be null.
  wtestresult[i] <- wilcox.test(table_only_columns[i-1, columns[[1]]],
    table_only_columns[i-1, columns[[2]]], mu=0)$p.value;
  if (is.na(wtestresult[i]))
    wtestresult[i] = 1.0;
}

# this defines if the p-value returned for each row is significant
wtestsignificant[wtestresult <= 0.05] <- "+"
wtestsignificant[wtestresult > 0.05] <- ""


# create two extra rows on the read_util$table, one for p-values and other
# for siginificance
#TODO: ou colocar perto da intensidade que se refere ou na 3ª coluna
read_util$table[paste0("Wilcox.test.result.", read_util$code)] <- NA;
read_util$table[paste0("Wilcox.test.result.", read_util$code)] <- wtestresult;
read_util$table[paste0("Wilcox.test.significant.", read_util$code)] <- NA;
read_util$table[paste0("Wilcox.test.significant.", read_util$code)] <- wtestsignificant;

# write out the read_util$table
writeout(options$outputfile_name, read_util$table);
