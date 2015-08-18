#!/usr/bin/env Rscript

# t-test.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: April 2015
#
# Required packages to work: (getopt", "gtools")
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.
require('gtools', quietly=TRUE);
require('getopt', quietly=TRUE);
#include and execute the read util script
library('../r_utils/read_util.R');
library('../r_utils/write_util.R');

#define de options input that the read_util$code will have
opt = matrix(c(
    'inputfile_name', 'i', 1, 'character',
    'type', 't', 1, 'character',
    'outputfile_name', 'o', 1, 'character'
),byrow=TRUE, ncol=4);

# parse de input
options = getopt(opt);

read_util <- read_function(options);

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
ttestresult <- c("");
ttestsignificant <- c("");
if (length(read_util$diff_cat) < 2) {
  print(sprintf("Can't calculate t-test. There is only one category for %s collumns", read_util$code));
  q(1,save="no");
}

for (i in seq(2, nrow(table_only_columns)+1)) {
  # the t-test arguments are the control values vector, the treatment values vector
  # and some extra arguments. var.equal says it's a student t-test with stardard
  # deviations assumed equal. mu=0 sets the hipothesis to be null.
  ttestresult[i] <- t.test(table_only_columns[i-1, columns[[1]]],
    table_only_columns[i-1, columns[[2]]], var.equal=TRUE, mu=0)$p.value;
  if (is.na(ttestresult[i]))
    ttestresult[i] = 1.0
}

# this defines if the p-value returned for each row is significant
ttestsignificant[ttestresult <= 0.05] <- "+"
ttestsignificant[ttestresult > 0.05] <- ""


# create two extra rows on the read_util$table, one for p-values and other
# for siginificance
#TODO: ou colocar perto da intensidade que se refere ou na 3ª coluna
read_util$table[paste0("T.test.result.", read_util$code)] <- NA;
read_util$table[paste0("T.test.result.", read_util$code)] <- ttestresult;
read_util$table[paste0("T.test.significant.", read_util$code)] <- NA;
read_util$table[paste0("T.test.significant.", read_util$code)] <- ttestsignificant;




# write out the read_util$table
writeout(options$outputfile_name, read_util$table);
