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

#include and execute the read util script
source('../r_utils/read_util.R');

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
table_only_columns <- table[aux]

# this loop computes the ttest result for each row
# and adds it to a vector
i <- 1;
ttestresult <- c();
ttestsignificant <- c();
if (length(different_categories) < 2) {
  print(sprintf("Can't calculate t-test. There is only one category for %s collumns", code));
  q(1,save="no");
}

for (i in seq(1, nrow(table_only_columns))) {
  # the t-test arguments are the control values vector, the treatment values vector
  # and some extra arguments. var.equal says it's a student t-test with stardard
  # deviations assumed equal. mu=0 sets the hipothesis to be null.
  ttestresult[i] <- t.test(table_only_columns[i, columns[[1]]],
    table_only_columns[i, columns[[2]]], var.equal=TRUE, mu=0)$p.value;
  if (is.na(ttestresult[i]))
    ttestresult[i] = 1.0
}

# this defines if the p-value returned for each row is significant
ttestsignificant[ttestresult <= 0.05] <- "+"
ttestsignificant[ttestresult > 0.05] <- ""


# create two extra rows on the table, one for p-values and other
# for siginificance
#TODO: ou colocar perto da intensidade que se refere ou na 3ª coluna
table[paste0("T.test.result.", code)] <- NA;
table[paste0("T.test.result.", code)] <- ttestresult;
table[paste0("T.test.significant.", code)] <- NA;
table[paste0("T.test.significant.", code)] <- ttestsignificant;




# write out the table
source('../r_utils/write_util.R');
