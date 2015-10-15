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
opt = matrix(c(
    'inputfile_name', 'i', 1, 'character',
    'type', 't', 1, 'character',
    'outputfile_name', 'o', 1, 'character'
),byrow=TRUE, ncol=4);

# parse de input
options = getopt(opt);

table <- read.delim(options$inputfile_name, header=TRUE, fill=TRUE);

if (options$type == "lfqlog2") {
  regexpr <- "^[[:digit:]]+.*([^[:digit:]]+)[.]lfq[.]intensity$";
  code <- "LFQ";
} else if (options$type == "intensity") {
  regexpr <- "^[[:digit:]]+.*([^[:digit:]]+)[.]intensity$";
  code <- "INT";
} else {
  regexpr <- "^[[:digit:]]+.*([^[:digit:]]+)[.]speccount$";
  code <- "MS";
}
if (!(TRUE %in% grepl(regexpr, colnames(table)))) {
  write(sprintf("Error: No columns of type %s in input table", code), stderr());
  q(status=3);
}
#define de options input that the code will have

# define the columns that will be taken in account for the anova
columns_names <- grep(regexpr, colnames(table), value=TRUE);

# here I extract the different experiment names in an array for easier
# manipulation, ordering them
experiment_names <- mixedsort(gsub("^([[:digit:]]+.*[^[:digit:]]+)[.].*$", "\\1",
                                    columns_names));

# extract from the experiment names all the different categories in the table
different_categories <- unique(gsub("^[[:digit:]]+.*([^[:digit:]]+).*$", "\\1",
                                    experiment_names));

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

# this loop computes the ttest result for each row
# and adds it to a vector
i <- 1;
ttestresult <- c("");
ttestsignificant <- c("");
if (length(different_categories) < 2) {
  print(sprintf("Can't calculate t-test. There is only one category for %s collumns", code));
  q(1,save="no");
}

for (i in seq(1, nrow(table_only_columns))) {
  # the t-test arguments are the control values vector, the treatment values vector
  # and some extra arguments. var.equal says it's a student t-test with stardard
  # deviations assumed equal. mu=0 sets the hipothesis to be null.
  ttestresult[i+1] <- t.test(as.numeric(table_only_columns[i, columns[[1]]]), as.numeric(table_only_columns[i, columns[[2]]]), var.equal=TRUE, mu=0)$p.value;
  if (is.na(ttestresult[i+1]))
    ttestresult[i+1] = 1.0
}

# this defines if the p-value returned for each row is significant
ttestsignificant[as.numeric(ttestresult) <= 0.05] <- "+"
ttestsignificant[as.numeric(ttestresult) > 0.05] <- ""
ttestsignificant[1] <- "";

# create two extra rows on the table, one for p-values and other
# for siginificance
table[paste0("T.test.result.", code)] <- NA;
table[paste0("T.test.result.", code)] <- ttestresult;
table[paste0("T.test.significant.", code)] <- NA;
table[paste0("T.test.significant.", code)] <- ttestsignificant;




# write out the table
output_handler <- file(options$outputfile_name, "w")
write.table(table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
