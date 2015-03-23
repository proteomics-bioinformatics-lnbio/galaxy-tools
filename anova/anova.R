#!/usr/bin/env Rscript

# anova.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: November 2014
#
# Required packages to work: (getopt", )
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.

require('getopt', quietly=TRUE);

#define de options input that the code will have
opt = matrix(c(
    'inputfile_name', 'i', 1, 'character',
    'type', 't', 1, 'character',
    'outputfile_name', 'o', 1, 'character'
),byrow=TRUE, ncol=4);

# parse de input
options = getopt(opt);

# reads the table from input
table <- read.delim(options$inputfile_name, header=TRUE, fill=TRUE);

# get the defined regex from the requested type
if (options$type == "lfqlog2") {
  regexpr <- "LFQ[.]intensity[.]([^[:digit:]]+)[[:digit:]]+";
  code <- "LFQ";
} else if (options$type == "intensity") {
  regexpr <- "Intensity[.]([^[:digit:]]+)[[:digit:]]+";
  code <- "INT";
} else {
  regexpr <- "MS[.]MS[.]Count[.]([^[:digit:]]+)[[:digit:]]+";
  code <- "MS";
}

# define the columns that will be taken in account for the t-test
columns_names <- grep(regexpr, colnames(table), value=TRUE);

#TODO make it don't depende on hard coded "C" and "T"
# two samples: control and treatment
control_columns <- columns_names[gsub(regexpr, "\\1", columns_names) == "C"]
treatment_columns <- columns_names[gsub(regexpr, "\\1", columns_names) == "T"]

# this is a filtered table to help with calculations
table_only_columns <- table[c(control_columns, treatment_columns)]

# this loop computes the ttest result for each row
# and adds it to a vector
i <- 1;
anovaresult <- c();
anovasignificant <- c();
for (i in seq(1, nrow(table_only_columns))) {
  # the t-test arguments are the control values vector, the treatment values vector
  # and some extra arguments. var.equal says it's a student t-test with stardard
  # deviations assumed equal. mu=0 sets the hipothesis to be null.
  x <- c(table_only_columns[i, control_columns],
    table_only_columns[i, treatment_columns], recursive=TRUE);
  y <- factor(rep(c("C", "T"),
    c(length(table_only_columns[i, control_columns]),
    length(table_only_columns[i, treatment_columns]))),)
  anovaresult[i] <- oneway.test(x~y, var.equal=TRUE)$p.value;
  if (is.na(anovaresult[i]))
    anovaresult[i] = 1.0
}

# this defines if the p-value returned for each row is significant
anovasignificant[anovaresult <= 0.05] <- "+"
anovasignificant[anovaresult > 0.05] <- ""


# create two extra rows on the table, one for p-values and other
# for siginificance
#TODO: ou colocar perto da intensidade que se refere ou na 3ª coluna
table[paste0("ANOVA.result.", code)] <- NA;
table[paste0("ANOVA.result.", code)] <- anovaresult;
table[paste0("ANOVA.significant.", code)] <- NA;
table[paste0("ANOVA.significant.", code)] <- anovasignificant;




# write out the table
output_handler <- file(options$outputfile_name, "w")
write.table(table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
