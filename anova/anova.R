#!/usr/bin/env Rscript

# anova.R
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
source('../R_util/read-utils.R');

args <- get_cmd_options(TRUE);
#define de options input that the code will have

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

# this loop computes the ttest result for each row
# and adds it to a vector
i <- 2;
anovaresult <- c("");
anovasignificant <- c("");
for (i in seq(2, nrow(table_only_columns)+1)) {
  # i make two lists that are going to be the arguments for the anova function
  # the oneway.test. the first list is the all data
  # the second list is a correlation to indentify of which categoty each data
  j<-1;
  x<-c();
  aux1 <- c();
  for (j in 1:length(different_categories)) {
    x <- c(x, table_only_columns[i-1, columns[[j]]], recursive=TRUE)
    aux1 <- c(aux1, length(table_only_columns[i-1, columns[[j]]]))
  }
  y <- factor(rep(different_categories,
    aux1),)
  # i get the p-value for the test aplied on the current row.
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
args$table[paste0("ANOVA.result.", args$code)] <- NA;
args$table[paste0("ANOVA.result.", args$code)] <- anovaresult;
args$table[paste0("ANOVA.significant.", args$code)] <- NA;
args$table[paste0("ANOVA.significant.", args$code)] <- anovasignificant;




# write out the table
output_handler <- file(args$options$outputfile_name, "w")
write.table(args$table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
