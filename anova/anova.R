#!/usr/bin/env Rscript

# anova.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: October 2015
#
# Required packages to work: ("getopt", "gtools")
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.

require('gtools', quietly=TRUE);
require('getopt', quietly=TRUE);
opt <- matrix(c(
    'inputfile_name', 'i', 1, 'character',
    'type', 't', 1, 'character',
    'outputfile_name', 'o', 1, 'character'
),byrow=TRUE, ncol=4);

options <- getopt(opt);

table <- read.delim(options$inputfile_name, header=TRUE, fill=TRUE);

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
  write(sprintf("Error: No columns of type %s in input table", code), stderr());
  q(status=3);
}
#define de options input that the code will have

# define the columns that will be taken in account for the anova
columns_names <- grep(regexpr, colnames(table), value=TRUE);

# here I extract the different experiment names in an array for easier
# manipulation, ordering them
experiment_names <- mixedsort(gsub("([^[:digit:]]+.*[[:digit:]]+)[.].*", "\\1",
                                    columns_names));

# extract from the experiment names all the different categories in the table
different_categories <- unique(gsub("([^[:digit:]]+).*", "\\1",
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

# this loop computes the anova result for each row
# and adds it to a vector
i <- 1;
anovaresult <- c("");
anovasignificant <- c("");
for (i in seq(1, nrow(table_only_columns))) {
  # i make two lists that are going to be the arguments for the anova function
  # the oneway.test. the first list is the all data
  # the second list is a correlation to indentify of which categoty each data
  j<-1;
  x<-c();
  aux1 <- c();
  for (j in 1:length(different_categories)) {
    x <- c(x, table_only_columns[i, columns[[j]]], recursive=TRUE)
    aux1 <- c(aux1, length(table_only_columns[i, columns[[j]]]))
  }
  y <- factor(rep(different_categories,
    aux1),)
  # i get the p-value for the test aplied on the current row.
  anovaresult[i+1] <- oneway.test(x~y, var.equal=TRUE)$p.value;
  if (is.na(anovaresult[i+1]))
    anovaresult[i+1] = 1.0
}

# this defines if the p-value returned for each row is significant
anovasignificant[as.numeric(anovaresult) <= 0.05] <- "+"
anovasignificant[as.numeric(anovaresult) > 0.05] <- ""
anovasignificant[1] <- "";


# create two extra rows on the table, one for p-values and other
# for siginificance
table[paste0("ANOVA.result.", code)] <- NA;
table[paste0("ANOVA.result.", code)] <- anovaresult;
table[paste0("ANOVA.significant.", code)] <- NA;
table[paste0("ANOVA.significant.", code)] <- anovasignificant;


# write out the table
output_handler <- file(options$outputfile_name, "w")
write.table(table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
