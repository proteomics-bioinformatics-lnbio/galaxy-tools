#!/usr/bin/Rscript

# filter_table_by_log_LFQ_intensity.R
# AUTHOR: Daniel Travieso (modified from script by MY Friend)
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: September 2015
#
# Required packages to work: ("gtools", "getopt", )
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.

# DECLARATION OF FUNCTIONS USED IN THE SCRIPT

#' Gets the max number of NA logs permitted by the filter
#' Which is the number of experiment the category that
#' has the least number of experiments, divided by 2
#' technically: floor(min(num_experiments_per_categorie)/2)
#'
#' @param experiments The array of experiment names
#' @param categories The array of categories names
#' @return truncated half of the minimum of experiments per categories
#' @examples
#' get_max_log_NA(c("C1", "C2", "T1"), c("C", "T"))
#' returns: floor(1/2) that is: 0;
get_max_log_NA <- function(experiments, categories) {
    i<- 1
    all_sums <- array()
    for (category in categories) {
        all_sums[i] <- sum(grepl(category, experiments)==TRUE)
        i<-i+1
    }
    return(as.integer(min(all_sums)/2))
}

# HERE IS THE START OF THE CODE ITSELF
# package gtools to use "mixedsort" function
require("gtools", quietly=TRUE);
require('getopt', quietly=TRUE);

opt <- matrix(c(
    'inputfile_name', 'i', 1, 'character',
    'outputfile_name', 'o', 1, 'character'
    ),byrow=TRUE, ncol=4);

options <- getopt(opt);

table <- read.delim(options$inputfile_name, header=TRUE, fill=TRUE, stringsAsFactors=FALSE);

if (!(TRUE %in% grepl("[[:digit:]]+.*[^[:digit:]]+[.]lfq[.]intensity$", colnames(table)))) {
  write("Error: No columns of type lfq.intensity in input table", stderr());
  q(status=3);
}
# this stores an array for the collumn names that has a pattern like
# "LFQ.intensity.[1 or more non numbers][1 or more numbers]"
lfq_column_names <- grep("[[:digit:]]+.*[^[:digit:]]+[.]lfq[.]intensity$",
                            colnames(table), value=TRUE);
                           }
# here I extract the different experiment names in an array for easier
# manipulation, ordering them
experiment_names <- mixedsort(gsub("([[:digit:]]+.*[^[:digit:]]+)[.].*$", "\\1",
                                    lfq_column_names));
# extract from the experiment names all the different categories in the table
different_categories <- unique(gsub("[[:digit:]]+.*([^[:digit:]]+).*$", "\\1",
                                    experiment_names));
# get the maximum of NA's permitted by the filter
max_NAs_allowed <- get_max_log_NA(experiment_names, different_categories);

#vector that contains all excluded row numbers
excluded_rows <- c();
#iterate in rows
for (row in seq(2, nrow(table))) {
  #nacount is individual for each category
  for (category in different_categories) {
    nacount <- 0;
    for(exp in experiment_names) {
      #current col name
      col <- paste0(exp, ".lfq.intensity");
      if (gsub("[[:digit:]]+.*([^[:digit:]]+).*$", "\\1", exp) == category) {
        #the log fo the number on the current cell
        logresult <- log(as.numeric(as.character(table[row, col])), base=2);
        #if result was -Inf, the log is NA
        if (logresult == "-Inf") {
          logresult<-"NA";
          nacount <- nacount+1;
        }
        table[row, col] <- logresult;
      }
    }
    #check nacount if surpasses limit
    if (nacount >= max_NAs_allowed) {
      excluded_rows <- c(excluded_rows, row);
    }
  }
}

#sort an make unique the values on the excluded rows
excluded_rows <- mixedsort(unique(excluded_rows))

#included rows are all non excluded
included_rows <- seq(2, nrow(table))[seq(2, nrow(table)) %in% excluded_rows == FALSE];

#refactor table to contain only filtered rows
table <- table[c(1, included_rows),];

# write the final table in a .txt file separated by tabs, just as the input
# file was!
output_handler <- file(options$outputfile_name, "w")

write.table(table, file=output_handler, sep="\t", row.names=FALSE);

close(output_handler)
