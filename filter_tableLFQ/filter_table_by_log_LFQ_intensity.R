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

#' Replaces every occurrence of "-Inf", which is the mathematical
#' result of log2(0), with the programming term NA.
#'
#' @param table The table that will have it's terms replaced
#' @return replaced_table The result of the replacing
#' @examples
#' replace_with_NA(c(c("-Inf", 0, "-Inf"), c(1, 2, 1.34), c(2.666, "-Inf", 3)))
#' returns:
#'      [,1]  [,2]  [,3]
#'[1,]   NA    0   NA
#'[2,]    1     2   1.34
#'[3,]  2.666  NA  3
replace_with_NA <- function(table) {
    replaced_table <- table
    for(row in seq(1, nrow(replaced_table))) {
        range <- grep("LFQ[.]intensity[.][^[:digit:]]+[[:digit:]]+", names(replaced_table))
        for (collumn in seq(1, ncol(replaced_table))) {
            if (collumn %in% range) {
                if (grepl("-Inf", replaced_table[row, collumn])==TRUE) {
                    replaced_table[row, collumn] <- "NA"
                }
            }
        }
    }
    return(replaced_table)
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

# this stores an array for the collumn names that has a pattern like
# "LFQ.intensity.[1 or more non numbers][1 or more numbers]"
lfq_column_names <- grep("LFQ[.]intensity[.][^[:digit:]]+[[:digit:]]+",
                            colnames(table), value=TRUE);
# here I extract the different experiment names in an array for easier
# manipulation, ordering them
experiment_names <- mixedsort(gsub(".*[.]([^[:digit:]]+[[:digit:]]+).*", "\\1",
                                    lfq_column_names));
# extract from the experiment names all the different categories in the table
different_categories <- unique(gsub("([^[:digit:]]+).*", "\\1",
                                    experiment_names));
# get the maximum of NA's permitted by the filter
max_NAs_allowed <- get_max_log_NA(experiment_names, different_categories);
#TODO NEED TO SET ALL LFQ COLUMNS AS NUMERIC FOR LOG TO WORK
# calculate the log on base 2 of the LFQ.intensity column of the input table
first_row <- table[1,];
number_of_rows <- nrow(table[-1, lfq_column_names]);
table_log_LFQ_intensity <- lapply(table[-1,lfq_column_names], function(x) { log(as.numeric(as.character(x)), base=2)})
# count the number of NA logs in the log table, to get
# which rows are going to be filtered out.
# initial definitions
table_log_LFQ_intensity <- data.frame(matrix(unlist(table_log_LFQ_intensity), nrow=number_of_rows, byrow=TRUE), stringsAsFactors=FALSE)
colnames(table_log_LFQ_intensity) <- lfq_column_names;

i<-1;
counter <- c();
# filter first for each category that presents more than the max
# allowed number of NA's logs
for (category in different_categories) {
    # iterate in the rows of this table
    range <- grep(paste0("LFQ[.]intensity[.]", category, "[[:digit:]]+.*"),colnames(table_log_LFQ_intensity), value=TRUE);
    counter.tmp <- as.vector(by(table_log_LFQ_intensity[,range], seq(1, nrow(table_log_LFQ_intensity)),  function(x) { sum(grepl("-Inf", x))>=max_NAs_allowed }));
    if (length(counter) == nrow(table_log_LFQ_intensity)) {
        counter <- counter.tmp & counter;
    } else {
        counter <- counter.tmp;
    }
}
# sort the array of excluded rows for better organization

# the array of rows that are included, that is, that does not have more
# than the max number of NAs allowed, have every row that belongs to the
# intervall between 1 and the number of rows in the table, except the ones
# that are in the excluded_rows array
included_rows <- !counter;

# the filtered table is the one that has only the included rows
included_rows_right <- c(TRUE, included_rows);
filtered_table <- table[included_rows_right,];
number_of_rows <- nrow(filtered_table);
# count the number of NA logs in the log table, to get
# which rows are going to be filtered out.
# initial definitions
filtered_table[-1,lfq_column_names] <- table_log_LFQ_intensity[included_rows,];

print(included_rows_right);
print(seq(1, length(included_rows_right)));
print(seq(1, length(included_rows_right))[included_rows_right])

# and the "-Inf" replaced with NA
filtered_table <- replace_with_NA(filtered_table);
# write the final table in a .txt file separated by tabs, just as the input
# file was!
output_handler <- file(options$outputfile_name, "w")

write.table(filtered_table, file=output_handler, sep="\t", row.names=FALSE);

close(output_handler)
