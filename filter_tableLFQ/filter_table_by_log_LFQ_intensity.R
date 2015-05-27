#!/usr/bin/Rscript

# filter_table_by_log_LFQ_intensity.R
# AUTHOR: Daniel Travieso (modified from script by MY Friend)
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: November 2014
#
# Required packages to work: ("gtools", "getopt", )
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.

# DECLARATION OF FUNCTIONS USED IN THE SCRIPT

#' Gets the max number of NaN logs permitted by the filter
#' Which is the number of experiment the categorie that
#' has the least number of experiments, divided by 2
#' technically: floor(min(num_experiments_per_categorie)/2)
#'
#' @param experiments The array of experiment names
#' @param categories The array of categories names
#' @return truncated half of the minimum of experiments per categories
#' @examples
#' get_max_log_NaN(c("C1", "C2", "T1"), c("C", "T"))
#' returns: floor(1/2) that is: 0;
get_max_log_NaN <- function(experiments, categories) {
    i<- 1
    all_sums <- array()
    for (categorie in categories) {
        all_sums[i] <- sum(grepl(categorie, experiments)==TRUE)
        i<-i+1
    }
    return(as.integer(min(all_sums)/2))
}

#' Replaces every occurrence of "-Inf", which is the mathematical
#' result of log2(0), with the programming term NaN.
#'
#' @param table The table that will have it's terms replaced
#' @return replaced_table The result of the replacing
#' @examples
#' replace_with_NaN(c(c("-Inf", 0, "-Inf"), c(1, 2, 1.34), c(2.666, "-Inf", 3)))
#' returns:
#'      [,1]  [,2]  [,3]
#'[1,]   NaN    0   NaN
#'[2,]    1     2   1.34
#'[3,]  2.666  NaN  3
replace_with_NaN <- function(table) {
    replaced_table <- table
    for(row in seq(1, nrow(replaced_table))) {
        range <- grep("LFQ[.]intensity[.][^[:digit:]]+[[:digit:]]+", names(replaced_table))
        for (collumn in seq(1, ncol(replaced_table))) {
            if (collumn %in% range) {
                if (grepl("-Inf", replaced_table[row, collumn])==TRUE) {
                    replaced_table[row, collumn] <- "NaN"
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
source('../R_util/read-utils.R', quietly=TRUE);

args <- get_cmd_options(FALSE);
# this stores an array for the collumn names that has a pattern like
# "LFQ.intensity.[1 or more non numbers][1 or more numbers]"
lfq_collumns_names <- grep("LFQ[.]intensity[.][^[:digit:]]+[[:digit:]]+",
                            colnames(args$table), value=TRUE);

# here I extract the different experiment names in an array for easier
# manipulation, ordering them
experiment_names <- mixedsort(gsub(".*[.]([^[:digit:]]+[[:digit:]]+).*", "\\1",
                                    lfq_collumns_names));

# extract from the experiment names all the different categories in the table
different_categories <- unique(gsub("([^[:digit:]]+).*", "\\1",
                                    experiment_names));

# get the maximum of NaN's permitted by the filter
max_NaNs_allowed <- get_max_log_NaN(experiment_names, different_categories);

# calculate the log on base 2 of the LFQ.intensity collumns of the input table
table_log_LFQ_intensity <- log(args$table[lfq_collumns_names], base=2);

# count the number of NaN logs in the log table, to get
# which rows are going to be filtered out.

# initial definitions
i<-1;
excluded_rows <- array();

# filter first for each categorie that presents more than the max
# allowed number of NaN's logs
for (categorie in different_categories) {
    # iterate in the rows of this table
    for(row in seq(1, nrow(table_log_LFQ_intensity))) {
        # reset the counter for this row
        NaN_counter <- 0;

        # here this range defines which categorie or group of
        # lfq intensities will be checked if there is or not
        # a number greater than the max allowed NaNs on the logs.
        range <- grep(paste0("LFQ[.]intensity[.]", categorie, "[[:digit:]]+.*"),
                        names(table_log_LFQ_intensity));

        # iterate in the collumns of the table
        for (collumn in seq(1, ncol(table_log_LFQ_intensity))) {
            # check if current collumn corresponds the the range of
            # log LFQ values being investigated
            if (collumn %in% range) {
                # if the log result was "-Inf", it is a NaN, so the counter gets
                # incremented
                if (grepl("-Inf", table_log_LFQ_intensity[row, collumn])) {
                    NaN_counter <- NaN_counter + 1;
                }
            }

            # to reserve some time, check if the counter already broke the
            # maximum allowed, if so, stop the counting
            if(NaN_counter >= max_NaNs_allowed) {
                break;
            }
        }

        # after all the counting, for this row, check if the counter
        # broke the max, and if so insert the corresponding row on the
        # array of excluded rows.. Also check if it's not there already
        if(NaN_counter >= max_NaNs_allowed && !(row %in% excluded_rows)) {
            excluded_rows[i] <- row;
            i<- i+1;
        }
    }
}

# sort the array of excluded rows for better organization
excluded_rows <- mixedsort(excluded_rows);

# the array of rows that are included, that is, that does not have more
# than the max number of NaNs allowed, have every row that belongs to the
# intervall between 1 and the number of rows in the table, except the ones
# that are in the excluded_rows array
included_rows <- seq(1, nrow(args$table))[seq(1, nrow(args$table)) %in% excluded_rows == FALSE];

# the filtered table is the one that has only the included rows
filtered_table <- args$table[included_rows,];
# with the log applied to it's LFQ collumns
filtered_table[lfq_collumns_names] <- log(filtered_table[lfq_collumns_names], base=2);
# and the "-Inf" replaced with NaN
filtered_table <- replace_with_NaN(filtered_table);
# write the final table in a .txt file separated by tabs, just as the input
# file was!

output_handler <- file(args$options$outputfile_name, "w")

write.table(filtered_table, file=output_handler, sep="\t", row.names=FALSE);

close(output_handler)
