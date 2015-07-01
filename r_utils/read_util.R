# auxiliar script to help with the read of all R scripts
read_function <- function(options) {

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
    if (!(TRUE %in% grepl(regexpr, colnames(table)))) {
      print (sprintf("Error: No columns of type %s in input table", code));
      q(1,save="no");
    }

    # define the columns that will be taken in account for the t-test
    columns_names <- grep(regexpr, colnames(table), value=TRUE);

    # here I extract the different experiment names in an array for easier
    # manipulation, ordering them
    experiment_names <- mixedsort(gsub(".*[.]([^[:digit:]]+[[:digit:]]+).*", "\\1",
                                        columns_names));

    # extract from the experiment names all the different categories in the table
    different_categories <- unique(gsub("([^[:digit:]]+).*", "\\1",
                                        experiment_names));

    read_list <- list(table=table, regex=regexpr, code=code, col_names=columns_names, ex_names=experiment_names, diff_cat=different_categories);

    return(read_list);
}
