#!/usr/bin/Rscript
# Import some required libraries
library('getopt');

# Make an option specification, with the following arguments:
# 1. long flag
# 2. short flag
# 3. argument type: 0: No argument, 1: required, 2: optional
# 4. target data type
option_specification = matrix(c(
    'input', 'f', 1, 'character',
    'output', 'g', 1, 'character'
), byrow=TRUE, ncol=4);

# Parse options
options = getopt(option_specification);

if(!is.null(options$input) && !is.null(options$output)){

    # inputname = "proteins_ranking_by_Beta-Binomial.csv"
    # outputname = "proteins_ranking_by_Beta-Binomial.data"

    csv = read.csv(file=options$input, header=TRUE)

    csv[1:3,]
    size = nrow(csv)
    n_attributes = ncol(csv) -2
    attributes_names = colnames(csv)

    cat("DY",file=options$output)
    cat("\n",file=options$output,append=TRUE)

    cat(size,file=options$output,append=TRUE)
    cat("\n",file=options$output,append=TRUE)

    cat(n_attributes,file=options$output,append=TRUE)
    cat("\n",file=options$output,append=TRUE)

    cat(attributes_names[2:(n_attributes+1)],file=options$output,append=TRUE,sep=";")
    cat("\n",file=options$output,append=TRUE)


    for(i in 1:nrow(csv)){
        cat(as.character(csv[i,1]),file=options$output,append=TRUE)  
        
        for(j in 2:ncol(csv)){
            cat(";",file=options$output,append=TRUE)
            cat(as.character(csv[i,j]), file=options$output, append=TRUE)
            
        }
        
        cat("\n",file=options$output,append=TRUE)
    }


}
