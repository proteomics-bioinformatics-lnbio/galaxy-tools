#!/usr/bin/env Rscript
# heatmap.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: November 2014
#
# Required packages to work: (getopt, )
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.
require(getopt, quietly=TRUE)
# Specification of how arguments will work
spec = matrix(c(
    'input', 'i', 1, "character",
    'type', 't', 1, "character",
    'output', 'o', 1, "character"
    ), byrow = TRUE, ncol=4);
opt <- getopt(spec);

# sets the type the user wants for the heatmap
if (!is.null(opt$type)) {
    type = opt$type;
    # sees if the type is known if not, halt
    if (type != "lfqlog2" && type != "intensity"
        && type != "mscount") {
        write("Error, unknown type", stderr());
        q(status=3);
    }
} else {
    # default type
    type = "lfqlog2";
}

table <- read.delim(file=opt$input, header=TRUE, fill=TRUE);

if (type == "lfqlog2") {
    # this stores an array for the collumn names that has a pattern like
    # "LFQ.intensity.[1 or more non numbers][1 or more numbers]"
    lfq_collumns_names <- grep("[^[:digit:]]+[[:digit:]]+[.]lfq[.]intensity",
                               colnames(table), value=TRUE);
    table <- scale(data.matrix(table[2:nrow(table), lfq_collumns_names]));
} else if (type == "intensity") {
    # this stores an array for the collumn names with a pattern as
    # "Intensity.[1 or more non numbers][1 or more numbers]"
    intensitynames<- grep("[^[:digit:]]+[[:digit:]]+[.]intensity",
                               colnames(table), value=TRUE);
    table <- scale(data.matrix(table[2:nrow(table), intensitynames]));
} else {
    # this stores an array for the collumn names with a pattern as
    # "MS.MS.Count.[1 or more non numbers][1 or more numbers]"
    mscountnames<- grep("[^[:digit:]]+[[:digit:]]+[.]speccount",
                               colnames(table), value=TRUE);
    table <- scale(data.matrix(table[2:nrow(table), mscountnames]));
}

pdf(file=opt$output,
    pointsize = 6);
invisible(heatmap(table, Rowv=NA, Colv=NA, na.rm=TRUE, keep.dendro=FALSE, verbose=FALSE));
dev.off();
q(status=0);
