#!/usr/bin/env Rscript
# newScatterplot.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: October 2015
#
# Required packages to work: (ggplot2, GGally, getopt, )
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.
require(ggplot2, quietly=TRUE)
require(GGally, quietly=TRUE)
require(getopt, quietly=TRUE)

spec = matrix( c(
    'input', 'i', 1, "character",
    'allvsall', 'a', 1, "character",
    'onevsone', 'o', 1, "character"
), byrow=TRUE, ncol=4);

opt = getopt(spec);
# Load datasets
df <- read.delim(opt$input, header=TRUE, fill=TRUE)
if (!(TRUE %in% grepl("[[:digit:]]+.*[^[:digit:]]+[.]lfq[.]intensity$", colnames(table)))) {
  write("Error: No columns of type lfq.intensity in input table", stderr());
  q(status=3);
}
cols <- grep("[[:digit:]]+.*[^[:digit:]]+[.]lfq[.]intensity$",
             colnames(df), value=TRUE)
newdf <- df[-1, cols]
lm_eqn <- function(df, col1, col2) {
    f <- paste(col1, "~", col2)
    m <- lm(f, data=df)
    eq <- paste("r=", format(summary(m)$r.squared, digit=3))
    as.character(as.expression(eq));
}

p <- ggpairs(newdf,
             # What to include above diagonal
             # list(continuous = "points") to mirror
             upper = list(continuous="smooth", params=c(colour="blue")),
             legends=T,

             # What to include below diagonal
             lower=list(continuous="smooth", params=c(colour="blue")),

             # What to include in the diagonal
             diag=list(axisLabels='show'),

             # How to label inner plots
             # internal, none, show
             axisLabels = "none",
             # Other aes() parameters
             title = "Scatterplot Proteomics", verbose=FALSE, na.rm=TRUE
)

for (i in 1:length(cols)) {
    # Address only the diagonal elements
    # Get plot out of matrix
    inner <- getPlot(p, i, i);
    # Add any ggplot2 settings you want
    inner <- inner + theme(panel.grid = element_blank()) + theme(axis.text.x = element_blank())
    # Put it back into the matrix
    p <- putPlot(p, inner, i, i)

    for (j in 1:length(cols)){
        if((i==1 & j==1)){
            inner <- getPlot(p, i, j)
            inner <- inner + theme(legend.position=c(length(cols)-0.25,0.50))
            p <- putPlot(p, inner, i, j)
        }
        else{
            inner <- getPlot(p, i, j)
            inner <- inner + theme(legend.position="none")
            p <- putPlot(p, inner, i, j)
        }
    }
}
pdf(file=opt$allvsall, height=50, width=50, pointsize=12)
suppressWarnings(print(p));
suppressMessages(dev.off());
pdf(file=opt$onevsone)
for (col1 in cols) {
    for (col2 in cols) {
        if (col1 != col2) {
            p <- ggplot(data=df, aes_string(x=col1, y=col2), na.rm=TRUE) +
                geom_point(colour="blue", size=1) + geom_text(aes(x=25, y=35,
                                                                    label=lm_eqn(df, col1, col2))) +
                geom_smooth(method=lm, se=FALSE);
            suppressWarnings(print(p));
        }
    }
}
suppressMessages(dev.off());
