#!/usr/bin/Rscript
# Import some required libraries
library('getopt')
library('scatterplot3d');
 
# Make an option specification, with the following arguments:
# 1. long flag
# 2. short flag
# 3. argument type: 0: No argument, 1: required, 2: optional
# 4. target data type
option_specification = matrix(c(
  'infilename', 'i', 2, 'character',
  'pdffile', 'f', 2, 'character',
  'pngfile', 'g', 2, 'character'
), byrow=TRUE, ncol=4);
 
# Parse options
options = getopt(option_specification);
 
# Create some simple test data

table = read.table("Correlation_CID_2scatter_18072014.txt", sep="\t")
GradientTime = as.integer(table[,1])
scatterplot3d ("table")
 
# Produce PDF file
if (!is.null(options$pdffile)) {
  pdf(options$pdffile);
  plot(x, y);
  dev.off();
}
 
# Produce PNG file
if (!is.null(options$pngfile)) {
  png(options$pngfile);
  plot(x, y);
  dev.off();
}
