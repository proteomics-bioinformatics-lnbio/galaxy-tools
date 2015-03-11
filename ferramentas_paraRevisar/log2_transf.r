#!/usr/bin/Rscript
# Import some required libraries
library('getopt');
library('scatterplot3d');
 
# Make an option specification, with the following arguments:
# 1. long flag
# 2. short flag
# 3. argument type: 0: No argument, 1: required, 2: optional
# 4. target data type
option_specification = matrix(c(
  'infilename', 'i', 1, 'character',
  'outfilename', 'o', 2, 'logical',
), byrow=TRUE, ncol=4);
 
# Parse options
options = getopt(option_specification);

normalization = function(infilename, head) {

## Read the data table and create a vector with the columns of interest
table = read.table(infilename, header=TRUE, sep="\t", dec=".")

head_input = colnames(table)
coluns = grep(head, head_input, perl=TRUE, value=TRUE)

data = as.matrix(table[,coluns])

## Makes the log2 transformation
l2 = log2(data)
l2 = signif(l2, 7)

	col_names = c("Intensity.txt")

	nor_list = list(data, l2)

	# Creates a output files with the data Log2 transformed columns #
	for (i in 1:7) {
		outfilename = file(col_names[i])
		write.table(nor_list[i], outfilename, dec=".", sep="\t", row.names=FALSE, col.names=TRUE)
	}


}


