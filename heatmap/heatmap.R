#########################################################
### A) Installing and loading required packages
#########################################################

## if (!require("gplots")) {
##    install.packages("gplots", dependencies = TRUE)
##    library(gplots)
##    }
## if (!require("RColorBrewer")) {
##    install.packages("RColorBrewer", dependencies = TRUE)
##    library(RColorBrewer)
##    }

require(getopt, quietly=TRUE)
require("gplots", quietly=TRUE)
## require(RColorBrewer, quietly=TRUE)
#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

spec = matrix(c(
    'input', 'i', 1, "character",
    'output', 'o', 1, "character"
    ), byrow = TRUE, ncol=4);
opt <- getopt(spec);


data <- read.csv(opt$input, sep="\t", comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-n into a matrix
rownames(mat_data) <- rnames                  # assign row names 


#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# creates a 5 x 5 inch image
pdf(paste(opt$output,".pdf",sep=""),    # create PNG for the heat map        
  pointsize = 8)        # smaller font size

row_distance = dist(mat_data, method = "euclidean")
row_cluster = hclust(row_distance, method = "ward")
col_distance = dist(t(mat_data), method = "euclidean")
col_cluster = hclust(col_distance, method = "ward")

## heatmap(mat_data,
##         Rowv = as.dendrogram(row_cluster),
##         Colv = as.dendrogram(col_cluster),
##         col = my_palette,
##         scale="row",
##         margins=c(10,10),
##         keep.dendro=TRUE
##         )


heatmap.2(mat_data,
  cellnote = mat_data,  # same data set for cell labels
  main = "Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier 
  ## breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="row",     # only draw a row dendrogram
  scale="row",
  Rowv = as.dendrogram(row_cluster),
  Colv = as.dendrogram(col_cluster))


dev.off()               # close the pdf device
