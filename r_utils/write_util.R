writeout <- function (filename, table) {
    output_handler <- file(filename, "w")
    write.table(table, file=output_handler, sep="\t", row.names=FALSE);
    close(output_handler)
}
