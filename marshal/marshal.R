#!/usr/bin/env Rscript

# marshal.R
# AUTHOR: Daniel Travieso
# E-mail: danielgtravieso@gmail.com
# LAST REVISED: August 2015
#
# Required packages to work: ("getopt")
# Laboratory of Mass Spectrometry at Brazilian Biosciences National Laboratory
# http://lnbio.cnpem.br/
# Copyright CC BY-NC-SA (c) 2014  Brazilian Center for Research in Energy and Materials
# All rights reserved.
require('getopt', quietly=TRUE);
require('RMySQL', quietly=TRUE);

#define de options input that the code will have
opt = matrix(c(
    'inputfile_name', 'i', 1, 'character',
    'tax_id', 't', 1, 'character',
    'outputfile_name', 'o', 1, 'character'
),byrow=TRUE, ncol=4);

# parse de input
options = getopt(opt);

table <- read.delim(options$inputfile_name, header=TRUE, fill=TRUE);

#Database connections and queries
db.connection <- dbConnect(RMySQL::MySQL(), user='galaxy', host='localhost', dbname='conversionMarcelo', password='123456', unix.sock='/tmp/mysql.sock');

# the '?' will be replaced in the query
db.sql.synonym <- "SELECT synonyms FROM Synonyms2Uniprot WHERE uniprot = %s";
db.sql.uniprot <- "SELECT uniprot FROM Synonyms2Uniprot WHERE synonyms = %s";
db.sql.all <- "SELECT * FROM Synonyms2Uniprot WHERE synonyms = %s";


#Definition of all regular expressions to be used
regex.intensity <- "^Intensity[.]([^[:digit:]]+)[[:digit:]]+$";
regex.lfqintensity <- "^LFQ[.]intensity[.]([^[:digit:]]+)[[:digit:]]+$";
regex.spectral <- "^MS[.]MS[.]Count[.]([^[:digit:]]+)[[:digit:]]+$";
regex.proteinIDs <- "^Protein[.]IDs$";
regex.id.uniprot.1 <- "^[OPQ][[:digit:]][[:upper:][:digit:]]{3}[[:digit:]]";
regex.id.uniprot.2 <- "^[A-NR-Z][[:digit:]]([[:upper:]][[:upper:][:digit:]]{2}[[:digit:]]){1,2}";
regex.id.ipi <- "^IPI[[:digit:]]+$";
regex.id.tair <- "^A[Tt][1-5MC][:alnum:][[:digit:]]+$";
regex.id.ensembl <- "^ENS[[:alnum:]]+$";
regex.id.refseq <- "^AC_|^N[CGTWSZMR]_|^X[MR]_|^[ANYXZ]P_";
regex.id.contaminant_reversed <- "^CON__|REV__";

#get the names for the columns and separate them for better use
column_names.intensity <- grep(regex.intensity, colnames(table), value=TRUE);
column_names.lfqintensity <- grep(regex.lfqintensity, colnames(table), value=TRUE);
column_names.spectral <- grep(regex.spectral, colnames(table), value=TRUE);
column_names.proteinIDs <- grep(regex.proteinIDs, colnames(table), value=TRUE);
column_names.uniprot_conversion <- "Uniprot.Conversion";

#create names of categories for extra row to be included
categories.intensity <- gsub(regex.intensity, '\\1.intensity', column_names.intensity);
categories.lfqintensity <- gsub(regex.lfqintensity, '\\1.lfq.intensity.', column_names.lfqintensity);
categories.spectral <- gsub(regex.spectral, '\\1.speccount', column_names.spectral);

#add a blank column
table <- cbind(table[,c(column_names.proteinIDs)], rep("", nrow(table)), table[,c(column_names.intensity, column_names.lfqintensity, column_names.spectral)])
#rename the column names
colnames(table) <- c(column_names.proteinIDs,column_names.uniprot_conversion,column_names.intensity, column_names.lfqintensity, column_names.spectral)
#set the protein IDs column as character
table[,column_names.proteinIDs] <- sapply(table[,column_names.proteinIDs], as.character)
#add a blank row
table <- rbind(table[0,], c("", "", categories.intensity, categories.lfqintensity, categories.spectral), table[seq(1, nrow(table)),])

cell.tax <- options$tax_id;
#finding the id type
for (row in seq(1, nrow(table))) {
    cell.row <- row;
    cell.value <- strsplit(table[cell.row, column_names.proteinIDs], ';')[[1]];
	cell.id <- cell.value[1];
    for (id_code in cell.value) {
        # remove contaminant or reversed protein ids from search. get the first valid id
        if(!grepl(regex.id.contaminant_reversed, id_code)) {
            cell.id <- id_code;
        }
    }
    # discover the type of id
    if (grepl(regex.id.uniprot.1, cell.id) || grepl(regex.id.uniprot.2, cell.id)) {
        cell.hash <- 'uniprot';
    } else if(grepl(regex.id.ipi, cell.id)) {
        cell.hash <- 'ipi';
    } else if(grepl(regex.id.tair, cell.id)) {
        cell.hash <- 'tair';
    } else if(grepl(regex.id.ensembl, cell.id)) {
        cell.hash <- 'ensembl';
    } else if(grepl(regex.id.refseq, cell.id)) {
        cell.hash <- 'refseq'
    } else {
        cell.hash <- 'genesymbol'
    }
    # write the id in the first row of the new table
    table[cell.row, column_names.proteinIDs] <- cell.id;
    # if the id is not uniprot type, get the correlant uniprot for that id
    if (cell.hash != 'uniprot') {
        # if is genesymbol use the database to search for a uniprot with same tax number
        if (cell.hash == 'genesymbol') {
			db.select.all <- dbEscapeStrings(db.connection, sprintf(db.sql.all, cell.id));
			print(db.select.all);
            db.select.results <- fetch(dbSendQuery(db.connection, db.select.all), n=-1);
            for (row in seq(1, nrow(db.select.results))) {
                if (db.select.results[row, 4] == cell.tax) {
                    cell.id.uniprot <- db.select.results[row, 2];
                    if (db.select.results[row, 5] == "YES") {
                        break;
                    }
                }
            }
            # if is not genesymbol, query for all ids in the table, and get the result uniprot
        } else {
            db.select.uniprot <- dbEscapeStrings(db.connection, sprintf(db.sql.uniprot, cell.id));
			print(db.select.uniprot);
            db.select.results <- fetch(dbSendQuery(db.connection, db.select.uniprot), n=-1);
            for (row in seq(1, nrow(db.select.results))) {
                if (grepl(db.select.results[row], regex.id.uniprot.1) == TRUE ||
                grepl(db.select.results[row], regex.id.uniprot.2) == TRUE) {
                    cell.id.uniprot <- db.select.results[row];
                    break;
                }
            }
        }
    } else {
        # if the id is already uniprot, do nothing and store the uniprot id value
        cell.id.uniprot <- cell.id;
        cell.id <- "";
        cell.id.possibleid <- c();
        cell.id.uniprot <- sub('-[[:digit:]]', '', cell.id.uniprot);
		db.select.synonym <- dbEscapeStrings(db.connection, sprintf(db.sql.synonym, cell.id.uniprot));
		print(db.select.synonym);
        db.select.results <- fetch(dbSendQuery(db.connection, db.select.synonym), n=-1);
        for (row in seq(1, nrow(db.select.results))) {
            cell.id.possibleid <- c(cell.id.possibleid, db.select.results[row]);
            if (grepl('^IPI|^ENS|^A[Tt]',db.select.results[row])) {
                cell.id <- db.select.results[row];
                break;
            }
        }
        if (cell.id == "") {
            if (cell.id.possibleid) {
                cell.id <- sort(cell.id.possibleid)[1]
            }
        }
    }
    # write the uniprot conversion cell
    table[cell.row, column_names.uniprot_conversion] <- paste0(cell.id, "_", cell.id.uniprot);
}

#write out the file
output_handler <- file(options$outputfile_name, "w")
write.table(table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
