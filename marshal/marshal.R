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
require('methods', quietly=TRUE);
require('getopt', quietly=TRUE);
require('RMySQL', quietly=TRUE);

#define de options input that the code will have
opt = matrix(c(
    'inputfile_name', 'i', 1, 'character',
    'tax_id', 't', 1, 'character',
    'keepcon', 'k', 1, 'character',
    'outputfile_name', 'o', 1, 'character'
),byrow=TRUE, ncol=4);

# parse de input
options = getopt(opt);

table <- read.delim(options$inputfile_name, header=TRUE, fill=TRUE);
db.socket <- '/var/run/mysqld/mysqld.sock';
#Database connections and queries
db.connection <- dbConnect(RMySQL::MySQL(), user='galaxy', host='localhost', dbname='conversionMarcelo', password='123456', unix.sock=db.socket);

# the '?' will be replaced in the query
db.sql.synonym <- "SELECT synonyms FROM Synonyms2Uniprot WHERE uniprot = ";
db.sql.uniprot <- "SELECT uniprot FROM Synonyms2Uniprot WHERE synonyms = ";
db.sql.all <- "SELECT * FROM Synonyms2Uniprot WHERE synonyms = ";


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
regex.id.contaminant_reversed <- "^(CON__|REV__)";
regex.id.contaminant <- "^CON__";

#get the names for the columns and separate them for better use
column_names.intensity <- grep(regex.intensity, colnames(table), value=TRUE);
column_names.lfqintensity <- grep(regex.lfqintensity, colnames(table), value=TRUE);
column_names.spectral <- grep(regex.spectral, colnames(table), value=TRUE);
column_names.proteinIDs <- grep(regex.proteinIDs, colnames(table), value=TRUE);
column_names.uniprot_conversion <- "Uniprot.Conversion.ID";

#create names of categories for extra row to be included
categories.intensity <- gsub(regex.intensity, '\\1.intensity', column_names.intensity);
categories.lfqintensity <- gsub(regex.lfqintensity, '\\1.lfq.intensity.', column_names.lfqintensity);
categories.spectral <- gsub(regex.spectral, '\\1.speccount', column_names.spectral);

#add a blank column
table <- cbind(table[,c(column_names.proteinIDs)], rep("", nrow(table)), table[,c(column_names.intensity, column_names.lfqintensity, column_names.spectral)])
#rename the column names
colnames(table) <- c(column_names.proteinIDs,column_names.uniprot_conversion,column_names.intensity, column_names.lfqintensity, column_names.spectral)
#set the protein IDs column as character
table[,column_names.proteinIDs] <- sapply(table[,column_names.proteinIDs], as.character);
#add a blank row
table <- rbind(table[0,], c("", "", categories.intensity, categories.lfqintensity, categories.spectral), table[seq(1, nrow(table)),]);
table[,column_names.uniprot_conversion] <- sapply(table[,column_names.uniprot_conversion], as.character);
cell.tax <- options$tax_id;
#finding the id type
rows.excluded <- c();
for (row in seq(2, nrow(table))) {
    row.deleted <- FALSE;
    cell.id.notfound <- FALSE;
    cell.row <- row;
    cell.value <- strsplit(table[cell.row, column_names.proteinIDs], ';')[[1]];
                                        #print("Current table row:");
                                        #print(cell.value);
    cell.id <- "";
    for (id_code in cell.value) {
                                        # remove contaminant or reversed protein ids from search. get the first valid id
        if(!is.na(id_code) && !grepl(regex.id.contaminant_reversed, id_code)) {
            cell.id <- id_code;
            break;
        }
    }
    if (cell.id == "") {

        cell.id.notfound <- TRUE;
        # NOT FOUND ANY RELEVANT IDs (All contaminat or reversed)
        if (options$keepcon == "Yes") {
            cell.id <- grep(regex.id.contaminant, cell.value, value=TRUE)[1];
            if (is.na(cell.id)) {
                rows.excluded <- c(rows.excluded, row);
                row.deleted <- TRUE;
            } else {
                cell.id <- gsub(regex.id.contaminant, '', cell.id);
            }
        } else {
            rows.excluded <- c(rows.excluded, row);
            row.deleted <- TRUE;
        }
    }
    if (!row.deleted) {
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
        # if the id is not uniprot type, get the correlant uniprot for that id
        if (cell.hash != 'uniprot') {
            # if is genesymbol use the database to search for a uniprot with same tax number
            if (cell.hash == 'genesymbol') {
                cell.id.escapeChars <- dbEscapeStrings(db.connection, cell.id);
                db.select.all <- paste(db.sql.all, "'", cell.id.escapeChars, "';", sep='');
                db.select.results <- dbGetQuery(db.connection, db.select.all);

                for (rrow in seq(1, nrow(db.select.results))) {
                    db.select.results.row <- db.select.results[row,];

                    if (db.select.results.row[[4]] == cell.tax) {
                        cell.id.uniprot <- db.select.results.row[[2]];
                        if (db.select.results.row[[5]] == "YES") {
                            break;
                        }
                    }
                }
                # if is not genesymbol, query for all ids in the table, and get the result uniprot
            } else {
                cell.id.escapeChars <- dbEscapeStrings(db.connection, cell.id);
                db.select.uniprot <- paste(db.sql.uniprot, "'", cell.id.escapeChars, "';", sep='');
                db.select.results <- dbGetQuery(db.connection, db.select.uniprot);
                db.select.results.ids <- as.vector(db.select.results[[1]]);
                for (res in db.select.results.ids) {
                    if (grepl(res, regex.id.uniprot.1) == TRUE ||
                        grepl(res, regex.id.uniprot.2) == TRUE) {
                        cell.id.uniprot <- res;
                        break;
                    }
                }
            }
            # write the uniprot conversion cell
            if (cell.id.notfound) {
                cell.id <- paste0("CON__", cell.id);
            }
            table[cell.row, column_names.proteinIDs] <- cell.id;
            table[cell.row, column_names.uniprot_conversion] <- paste0(cell.id, "_", cell.id.uniprot);
        } else {
            # if the id is already uniprot, do nothing and store the uniprot id value
            cell.id.uniprot <- cell.id;
            cell.id <- "";
            cell.id.possibleid <- c();
            cell.id.uniprot <- sub('-[[:digit:]]', '', cell.id.uniprot);
            cell.id.escapeChars <- dbEscapeStrings(db.connection, cell.id.uniprot);
            db.select.synonym <- paste(db.sql.synonym, "'", cell.id.escapeChars, "';", sep='');
            db.select.results <- dbGetQuery(db.connection, db.select.synonym);
            db.select.results.ids <- as.vector(db.select.results[[1]]);
            for (res in db.select.results.ids) {
                cell.id.possibleid <- c(cell.id.possibleid, res);
                if (!grepl('^IPI|^ENS|^A[Tt]',res)) {
                    cell.id <- res;
                    break;
                }
            }
            if (cell.id == "") {
                if (!is.null(cell.id.possibleid)) {
                    cell.id <- sort(cell.id.possibleid)[1]
                } else {
                    cell.id <- cell.id.uniprot;
                }
            }
            # write the uniprot conversion cell
            if (cell.id.notfound) {
                table[cell.row, column_names.proteinIDs] <- paste0("CON__", cell.id.uniprot);
                table[cell.row, column_names.uniprot_conversion] <- paste0("CON__", cell.id, "_", cell.id.uniprot);
            } else {
                table[cell.row, column_names.proteinIDs] <- cell.id.uniprot;
                table[cell.row, column_names.uniprot_conversion] <- paste0(cell.id, "_", cell.id.uniprot);
            }
        }
    }
}
table[-1, -c(1, 2)] <- apply(table[-c(1),-c(1, 2)], 2, as.numeric);
table <- table[-rows.excluded,];
#dbDisconnect(db.connection);
#write out the file
output_handler <- file(options$outputfile_name, "w")
write.table(table, file=output_handler, sep="\t", row.names=FALSE);
close(output_handler)
