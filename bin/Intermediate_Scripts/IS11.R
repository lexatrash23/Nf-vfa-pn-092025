#!/usr/bin/env Rscript
arg <- commandArgs(trailingOnly = TRUE)
#read in command line arguments
# read in sample_name
Sample <- arg[1]
#read in blastx0 file
blastx <- arg[2]
blastx_file <- readLines(blastx)

#obtain database information from blastx file
blastx_db_line <- grep("^Database:", blastx_file)
blastx_db_info <- blastx_file[blastx_db_line:(blastx_db_line +1)]
writeLines(blastx_db_info, "db_info_blastx.txt")


#creating query chunks
blastx_result_lines <- blastx_file[(blastx_db_line + 2):length(blastx_file)]
blastx_query_starts <- grep ("^Query=", blastx_result_lines)
blastx_query_ends <- grep ("^Effective search space used:", blastx_result_lines)
blastx_query_chunks <- mapply(function(start,end) {
    blastx_result_lines[start:end]
}, blastx_query_starts, blastx_query_ends, SIMPLIFY=FALSE)

#save chunks
saveRDS(blastx_query_chunks, "query_chunks_blastx.rds", compression = "xz")

#read in blastp0 file
blastp <- arg[3]
blastp_file <- readLines(blastp)

#obtain database information from blastx file #this is going to be the same as the blastp file because it's the same database
blastp_db_line <- grep("^Database:", blastp_file)
blastp_db_info <- blastp_file[blastp_db_line:(blastp_db_line +1)]
writeLines(blastp_db_info,"db_info_blastp.txt" )

#creating query chunks. Each chunk starts with Query=  and ends with effective space used..
blastp_result_lines <- blastp_file[(blastp_db_line + 2):length(blastp_file)]
blastp_query_starts <- grep ("^Query=", blastp_result_lines)
blastp_query_ends <- grep ("^Effective search space used:", blastp_result_lines)
blastp_query_chunks <- mapply(function(start,end) {
    blastp_result_lines[start:end]
}, blastp_query_starts, blastp_query_ends, SIMPLIFY=FALSE)

#save chunks
saveRDS(blastp_query_chunks, "query_chunks_blastp.rds", compression = "xz" )


blastn <- arg[4]

#read in blastn0 file (IF PRESENT)

if (!is.na(blastn)) {
blastn_file <- readLines(blastn)


#obtain database information from blastx file #this is going to be the same as the blastp file because it's the same database
blastn_db_line <- grep("^Database:", blastn_file)
blastn_db_info <- blastn_file[blastn_db_line:(blastn_db_line +1)]
writeLines(blastn_db_info, "db_info_blastn.txt")

#creating query chunks
blastn_result_lines <- blastn_file[(blastn_db_line + 2):length(blastn_file)]
blastn_query_starts <- grep ("^Query=", blastn_result_lines)
blastn_query_ends <- grep ("^Effective search space used:", blastn_result_lines)
blastn_query_chunks <- mapply(function(start,end) {
    blastn_result_lines[start:end]
}, blastn_query_starts, blastn_query_ends, SIMPLIFY=FALSE)



#save chunks
saveRDS(blastn_query_chunks, "query_chunks_blastn.rds", compression = "xz")
} else {
    print("No blastn file provided.")
}
