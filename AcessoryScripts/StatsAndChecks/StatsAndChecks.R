library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1] 
sample <- args[2] 

#how many complete ORFs pep
completeORF_dir <- paste0(dir, "Pipelines/Venomflow/ORFprediction/Combined/Complete/")
completeORFfasta <- list.files(path = completeORF_dir, pattern = "\\.combine\\.complete\\.cleaned\\.pep$", full.names = TRUE)
headers <- grep("^>", readLines(completeORFfasta), value = TRUE)
num_completesequences <- length(headers) #A

#how many secreted ORFs pep
secretedORFs_dir1 <- paste0(dir, "Pipelines/Venomflow/Secreted/Mature/Combined/")
secretedORF_file1 <- list.files(path = secretedORFs_dir1, pattern = "\\.fasta$", full.names = TRUE)
secretedORFs_dir2 <- paste0(dir, "Pipelines/Venomflow/Secreted/Mature/Signalp/")
secretedORF_file2 <- list.files(path = secretedORFs_dir2, pattern = "\\_mature\\.fasta$", full.names = TRUE)

if (length(secretedORF_file1) != 0) {
  secretedORFfile <- secretedORF_file1
} else {
  secretedORFfile <- secretedORF_file2
}
headers2 <- grep("^>", readLines(secretedORFfile), value = TRUE)
num_secretedsequences <- length(headers2) #B

#how many uniquesequences transdf 
transdfdir <- paste0(dir, "Pipelines/Analysis/IntermediateFiles/Dataframes/ORFs/")
transdffile <- list.files(path = transdfdir, pattern = "transdf\\.csv$", full.names = TRUE)
transdf <- read.csv(transdffile, header = TRUE)
transdfunique <- transdf %>%
  distinct(Transdecoder_ID, .keep_all = TRUE)
number_sequences_transdffile <- nrow(transdfunique) #C should equal A 1
filteredout <- transdfunique %>%
  filter(SP == "SP" & TMHMM == "TRUE")
number_filteredout <- nrow(filteredout) #D
#how mnay rows transdf_distinct
transdfdistinctfile <- list.files(path = transdfdir, pattern = "distinct\\.csv$", full.names = TRUE)
transdfdistinct <- read.csv(transdfdistinctfile, header = TRUE)
transdfdistinct <- transdfdistinct %>%
  distinct(Transdecoder_ID, .keep_all = TRUE)
number_sequences_transdfdistinct <- nrow(transdfdistinct) #E E+D = B 2

#unfiltered how many rows 
unfiltereddir <- paste0(dir, "Pipelines/Analysis/Overview/Dataframes/Unannotated/")
unfilteredfile <- list.files(path = unfiltereddir, pattern = "unfiltered\\.csv$", full.names = TRUE)
unfiltered <- read.csv(unfilteredfile, header = TRUE)
number_sequences_unfiltered <- nrow(unfiltered) #F equal E 3

laxfile <- list.files(path = unfiltereddir, pattern = "lax\\.csv$", full.names = TRUE)
lax <- read.csv(laxfile, header = TRUE)
number_sequences_totallax <- nrow(lax) #G
number_sequences_laxonly <- nrow(lax[lax$Filter == "Lax",]) #H
number_sequences_strictonly <- nrow(lax[lax$Filter == "Strict",]) #I H+I = G 4

strictfile <- list.files(path = unfiltereddir, pattern = "strict\\.csv$", full.names = TRUE)
strict <- read.csv(strictfile, header = TRUE)
number_sequences_totalstrict <- nrow(strict) #j EQUAL I 5

#annoatated all how many rows 
annotateddir <-paste0(dir, "Pipelines/Analysis/Overview/Dataframes/Annotated/")
anntoatedselect <- list.files(path = annotateddir, pattern = "Select_Annotated_df\\.csv$", full.names = TRUE)
anntoatedselect <- read.csv(anntoatedselect, header = TRUE)
putativetoxinstotal <-  nrow(anntoatedselect) #K = L+M +N 6

number_sequences_cat1 <- nrow(anntoatedselect[anntoatedselect$Category == "Category 1",]) #L
number_sequences_cat2 <- nrow(anntoatedselect[anntoatedselect$Category == "Category 2",]) #M
number_sequences_cat3 <- nrow(anntoatedselect[anntoatedselect$Category == "Category 3",]) #N



anntoatedall <-list.files(path = annotateddir, pattern = "all_Annotated_df\\.csv$", full.names = TRUE)
anntoatedall <- read.csv(anntoatedall, header = TRUE)
annotataedtotal <-  nrow(anntoatedall) #o should equal G 7

number_sequences_cat1_2 <- nrow(anntoatedall[anntoatedall$Category == "Category 1",]) #P = L 8
number_sequences_cat2_2 <- nrow(anntoatedall[anntoatedall$Category == "Category 2",]) #Q =M 9
number_sequences_cat3_2 <- nrow(anntoatedall[anntoatedall$Category == "Category 3",]) #R =N 10 
number_sequences_notox <- nrow(anntoatedall[anntoatedall$Category == "Putative non-toxins",]) #S
number_sequences_lowgenome <- nrow(anntoatedall[anntoatedall$Category == "Lower genome support transcripts",]) #T
#p q r s t = O 11

#html table check 
#completeorfs



#Checks 
Check1 <- number_sequences_transdffile == num_completesequences
Check2 <- number_sequences_transdfdistinct + number_filteredout == num_secretedsequences
Check3 <- number_sequences_unfiltered == number_sequences_transdfdistinct
Check4 <-number_sequences_strictonly+ number_sequences_laxonly == number_sequences_totallax
Check5 <-number_sequences_totalstrict == number_sequences_strictonly
Check6 <-number_sequences_cat1 + number_sequences_cat2 + number_sequences_cat3 == putativetoxinstotal
Check7 <- annotataedtotal ==number_sequences_totallax
Check8 <- number_sequences_cat1_2 == number_sequences_cat1
Check9 <- number_sequences_cat2_2 == number_sequences_cat2
Check10 <- number_sequences_cat3_2 == number_sequences_cat3
Check11 <- number_sequences_cat1_2+number_sequences_cat2_2+number_sequences_cat3_2 + number_sequences_notox + number_sequences_lowgenome ==annotataedtotal

All_Checks_True <- Check1 & Check2 & Check3 & Check4 & Check5 & 
  Check6 & Check7 & Check8 & Check9 & Check10 & Check11

CheckPass <- ifelse(All_Checks_True, "PASS", "FAIL")
CheckPass[is.na(CheckPass)] <- "FAIL"


final_df <- data.frame(
  sample,
  num_completesequences,
  number_sequences_transdfdistinct,
  number_sequences_totallax,
  number_sequences_strictonly,
  number_sequences_notox,
  number_sequences_lowgenome,
  number_sequences_cat1,
  number_sequences_cat2,
  number_sequences_cat3,
  putativetoxinstotal,
  CheckPass,
  stringsAsFactors = FALSE
)

colnames(final_df) <- c("Sample", "Total ORFs", "Total Secreted", "Lax Filter", "Strict Filter", 
                        "Putative Non-Toxins", "LGSS", "Cat 1", "Cat 2", "Cat 3", 
                        "Annotated Toxins Total", "CheckPass")

write.csv(final_df, paste0 (dir,sample, "_Overview_stats.csv"), row.names = FALSE)
