## ---------------------------
##
## Script name: Main.R
##
## Purpose of script: 
##  Solutions to the 2nd exercise of the 1st Assignemt for the BNP54 class
##
## Author: Konstantinos Perperidis
##
## Date Created: 2024-03-24
##
## Email: k.a.perperidis@gmail.com
##        std531050@ac.eap.gr
##
## ---------------------------
library("stringi");
library("scales");

source("CalcGCContent.R");

# Part 1 - Read file contents as a character sequence into a variable.
dnaSequence <- paste(readLines("dna_sequence.txt"))[1];

# Part 2 - Find the GC-Content inside the DNA sequence.
print(sprintf("The total G-C Content for the given sequence equals to: %s", percent(calcGCContent(dnaSequence), accuracy = 0.2)), quote = FALSE);

# Find the GC-Content in batches of 500 and the amount of GC binucleotides within each batch.

# Generate a vector of vectors containing binucleotides extracted from the sequence.
binucleotides <- regmatches(dnaSequence, gregexpr(".{2}", dnaSequence))[[1]];

# Break the vector of 2.500 bi-nucleotides evenly into sub vectors w. 250 base pairs each.
chunkSize <- 250;
binucleotidesVectors <- split(binucleotides, ceiling(seq_along(binucleotides) / chunkSize));

# Consolidate all bi-nucleotides into a data frame
binucleotide.abs.dataFrame <- as.data.frame(rbind(
      table(binucleotidesVectors[[1]]), 
      table(binucleotidesVectors[[2]]),
      table(binucleotidesVectors[[3]]),
      table(binucleotidesVectors[[4]]),
      table(binucleotidesVectors[[5]]),
      table(binucleotidesVectors[[6]]),
      table(binucleotidesVectors[[7]]),
      table(binucleotidesVectors[[8]]),
      table(binucleotidesVectors[[9]]),
      table(binucleotidesVectors[[10]])));

# Iterate over all rows and print out the results for the G-C Content
for (index in 1:10) {
  currentSubSeq <- paste(binucleotidesVectors[[index]], collapse = '');
  currentGCContent <- calcGCContent(currentSubSeq);
  print(sprintf("============= Batch no. %s =============", index), quote = FALSE);
  print(sprintf("G-C Content: %s", percent(currentGCContent, accuracy = 0.2)), quote = FALSE);
  print(sprintf("GC-Binucleotide Amount: %s", binucleotide.abs.dataFrame$GC[index]), quote = FALSE);
}

