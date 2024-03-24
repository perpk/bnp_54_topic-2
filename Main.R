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

# Part 1 - Read file contents as a character sequence into a variable.
dnaSequence <- paste(readLines("dna_sequence.txt"))[1];

# Part 2 - Find the GC-Content inside the DNA sequence.

## 2.1 Find the total amount of G-C pairs.
totalContentInSeq <- stringi::stri_count_fixed(dnaSequence, "GC");
print(sprintf("The total amount of G-C pairs in the sequence: %s", totalContentInSeq), quote=FALSE);

## 2.2 Find the GC-Content in batches of 500 (totals in 10 batches for the given 5000 base long sequnece).
# Split to bi-nucleotides and vectorize
binucleotides <- regmatches(dnaSequence, gregexpr(".{2}", dnaSequence))[[1]];

# Break the vector of 2.500 bi-nucleotides into equidistant sub vectors w. 250 base pairs each.
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

binucleotide.pct.dataFrame <- as.data.frame(rbind(
    prop.table(table(binucleotidesVectors[[1]])), 
    prop.table(table(binucleotidesVectors[[2]])),
    prop.table(table(binucleotidesVectors[[3]])),
    prop.table(table(binucleotidesVectors[[4]])),
    prop.table(table(binucleotidesVectors[[5]])),
    prop.table(table(binucleotidesVectors[[6]])),
    prop.table(table(binucleotidesVectors[[7]])),
    prop.table(table(binucleotidesVectors[[8]])),
    prop.table(table(binucleotidesVectors[[9]])),
    prop.table(table(binucleotidesVectors[[10]]))));

# Iterate over all rows and print out the results for the G-C Content
for (index in 1:10) {
  outContent <- sprintf("The G-C Content for batch-no. %s: %s", index, percent(binucleotide.pct.dataFrame$GC[index], accuracy=0.1));
  outAmount <- sprintf("The G-C Amount for batch-no. %s: %s", index, binucleotide.abs.dataFrame$GC[index]);
  print(outContent, quote=FALSE);
  print(outAmount, quote=FALSE);
  print("---------------------------------------------", quote=FALSE);
}

