## ---------------------------
##
## Script name: CalcGCContent.R
##
## Purpose of script:
##  Calculate the GC Content for a given sequence of bases.
##
## Author: Konstantinos Perperidis
##
## Date Created: 2024-03-25
##
## Email: k.a.perperidis@gmail.com
##        std531050@ac.eap.gr
##
## ---------------------------
library("stringi");

calcGCContent <- function(seq) {
  seq.trimmed <- stri_trim(seq);
  guanine.count <- stri_count_fixed(seq.trimmed, "G");
  cytosine.count <- stri_count_fixed(seq.trimmed, "C");
  total.gc.count <- guanine.count + cytosine.count;
  total.base.count <- stri_length(seq.trimmed);
  return (total.gc.count / total.base.count) * 100;
}