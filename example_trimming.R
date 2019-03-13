
library(ShortRead)

fq <- readFastq("fq_example.fastq")
bed  <- read.table("bed_example.txt", header = TRUE)
bed <- as(bed, "GRanges")

adapter <-  DNAString("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG")

# summary of sensibility, specificity, positive predictive value, negative predictive value, Matthews correlation 
summarize_extern <- function(out, bed, readlen) {
  
  TP <- sum((width(out) == width(bed)) & width(bed) < readlen ) / sum(width(bed) < readlen)
  TN <- sum((width(out) == width(bed)) & width(bed) == readlen) / sum(width(bed) == readlen)
  FP <- sum(width(out) < width(bed)) / length(bed)
  FN <- sum(width(out) > width(bed) & width(bed) < readlen ) / sum(width(bed) < readlen)
  
  outlist <- list(SEN = TP/(TP + FN), 
                  SPC = TN / (FP + TN),
                  PPV =  TP / (TP + FP), 
                  NVP = TN/(TN+FN), 
                  MCC = ((TP * TN) - (FP  * FN))/ sqrt((TP + FN) * (TP + FN) * (TN + FP) + (TN + FN)))
  
  outlist
  
}

  outlist <- list()
  for(i in seq(1, 50, 1)) {
    out <- trimLRPatterns(subject=fq, 
                          Rpattern = adapter,
                          with.Rindels = FALSE, 
                          max.Rmismatch = i)
    
    outlist[[i]] <-summarize_extern(out, bed, 150)
  }
  
  outlist <- do.call("rbind", outlist)
  
  outlist
  