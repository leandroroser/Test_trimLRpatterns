library(ShortRead)
library(girafe)

fq <- readFastq("fq_example.fastq")
bed  <- read.table("bed_example.txt", header = TRUE)
bed <- as(bed, "GRanges")

adapter <-  DNAString("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG")

# summary of sensitivity, specificity, positive predictive value, negative predictive value, Matthews correlation 
summarize_extern <- function(out, bed, readlen) {
  
  TP <- sum((width(out) == width(bed)) & width(bed) < readlen ) / sum(width(bed) < readlen)
  TN <- sum((width(out) == width(bed)) & width(bed) == readlen) / sum(width(bed) == readlen)
  FP <- sum(width(out) < width(bed)) / length(bed)
  FN <- sum(width(out) > width(bed) & width(bed) < readlen ) / sum(width(bed) < readlen)
  
  outlist <- list(SEN = TP /(TP + FN), 
                  SPC = TN / (FP + TN),
                  PPV =  TP / (TP + FP), 
                  NVP = TN / (TN + FN), 
                  MCC = ((TP * TN) - (FP  * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  
  outlist
  
}

outlist_indels <- outlist_not_indels <- outlist_girafe <- list()
for(i in seq(1, 100, 1)) {
  out_indels <- trimLRPatterns(subject = fq, 
                               Rpattern = adapter,
                               with.Rindels = TRUE,
                               max.Rmismatch = i/100)
  
  out_not_indels <- trimLRPatterns(subject = fq, 
                                   Rpattern = adapter,
                                   with.Rindels = FALSE,
                                   max.Rmismatch = i/100)
  
  outlist_indels[[i]] <- summarize_extern(out_indels, bed, 150)
  outlist_not_indels[[i]] <- summarize_extern(out_not_indels, bed, 150)
  outlist_girafe[[i]] <- summarize_extern(out_girafe, bed, 150)
}

out_girafe <- trimAdapter(fq, adapter)
outlist_girafe <- summarize_extern(out_girafe, bed, 150)

outlist_indels <- do.call("rbind", outlist_indels)
outlist_not_indels <- do.call("rbind", outlist_not_indels)

outlist_indels
outlist_not_indels
outlist_girafe

