library(ape)
## Gene bank accession numbers taken from http://www.jcsantosresearch.org/Class_2014_Spring_Comparative/pdf/week_2/Jan_13_15_2015_GenBank_part_2.pdf
lizards_accession_numbers <- c("JF806202", "HM161150", "FJ356743", "JF806205", 
                               "JQ073190", "GU457971", "FJ356741", "JF806207",
                               "JF806210", "AY662592", "AY662591", "FJ356748",       
                               "JN112660", "AY662594", "JN112661", "HQ876437", 
                               "HQ876434", "AY662590", "FJ356740", "JF806214", 
                               "JQ073188", "FJ356749", "JQ073189", "JF806216", 
                               "AY662598", "JN112653", "JF806204", "FJ356747", 
                               "FJ356744", "HQ876440", "JN112651", "JF806215",
                               "JF806209") 
lizards_sequences<-ape::read.GenBank(lizards_accession_numbers)
print(lizards_sequences)
ape::write.dna(lizards_sequences, file ="lizard_seqs.fasta", format = "fasta", colsep = "")


