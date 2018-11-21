## Lab 2

str(lizards_sequences)
label <- attributes(lizards_sequences)
desc <- names(lizards_sequences)
print(lizards_sequences)

### Question 1
## Part 1.1

artificialDNA_simulator <- function(length_seq,composition)
{
  nucleotides <- c("a","c","g","t")
  seq <- sample(nucleotides,length_seq,rep=TRUE,prob = composition)
  
}

seq_len <- as.numeric(lengths(lizards_sequences))
list_seq <- list()

for (i in 1:length(lizards_sequences))
{
  compos <- base.freq(lizards_sequences[i])
  len_seq <- seq_len[i]
  sim_seq <- artificialDNA_simulator(len_seq,compos)
  list_seq[i] <- list(sim_seq)
  
}

simulated_seq <- as.DNAbin(list_seq)
names(simulated_seq) <- desc

base.freq(lizards_sequences[1])
base.freq(simulated_seq[1])

base.freq(lizards_sequences[33])
base.freq(simulated_seq[33])


## 1.2

