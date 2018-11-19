## Lab 2
View(lizards_sequences)

## Question 1: DNA sequence acquisition and simulation
len_list <- lengths(lizards_sequences)

generate_nt <- function(uniform_limits){
  #print(uniform_limits)
  u <- runif(1,0,1)
  nk <- c("a","c","g","t")
  ind <- which(u<limits)[1] -1
  result <- nk[ind]
  return(result)
}

new_seq_list <- list()
new_seq <- c()
for (i in 1:length(lizards_sequences)){
  # extract the base composition of each sequence
  comp <- ape::base.freq(lizards_sequences[i])
  # calculate the limits for the uniform simulation method
  limits <- c(0,comp[1],comp[1]+comp[2],comp[1]+comp[2]+comp[3],1)
  limits_matrix <- matrix(rep.int(limits,len_list[i]), nrow = len_list[i],byrow = TRUE)
  # simulate the sequence of the same composition
  new_seq <- apply(limits_matrix,1,generate_nt)
  new_seq_list[i] <- list(new_seq)
}
# convert result to the DNAbin class
simulated_seq <- as.DNAbin(new_seq_list)

ape::base.freq(lizards_sequences[1])
ape::base.freq(simulated_seq[1])

ape::base.freq(lizards_sequences[29])
ape::base.freq(simulated_seq[29])
