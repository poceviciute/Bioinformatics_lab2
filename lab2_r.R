## Lab 2
setwd("C:/Users/Milda/Documents/Documents/Master degree/Bioinformatics/Bioinformatics_lab2")
#install.packages("ape")
library(ape)
lizards_sequences<-ape::read.FASTA("lizard_seqs.fasta")
#View(lizards_sequences)
#class(lizards_sequences)

# Question 1: DNA sequence acquisition and simulation

## Q1.1

len_list <- lengths(lizards_sequences)

new_seq_list <- list()
new_seq <- c()
for (i in 1:length(lizards_sequences)){
    # extract the base composition of each sequence
    freqx <- ape::base.freq(lizards_sequences[i])
    new_seq <- sample(c("a","c","g","t"),len_list[i],rep=TRUE,prob=freqx)
    new_seq_list[i] <- list(new_seq)
}

# convert result to the DNAbin class
simulated_seq1 <- as.DNAbin(new_seq_list)
names(simulated_seq1) <- names(lizards_sequences)

# Print out the base compositions
for (i in 1:length(lizards_sequences)){
print(ape::base.freq(lizards_sequences[i]))
print(ape::base.freq(simulated_seq1[i]))
}


## Q1.2

new_seq_list <- list()
new_seq <- c()
for (i in 1:length(lizards_sequences)){
    # extract the base composition of each sequence
    freqx <- ape::base.freq(lizards_sequences[i])
    new_seq <- sample(c("a","c","g","t"),len_list[i],rep=TRUE,prob=freqx)
    new_seq_list[i] <- list(new_seq)
}

# convert result to the DNAbin class
simulated_seq2 <- as.DNAbin(new_seq_list)
names(simulated_seq2) <- paste0(names(lizards_sequences),"_2")
lengths(simulated_seq2)
lengths(lizards_sequences)

# simulate a tree with 33 nodes

tree <- rtree(n = 33, rooted = TRUE, tip.label = names(lizards_sequences) )
plot(tree, edge.width = 2)

# Simulate sequence from the tree comperable to the lizard sequences

library(phangorn)
#install.packages("phangorn")
data <- c()
simulated_seq_tree <- list()

# the transition matrix is based on the longest sequence
    index <- which.max(len_list)
    #DNA_list <- lizards_sequences[index]$AY662592
    DNA_list <- unlist(as.character(lizards_sequences),use.names = FALSE)
    DNA_unique <- c("a","c","g","t")
    matrix <-  matrix(0,
                      ncol = length(DNA_unique),
                      nrow = length(DNA_unique))
    
    for (i in 1:(length(DNA_list) - 1)) {
        index_of_i <- DNA_unique == DNA_list[i]
        index_of_i_plus_1 <- DNA_unique == DNA_list[i + 1]
        matrix[index_of_i, index_of_i_plus_1] = matrix[index_of_i, index_of_i_plus_1] + 1
    }
    
    Q <- matrix / rowSums(matrix)
    Q_lower <- Q[lower.tri(Q)]



    size <- round(mean(len_list))
    bf <-  ape::base.freq(lizards_sequences)
    print(bf)
    # generate the data from the tree

    simulated_seq_tree <-  phangorn::simSeq(
            tree,
            l = size,
            type = "DNA",
            bf = bf,
            Q = Q_lower)
    
simulated_seq_tree <- as.DNAbin(simulated_seq_tree)
ape::base.freq(lizards_sequences)
ape::base.freq(simulated_seq_tree)
ape::write.dna(simulated_seq_tree, file ="lizard_sim_seqs.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)

# Question 2: Sequence analysis

## Q2.1

# Print out the individual base compositions
for (i in 1:length(lizards_sequences)){
    print(ape::base.freq(lizards_sequences[i]))
    print(ape::base.freq(simulated_seq1[i]))
    print(ape::base.freq(simulated_seq_tree[i]))
}
# CG content
comp_true <- ape::base.freq(lizards_sequences)
cg_true <- sum(comp_true[2:3])/sum(comp_true)
comp_s1 <- ape::base.freq(simulated_seq1)
cg_s1 <- sum(comp_s1[2:3])/sum(comp_s1)
comp_s2 <- ape::base.freq(simulated_seq_tree)
cg_s2 <- sum(comp_s2[2:3])/sum(comp_s2)

# A,C,G,T content

print(ape::base.freq(lizards_sequences))
print(ape::base.freq(simulated_seq1))
print(ape::base.freq(simulated_seq_tree))

# First locate the first start codon, and from there look for a stop codon

lizard <- as.character(lizards_sequences)
sim1 <- as.character(simulated_seq1)
sim2 <- as.character(simulated_seq2)

find_start <- function(seq){
    start_cod <- c("a","t","c")
    limit <- length(seq)
    # checking from the 1st possition
    j = 1
    result1 <- c(999999)
    while (j < limit-2){ 
        acid <- seq[j:(j+2)]
        #print(acid)
        if (all(acid == start_cod)){
            result1[1] <- j
            break
        }
        j= j+3
    }  
    # checking from the 2nd possition
    j = 2
    result2 <- c(999999)
    while (j < limit-2){ 
        acid <- seq[j:(j+2)]
        #print(acid)
        if (all(acid == start_cod)){
            result2[1] <- j
            break
        }
        j= j+3
    } 
    # checking from the 3rd possition
    j = 3
    result3 <- c(999999)
    while (j < limit-2){ 
        acid <- seq[j:(j+2)]
        #print(acid)
        if (all(acid == start_cod)){
            result3[1] <- j
            break
        }
        j= j+3
    }
    js <- c(1,2,3)
    results <- c(result1,result2,result3)
    indx <- which.min(results)
    final_result <- list(pos = results[indx], start = js[indx])
    return(final_result)
}



find_stops <- function(seq){
    #print(seq)
    start_pos <- find_start(seq)$pos
    if (start_pos == 999999){start_pos = c(1,2,3)}
    j = start_pos
    limit  <- length(seq)
    stop_co1 <-c("t","a","a")
    stop_co2 <-c("t","a","g")
    stop_co3 <-c("t","g","a")
    if (length(j)==1){
    i = 1
    result <- c()
    while (j < limit-2){ 
        acid <- seq[j:(j+2)]
        #print(acid)
        if (all(acid == stop_co1) || all(acid ==stop_co2) || all(acid == stop_co3)){
            result[i] <- j
            i <- i+1
        }
        j= j+3
    }
    }else{
        #print(j)
        result_final <- list()
        for (k in 1:length(j)){
            result2 <- c()
            i = 1
            jj = k
        while (jj < limit-2){ 
            acid <- seq[jj:(jj+2)]
            #print(acid)
            if (all(acid == stop_co1) || all(acid ==stop_co2) || all(acid == stop_co3)){
                result2[i] <- jj
                i <- i+1
            }
            jj= jj+3
        }
         #   print(result2)
            result_final[k] <- list(result2)
        }
     result <- result_final[which.max(lengths(result_final))]
    }
    return( list(stop_loc= result))
}

find_start(lizard[2]$JF806202)
find_stops(lizard[2]$JF806202)
lizard_stops <- lapply(lizard,find_stops)
sim1_stops <- lapply(sim1,find_stops)
sim2_stops <- lapply(sim2,find_stops)

no_stop <- function(set) {
    no_stop_codon <- c()
    k = 1
    for (i in 1:33) {
        if (is.null(set[[i]]$stop_loc)) {
            no_stop_codon[k] = i
            k = k + 1
        }
    }
    return(no_stop_codon)
}

is.null(lizard_stops[[32]]$stop_loc)
no_stop(lizard_stops)
no_stop(sim1_stops)
no_stop(sim2_stops)

## Q2.2
#install.packages("markovchain")
library(markovchain)
#mcX <- markovchainFit(lizard[2])
#mcX 
dim(mcX$estimate)
which(lizard[1] == "y")

flat_lizard <- unlist(lizard,use.names = FALSE)
mc_lizard <-  markovchainFit(flat_lizard)

# the original data contains some "weird nucleotides" - what to do with them?
which(flat_lizard == "y")
which(flat_lizard == "m")
which(flat_lizard == "r")
which(flat_lizard == "s")

flat_sim1 <- unlist(sim1,use.names = FALSE)
mc_sim1 <-  markovchainFit(flat_sim1)

flat_sim2 <- unlist(sim2,use.names = FALSE)
mc_sim2 <-  markovchainFit(flat_sim2)

mc_lizard$estimate
mc_sim2$estimate
Q

## Q2.3
library(ape)
ape::write.dna(flat_lizard, file ="lizard_flat.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)

ape::write.dna(flat_sim1, file ="sim1_flat.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)

ape::write.dna(flat_sim2, file ="sim2_flat.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)

#ape::dist.dna(as.DNAbin(lizard[1]),as.DNAbin(sim1[1]),model="F81")
#ape::clustal()
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("msa", version = "3.8")
# library(msa)
filepath <- system.file("lizard_flat.fasta", "sim1_flat.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)
msa(mySeqs)
