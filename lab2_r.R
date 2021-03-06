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

ape::write.dna(simulated_seq1, file ="sim1_seq.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)

## Q1.2

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
#ape::write.dna(simulated_seq_tree, file ="lizard_sim_seqs.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)
simulated_seq_tree <- ape::read.FASTA("lizard_sim_seqs.fasta")
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
sim2 <- as.character(simulated_seq_tree)

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
#dim(mcX$estimate)
#which(lizard[1] == "y")

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
# https://www.biostars.org/p/179953/
library(msa)
library(seqinr)
library(lattice)
#install.packages("seqinr")

aligment_lizard <- msa("lizard_seqs.fasta", type="dna")
aligment_sim1 <- msa("sim1_seq.fasta", type="dna")
aligment_sim2 <- msa("lizard_sim_seqs.fasta", type="dna")

aligment_lizard2 <- msaConvert(aligment_lizard, type="seqinr::alignment")
dist_lizard <- seqinr::dist.alignment(aligment_lizard2,matrix="identity")

aligment_sim1_2 <- msaConvert(aligment_sim1, type="seqinr::alignment")
dist_sim1 <- seqinr::dist.alignment(aligment_sim1_2,matrix="identity")

aligment_sim2_2 <- msaConvert(aligment_sim2, type="seqinr::alignment")
dist_sim2 <- seqinr::dist.alignment(aligment_sim2_2,matrix="identity")

# the seqinr::dist.alignment computes the square root distance
dm_lizard <- as.matrix(dist_lizard)
dm_sim1 <- as.matrix(dist_sim1)
dm_sim2 <- as.matrix(dist_sim2)

# heatmaps
new.palette=colorRampPalette(c("black","red","yellow","white"),space="rgb")
levelplot(dm_lizard[1:ncol(dm_lizard),ncol(dm_lizard):1],col.regions=new.palette(20))
levelplot(dm_sim1[1:ncol(dm_sim1),ncol(dm_sim1):1],col.regions=new.palette(20))
levelplot(dm_sim2[1:ncol(dm_sim2),ncol(dm_sim2):1],col.regions=new.palette(20))

# Question 3: Phylogeny reconstruction
## Q 3.1
nj_tree <- function(x){
    tr <- nj(x)
    tr$tip.label <- names(lizards_sequences)
    return(tr)
} 


hemoTree_lizard <- nj_tree(dist_lizard)



hemoTree_sim1 <- nj_tree(dm_sim1)


hemoTree_sim2 <- nj_tree(dm_sim2)

# Perform a phylogenetic bootstrap analysis


boot_lizard <-  boot.phylo(hemoTree_lizard, dm_lizard, FUN=nj_tree, quiet = TRUE)
boot_s1     <- boot.phylo(hemoTree_sim1, dm_sim1, FUN=nj_tree, quiet = TRUE)
boot_s2     <- boot.phylo(hemoTree_sim2, dm_sim2, FUN=nj_tree, quiet = TRUE)

#it is the number of bootstraped trees that have the same breaks as your main input tree


## Q3.2

plot(hemoTree_lizard, main="Phylogenetic Tree of Lizard")
plot(hemoTree_sim1, main="Phylogenetic Tree of Lizard DNA Sim1")
plot(hemoTree_sim2, main="Phylogenetic Tree of Lizard DNA Sim2")

comp1 <- ape::comparePhylo(hemoTree_lizard,hemoTree_sim1)
comp2 <- ape::comparePhylo(hemoTree_lizard,hemoTree_sim2)
comp12 <- ape::comparePhylo(hemoTree_sim1,hemoTree_sim2)
# An ultrametric tree - a tree where all the path-lengths from the root to the tips are equal








