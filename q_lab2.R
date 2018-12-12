## Lab 2
#rm(list = ls())
str(lizards_sequences)
library(ape)
library(seqinr)

lizard_seq_seqinr_format <- read.fasta(file = "lizard_seqs.fasta", seqtype = "DNA",
                                       as.string = TRUE, forceDNAtolower = FALSE)

### Question 1
## Part 1.1
lizards_sequences<-read.FASTA("lizard_seqs.fasta")
desc <- names(lizards_sequences)
print(lizards_sequences)

seq_len <- as.numeric(lengths(lizards_sequences))
list_seq <- list()


artificialDNA_simulator <- function()
{
  nucleotides <- c("a","c","g","t")
  
  for (i in 1:length(lizards_sequences))
  {
    composition <- base.freq(lizards_sequences[i])
    sim_seq <- sample(nucleotides,seq_len[i],rep=TRUE,prob = composition)
    list_seq[i] <- list(sim_seq)
  }
  return(list_seq)
}

simulated_seq <- artificialDNA_simulator()
simulated_seq <- as.DNAbin(simulated_seq)
names(simulated_seq) <- desc

write.dna(simulated_seq, file ="simulated_seqs.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)


## 1.2

library(TreeSim)
## simulating artificial seq2

## simulating phylogenetic tree using TreeSim
phyltree<- TreeSim::sim.bd.taxa(n=33, numbsim=1, lambda=1, mu=0)[[1]]
phyltree$tip.label <- desc

##plotting tree
plot(phyltree,type = "phylogram",show.tip.label = TRUE)

## simulating from tree
library(phangorn)
seq_len <- as.numeric(lengths(lizards_sequences))

#transition_matrix <- matrix(0,nrow = 4,ncol = 4)

#equillibrium_freq <- rep(0.25,4)

index <- which.max(seq_len)
#DNA_list <- lizards_sequences[index]$AY662592
DNA_list <- unlist(as.character(lizards_sequences),use.names = FALSE)
DNA_unique <- c("a","c","g","t")
rate_matrix <-  matrix(0,
                  ncol = length(DNA_unique),
                  nrow = length(DNA_unique))

for (i in 1:(length(DNA_list) - 1)) {
  index_of_i <- DNA_unique == DNA_list[i]
  index_of_i_plus_1 <- DNA_unique == DNA_list[i + 1]
  rate_matrix[index_of_i, index_of_i_plus_1] = rate_matrix[index_of_i, index_of_i_plus_1] + 1
}

Q <- rate_matrix/rowSums(rate_matrix)

simulated_tree <- simSeq(phyltree, l = round(mean(seq_len)), Q = Q, type = "DNA", bf = base.freq(lizards_sequences))


simulated_tree2 <- lapply(simulated_tree,function(seq){
  #seq = a
  states = c("a","c","g","t") 
  seq[which(seq==1)] <- states[1]
  seq[which(seq==2)] <- states[2]
  seq[which(seq==3)] <- states[3]
  seq[which(seq==4)] <- states[4]
  return(seq)
})

simulated_seq2 <- as.DNAbin(simulated_tree2)

# Write the simulated sequence to a fasta file
#write.dna(simulated_seq2, file ="lizard_sim_seqs.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)


cat("Base composition of the data from GenBank:")
cat("\n")
ape::base.freq(lizards_sequences)
cat("\n")
cat("Base composition of the simulated sequences Q1.1:")
cat("\n")
ape::base.freq(simulated_seq)
cat("\n")
cat("Base composition of the simulated sequences Q1.2:")
cat("\n")
ape::base.freq(simulated_seq2)


### using TreeSim and complete Q the result is almost same


#### 2.1

for (i in 1:length(lizards_sequences)){
  cat(paste0("Individual base composition of sequence ",i))
  cat("\n")
  cat("The data from the GenBank:")
  cat("\n")
  print(base.freq(lizards_sequences[i]))
  cat("The simulated data in Q1.1:")
  cat("\n")
  print(base.freq(simulated_seq[i]))
  cat("The simulated data in Q1.2:")
  cat("\n")
  print(base.freq(simulated_seq2[i]))
  cat("\n")
}


# GC content

comp_true <- base.freq(lizards_sequences)
cg_true <- sum(comp_true[2:3])/sum(comp_true)
comp_s1 <- base.freq(simulated_seq)
cg_s1 <- sum(comp_s1[2:3])/sum(comp_s1)
comp_s2 <- ape::base.freq(simulated_seq2)
cg_s2 <- sum(comp_s2[2:3])/sum(comp_s2)

cat("In the data from the GenBank:")
cat("\n")
print(cg_true)
cat("In the simulated data in Q1.1:")
cat("\n")
print(cg_s1)
cat("In the simulated data in Q1.2:")
cat("\n")
print(cg_s2)
cat("\n")

# GC content of sequnce
# gc_content = data.frame()
# lizard_seq_seqinr_format <- read.fasta(file = "lizard_seqs.fasta", seqtype = "DNA",
#                                        as.string = TRUE, forceDNAtolower = FALSE)
# for (i in 1:length(lizard_seq_data))
# {
#   #string manuplution
#   original_seq = unlist(strsplit(lizard_seq_seqinr_format[[i]], split=""))
#   original_seq = original_seq[original_seq != " "]
#   
#   
#  gc_content = rbind(gc_content,c(
#                     GC(x),
#                       GC(unlist(as.character(simulated_seq_1[[i]]),use.names = FALSE)),
#                       GC(unlist(as.character(simulated_seq_2[[i]]),use.names = FALSE)))
#                     )
#                     
# }
# colnames(gc_content) = c("original","sim1","sim2")



#AT Content

comp_true <- base.freq(lizards_sequences)
cg_true <- sum(comp_true[c(1,4)])/sum(comp_true)
comp_s1 <- base.freq(simulated_seq)
cg_s1 <- sum(comp_s1[c(1,4)])/sum(comp_s1)
comp_s2 <- ape::base.freq(simulated_seq2)
cg_s2 <- sum(comp_s2[c(1,4)])/sum(comp_s2)

cat("In the data from the GenBank:")
cat("\n")
print(cg_true)
cat("In the simulated data in Q1.1:")
cat("\n")
print(cg_s1)
cat("In the simulated data in Q1.2:")
cat("\n")
print(cg_s2)
cat("\n")


#smiluated sequence stop codon search

amino_lizard_seq <- read.fasta("lizard_amino_seq.fasta")
amino_sim1_seq <- read.fasta("sim_aminio_seqs_1.fasta")
amino_sim2_seq <- read.fasta("sim_aminio_seqs_2.fasta")

stop_code_count <- function(amino_acid_seq){
  
  stop_count = c()
  for (i in 1:length(amino_acid_seq)) {
    
    sequence = amino_acid_seq[[i]]
    stop <- sequence[which(sequence == "*")]
    stop_count[i] = length(stop)
  }
  return(stop_count)
}


data.frame(original_seq= stop_code_count(amino_lizard_seq),
           sim1 = stop_code_count(amino_sim1_seq),
           sim2 = stop_code_count(amino_sim2_seq)
)


### For some sequences the stop codons more in true in sequence than in simulated sequences but there are some sequences for which
# the stop codon does not appear at all in the original sequence but it appears many times in the simulated sequences.
# It could be possible because the nucleotides forming the stop codon was simulated more in simulated seq.

### 2.2
library(markovchain)

lizard <- as.character(lizards_sequences)
sim1 <- as.character(simulated_seq)
sim2 <- as.character(simulated_seq2)
flat_lizard <- unlist(lizard,use.names = FALSE)
mc_lizard <-  markovchainFit(flat_lizard,method = "bootstrap", nboot = 3, name = "Bootstrap Mc")
# the original data contains some "weird nucleotides" - what to do with them?
which(flat_lizard == "y")
which(flat_lizard == "m")
which(flat_lizard == "r")
which(flat_lizard == "s")
flat_sim1 <- unlist(sim1,use.names = FALSE)
mc_sim1 <-  markovchainFit(flat_sim1,method = "bootstrap", nboot = 3, name = "Bootstrap Mc")
flat_sim2 <- unlist(sim2,use.names = FALSE)
mc_sim2 <-  markovchainFit(flat_sim2,method = "bootstrap", nboot = 3, name = "Bootstrap Mc")

# comparison of the estimated transition matrices
#mc_lizard$estimate
#flat_sim1$estimate
#mc_sim2$estimate
#Q

##We see that the original data set from GenBank contains a small amount of additional "nucleotides": they are called "y", "m", and "r".
# These are the nucleotides which we are not sure about where Y indicates unsurity in A or G while R shows ambiguity in C or T and M appears when there is no surity at all.

###2.3

#In this part the three DNA collections of the sequences are aligned using the _msa_ package. Furthermore the distances between the alignments are calculated. They are used to produce the heatmaps
library(msa)
library(lattice)

aligment_lizard <- msa("lizard_seqs.fasta", type="dna")
aligment_sim1 <- msa("sim1_seq.fasta", type="dna")
aligment_sim2 <- msa("lizard_sim_seqs.fasta", type="dna")
aligment_lizard2 <- msaConvert(aligment_lizard, type="seqinr::alignment")
dist_lizard <- dist.alignment(aligment_lizard2,matrix="identity")
aligment_sim1_2 <- msaConvert(aligment_sim1, type="seqinr::alignment")
dist_sim1 <- dist.alignment(aligment_sim1_2,matrix="identity")
aligment_sim2_2 <- msaConvert(aligment_sim2, type="seqinr::alignment")
dist_sim2 <- seqinr::dist.alignment(aligment_sim2_2,matrix="identity")

# the seqinr::dist.alignment computes the square root distance
dm_lizard <- as.matrix(dist_lizard)
dm_sim1 <- as.matrix(dist_sim1)
dm_sim2 <- as.matrix(dist_sim2)

new.palette=colorRampPalette(c("black","red","yellow","white"),space="rgb")
levelplot(dm_lizard[1:ncol(dm_lizard),ncol(dm_lizard):1],col.regions=new.palette(20))
levelplot(dm_sim1[1:ncol(dm_sim1),ncol(dm_sim1):1],col.regions=new.palette(20))
levelplot(dm_sim2[1:ncol(dm_sim2),ncol(dm_sim2):1],col.regions=new.palette(20))


##3.1

upgma_tree <- function(x){
  tree <- upgma(x)
  tree$tip.label <- names(lizards_sequences)
  return(tree)
} 

Tree_lizard <- upgma_tree(dist_lizard)
Tree_sim1 <- upgma_tree(dm_sim1)
Tree_sim2 <- upgma_tree(dm_sim2)

plot(Tree_lizard, main = "Expected Tree")
plot(Tree_sim1, main = "Simulation 1 Tree")
plot(Tree_sim2, main= "Simulation 2 Tree")

# Perform a phylogenetic bootstrap analysis
boot_lizard <-  boot.phylo(Tree_lizard, dm_lizard, FUN=upgma_tree, trees = TRUE)$trees
lizard_clade <- prop.clades(Tree_lizard,boot_lizard, rooted = TRUE)

layout(1)
par(mar = rep(2, 4))
plot(Tree_lizard, main = "Individual Clade Support Values")

nodelabels(lizard_clade)
legend("bottomright", legend = c("Clades"), pch = 22,
       pt.bg = c("lightblue"), pt.cex = 2.5)


boot_sim1 <- boot.phylo(Tree_sim1, dm_sim1, FUN=upgma_tree, trees = TRUE)$trees
sim1_clade <- prop.clades(Tree_sim1,boot_sim1, rooted = TRUE)

plot(Tree_lizard, main = "Individual Clade Support Values - Sim1")

nodelabels(sim1_clade)
legend("bottomright", legend = c("Clades"), pch = 22,
       pt.bg = c("lightblue"), pt.cex = 2.5)


boot_sim2 <- boot.phylo(Tree_sim2, dm_sim2, FUN=upgma_tree, trees = TRUE)$trees
sim2_clade <- prop.clades(Tree_sim2,boot_sim2, rooted = TRUE)

plot(Tree_lizard, main = "Individual Clade Support Values - Sim2")

nodelabels(sim2_clade)
legend("bottomright", legend = c("Clades"), pch = 22,
       pt.bg = c("lightblue"), pt.cex = 2.5)


### 3.2

path_diff_sim1_lizard <- path.dist(Tree_sim1,Tree_lizard)
path_diff_sim2_lizard <- path.dist(Tree_sim2,Tree_lizard)
path_diff_sim1_sim2 <- path.dist(Tree_sim1,Tree_sim2)


felenstein_diff_sim1_lizard <- KF.dist(Tree_sim1,Tree_lizard)
felenstein_diff_sim2_lizard <- KF.dist(Tree_sim2,Tree_lizard)
felenstein_diff_sim1_sim2 <- KF.dist(Tree_sim1,Tree_sim2)


comp1 <- ape::comparePhylo(Tree_lizard,Tree_sim1)
comp2 <- ape::comparePhylo(Tree_lizard,Tree_sim2)
comp12 <- ape::comparePhylo(Tree_sim1,Tree_sim2)
