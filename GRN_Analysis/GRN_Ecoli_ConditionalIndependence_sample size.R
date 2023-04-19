#######################################################
# Conditional independencies of a GRN vs sample size  # 
#######################################################

# Load packages
library(dagitty)
library(lavaan)
library(sna)
library(DOT)
library(igraph)
library(pcalg)
library(SEMgraph)
library(reticulate)
library(dplyr)
library(ggdag)
library(rgbif)
library(readr)

#############
# Functions #
#############

# Function to create a agrapgh from interaction file 

path <- function (GRN){
  relationships <- list()# open the file
  paths <- list()
  # split the line into child and parents
  for (line in 1:nrow(GRN)) {
    line_elements <- strsplit(GRN[line,], ",") [[1]]
    #child <- line_elements[1]
    child <- as.numeric(line_elements[1])
    child <- GeneIDs[child+1]
    Num.regs <- as.numeric(line_elements[2])
    parents <- strsplit(line_elements[3:(2+Num.regs)]," ")
    parents <-lapply(parents,as.numeric)
    parents <- unlist(parents)
    parents <- GeneIDs[parents +1]
    for (parent in parents) {
      if (!(parent %in% names(relationships))) {
        relationships[[parent]] <- list()
      }
      relationships[[parent]] <- c(relationships[[parent]], child)
    }
  }
  #relationships # each child with their parents
  
  # create the paths
  paths <- list()  
  for (parent in names(relationships)) {
    for (child in relationships[[parent]]) {
      path <- child
      current <- child
      while (current %in% names(relationships)) {
        current <- relationships[[current]][[1]]
        path <- paste0(path, "->" ,current)
      }
      path <- paste0(parent,"->",path)
      paths <- c(paths, path)
    }
  }
  return (unlist(paths))
} 


# Function to extract simulated data for each cell type 
#get.data.celltype <- function(df, bins, celltype){
#  if (celltype == bins){
#    temp <- df %>%
#      slice(which(row_number() %% bins == 0))
#  }
#  else {
#    temp <- df %>%
#      slice(which(row_number() %% bins == celltype))
#  }
#  return(temp)
#}

# Function to extract simulated data for each cell type 
get.data.celltype <- function(df, bins, celltype){
  n.cells <- nrow(df)/bins
  temp <- split( df , cut(1:nrow(df), bins) )[[celltype]]
  if(nrow(temp)==n.cells){
    return(temp)
  }
  else{
    return("number of simulated cells doesn't match with the extracted data")
  }
}


# Function for conditional independence test (compute independencies with at most 3 conditioning variables)
con.ind <- function(df, bins, celltype, g){
  d <- get.data.celltype(df, bins, celltype)
  d <- data.frame(scale(d))
  imp <- Filter(function(x) length(x$Z)<4, impliedConditionalIndependencies(g))
  CI <- localTests( g, d, "cis", R=100, tests=imp, loess.pars=list(span=0.6) )
  return (CI)
}

###########################################################
### Load E.coli GRN graph, gene IDs, and simulated data ###
###########################################################

# Load dot file and GRN interaction file for the GRN grapgh
dot <- read.dot("C:/Projects/GitHub/Sergio_NM/GNW_sampled_GRNs/Ecoli_100_net1.dot") # GRN graph
GRN <- read.table("C:/Projects/GitHub/Sergio_NM/data_sets/De-noised_100G_9T_300cPerT_4_DS1/Interaction_cID_4.txt", quote="\"", comment.char="")

# Load simulated gene expressions
# GRN_out <- read.csv("~/Missing data/SERGIO/Demo/GeneExpression_SingleCell.csv") # Simulated single cell data
GRN_out <- read.csv("C:/Projects/GitHub/Sergio_NM/data_sets/De-noised_100G_9T_300cPerT_4_DS1/simulated_noNoise_0.csv") # Existing data

# Load gene IDs
#GeneIDs <- read.delim("~/Missing data/GRN_Analysis/GeneIDs.txt", header=FALSE) # Gene IDs
#GeneIDs <- gsub(';','',GeneIDs[,2])
GeneIDs <- read.delim("C:/Projects/GitHub/Sergio_NM/GNW_sampled_GRNs/Ecoli_100_net1.dot", header=FALSE)
GeneIDs <-GeneIDs$V2[-c(1,2)][1:100]
GeneIDs <- gsub(';','',GeneIDs)
gene.geneID <- data.frame(cbind(Gene=GeneIDs, ID=0:99))
gene.geneID$ID <- as.numeric(gene.geneID$ID)
gene.geneID

# Master regulators
gene.geneID[gene.geneID$ID %in% c(67,17,56,44,62,93,84),]

# Load data from Jeremy's DS.daggity file
grn.dagitty <- read.delim("C:/Projects/GitHub/Sergio_NM/GRN_Analysis/DS1.dagitty", header=FALSE)


############################# 
# Summary of simulated data #  
#############################

summary(GRN_out)
GRN_out <- t(GRN_out[,-1])
colnames(GRN_out) <- GeneIDs
dim(GRN_out)
head(GRN_out)

# Pearson's correlation plot (Linear correlation)
res <- cor(GRN_out)
corrplot::corrplot(res,type = "upper", tl.col = "black", tl.cex = 0.5)

#################
# The GRN graph #  
#################

# Load from .dot file (This turned out to be inconsistent with the data in SERGIO)
#g.graph <- graph.adjacency(dot)
#plot(g.graph)
#g <- graph2dagitty(g.graph)
#g <- graphLayout(g)
#plot(g)

# Create the GRN graph from interaction file
unlisted <- path(GRN)
unlisted[duplicated(unlisted)] # check duplicated paths
write.table(unlisted,"paths.txt",col.names = FALSE, row.names = FALSE)
path_strings <- paste(unlisted,collapse=' ')
g <- dagitty(paste0("dag {",path_strings,"}"))
g <- dagitty(paste0('dag { hybC [pos="0,1"] rutR [pos="1,1"]',path_strings,'}'))
plot(g)

# Create the GRN graph from Jeremy's daggity file
path_strings <- paste(grn.dagitty$V1,collapse=' ')
g <- dagitty(path_strings )
plot(g)

# Conditional independencies
# For each pair of non-adjacent nodes in this graph, the set of variables that d-separates that pair.
# i.e., for each non-adjacent pair of variables, all minimal sets that we can condition on to render that pair independent.
true.ind <- impliedConditionalIndependencies(g, type = "missing.edge", max.results = Inf)
true.ind

# Draw paths from icd to fadI
#paths(g, "icd", "fadI" )$paths

# Draw directed paths from icd to fadI
#paths(g, "potI", "nemA", directed=TRUE)$paths


##########################
# Non-linear correlation #
##########################

df <- data.frame(GRN_out) # 100 genes with 300 cells with first 9 cell types

##############################################################################
# For cell type 1 (compute independencies with at most 3 conditioning variables)
set.seed(1000)
CI <- con.ind(df, 9, 1, g) # Conditional independence test using data from cell type 1
plotLocalTestResults(head(CI,10))
plotLocalTestResults(head(CI[order(CI$estimate),],20))
ind <- CI[(CI$`2.5%` <0 & CI$`97.5%`>0),]
dep <- CI[!(CI$`2.5%` <0 & CI$`97.5%`>0),]

# Plot the conditional independencies that we were not able to prove to be dependent (we cannot conclude them to be independent: do not reject null hypothesis that beta=0)
plotLocalTestResults(ind)

# Plot the conditional independencies that we prove to be dependent
plotLocalTestResults(dep)

# When sample size increases the probability of detecting dependencies increases. 
p <- nrow(dep)/(nrow(dep)+nrow(ind))
p

# Investigate on dependence from the local test
df.celltype1 <- data.frame(scale(get.data.celltype(df,9,1))) # Cell type 1
plot(df.celltype1$argI,df.celltype1$mgtA)
fit.lm <- lm(argI ~ mgtA +betI+ddpC , data=df.celltype1 ) # Look at linear relationship
confint(fit.lm) # Result: argI vs mgtA: 95% CIs for coefficients do not contain zero. The conditional independence is not verified.

##############################################################################
# For cell type 2 (compute independencies with at most 3 conditioning variables)

# Compute independencies with at most 3 conditioning variables
set.seed(1000)
CI <- con.ind(df, 9, 2, g) # Conditional independence test using data from cell type 2
plotLocalTestResults(head(CI,10))
plotLocalTestResults(head(CI[order(CI$estimate),],100))
ind <- CI[(CI$`2.5%` <0 & CI$`97.5%`>0),]
dep <- CI[!(CI$`2.5%` <0 & CI$`97.5%`>0),]

# Plot the conditional independencies that we were not able to prove to be dependent (we cannot conclude them to be independent: do not reject null hypothesis that beta=0)
plotLocalTestResults (ind)

# Plot the conditional independencies that we prove to be dependent
plotLocalTestResults (dep)

# When sample size increases the probability of detecting dependencies increases. 
p <- nrow(dep)/(nrow(dep)+nrow(ind))
p

##############################################################################
# For cell type 3 (compute independencies with at most 3 conditioning variables)

# Compute independencies with at most 3 conditioning variables
set.seed(1000)
CI <- con.ind(df, 9, 3, g) # Conditional independence test using data from cell type 3
plotLocalTestResults (head(CI,50))
plotLocalTestResults(head(CI[order(CI$estimate),],20))
ind <- CI[(CI$`2.5%` <0 & CI$`97.5%`>0),]
dep <- CI[!(CI$`2.5%` <0 & CI$`97.5%`>0),]

# Plot the conditional independencies that we were not able to prove to be dependent (we cannot conclude them to be independent: do not reject null hypothesis that beta=0)
plotLocalTestResults(ind)

# Plot the conditional independencies that we prove to be dependent
plotLocalTestResults(dep)

# When sample size increases the probability of detecting dependencies increases. 
p <- nrow(dep)/(nrow(dep)+nrow(ind))
p

##############################################################################
# For cell type 4 (compute independencies with at most 3 conditioning variables)

# Compute independencies with at most 3 conditioning variables
set.seed(1000)
CI <- con.ind(df, 9, 4, g) # Conditional independence test using data from cell type 3
plotLocalTestResults (head(CI,50))
plotLocalTestResults(head(CI[order(CI$estimate),],20))
ind <- CI[(CI$`2.5%` <0 & CI$`97.5%`>0),]
dep <- CI[!(CI$`2.5%` <0 & CI$`97.5%`>0),]

# Plot the conditional independencies that we were not able to prove to be dependent (we cannot conclude them to be independent: do not reject null hypothesis that beta=0)
plotLocalTestResults(ind)

# Plot the conditional independencies that we prove to be dependent
plotLocalTestResults(dep)

# When sample size increases the probability of detecting dependencies increases. 
p <- nrow(dep)/(nrow(dep)+nrow(ind))
p

##############################################################################
# For cell type 5 (compute independencies with at most 3 conditioning variables)

# Compute independencies with at most 3 conditioning variables
set.seed(1000)
CI <- con.ind(df, 9, 5, g) # Conditional independence test using data from cell type 3
plotLocalTestResults (head(CI,50))
plotLocalTestResults(head(CI[order(CI$estimate),],20))
ind <- CI[(CI$`2.5%` <0 & CI$`97.5%`>0),]
dep <- CI[!(CI$`2.5%` <0 & CI$`97.5%`>0),]

# Plot the conditional independencies that we were not able to prove to be dependent (we cannot conclude them to be independent: do not reject null hypothesis that beta=0)
plotLocalTestResults(ind)

# Plot the conditional independencies that we prove to be dependent
plotLocalTestResults(dep)

# When sample size increases the probability of detecting dependencies increases. 
p <- nrow(dep)/(nrow(dep)+nrow(ind))
p

##############################################################################
# For cell type 6 (compute independencies with at most 3 conditioning variables)

# Compute independencies with at most 3 conditioning variables
set.seed(1000)
CI <- con.ind(df, 9, 6, g) # Conditional independence test using data from cell type 3
plotLocalTestResults (head(CI,50))
plotLocalTestResults(head(CI[order(CI$estimate),],20))
ind <- CI[(CI$`2.5%` <0 & CI$`97.5%`>0),]
dep <- CI[!(CI$`2.5%` <0 & CI$`97.5%`>0),]

# Plot the conditional independencies that we were not able to prove to be dependent (we cannot conclude them to be independent: do not reject null hypothesis that beta=0)
plotLocalTestResults(ind)

# Plot the conditional independencies that we prove to be dependent
plotLocalTestResults(dep)

# When sample size increases the probability of detecting dependencies increases. 
p <- nrow(dep)/(nrow(dep)+nrow(ind))
p

##############################################################################
# For cell type 7 (compute independencies with at most 3 conditioning variables)

# Compute independencies with at most 3 conditioning variables
set.seed(1000)
CI <- con.ind(df, 9, 7, g) # Conditional independence test using data from cell type 3
plotLocalTestResults (head(CI,50))
plotLocalTestResults(head(CI[order(CI$estimate),],20))
ind <- CI[(CI$`2.5%` <0 & CI$`97.5%`>0),]
dep <- CI[!(CI$`2.5%` <0 & CI$`97.5%`>0),]

# Plot the conditional independencies that we were not able to prove to be dependent (we cannot conclude them to be independent: do not reject null hypothesis that beta=0)
plotLocalTestResults(ind)

# Plot the conditional independencies that we prove to be dependent
plotLocalTestResults(dep)

# When sample size increases the probability of detecting dependencies increases. 
p <- nrow(dep)/(nrow(dep)+nrow(ind))
p

##############################################################################
# For cell type 8 (compute independencies with at most 3 conditioning variables)

# Compute independencies with at most 3 conditioning variables
set.seed(1000)
CI <- con.ind(df, 9, 8, g) # Conditional independence test using data from cell type 3
plotLocalTestResults (head(CI,50))
plotLocalTestResults(head(CI[order(CI$estimate),],20))
ind <- CI[(CI$`2.5%` <0 & CI$`97.5%`>0),]
dep <- CI[!(CI$`2.5%` <0 & CI$`97.5%`>0),]

# Plot the conditional independencies that we were not able to prove to be dependent (we cannot conclude them to be independent: do not reject null hypothesis that beta=0)
plotLocalTestResults(ind)

# Plot the conditional independencies that we prove to be dependent
plotLocalTestResults(dep)

# When sample size increases the probability of detecting dependencies increases. 
p <- nrow(dep)/(nrow(dep)+nrow(ind))
p

##############################################################################
# For cell type 9 (compute independencies with at most 3 conditioning variables)

# Compute independencies with at most 3 conditioning variables
set.seed(1000)
CI <- con.ind(df, 9, 9, g) # Conditional independence test using data from cell type 3
plotLocalTestResults (head(CI,50))
plotLocalTestResults(head(CI[order(CI$estimate),],20))
ind <- CI[(CI$`2.5%` <0 & CI$`97.5%`>0),]
dep <- CI[!(CI$`2.5%` <0 & CI$`97.5%`>0),]

# Plot the conditional independencies that we were not able to prove to be dependent (we cannot conclude them to be independent: do not reject null hypothesis that beta=0)
plotLocalTestResults(ind)

# Plot the conditional independencies that we prove to be dependent
plotLocalTestResults(dep)

# When sample size increases the probability of detecting dependencies increases. 
p <- nrow(dep)/(nrow(dep)+nrow(ind))
p
