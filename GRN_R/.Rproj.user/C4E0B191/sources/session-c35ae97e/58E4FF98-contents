
# From Jeremy
text_to_paths <- function(file_path){
  # create an empty list to store paths
  paths <- list()# create a dictionary to store the parent-child relationships
  relationships <- list()# open the file
  file_content <- readLines(file_path)
  # iterate over the lines in the file
  for (line in file_content) {
    # split the line into child and parents
    line_elements <- strsplit(line, ":")[[1]]
    child <- line_elements[1]
    parents <- strsplit(line_elements[2], " ")[[1]]
    # store the parent-child relationships in the dictionary
    for (parent in parents) {
      if (!(parent %in% names(relationships))) {
        relationships[[parent]] <- list()
      }
      relationships[[parent]] <- c(relationships[[parent]], child)
    }
  }
  # create the paths
  for (parent in names(relationships)) {
    for (child in relationships[[parent]]) {
      path <- child
      current <- child
      while (current %in% names(relationships)) {
        current <- relationships[[current]][[1]]
        path <- paste0(current,path)
        paths <- c(paths, path)
      }
    }
  }
}

file_path = "graph.txt"
paths <- text_to_paths(file_path)

###############
### My code ###
###############

# Load data
#GRN <- read.table("~/Missing data/SERGIO/Demo/differentiation_input_GRN.txt", quote="\"", comment.char="")
GRN <- read.table("~/Missing data/SERGIO/Demo/single_cell_GRN.txt", quote="\"", comment.char="")
GRN_out <- read.csv("~/Missing data/SERGIO/Demo/GeneExpression.csv")

#GRN <- data.frame(rbind(c("bafB,1,bifA"),  c("bifA,1,exbD"), c("expD,1,fur"), c("bupC,2,bafB,bupB"),c("bupB,1,fofH"),c("fofH,1,fur")))

###############
### Input #####
###############

relationships <- list()# open the file
paths <- list()
# split the line into child and parents
for (line in 1:nrow(GRN)) {
  line_elements <- strsplit(GRN[line,], ",") [[1]]
  #child <- as.character(as.numeric(line_elements[1]))
  child <- line_elements[1]
  child <- paste0("Gene",child)
  Num.regs <- as.numeric(line_elements[2])
  parents <- strsplit(line_elements[3:(2+Num.regs)]," ")
  #parents <-lapply(parents,as.numeric)
  parents <-lapply(parents,as.character)
  for (parent in parents) {
    if (!(parent %in% names(relationships))) {
      relationships[[parent]] <- list()
    }
    relationships[[parent]] <- c(relationships[[parent]], child)
    #relationships[[parent]] <-lapply(relationships[[parent]],as.numeric)
  }
}
relationships

# create the paths
paths <- list()  
for (parent in names(relationships)) {
  for (child in relationships[[parent]]) {
    path <- child
    current <- child
    while (current %in% names(relationships)) {
      print(current[[1]])
      current <- relationships[[current]][[1]]
      print(current[[1]])
      path <- paste0(current,"-" ,path)
      print("itr")
    }
    path <- paste0(path,"-",parent)
    print("itrchild")
    paths <- c(paths, path)
  }
}
paths

# Eliminate duplicate paths
eliminate <- c()
for (i in 1:length(paths)){
  if(sum(grepl(paths[i], paths[-i], fixed = TRUE))>0 ){
    eliminate <- c(eliminate, i)
  }
}
final_paths <- paths[-eliminate]
final_paths # from children.... to parents...


###############
### Output ####
###############

summary(GRN_out)
dim(GRN_out)
res <- cor(t(GRN_out[,-1]))
corrplot::corrplot(res,type = "upper", tl.col = "black", tl.cex = 0.5)


######################################################
# Daggity network to obtain conditional independence # 
######################################################

library(dagitty)
library(lavaan)
library(sna)
library(DOT)
library(Rgraphviz)
library(igraph)
library(pcalg)
library(SEMgraph)

dot <- read.dot("~/Missing data/SERGIO/GNW_sampled_GRNs/Ecoli_100_net1.dot")
g.graph <- graph.adjacency(dot)
plot(g.graph)
g <- graph2dagitty(g.graph)
g <- graphLayout(g)
plot(g)

# Conditional independencies
# For each pair of non-adjacent nodes in this graph, the set of variables that d-separates that pair.
# i.e., for each non-adjacent pair of variables, all minimal sets that we can condition on to render that pair independent.
impliedConditionalIndependencies(g, type = "missing.edge", max.results = Inf)

# Draw paths from icd to fadI
paths(g, "icd", "fadI" )$paths

# Draw directed paths from icd to fadI
paths(g, "potI", "nemA", directed=TRUE)$paths


#dot("digraph {A -> B;}")
#dot(dotfile)
#g <- dagitty( "dag{ x -> m -> y }" )
#impliedConditionalIndependencies( g ) # one
#latents( g ) <- c("m")
#impliedConditionalIndependencies( g ) # none


###########################
### A simulation example ##
###########################

N <- 10000 # sample size
Ux <- rnorm( N )
Uy <- rnorm( N )
Uz <- rnorm( N )
X <- Ux
Y <- 1/3*X + Uy
Z <- 1/16*Y + Uz
d <- data.frame(X=X,Y=Y,Z=Z)

g <- dagitty("dag {
    Ux -> X -> Y -> Z <- Uz
    Uy -> Y
}")
coordinates(g) <- list(
  x=c(Ux=1,Uy=2,Uz=3,X=1,Y=2,Z=3),
  y=c(Ux=1,Uy=1,Uz=1,X=0,Y=0,Z=0) )
plot(g)

# For each pair of non-adjacent nodes in this graph, the set of variables that d-separates that pair.
# i.e., for each non-adjacent pair of variables, all minimal sets that we can condition on to render that pair independent.
impliedConditionalIndependencies(g, type = "missing.edge", max.results = Inf)

# For each pair of non-adjacent nodes in the graph, determine whether they are independent conditional on all other variables.
pairs <- combn( names(g), 2 )
apply( pairs, 2, function(x){
  all.other.variables <- setdiff( names(g), x )
  if( dseparated(g, x[1], x[2], all.other.variables ) ){
    message( x[1]," and ",x[2]," are independent given ", 
             paste( all.other.variables, collapse=",") )
  }
} )



###########################
### A simulation example ##
###########################

g<- dagitty("dag {
    X -> R -> S -> T <- U <- V -> Y
    T -> P
}")
g<-graphLayout(g)
plot(g)

# For each pair of non-adjacent nodes in this graph, the set of variables that d-separates that pair.
# i.e., for each non-adjacent pair of variables, all minimal sets that we can condition on to render that pair independent.
impliedConditionalIndependencies(g, type = "missing.edge", max.results = Inf)

# Independence of X and Y  given {R,V}
paths(g, "X", "Y", c("R","V") )

# all pairs of variables in the graph that are independent conditional on the set Z={R,V}
pairs <- combn( names(g), 2 )
apply( pairs, 2, function(x){
  p <- paths(g, x[1], x[2], c("R","V") )
  if( !p$open ){
    message( x[1]," and ",x[2]," are independent given {R,V}" )
  } else {
    message( x[1]," and ",x[2]," are possibly dependent given {R,V}" )
  }
} )

# For each pair of non-adjacent nodes in the graph, determine whether they are independent conditional on all other variables.
apply( pairs, 2, function(x){
  all.other.variables <- setdiff( names(g), x )
  if( dseparated(g, x[1], x[2], all.other.variables ) ){
    message( x[1]," and ",x[2]," are independent given ", 
             paste( all.other.variables, collapse=",") )
  }
} )

# Fitting models with fixed effects (no latent variables)
# Y=a+bX+cR+dS+eT+fP
predictors <- c("X","R","S","T","P")

# which of the five variables X,R,S,T,P are d-separated from Y given the same five variables
dseparated( g, "Y", list(), predictors)
intersect(predictors, dseparated( g, "Y", list(), predictors)) #Result: X,R,P

# Simulate data
d <- simulateSEM( g, .7, .7, N=10000 )
head(d)
fit.lm <- lm( Y ~ X + R + S + T + P, data=d )
confint(fit.lm) #Result: X,R,P coeficients are zero. The conditional independence is verified
