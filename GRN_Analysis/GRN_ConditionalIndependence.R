###############################################################
# Conditional independencies of a GRN - Using daggity network # 
###############################################################

# Load packages
library(dagitty)
library(lavaan)
library(sna)
library(DOT)
library(igraph)
library(pcalg)
library(SEMgraph)

###########################
### Simulation example 1 ##
###########################

g<- dagitty("dag {
    X -> R -> S -> Z <- U <- V -> Y
    Z -> P
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
predictors <- c("X","R","S","Z","P")

# Which of the five variables X,R,S,T,P are d-separated from Y given the same five variables
dseparated( g, "Y", list(), predictors)
intersect(predictors, dseparated( g, "Y", list(), predictors)) #Result: X,R,P

# Simulate data
d <- simulateSEM( g, .7, .7, N=10000 )
head(d)
fit.lm <- lm(Y ~ X + R + S + Z + P, data=d )
confint(fit.lm) # Result: X,R,P 95% CIs for coefficients contain zero. The conditional independence is verified.

fit.lm <- lm(Y ~ X +Z , data=d )
confint(fit.lm) # Result: X,R,P 95% CIs for coefficients contain zero. The conditional independence is verified.


###########################
### Simulation example 2 ##
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


##################
### E.coli GRN ###
##################

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
