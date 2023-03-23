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


###################################################
# Implement a function for conditional dependence #
###################################################

con.ind <- function (Y, X, R){
  fit.lm <- lm(Y ~ X +R)
  CI <- as.data.frame(confint(fit.lm)) 
  CI$Ind <- with(CI, ifelse(CI$`2.5 %`<0 & CI$`97.5 %`<0 , "Negative conditional dependence",
                            ifelse(CI$`2.5 %`>0 & CI$`97.5 %`>0 , "Positive conditional dependence",
                                   "Not enough evidence for conditional dependence")))
  return(CI["X",])
}
con.ind(d$Y, d$X, d$R)



con.ind <- function (y_seq, x_seq, df){
  dtemp.all <- t(combn(y_seq, 2))
  for(i in 1:nrow(dtemp.all)){
    dtemp <- dtemp.all[i,]
    fit.lm <- lm(df[,dtemp[1]] ~ df[,dtemp[2]] + df[,x_seq])
    CI <- as.data.frame(confint(fit.lm)) 
    CI$Ind <- with(CI, ifelse(CI$`2.5 %`<0 & CI$`97.5 %`<0 , paste0("Negative conditional dependence given"," ",x_seq),
                              ifelse(CI$`2.5 %`>0 & CI$`97.5 %`>0 , paste0("Positive conditional dependence given"," ", x_seq),
                                     paste0("Not enough evidence for conditional dependence given", " ", x_seq))))
    
    print(dtemp)
    #print(CI)
    print(CI["df[, dtemp[2]]","Ind"])
    
  }
  #return(CI["d[, dtemp[2]]","Ind"])
}

y_seq <- c("X","Y", "P", "R","S") 
x_seq <- "Z" #given
con.ind(y_seq, x_seq, df=d)

fit.lm <- lm(Y ~ X , data=d )
confint(fit.lm) 

pairs <- combn( names(g), 2 )
apply( pairs, 2, function(x){
  p <- paths(g, x[1], x[2], "Z" )
  if( !p$open ){
    message( x[1]," and ",x[2]," are independent given Z" )
  } else {
    message( x[1]," and ",x[2]," are possibly dependent given Z" )
  }
} )

################################################################################
# For presentation                                                             #
# (Identification of causal queries using observational and intervention data) #
################################################################################

library(truncnorm)

set.seed(100)
x <- rtruncnorm(10,1,10)
z <- 2*x + rnorm(10,0,1)
s <- rnorm(10,0,1)
beta <- c(1,-1,2)
y <- beta[1] + beta[2]*x + beta[3]*z +s

df <- data.frame(nemA=y, nemR=x, rutA=z)
df <- sapply(df, round, 3)
View(df)

fit.lm1 <- lm(df$nemA ~ df$nemR) 
fit.lm2 <- lm(df$nemA ~ df$nemR + df$rutA) 
confint(fit.lm1)
confint(fit.lm2)
fit.lm1 
fit.lm2











