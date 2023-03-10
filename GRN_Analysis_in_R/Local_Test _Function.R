# Simulate full mediation model with measurement error of M1
set.seed(123)
d <- simulateSEM("dag{X->{U1 M2}->Y U1->M1}",.6,.6)

# Postulate and test full mediation model without measurement error
plotLocalTestResults(localTests( "dag{ X -> {M1 M2} -> Y }", d, "cis" ))

# Simulate data from example SEM
g <- getExample("Polzer")
d <- simulateSEM(g,.1,.1)

# Compute independencies with at most 3 conditioning variables
imp <- Filter(function(x) length(x$Z)<4, impliedConditionalIndependencies(g))
plotLocalTestResults(localTests( g, d, "cis.loess", R=100, tests=imp, loess.pars=list(span=0.6) ))