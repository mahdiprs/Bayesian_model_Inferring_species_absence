
# Calculate the moments and PMF for discrete distributions from their PGF
# Support is for non-negative integers
# 

library(Deriv)

calc_moments<- function(pgf, pars) {
# calculate first two moments (i.e mean, variance)
 require(Deriv)
    d1<- Deriv(pgf, "s")
    d2<- Deriv(pgf, "s", nderiv=2)
    ex<- d1(pars, s=1)
    vx<- d2(pars,s=1) + ex - ex^2
    list(Ex = ex, Vx = vx)
}

calc_pmf<- function(pgf, pars, support) {
# Given PGF, calculate PMF for the support
  require(Deriv)
  probs<- rep(NA, length(support))
  probs[1]<- pgf(pars, 0)
  for(i in 2:length(support)) {
    dn<- Deriv(pgf, "s", nderiv=support[i])
    lprob<- log(dn(pars, 0)) - lfactorial(support[i])
    probs[i]<- exp(lprob)
  }
  data.frame(n=support,p=probs)
}

#---------------------------------------------
# Poisson
pgf <- function(pars,s) {exp(pars * (s - 1))}

calc_moments(pgf, 12)

probs<- calc_pmf(pgf, 12, 0:30) # PMF on 0:30

all.equal(probs$p, dpois(0:30,12)) # should be TRUE

win.graph(8,8)
plot(probs$n,probs$p, type="l")

#---------------------------------------------
# Binomial
pgf <- function (pars, s)  {(pars[1]*s + (1-pars[1]))^pars[2]}

calc_moments(pgf, c(0.25, 20))

probs<- calc_pmf(pgf, c(0.25,20), 0:20) # p = 0.25, N=20

all.equal(probs$p, dbinom(0:20,20,0.25))

win.graph(8,8)
plot(probs$n,probs$p, type="l")

#-----------------------------------------------
# Inference without population growth
# Equation 6

pgf<- function(pars, s) {
  n<- length(pars)
  p<- pars[1]
  delta<- pars[2:n]
  (1-p + p*prod(delta)*s)/(1-p + p*prod(delta))
}

# one survey
calc_moments(pgf, c(0.5, 0.2))  # p = 0.5, delta_bar = 0.2 

# two surveys
calc_moments(pgf, c(0.5, 0.2, 0.2))

# three surveys
calc_moments(pgf, c(0.5, 0.2, 0.2, 0.2))

# P(X=0|D)

pgf(c(0.5, 0.2), s=0) # one survey
pgf(c(0.5, 0.2, 0.2), s=0) # two surveys
pgf(c(0.5, 0.2, 0.2, 0.2), s=0) # three surveys

probs<- calc_pmf(pgf, c(0.5, 0.2, 0.2, 0.2), 0:1)

#-----------------------------------------------
# Deterministic population growth

expon<- function(lam, k){
# recursive function to calculate deterministic growth after k surveys
  if(k == 1) return(lam)
  else return(lam^k + expon(lam,k-1))
}

pgf<- function(pars, s) {
  # pgf for probability X > 0
  (1- pars[1] + pars[1]*(1-pars[2]*pars[3]/pars[4])^expon(pars[5], pars[6]) * s)/
    (1 - pars[1] + pars[1]*(1-pars[2]*pars[3]/pars[4])^expon(pars[5],pars[6]))
} 


parms<- c(0.5, 0.2, 8, 10, 2, 1) # p, gamma, n, N, lambda, kappa


pgf(parms, s=0) # P(X=0|Dbar)

# one survey
calc_moments(pgf, parms) 

# two surveys
parms<- c(0.5, 0.2, 8, 10, 2, 2) 
calc_moments(pgf, parms)

# three surveys
parms<- c(0.5, 0.2, 8, 10, 2, 3) 
calc_moments(pgf, parms)


probs<- calc_pmf(pgf, parms, 0:1)


#----------------------------------------
# Stochastic growth 

# Case of single survey (kappa=1)

pgf<- function(pars, s) {
# single -ve survey
    (1-pars[1] + pars[1]*exp(pars[5]*((1-pars[2]*pars[3]/pars[4]) * s - 1)))/
      (1-pars[1] + pars[1]*exp(-pars[5]*(pars[2]*pars[3]/pars[4])))
}

parms<- c(0.5, 0.2, 8, 10, 2) #p, gamma, n, N, lambda

pgf(parms, s=0)

calc_moments(pgf, parms) 

probs<- calc_pmf(pgf, parms, 0:10)

win.graph(8,8)
plot(probs$n,probs$p, type="b")

# Now try for kappa > 1
# Since the PGF is recursive in s Need to use non-standard evaluation 
# so that you can still get symbolic derivatives

# Helper functions

make_num<- function(kappa) {
  if(kappa==1){
    num<- paste0("exp(lambda*((1-gamma*n/N)*s-1))")
  }
  else {
    num<- paste0("exp(lambda*((1-gamma*n/N)*",make_num(kappa-1),"-1))")
  }
  num
}


make_denom<- function(kappa) {
  if(kappa==1){
    denom<- paste0("exp(-lambda*gamma*n/N)")
  }
  else {
    denom<- paste0("exp(lambda*((1-gamma*n/N)*",make_denom(kappa-1),"-1))")
  }
  denom
}


pgf<- function(p,lambda,gamma,n,N,kappa,s, what=c("eval","moments","pmf"), support=NULL) {
  what<- match.arg(what)
  num<- make_num(kappa)
  denom<- make_denom(kappa)
  e<- str2expression(text=paste0("(1-p+p*",num,")/(","1-p+p*",denom,")"))
  if(what == "eval") return(eval(e))
  else if(what == "moments") {
    d1<- Deriv(e, "s")
    d2<- Deriv(e, "s", nderiv = 2)
    ex<- eval(d1)
    vx<- eval(d2) + ex - ex^2
    return(list(Ex=ex,Vx=vx))
  }
  else if(what == "pmf") {
    if(is.null(support)) stop("option 'pmf' requires integer support")
    probs<- rep(NA, length(support))
    probs[1]<- eval(e)
    for(i in 2:length(support)) {
      cat("doing pmf for ",support[i],"\n")
      dn<- Deriv(e, "s", nderiv=support[i])
      lprob<- log(eval(dn)) - lfactorial(support[i])
      probs[i]<- exp(lprob)
    }
    return(data.frame(n=support,p=probs))
  }
}

# single survey (check with above)
pgf(p=0.5, lambda=2, gamma=0.2, n=8, N=10, kappa=1, s=0, what="eval") #P(X = 0|Dbar)

# single survey (check with above)
pgf(p=0.5, lambda=2, gamma=0.2, n=8, N=10, kappa=1, s=1, what="moments")

#two surveys
pgf(p=0.5, lambda=2, gamma=0.2, n=8, N=10, kappa=2, s=1, what="moments")

#three surveys
pgf(p=0.5, lambda=2, gamma=0.2, n=8, N=10, kappa=3, s=1, what="moments")


# now try PMF (**note this struggles when kappa >= 3)

pk1<- pgf(p=0.5, lambda=2, gamma=0.2, n=8, N=10, kappa=1, s=0, what="pmf", support=0:6)
pk2<- pgf(p=0.5, lambda=2, gamma=0.2, n=8, N=10, kappa=2, s=0, what="pmf", support=0:6)
pk3<- pgf(p=0.5, lambda=2, gamma=0.2, n=8, N=10, kappa=3, s=0, what="pmf", support=0:6)


win.graph(8,8)
plot(pk1$n,pk1$p,type="b", pch=0, ylim=c(0,1))
lines(pk2$n,pk2$p, type="b", lty=2, pch = 1)
lines(pk3$n,pk3$p, type="b", lty=3, pch = 2)
lines(pk3$n,pk4$p, type="b", lty=4, pch = 3)



#============================================
# Alternative implementation using symengine
#============================================

library(symengine)


calc_moments<- function(pgf, parms) {
  # calculate first two moments (i.e mean, variance)
  require(symengine)
  e<- S(pgf)
  parms$s<- 1
  nam<- names(parms)
  args<- list()
  for(i in 1: length(parms)) {
    assign(nam[i],parms[[i]])
    args[[i]]<- S(eval(nam[i]))
  }
  args<- Vector(args)
  d1<- symengine::D(e, "s")
  d2<- symengine::D(d1, "s")
  fd1<- as.function(d1, args=args)
  fd2<- as.function(d2, args=args)
    ex<- do.call(fd1, parms)
    vx<- do.call(fd2, parms) + ex - ex^2
    return(list(Ex=ex,Vx=vx))
}

calc_pmf<- function(pgf, parms, support) {
  # Given PGF, calculate PMF for the support
  require(symengine)
  probs<- rep(NA, length(support))
  parms$s<- 0
  nam<- names(parms)
  args<- list()
  for(i in 1: length(parms)) {
    assign(nam[i],parms[[i]])
    args[[i]]<- S(eval(nam[i]))
  }
  args<- Vector(args)
  dn<- list()
  e<- symengine::S(pgf)
  dn[[1]]<- e
  fe<- as.function(e, args=args)
  probs[1]<- do.call(fe, parms)
  for(i in 2:length(support)) {
    dn[[i]]<- symengine::D(dn[[i-1]], "s")
    fe<- as.function(dn[[i]], args=args)
    ev<- do.call(fe, parms)
    lprob<- log(ev) - lfactorial(support[i])
    probs[i]<- exp(lprob)
  }
  data.frame(n=support,p=probs)
}

#---------------------------------------------
# Poisson
pgf <- "exp(lambda * (s - 1))"

parms<- list(lambda=12)

calc_moments(pgf, parms)

probs<- calc_pmf(pgf, parms, 0:30) # PMF on 0:30

all.equal(probs$p, dpois(0:30,12)) # should be TRUE

win.graph(8,8)
plot(probs$n,probs$p, type="l")

#---------------------------------------------
# Binomial
pgf<- "(p*s + (1-p))^N"

parms<- list(p=0.3, N=20)

calc_moments(pgf, parms)

probs<- calc_pmf(pgf, parms, 0:20) # p = 0.25, N=20

all.equal(probs$p, dbinom(0:20,20,0.3))

win.graph(8,8)
plot(probs$n,probs$p, type="l")


#-----------------------------------------------
# Inference without population growth
# Equation 6


pgf<- "(1-p + p*delta*s)/(1-p + p*delta)"


# one survey
parms<- list(p=0.5, delta = 0.2) # p = 0.5, delta_bar = 0.2 
calc_moments(pgf, parms)  

# two surveys
parms<- list(p=0.5, delta = prod(0.2,0.2)) # p = 0.5, delta_bar = 0.2 
calc_moments(pgf, parms)

# three surveys
parms<- list(p=0.5, delta = prod(0.2,0.2,0.2)) # p = 0.5, delta_bar = 0.2 
calc_moments(pgf, parms)


parms<- list(p=0.5, delta = 0.2) # p = 0.5, delta_bar = 0.2 
probs<- calc_pmf(pgf, parms, 0:1)

#----------------------------------------
# Inference with stochastic growth 
# Since the PGF is recursive in s Need to use non-standard evaluation 
# so that you can still get symbolic derivatives

# Helper functions


make_num<- function(kappa) {
  if(kappa==1){
    num<- paste0("exp(lambda*((1-gamma*n/N)*s-1))")
  }
  else {
    num<- paste0("exp(lambda*((1-gamma*n/N)*",make_num(kappa-1),"-1))")
  }
  num
}


make_denom<- function(kappa) {
  if(kappa==1){
    denom<- paste0("exp(-lambda*gamma*n/N)")
  }
  else {
    denom<- paste0("exp(lambda*((1-gamma*n/N)*",make_denom(kappa-1),"-1))")
  }
  denom
}

make_pgf<- function(kappa) {
  num<- make_num(kappa)
  denom<- make_denom(kappa)
  pgf<- paste0("(1-p+p*",num,")/(","1-p+p*",denom,")")
}
  


# moments for 1,2 and 3 surveys (i.e. kappa = 1, 2 or 3)


parms<-list(p=0.5, lambda=2, gamma=0.2, n=8, N=10)

# Moments

pgf<- make_pgf(kappa=1) 
calc_moments(pgf, parms) #P(X = 0|Dbar)

pgf<- make_pgf(kappa=2) 
calc_moments(pgf, parms) #P(X = 0|Dbar)

pgf<- make_pgf(kappa=3) 
calc_moments(pgf, parms) #P(X = 0|Dbar)



# now try PMF 
# symengine allows increased support compared with Deriv

pgf<- make_pgf(kappa=1) 
pk1<- calc_pmf(pgf, parms, support=0:15)

pgf<- make_pgf(kappa=2) 
pk2<- calc_pmf(pgf, parms, support=0:15)

pgf<- make_pgf(kappa=3) 
pk3<- calc_pmf(pgf, parms, support=0:15)


win.graph(8,8)
plot(pk1$n,pk1$p,type="b", pch=0, ylim=c(0,1))
lines(pk2$n,pk2$p, type="b", lty=2, pch = 1)
lines(pk3$n,pk3$p, type="b", lty=3, pch = 2)
