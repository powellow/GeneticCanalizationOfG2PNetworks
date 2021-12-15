# defining functions

#### calculating coefficients  ----
### (i.e. amount produced relative to amount that's genetically possible)

fun.proportion <- function(genetic.values,genetic.effects){
  maximum.value <- sum(2*genetic.effects)
  proportion <- genetic.values/maximum.value
  return(proportion)
}

# attention! it's 1 - relative value
fun.inv.proportion <- function(genetic.values,genetic.effects){
  maximum.value <- sum(2*genetic.effects)
  inv.proportion <- 1 - genetic.values/maximum.value
  return(inv.proportion)
}

#### Functions that describe the relationships ----

fun.calc.edge <- function(focal.node,interacting.node,relationship,coef,sign){
  ifelse(sign=="+",linear.output <- focal.node + coef*interacting.node,linear.output <- focal.node - coef*interacting.node)
  if(relationship=="linear"){
    return(linear.output)
  }
  if(relationship=="logistic"){
    logistic.output <- exp(linear.output)/(1 + exp(linear.output))
    return(logistic.output)
  }
}

#### Propensity to Branch ----
#### Final Ratio of Network Components

fun.g.val.ratio  = function(g.val.suc, g.val.sl){

  (g.val.suc/g.val.sl) + aux.cyt.coeff * (g.val.suc/g.val.sl)

}


#### Physiological Relationships for Network Nodes ----

fun.delta.cyt <- function(a,b,c,d,k,auxin,sucrose,cytokinin,Nutrients="No"){
  aux.repressed.cyto <- (c/(1 + (b*auxin)))
  sucr.stimulated.cyto <- (a*((sucrose^2)/(k + (sucrose^2))))
  cyto.degredatrion.cyto <- d
  delta.cyt <- (aux.repressed.cyto + sucr.stimulated.cyto)/cyto.degredatrion.cyto
  if(Nutrients=="Nitrogen"|Nutrients=="Nitrogen&Phosphorus"){
    nitrogen.effect <- 10
    delta.cyt <- delta.cyt*nitrogen.effect
  }
  return(delta.cyt)
}

fun.delta.strigo <- function(a,c,d,k,auxin,strigolactone,Nutrients="No"){
  base.synthesis.rate <- c
  aux.stimulated.strigo <- (a*((auxin^2)/(k + (auxin^2))))
  strigo.degredatrion.strigo <- d
  delta.strigo <- (base.synthesis.rate + aux.stimulated.strigo)/strigo.degredatrion.strigo
  if(Nutrients=="Phosphorus"){
    phosphorus.effect <- 50
    delta.strigo <- delta.strigo*(1/phosphorus.effect)
  }
  if(Nutrients=="Nitrogen&Phosphorus"){
    phosphorus.effect <- 50
    nitrogen.effect <- 10
    delta.strigo <- delta.strigo*(1/phosphorus.effect)*nitrogen.effect
  }
  return(delta.strigo)
}

fun.delta.suc <- function(a,sucrose,cytokinin){
  cyto.stimulated.suc <- a*cytokinin
  #cyto.stimulated.suc <- 0
  delta.suc <- sucrose + cyto.stimulated.suc
  return(delta.suc)
}

fun.delta.signal <- function(a.I.3,a.I.4,u.I.1,u.I.2,c.I,d.I,k.I,strigolactone,sucrose,cytokinin,signal){
  base.synthesis.rate <- c.I
  sucrose.striga.coeff <- (u.I.1 + (u.I.2*(sucrose^2)))
  strigo.response <- (a.I.3*((strigolactone^2)/(1 + (sucrose.striga.coeff*(strigolactone^2))))) #Sucrose might have to be on absolute scale...
  cyto.response <- (a.I.4*(1/(1 + (k.I*(cytokinin^2)))))
  signal.degredation.signal <- d.I
  delta.signal <- (base.synthesis.rate + strigo.response + cyto.response)/signal.degredation.signal
  return(delta.signal)
}

fun.time.bud.threshold <- function(intercept,sensitivity,threshold,signal){
  time <- ifelse(signal<threshold,intercept + (sensitivity*signal),8.3)
  time[time<0] <- 0.1
  return(as.numeric(time))
}

fun.time.bud <- function(intercept,sensitivity,threshold,signal){
  time <- intercept + (sensitivity*signal)
  return(as.numeric(time))
}

fun.bud.growth <- function(signal,threshold){
  bud.growth <- ifelse(signal<threshold,"yes","no")
  return(bud.growth)
}

fun.bud.lambda <- function(signal,min.lambda,max.lambda){
  min.signal <- min(signal[!is.na(signal)])
  max.signal <- max(signal[!is.na(signal)])
  lambdas <- min.lambda + (((signal - min.signal)*(max.lambda-min.lambda))/(max.signal-min.signal))
  return(lambdas)
}

fun.num.tiller <- function(nInd,tiller.signal){
  num.tillers <- rpois(nInd,tiller.signal) #This is adding stochasticity
  #We need to remember this when we want to specficy/calculate heritability of traits etc...
  return(num.tillers)
}

fun.alpha <- function(auxin,sucrose,max.sucrose){
  main_term = (auxin + 1)/(sucrose + 0.2)
  correction = 1 - (0.15/(sucrose + 0.2))
  alpha = main_term*correction
  return(alpha)
}

### Define Model Parameters

## Cytokinin

# Maximum induction of CK synthesis by Sucrose (mol per second)
a.CK = 0.25
# Strength of CK synthesis inhibition by Auxin (per mol)
b.CK = 0.96
# Base synthesis rate of CK without sucrose and auxin (mol per second)
c.CK = 0.79
# CK degradation rate (per second)
d.CK = 0.99
# Parameter of the Hill function relating sucrose and CK synthesis (mol squared)
k.CK = 0.19

## Strigolactone

# Maximum induction of SL synthesis by auxin (mol per second)
a.SL = 24.89
# Base synthesis rate of SL without auxin (mol per second)
c.SL = 0.34
# SL degradation rate (per second)
d.SL = 0.86
# Parameter of the Hill function relating relating auxin and SL synthesis (mol squared)
k.SL = 294.58

## Sucrose

#percentage increase in Sucrose driven by Cytokinin
a.Sucr = 0

## Signal Integrator (BRC1)

# Parameter relating the production rate of I to SL and sucrose (mol^-1.second^-1)
a.I.3 = 5.64
# Parameter relating the production rate of I to CK (mol.second^-1)
a.I.4 = 287.53
# Base production rate of I (mol.second^-1)
c.I = 0.33
# Constant degradation rate of I (second^-1)
d.I = 0.99
# Strength of CK effect on I production (mol^2)
k.I = 1000
# Minimum inhibiting effect of sucrose on SL response (mol^-2)
u.I.1 = 4.8*(10^(-13))
# Strength of sucrose inhibiting effect on SL response (mol^4)
u.I.2 = 7.10

## Time of Bud Outgrowth

# Intercept of the linear relationship between T and I (day)
m.T.0 = -2.2
# Sensitivity of the time at which elongation starts to I (day.mol^-1)
m.T.1 = 3.5
# Threshold of I above which bud elongation is completely prevented (mol)
I.T.0 = 3
