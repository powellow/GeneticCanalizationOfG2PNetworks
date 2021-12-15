#### Sample Genetic Effects of Shoot Branching Newtork ------------------------------------------------------------------------

if(scenario == paste0("EqualEffects_",n_qtn_per_node,"QtlPerNode")){
	### Genetic Effects for Auxin
	aux.eff <- rep(1.25/n_qtn_per_node,n_qtn_per_node)
	### Genetic Effects for Cytokinins
	cyt.eff <- rep(0.55/n_qtn_per_node,n_qtn_per_node)
	### Genetic Effects for Sucrose
	suc.eff <- rep(1.25/n_qtn_per_node,n_qtn_per_node)
	### Genetic Effects for Strigolatones
	strigo.eff <- rep(0.5/n_qtn_per_node,n_qtn_per_node)
}

if(scenario == paste0("SampledEffects_",n_qtn_per_node,"QtlPerNode")){
	### Genetic Effects for Auxin
	aux.eff <- rnorm(n_qtn_per_node,1.25/n_qtn_per_node,0.01)
	### Genetic Effects for Cytokinins
	cyt.eff <- rnorm(n_qtn_per_node,0.55/n_qtn_per_node,0.01)
	### Genetic Effects for Sucrose
	suc.eff <- rnorm(n_qtn_per_node,1.25/n_qtn_per_node,0.01)
	### Genetic Effects for Strigolatones
	strigo.eff <- rnorm(n_qtn_per_node,0.5/n_qtn_per_node,0.01)
}

if(scenario == paste0("ConstrainedSampledEffects_",n_qtn_per_node,"QtlPerNode")){
	### Genetic Effects for Auxin
	{## Contrain Sum of Effects to Particular Value
	aux.eff <- unlist(purrr::flatten(RandVec(a=0.001, b=0.625, s=1.25, n=n_qtn_per_node, Seed=sample(1:100000, size = 1))))
	if(sum(aux.eff)>1.26|sum(aux.eff)<1.24){stop("Auxin Effects Not Sampled Correctly")}
	}
	### Genetic Effects for Cytokinins
	{## Contrain Sum of Effects to Particular Value
	cyt.eff <- unlist(purrr::flatten(RandVec(a=0.001, b=0.275, s=0.55, n=n_qtn_per_node, Seed=sample(1:100000, size = 1))))
	if(sum(cyt.eff)>0.551|sum(cyt.eff)<0.549){stop("Cytokinin Effects Not Sampled Correctly")}
	}
	### Genetic Effects for Sucrose
	{## Contrain Sum of Effects to Particular Value
	suc.eff <- unlist(purrr::flatten(RandVec(a=0.001, b=0.625, s=1.25, n=n_qtn_per_node, Seed=sample(1:100000, size = 1))))
	if(sum(suc.eff)>1.26|sum(suc.eff)<1.24){stop("Sucrose Effects Not Sampled Correctly")}
	}
	### Genetic Effects for Strigolatones
	{## Contrain Sum of Effects to Particular Value
	strigo.eff <- unlist(purrr::flatten(RandVec(a=0.001, b=0.25, s=0.5, n=n_qtn_per_node, Seed=sample(1:100000, size = 1))))
	if(sum(strigo.eff)>0.51|sum(strigo.eff)<0.49){stop("Strigolactone QTL Effects Not Sampled Correctly")}
	}
}

## Sample QTL indexes

eff.index <- sample(1:(n_qtn*n_chr),n_qtn*n_chr,replace=F) #Random QTN order
#eff.index <- seq(1,(n_qtn*n_chr),1) #Trait QTN each on independent chromosome
aux.eff.index <- eff.index[1:n_qtn_per_node]
cyt.eff.index <- eff.index[(n_qtn_per_node+1):(n_qtn_per_node*2)]
suc.eff.index <- eff.index[((n_qtn_per_node*2)+1):(n_qtn_per_node*3)]
strigo.eff.index <- eff.index[((n_qtn_per_node*3)+1):(n_qtn_per_node*4)]

#### Set Maximum Total Genetic Values of Shoot Branching Newtork ------------------------------------------------------------------------
#### Only Appropriate for sampling equal QTN effects
max.auxin <- 2.5
max.cytokinin <- 1.356
max.strigolactone <- 1.764
max.sucrose <- 2.5
max.time.to.bud <- 17.562 #Corrected by 0.481 to allow normalization as minimum value from network is -0.481
max.signal <- 6

min.auxin <- 0
min.cytokinin <- 0
min.strigolactone <- 0
min.sucrose <- 0
min.time.to.bud <- 0
min.signal <- 0.491

### Simulate our Reference Population of Genotypes (RPG) &  Causal Loci-------------------------------------------------------------------------
founderPop = quickHaplo(nInd=n_founders,
                        nChr=n_chr,
                        segSites=n_qtn+n_snp,
                        inbred=FALSE)

SP = SimParam$new(founderPop) #Store Simulation Parameters

### Create 'empty slots' for Trait Values ---
SP$addTraitA(n_qtn,mean=10,var=3)$setVarE(h2=1.0) #This step is required to allow us to extract the RPG genotypes (line 84).
# *We will not use this trait info in the simulation* #
SP$setTrackPed(TRUE)
SP$setTrackRec(TRUE)
### Create haplotypes of Extreme Parents ----

RPG = newPop(founderPop)


first_parent = matrix(rep(rep(0,(n_qtn*n_chr)),each = ploidy),nrow = ploidy)
second_parent = matrix(rep(rep(1,(n_qtn*n_chr)),each = ploidy),nrow = ploidy)
parents <- rbind(first_parent,second_parent)

## Format haplotypes and Genetic Map to import into AlphaSimR

start <- seq(1,n_seg_sites,n_qtn+n_snp)
stop <- seq(n_qtn+n_snp,n_seg_sites,n_qtn+n_snp)

for (each in 1:n_chr){
  assign(paste0("chr",each),matrix(as.integer(parents[,c(start[each]:stop[each])]),ncol=(n_qtn+n_snp)))
}

haplotypes <- list(chr1,chr2,chr3,chr4)
genMapRep = rep(list(seq(0,1,length.out=n_qtn+n_snp)),n_chr)

## Import Parental Haplotypes  into AlphaSimR
founderPop = newMapPop(genMap = genMapRep,haplotypes)
RPG = newPop(founderPop)

### Create Backcross (Reference (0) Parent)  QTL Allele Frequency = 0.25 ---

##Identify Extreme Parents
ref_hom <- RPG[1] # Parent Homzygous for 0 all
ref_hom_haplo <- pullSegSiteHaplo(ref_hom)
alt_hom <- RPG[2]  # Parent Homzygous for 1 allele

##Create Uniform F1
F1 <- hybridCross(ref_hom,alt_hom)
# F1_geno <- pullSegSiteGeno(F1)
# F1_haplo <- pullSegSiteHaplo(F1)

##Create segregating backcross
#BC1 <- randCross2(ref_hom,F1,nCrosses=1,nProgeny=1000)
#bc1_haplo = pullSegSiteHaplo(BC1)
#bc1_AF <- colMeans(bc1_haplo)
#mean(bc1_AF); var(bc1_AF)

##Create segregating F2
F2 <- self(F1,nProgeny = n.off.f2)

## Assign BC1 Population as Reference Population of Genotypes
# RPG <- BC1

## Assign F2 Population as Reference Population of Genotypes
RPG <- F2

#### Selection Cycles ------------------------------------------------------------------------------------------------

## Create empty dataframes to store results
RPG_Trait <- vector(mode = "list", length = n.cycles)
RPG_VarComp <- vector(mode = "list", length = n.cycles)

trait_results = data.frame(rep=rep(NA,1),cycle = rep(NA,1),trait = rep(NA,1), Metric = rep(NA,1), Mean=rep(NA,1),Variance=rep(NA,1))
branching.network <- data.frame(Index=seq(1,n.off.f2,1),Sucrose.GeneticValue=rep(NA,n.off.f2),Auxin.GeneticValue=rep(NA,n.off.f2),Cytokinin.GeneticValue=rep(NA,n.off.f2),Strigolactone.GeneticValue=rep(NA,n.off.f2))
RPG_Trait[[1]] <- branching.network

t.b.o = c()
freq.total = c()
sel.coef.total = c()
tillers.all = c()

for (cycle.no in 1:n.cycles) {
  #print(cycle.no)

  #### Calculate QTL Allele Frequencies ----

  g.mat = pullSegSiteGeno(RPG)
  #g.mat[,1:5]
  freq = colMeans(g.mat)/ploidy
  freq.total = cbind(freq.total, freq)
  sel.coef = s_hat_function(t(g.mat),n_qtl = nrow(g.mat),n_chr=1)
  sel.coef.total <- cbind(sel.coef.total,sel.coef$allele_sel_coefficient)

  #### Heatmap for Deep Learning Project ----
  if(cycle.no==n.cycles){
    grays = rgb(red = 0:255/255, blue = 0:255/255, green = 0:255/255)
    par(mar = c(0, 0, 0, 0))
    if(heritability_constant==TRUE){
    png(paste0(selection.or.drift,"/",paste(selection.or.drift,paste0("ConstantH2=",H2),paste0("rep",rep_count),sep="_"),".png"),width=224,height=224)
    }else{
      png(paste0(selection.or.drift,"/",paste(selection.or.drift,paste0("H2=",H2),paste0("rep",rep_count),sep="_"),".png"),width=224,height=224)
    }
    heatmap(g.mat,Rowv=NA,Colv=NA,col = grays,scale="none",labRow=NA, labCol=NA)
    dev.off()
  }

  #### Calculate Hormone/Nutrient Genetic Values ----
  branching.network$Auxin.AdditiveGeneticValue  = g.mat[,aux.eff.index] %*% aux.eff
  branching.network$Cytokinin.AdditiveGeneticValue = g.mat[,cyt.eff.index]  %*% cyt.eff
  branching.network$Sucrose.AdditiveGeneticValue = g.mat[,suc.eff.index] %*% suc.eff
  branching.network$Strigolactone.AdditiveGeneticValue = g.mat[,strigo.eff.index]  %*% strigo.eff


  #### Update Hormone Genetic & Phenotypic Values due to Signalling ----
  ## No Change in Auxin
  branching.network$Auxin.TotalGeneticValue <- branching.network$Auxin.AdditiveGeneticValue

  ## Calculate Change in CYT
  branching.network$Cytokinin.TotalGeneticValue <- fun.delta.cyt(a.CK,b.CK,branching.network$Cytokinin.AdditiveGeneticValue,d.CK,k.CK,
  	branching.network$Auxin.TotalGeneticValue,
  	branching.network$Sucrose.AdditiveGeneticValue,
  	branching.network$Cytokinin.AdditiveGeneticValue,
  	Nutrients)

  ## Calculate Change in SL
  branching.network$Strigolactone.TotalGeneticValue <- fun.delta.strigo(a.SL,branching.network$Strigolactone.AdditiveGeneticValue,d.SL,k.SL,
  	branching.network$Auxin.TotalGeneticValue,
  	branching.network$Strigolactone.AdditiveGeneticValue,
  	Nutrients)

  #### Update Nutrient Values ----

  ## Calculate Change in Sucrose
  branching.network$Sucrose.TotalGeneticValue <- fun.delta.suc(a.Sucr,branching.network$Sucrose.AdditiveGeneticValue,
  	branching.network$Cytokinin.TotalGeneticValue)

  #### Calculate Change in Signal Integrator ----
  branching.network$sig.coeff = 0

  branching.network$`Integrator Signal` <- fun.delta.signal(a.I.3,a.I.4,u.I.1,u.I.2,c.I,d.I,k.I,
  	branching.network$Strigolactone.TotalGeneticValue,
  	branching.network$Sucrose.TotalGeneticValue,
  	branching.network$Cytokinin.TotalGeneticValue,
  	branching.network$sig.coeff)

  #### Calculate Simple Control Variable (Alpha) ----
  max.sucrose = sum(suc.eff*ploidy) #calculate maximum possible genetic value of sucrose

  branching.network$`Alpha` <- fun.alpha(branching.network$Auxin.TotalGeneticValue,
  	branching.network$Sucrose.TotalGeneticValue,
  	max.sucrose)

  #### To Bud or Not to Bud ----

  branching.network$`Bud Outgrowth`  <- fun.bud.growth(signal =  branching.network$`Integrator Signal`,threshold = I.T.0)
  if(threshold=="Yes"){
  	branching.network$TimeToBudOutgrowth.TotalGeneticValue <- fun.time.bud.threshold(m.T.0,m.T.1,I.T.0,
    	branching.network$`Integrator Signal`)
  }
  if(threshold=="No"){
  	branching.network$TimeToBudOutgrowth.TotalGeneticValue  <- fun.time.bud(m.T.0,m.T.1,I.T.0,
    	branching.network$`Integrator Signal`) + 0.481
  }

  ### Specify Broad Sense Hetritability for 'Time to Bud Outgrowth'
  if(heritability_constant==TRUE){
    var.t2bf <- var(branching.network$Sucrose.TotalGeneticValue)
    var.p <- var.t2bf/as.numeric(H2)
    var.e <- var.p - var.t2bf
  }else{
    if(cycle.no==1){
      var.t2bf <- var(branching.network$Sucrose.TotalGeneticValue)
      var.p <- var.t2bf/as.numeric(H2)
      var.e <- var.p - var.t2bf
    }
  }

  #Add random error to generate phenotypes for selection
  branching.network$`Time to Bud Outgrowth` <- branching.network$TimeToBudOutgrowth.TotalGeneticValue + rnorm(n.off.f2,0,var.e)
  branching.network$`Sucrose` <- branching.network$Sucrose.TotalGeneticValue + rnorm(n.off.f2,0,var.e)

  #### Calculate Interaction Genetic Values (Total - Additive Genetic Values) ----

  #### Store Current Reference Population of Genotypes (RPG) Data ----

  RPG_Trait[[cycle.no]] <- branching.network

  #### Estimate Variance Componenets of Branching Network ----
  # if(cycle.no%in%cycles.of.interest){
  #   #unique_records <- distinct_at(branching.network, vars(contains("TotalGeneticValue"),target.trait))
  #   unique_records <- branching.network
  #
  #   y <- unique_records$`Time to Bud Outgrowth`
  #   unique_records$`Auxin` <- unique_records$Auxin.TotalGeneticValue
  #   unique_records$`Cytokinin` <- unique_records$Cytokinin.TotalGeneticValue
  #   unique_records$`Strigolactone` <- unique_records$Strigolactone.TotalGeneticValue
  #   unique_records$`Sucrose` <- unique_records$Sucrose.TotalGeneticValue
  #
  #   full_model <- asreml(fixed = y ~ 1, random = ~ `Auxin` + `Cytokinin` +
  #                          `Strigolactone` + `Sucrose` + `Auxin`:`Cytokinin` +
  #                          `Auxin`:`Strigolactone` +
  #                          `Cytokinin`:`Sucrose` +
  #                          `Strigolactone`:`Sucrose`,
  #                        data=unique_records,trace=F)
  #   # reduced_model <- asreml(fixed = y ~ 1, random = ~ `Auxin` + `Cytokinin` +
  #   #                        `Strigolactone` + `Sucrose`,
  #   #                      data=unique_records,maxit=20, trace=F)
  #
  #   RPG_VarComp[[cycle.no]] <- summary(full_model)$varcomp
  # }

  #### Calculate Trait Parameters Required for Figures ----

  # Calculate Relative Genetic Values
  branching.network$Normalised.Sucrose.AdditiveGeneticValue <- fun.proportion(branching.network$Sucrose.AdditiveGeneticValue,suc.eff)
  branching.network$Normalised.Auxin.AdditiveGeneticValue <- fun.proportion(branching.network$Auxin.AdditiveGeneticValue,aux.eff)
  branching.network$Normalised.Cytokinin.AdditiveGeneticValue <- fun.proportion(branching.network$Cytokinin.AdditiveGeneticValue,cyt.eff)
  branching.network$Normalised.Strigolactone.AdditiveGeneticValue <- fun.proportion(branching.network$Strigolactone.AdditiveGeneticValue,strigo.eff)

  branching.network$Normalised.Sucrose.TotalGeneticValue <- (branching.network$Sucrose.TotalGeneticValue - min.sucrose)/(max.sucrose - min.sucrose)
  branching.network$Normalised.Auxin.TotalGeneticValue <- (branching.network$Auxin.TotalGeneticValue - min.auxin)/(max.auxin - min.auxin)
  branching.network$Normalised.Cytokinin.TotalGeneticValue <- (branching.network$Cytokinin.TotalGeneticValue - min.cytokinin)/(max.cytokinin - min.cytokinin)
  branching.network$Normalised.Strigolactone.TotalGeneticValue <- (branching.network$Strigolactone.TotalGeneticValue - min.strigolactone)/(max.strigolactone - min.strigolactone)
  branching.network$Normalised.IntegratorSignal.TotalGeneticValue <- (branching.network$`Integrator Signal` - min.signal)/(max.signal - min.signal)
  branching.network$Normalised.TimeToBudOutgrowth.TotalGeneticValue <- (branching.network$TimeToBudOutgrowth.TotalGeneticValue - min.time.to.bud)/(max.time.to.bud - min.time.to.bud)

  # Relative Phenotyic Values
  branching.network$Normalised.TimeToBudOutgrowth <- (branching.network$`Time to Bud Outgrowth` - min.time.to.bud)/(max.time.to.bud - min.time.to.bud)

  gv_mean <- sapply(branching.network[,genetic_traits],function(x) mean(na.omit(x)))
  gv_variance <- sapply(branching.network[,genetic_traits],function(x) var(na.omit(x)))
  pheno_mean <- sapply(branching.network[,phenotypic_traits],function(x) mean(na.omit(x)))
  pheno_variance <- sapply(branching.network[,phenotypic_traits],function(x) var(na.omit(x)))
  cycle_genetic_trait_results <- cbind(rep = rep(rep_count,length(genetic_traits)),cycle = rep(cycle.no,length(genetic_traits)), trait = genetic_traits, Metric = rep("Genetic Value",length(genetic_traits)), Mean = gv_mean,Variance = gv_variance)
  cycle_phenotypic_trait_results <- cbind(rep = rep(rep_count,length(phenotypic_traits)),cycle = rep(cycle.no,length(phenotypic_traits)), trait = phenotypic_traits, Metric = rep("Phenotypic Value",length(genetic_traits)),Mean = pheno_mean,Variance = pheno_variance)
  cycle_trait_results <- rbind(cycle_genetic_trait_results,cycle_phenotypic_trait_results)
  trait_results = rbind(trait_results,cycle_trait_results)

  #### Select New Parents  ----

  if(selection.or.drift=="Selection"){
    if(selection.direction=="Min"){
      filtered.lines <-  branching.network
      # filtered.lines <-  branching.network %>% filter(Normalised.Auxin.TotalGeneticValue>=0.25,Normalised.Cytokinin.TotalGeneticValue>=0.25,Normalised.Strigolactone.TotalGeneticValue>=0.25) %>%
      #   filter(Normalised.Auxin.TotalGeneticValue<=0.80,Normalised.Cytokinin.TotalGeneticValue<=0.75,Normalised.Strigolactone.TotalGeneticValue<=0.75)
      selected.lines <- as.numeric(filtered.lines[order(filtered.lines[,target.trait],decreasing=F),"Index"][1:n.selected]) # Sort Individuals based on Time to Bud. Add 'decreasing = F' for min Time to Bud
    }
    if(selection.direction=="Max"){
      filtered.lines <-  branching.network
      # filtered.lines <-  branching.network %>% filter(Normalised.Auxin.TotalGeneticValue>=0.25,Normalised.Cytokinin.TotalGeneticValue>=0.25,Normalised.Strigolactone.TotalGeneticValue>=0.25) %>%
      #   filter(Normalised.Auxin.TotalGeneticValue<=0.80,Normalised.Cytokinin.TotalGeneticValue<=0.75,Normalised.Strigolactone.TotalGeneticValue<=0.75)
      selected.lines <- as.numeric(filtered.lines[order(filtered.lines[,target.trait],decreasing=T),"Index"][1:n.selected]) # Sort Individuals based on Time to Bud. Add 'decreasing = T' for max Time to Bud
    }
  }
  if(selection.or.drift=="Drift"|selection.or.drift=="drift"){
    selected.lines <- sample(as.numeric(noquote(rownames(branching.network))),n.selected) #Random Sample of Individuals
  }
  #selected.lines = sample(selected.lines) #randomise order before crossing

  #### Cross New Parents ----

  RPG = randCross(RPG[selected.lines],nCrosses = n.off.f2)

}

#### Allele Frequency Changes Over Selection Cycles at Causal QTL ----

head(freq.total)

# make long vector for the frequencies to plot in ggplot
long.vec  = c()
sel.coef = c()
cycle  = c()
trait_qtn_order = rep(NA,n_qtn*n_chr);
trait_qtn_order[aux.eff.index] <- "Auxin"; trait_qtn_order[cyt.eff.index] <- "Cytokinin";trait_qtn_order[suc.eff.index] <- "Sucrose";trait_qtn_order[strigo.eff.index] <- "Strigolactone"
gene.type = c()
qtl       = c()
effect.size  = rep(NA,n_qtn*n_chr);
effect.size [aux.eff.index] <- aux.eff; effect.size [cyt.eff.index] <- cyt.eff; effect.size[suc.eff.index] <- suc.eff; effect.size[strigo.eff.index] <- strigo.eff

### Correct for QTN index

for (i in 1:ncol(freq.total)){

  long.vec  = c(long.vec, freq.total[,i])
  sel.coef = c(sel.coef,sel.coef.total[,i])
  cycle  = c(cycle, rep(i,n_qtn*n_chr))
  gene.type = c(gene.type,trait_qtn_order)
  qtl       = seq(1,n_qtn*n_chr)
  effect.size = effect.size
}

QTL.info = data.frame (rep       = rep_count,
                    freq      = long.vec,
                    sel.coef  = sel.coef,
                    gen       = cycle,
                    gene.type = gene.type,
                    qtl       = qtl,
                    effect.size = effect.size)

#### Save Network Changes Over Selection Cycles ----

genetic_values <- trait_results %>%
  na.omit() %>%
  dplyr::filter(Metric=="Genetic Value") %>%
  mutate(across(Mean, as.numeric)) %>%
  mutate(across(cycle, as.integer))


phenotypic_values <- trait_results %>%
  na.omit() %>%
  dplyr::filter(Metric=="Phenotypic Value") %>%
  mutate(across(Mean, as.numeric)) %>%
  mutate(across(cycle, as.integer))

saveRDS(QTL.info,file = paste0(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),Nutrients,paste0("H2=",H2),sep="/"),"/",paste0("QtlData_Rep",rep_count,".rds")))
saveRDS(genetic_values,file = paste0(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),Nutrients,paste0("H2=",H2),sep="/"),"/",paste0("GeneticData_Rep",rep_count,".rds")))
saveRDS(phenotypic_values,file = paste0(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),Nutrients,paste0("H2=",H2),sep="/"),"/",paste0("PhenotypicData_Rep",rep_count,".rds")))
saveRDS(purrr::compact(RPG_VarComp),file = paste0(paste(cwd,"output",scenario,target.trait,selection.or.drift,selection.direction,paste0(threshold,"_threshold"),Nutrients,paste0("H2=",H2),sep="/"),"/",paste0("VarianceComponenstEstimates_Rep",rep_count,".rds")))
