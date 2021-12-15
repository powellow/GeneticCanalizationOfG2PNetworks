
### Estimate the values of two-locus genotype classes
# INPUT:
# data = tibble
# locus1 = column name of 1st locus
# locus2 = column name of 2nd locus
# trait_value =  column name of trait values to be supplied
# table = option to display output in a 3x3 table or in long/tidy format
# OUTPUT:
# 1. 3x3 matrix of two-locus epistatic values
# OR
# 2. 9x3 tibble (Primary QTL, Secondary QTL, trait value)

non_epistatic_values <- function(data,locus1,locus2,trait_value,table){
  if (table==TRUE){
    data %>%
      group_by_(locus1,locus2) %>%
      summarise(GV = mean(!! sym(trait_value),na.rm=T)) %>%
      spread(locus1,GV)
  } else {
    data %>%
      group_by_(locus1,locus2) %>%
      summarise(GV = mean(!! sym(trait_value),na.rm=T))
  }
}


### Simulate Digenic Epistasis for a Single Interaction
# INPUT:
# interaction_type = specify underlying epistatic gene action
# mean = average value of distribution from which epistasis values are sampled. Default is 0 i.e symmetrical epistasis
# variance = variance of distribution from which epistasis values are sampled
# dominance_degree = level of dominance underpinning epistasis for heterzygous genotypes
# supply_effects = option to supply the mean value of distribution from which epistasis values are sampled.i.e allows directed epistasis (average value != 0)
# OUTPUT:
# 3x3 matrix of two-locus epistatic values

epi_values <- function(interaction_type,mean,variance,dominance_degree,supply_effects){ #in the middle of integrating excel sheet nominated interactions with old function
  if(supply_effects==TRUE){
    sampled_epi_effect <- rnorm(1,mean,variance) #sample effect with mean=value in interaction matrix
  }
  else{
    sampled_epi_effect <- rnorm(1,0,variance) #sample effect with mean=0
  }
  epi_effects <- matrix(0,nrow=3,ncol=3) #create empty 3x3 Genotype matrix to store epistasis vlaues
  if (interaction_type =="AA"){
    # process to generate Additive-by-Additive Epistasis
    # sampled effect assigned to homozygote genotypes (AABB/aabb)
    # sampled effect used to calculate effects of other genotypes

    epi_effects[1,1] <- sampled_epi_effect
    epi_effects[1,3] <- -sampled_epi_effect
    epi_effects[3,1] <- -sampled_epi_effect
    epi_effects[3,3] <- sampled_epi_effect
  }
  if (interaction_type =="AD"){
    # process to generate Additive (Locus1)-by-Dominance (Locus2) Epistasis
    # sampled effect assigned to AABb genotype
    # negative of sampled effect assigned to AaBB genotype
    # dominance scaled sampled effect assigned to aa** genotypes
    # negative dominance scaled sampled effect assigned to AA** genotypes

    epi_effects[,1] <- (-dominance_degree)*sampled_epi_effect
    epi_effects[2,1] <- sampled_epi_effect
    epi_effects[,3] <- (dominance_degree)*sampled_epi_effect
    epi_effects[2,3] <- -sampled_epi_effect
  }
  if (interaction_type =="DA"){
    # process to generate Dominance (Locus1)-by-Additive (Locus2) Epistasis

    # sampled effect assigned to AaBB genotype
    # negative of sampled effect assigned to Aabb genotype
    # dominance scaled sampled effect assigned to **bb genotypes
    # negative dominance scaled sampled effect assigned to **BB genotypes

    epi_effects[1,] <- (-dominance_degree)*sampled_epi_effect
    epi_effects[1,2] <- sampled_epi_effect
    epi_effects[3,] <- (dominance_degree)*sampled_epi_effect
    epi_effects[3,2] <- -sampled_epi_effect
  }
  if (interaction_type =="DD"){
    # process to generate Dominance-by-Dominance Epistasis
    # sampled effect assigned to AaBb genotype
    # negative dominance scaled sampled effect assigned to **Bb genotypes
    # negative dominance scaled sampled effect assigned to Aa** genotypes
    # dominance^2 scaled sampled effect assigned to homozygous (AABB/AAbb/aaBB/aabb) genotypes

    epi_effects[1,1] <- (dominance_degree^2)*sampled_epi_effect
    epi_effects[1,3] <- (dominance_degree^2)*sampled_epi_effect
    epi_effects[,2] <- (-dominance_degree)*sampled_epi_effect
    epi_effects[2,] <- (-dominance_degree)*sampled_epi_effect
    epi_effects[2,2] <- sampled_epi_effect
    epi_effects[3,1] <-(dominance_degree^2)*sampled_epi_effect
    epi_effects[3,3] <- (dominance_degree^2)*sampled_epi_effect
  }
  if (interaction_type =="Random"){
    # process to a mix of sources of Epistasis
    # each genotype class has a epistasi value randomly sampled
    epi_effects <- matrix(rnorm(9,mean,variance),nrow=3,ncol=3)
  }
  return(epi_effects)
}


### Simulate Digenic Epistasis for a Multiple Interactions
# INPUT:
# interaction_matrix = interaction_matrix = matrix (QTL x QTL), Non-zero values indicate interaction.
# interaction_type = specify underlying epistatic gene action
# variance = variance of distribution from which epistasis values are sampled
# OUTPUT:
# list containing 3x3 matrices of two-locus epistatic values for genotype classes

per_locus_interaction_deviations <- function(interaction_matrix,interaction_type,variance){
  per_locus_interaction <- interaction_matrix %>%
    flatten %>%
    lapply(function(x) if(x!=0)  epi_values(interaction_type,mean=x,variance = variance,supply_effects=F) else x) %>%
    matrix(nrow=n_qtl*n_chr,ncol=n_qtl*n_chr)
  p_l_i_d <- per_locus_interaction[lapply(per_locus_interaction,length)>1]
  return(p_l_i_d)
}

### Identify & Summarise Interactions between QTL
# INPUT:
# interaction_matrix = matrix (QTL x QTL), Non-zero values indicate interaction.
# per_locus_interactions = list containing 3x3 matrices of two-locus epistatic values for genotype classes
# OUTPUT:
# tibble with columns specfiying;
# 1. Directionality of QTL-QTL interactions
# 2. Number of interactions per QTL
# 3. Epistatic Values of each QTL-QTL interaction (3x3 matrices)
# 4. Directionality of average epistatic values of each QTL-QTL interaction. i.e. "postive","neutral","negative".
# 4. Variance of each QTL-QTL interaction

summarise_interactions <- function(interaction_matrix,per_locus_interaction){
  interaction_dataframe <- which(interaction_matrix!=0,arr.ind=TRUE) %>%
    as.tibble() %>%
    mutate(position= which(interaction_matrix!=0)) %>%
    rename(Primary_QTL = col, Secondary_QTL = row) %>%
    relocate(Primary_QTL,Secondary_QTL) %>%
    mutate(interaction_deviations = per_locus_interaction[rank(position)]) %>%
    mutate(interaction_variance = unlist(lapply(interaction_deviations,function(x) sum(x^2)))) %>%
    group_by(Primary_QTL) %>%
    mutate(n_interactions = length(Secondary_QTL)) %>%
    mutate(average_interaction_sign = ifelse(mean(interaction_deviations[[1]])>0, "positive",
                                        ifelse(mean(interaction_deviations[[1]])<0, "negative", "neutral"))) %>%
    nest() %>%
    rename(interaction_info = data)
  return(interaction_dataframe)
}

### Create Genetic Architecure Object For Traits
# INPUT:
# trait = traint number (eg. 1, 2, etc)
# additive_effects = additive qtl effects for trait
# dominance_effects = dominance qtl effects for trait
# interaction_effects = epistatic effects for trait from qtl interactions
# append = Option to specfiy if you would like to add another trait to existing gen_architecture object
# gen_architecture = previous genetic architecture object (only called if append = T)
# OUTPUT:
# genetic architecture object (nrow = number of QTL, ncol = number of traits + 1[QTL ID's]) with nested columns for each trait
# nested columns for each trait store additive. dominace and interaction effects for each trait

gen_architecture <- function(trait,additive_effects,dominance_effects=NULL,interaction_effects=NULL,append=F,gen_architecture=NULL){
  Primary_QTL <- c(1:(length(additive_effects)))
  if(append==TRUE){
    if(missing(dominance_effects)==T & missing(interaction_effects)==T){
      gen_architecture  <- cbind(gen_architecture,additive_effects = additive_effects,dominance_effects = rep(0,length(Primary_QTL)),interaction_effects = rep(0,length(Primary_QTL))) %>% as.tibble()
    }
    if(missing(dominance_effects)==F & missing(interaction_effects)==T){
      gen_architecture  <- cbind(gen_architecture,additive_effects = additive_effects,dominance_effects = dominance_effects,interaction_effects = rep(0,length(Primary_QTL ))) %>% as.tibble()
    }
    if(missing(dominance_effects)==F & missing(interaction_effects)==F){
      gen_architecture  <- cbind(gen_architecture,additive_effects = additive_effects,dominance_effects = dominance_effects) %>% as.tibble() %>% left_join(interaction_effects,by="Primary_QTL")
    }
  } else{
    if(missing(dominance_effects)==T & missing(interaction_effects)==T){
      gen_architecture  <- cbind(Primary_QTL,additive_effects = additive_effects,dominance_effects = rep(0,length(Primary_QTL)),interaction_effects = rep(0,length(Primary_QTL))) %>% as.tibble() %>% rename(additive_effects = V2)
    }
    if(missing(dominance_effects)==F & missing(interaction_effects)==T){
      gen_architecture  <- cbind(Primary_QTL,additive_effects = additive_effects,dominance_effects = dominance_effects,interaction_effects = rep(0,length(Primary_QTL ))) %>% as.tibble()
    }
    if(missing(dominance_effects)==F & missing(interaction_effects)==F){
      gen_architecture  <- cbind(Primary_QTL,additive_effects = additive_effects,dominance_effects = dominance_effects) %>% as.tibble() %>% left_join(interaction_effects,by="Primary_QTL")
    }
  }
  gen_architecture <- gen_architecture %>%
    nest(-c(Primary_QTL,contains("trait"))) %>%
    rename_with(.fn = ~paste0("trait", trait), .cols = "data")
  return(gen_architecture)
}

### Calculate Epistatic Values for Individuals
# INPUT:
# interactions = dataframe including interaction info (output from 'identify_interactions' function)
# genotypes = matrix of genotypes for reference population of genotypes
# OUTPUT:
# numeric vector of epistatic values. Same order as the genotypes matrix

extract_EV <- function(interactions,genotypes){
  unnested_interactions <- interactions %>% unnest(cols=c(interaction_info))
  genotypes_at_interacting_loci <- apply(unnested_interactions[,1:2],MARGIN=1,function(x) select(genotypes,all_of(x)))
  first_level <- seq_along(genotypes_at_interacting_loci)
  second_level <- seq_along(genotypes_at_interacting_loci[[1]]$Primary_QTL)
  individual_EV <- lapply(map(first_level,function(x) map(second_level, function(y) unnested_interactions$interaction_deviations[[x]][genotypes_at_interacting_loci[[x]]$Primary_QTL[y]+1,genotypes_at_interacting_loci[[x]]$Secondary_QTL[y]+1])),unlist) %>% pmap(.,.f=sum) %>% unlist()
  return(individual_EV)
}

extract_EV_diag <- function(interactions,genotypes){
  unnested_interactions <- interactions %>% unnest(cols=c(interaction_info))
  genotypes_at_interacting_loci <- apply(unnested_interactions[,1:2],MARGIN=1,function(x) select(genotypes,all_of(x)))
  individual_EV <- map(seq_along(genotypes_at_interacting_loci),function(x) diag(unnested_interactions$interaction_deviations[[x]][genotypes_at_interacting_loci[[x]]$Primary_QTL+1,genotypes_at_interacting_loci[[x]]$Secondary_QTL+1])) %>% pmap(.,.f=sum)
  individual_EV <- lapply(individual_EV,unlist)
  return(individual_EV)
}

### Calculate Total Genotypic Values for Individuals
# INPUT:
# interactions = dataframe including interaction info (output from 'identify_interactions' function)
# genotypes = matrix of genotypes for reference population of genotypes
# non_epistatic_qtl_effects = numeric vector of qtl effects not including epistasis (currently only additive effects)
# OUTPUT:
# RPG tibble with IDs, genotypes, non-epistatic trait values, epistatic trait values & total genotypic values

calc_GV <- function(id,genotypes,non_epistatic_qtl_effects,dominance_present,dom_qtl_effects,interactions_present,interactions,effects,population,trait,append){
  if(append==TRUE){
    genotypes <- population %>%
      unnest(qtl_genotypes) %>%
      ungroup() %>%
      select(-id,-contains("trait"))
  }
  dom_genotypes <-  as.matrix(genotypes); dom_genotypes[which(dom_genotypes!=1)] <- 0
  AGV <- as.matrix(genotypes)%*%non_epistatic_qtl_effects
  if(dominance_present==FALSE){
    DGV <- as.numeric(rep(0,nrow(genotypes)))
  }
  else{
    DGV <- as.matrix(dom_genotypes)%*%dom_qtl_effects
  }
  NEV <- AGV + DGV
  if(interactions_present==FALSE){
    ind_EV <- as.numeric(rep(0,nrow(genotypes)))
  }
  else{
    ind_EV <- extract_EV(interactions,genotypes)
  }
  if(append==TRUE){
    if(effects==TRUE){
      ind_data <- cbind(population,additive_values=AGV,dominance_values=DGV,non_epistatic_values=NEV,epistatic_values=ind_EV) %>%
        mutate(genotypic_values=non_epistatic_values+epistatic_values) %>%
        # group_by(id) %>%
        nest(-c(id,qtl_genotypes,paste0("trait",c(1:(trait-1))))) %>%
        rename_with(.fn = ~paste0("trait", trait), .cols = "data")%>%
        ungroup()
    }
    else{
      ind_data <- cbind(population,additive_values=AGV,dominance_values=DGV,non_epistatic_values=NEV,epistatic_values=ind_EV) %>%
        mutate(estimated_genotypic_values=non_epistatic_values+epistatic_values) %>%
        # group_by(id) %>%
        nest(-c(id,qtl_genotypes)) %>%
        rename_with(.fn = ~paste0("trait", trait), .cols = "data")%>%
        ungroup()
    }
  }
  else{
    if(effects==TRUE){
      ind_data <- cbind(id,genotypes,additive_values=AGV,dominance_values=DGV,non_epistatic_values=NEV,epistatic_values=ind_EV) %>%
        mutate(genotypic_values=non_epistatic_values+epistatic_values) %>%
        # rownames_to_column(var = "id") %>%
        # group_by(id) %>%
        nest(-c(id,additive_values,dominance_values,non_epistatic_values,epistatic_values,genotypic_values)) %>%
        rename(qtl_genotypes = data) %>%
        nest(-c(id,qtl_genotypes)) %>%
        rename_with(.fn = ~paste0("trait", trait), .cols = "data")%>%
        ungroup()
    }
    else{
      ind_data <- cbind(id,genotypes,additive_values=AGV,dominance_values=DGV,non_epistatic_values=NEV,epistatic_values=ind_EV) %>%
        mutate(estimated_genotypic_values=non_epistatic_values+epistatic_values) %>%
        # rownames_to_column(var = "id") %>%
        # group_by(id) %>%
        nest(-c(id,additive_values,dominance_values,non_epistatic_values,epistatic_values,estimated_genotypic_values)) %>%
        rename(qtl_genotypes = data) %>%
        nest(-c(id,qtl_genotypes)) %>%
        rename_with(.fn = ~paste0("trait", trait), .cols = "data") %>%
        ungroup()
    }
  }
  ind_data <- ind_data %>% mutate(id = as.character(id))
  return(ind_data)
}

### Calculate Total Genotypic Values for End-Point Trait (Sum Component Traits)
## Currently only works for: Trait 3 = Trait 1 + Trait 2
# INPUT:
# RPG tibble for Trait 1
# RPG tibble for Trait 2
# OUTPUT:
# RPG tibble for Trait 3

calc_GV_EPT <- function(trait_data,component_trait1,component_trait2,end_point_trait){
  trait1_genetic_values <- trait_data[[paste0("trait",component_trait1)]] %>%
    simplify_all() %>%
    transpose() %>%
    set_names(c("additive_values","dominance_values","non_epistatic_values","epistatic_values","genotypic_values"))
  trait2_genetic_values <- trait_data[[paste0("trait",component_trait2)]] %>%
    simplify_all() %>%
    transpose() %>%
    set_names(c("additive_values","dominance_values","non_epistatic_values","epistatic_values","genotypic_values"))
  end_point_trait_data <- trait_data %>% add_column(additive_values = 0,
                                                    dominance_values = 0,
                                                    non_epistatic_values = 0,
                                                    epistatic_values = 0,
                                                    genotypic_values = trait1_genetic_values[["genotypic_values"]] + trait2_genetic_values[["genotypic_values"]]) %>%
    nest(-c(id,qtl_genotypes,paste0("trait",component_trait1),paste0("trait",component_trait2))) %>%
    rename_with(.fn = ~paste0("trait", end_point_trait), .cols = "data")
  return(end_point_trait_data)
}

### Calculate Total Genotypic Values for End-Point Trait (Multiply Component Traits)
## Currently only works for: Trait 3 = Trait 1 + Trait 2
# INPUT:
# RPG tibble for Trait 1
# RPG tibble for Trait 2
# OUTPUT:
# RPG tibble for Trait 3

calc_GV_EPT_Product <- function(trait_data,component_trait1,component_trait2,end_point_trait){
  trait1_genetic_values <- trait_data[[paste0("trait",component_trait1)]] %>%
    simplify_all() %>%
    transpose() %>%
    set_names(c("additive_values","dominance_values","non_epistatic_values","epistatic_values","genotypic_values"))
  trait2_genetic_values <- trait_data[[paste0("trait",component_trait2)]] %>%
    simplify_all() %>%
    transpose() %>%
    set_names(c("additive_values","dominance_values","non_epistatic_values","epistatic_values","genotypic_values"))
  end_point_trait_data <- trait_data %>% add_column(additive_values = 0,
                                                    dominance_values = 0,
                                                    non_epistatic_values = 0,
                                                    epistatic_values = 0,
                                                    genotypic_values = trait1_genetic_values[["genotypic_values"]] * trait2_genetic_values[["genotypic_values"]]) %>%
    nest(-c(id,qtl_genotypes,paste0("trait",component_trait1),paste0("trait",component_trait2))) %>%
    rename_with(.fn = ~paste0("trait", end_point_trait), .cols = "data")
  return(end_point_trait_data)
}

### Calculate Total Genotypic Values for End-Point Trait (Quadratic Effect of Trait 1 on Trait 2)
# INPUT:
# RPG tibble for Trait 1
# RPG tibble for Trait 2
# OUTPUT:
# RPG tibble for Trait 3

calc_GV_EPT_Quad <- function(trait_data,component_trait1,component_trait2,end_point_trait){
  trait1_genetic_values <- trait_data[[paste0("trait",component_trait1)]] %>%
    simplify_all() %>%
    transpose() %>%
    set_names(c("additive_values","dominance_values","non_epistatic_values","epistatic_values","genotypic_values"))
  trait2_genetic_values <- trait_data[[paste0("trait",component_trait2)]] %>%
    simplify_all() %>%
    transpose() %>%
    set_names(c("additive_values","dominance_values","non_epistatic_values","epistatic_values","genotypic_values"))
  end_point_trait_data <- trait_data %>% add_column(additive_values = 0,
                                                    dominance_values = 0,
                                                    non_epistatic_values = 0,
                                                    epistatic_values = 0,
                                                    genotypic_values = ((trait2_genetic_values[["additive_values"]]) - (trait1_genetic_values[["additive_values"]] - ((-0.2117*(trait1_genetic_values[["additive_values"]]^2)) + (4.2341*(trait1_genetic_values[["additive_values"]])))))) %>%
    nest(-c(id,qtl_genotypes,paste0("trait",component_trait1),paste0("trait",component_trait2))) %>%
    rename_with(.fn = ~paste0("trait", end_point_trait), .cols = "data")
  return(end_point_trait_data)
}

### Calculate Trait Values with Trait Level Interactions
# INPUT:
# RPG tibble for Trait 1
# RPG tibble for Trait 2
# OUTPUT:
# RPG tibble for Trait 3

calc_GV_Trait_Interaction <- function(trait_data,trait1,trait2,interaction){
  trait1_genetic_values <- trait_data[[paste0("trait",trait1)]] %>%
    simplify_all() %>%
    transpose() %>%
    set_names(c("additive_values","dominance_values","non_epistatic_values","epistatic_values","genotypic_values"))
  trait2_genetic_values <- trait_data[[paste0("trait",trait2)]] %>%
    simplify_all() %>%
    transpose() %>%
    set_names(c("additive_values","dominance_values","non_epistatic_values","epistatic_values","genotypic_values"))
  interaction_trait_data <- trait_data %>% unnest(paste0("trait",trait2)) %>%
    mutate(genotypic_values := !!interaction) %>%
    nest(-c(id,qtl_genotypes,contains("trait"))) %>%
    rename_with(.fn = ~paste0("trait", trait2), .cols = "data")
  return(interaction_trait_data)
}


### Calculate Quant Gen Parameters with Least Squares
# INPUT:
# RPG = output from 'calc_GV' function
# OUTPUT:
# tibble (nrow=nQTL) with per QTL records: QTL ID's, allele substiution effects, estimated & true additive genetic variance, RPG genotypes, RPG non-epistatic trait values, RPG epistatic trait values & RPG total genotypic values

GLS_Estimates <- function(RPG,trait){
  parameter_estimates <- RPG %>%
    ungroup() %>%
    select(qtl_genotypes,paste0("trait",trait)) %>%
    unnest(cols = c(qtl_genotypes,paste0("trait",trait))) %>%
    gather(QTL, value, -c(additive_values,dominance_values,non_epistatic_values,epistatic_values,genotypic_values)) %>%
    nest(-QTL) %>%
    mutate(fit = map(data, ~ lm(value ~ genotypic_values,data=.x)), # estimate allele subst effect
           allele_freq =  as.numeric(map(data, ~ (mean(.$value)/2))), #estimate qtl allele frequency
           tidied = map(fit, tidy)) %>%
    unnest(tidied) %>%
    rename(RPG_data = data,trait_value_supplied = term) %>%
    filter(trait_value_supplied=="genotypic_values") %>%
    rename(allele_subst_effect = estimate) %>%
    mutate(add_genic_var = (2*allele_freq*(1-allele_freq)*(allele_subst_effect^2))) # estimate additive genic variance from est. allele subst. effect
  return(parameter_estimates)
}

### Assign Trait Value to ASR RPG Object
# INPUT:
# RPG_Trait_Data = output from 'calc_GV' function
# RPG = AlphaSimR population object
# trait = trait number to be selected on

assign_values <- function(trait_object,RPG,trait){
  RPG@gv <- matrix(transpose(trait_object[[paste0("trait",trait)]])[["genotypic_values"]])
  return(RPG)
}

### Select Top Individuals
# INPUT:
# trait_object =  RPG_Trait_DAta (output from 'calc_GV' function)
# nInd = Number of individuals to select
# selection_criteria = trait value to select on. i.e additive or genotypic value
# trait = trait number to be selected on

select_parents <- function(trait_object,nInd,selection_criteria,trait){
  selection_criteria <- as.name(selection_criteria)
  RPG_parents <- trait_object %>% unnest(paste0("trait",trait)) %>% ungroup()
  RPG_parents <- RPG_parents %>% arrange(desc(!!selection_criteria)) %>%
    slice(1:nInd) %>%   #ignores 'tied' individuals
    nest(-c(id,qtl_genotypes,contains("trait"))) %>%
    rename_with(.fn = ~paste0("trait", trait), .cols = "data")
  return(RPG_parents)
}

### Truncate ASR RPG Object to match selected individuals
# INPUT:
# selected_population = output from 'selec_parents' function
# RPG = AlphaSimR population object
select_ASR <- function(selected_population,ASR_RPG){
  ASR_RPG_parents <- ASR_RPG[which(ASR_RPG@id%in%selected_population$id)]
  return(ASR_RPG_parents)
}

### Estimate Breeding Values
# INPUT:
# RPG_Trait_Data = output from 'calc_GV' function
# Allele Substitution Effects = output from GLS_Estimates()
# Trait of Interest
# OUTPUT:
# Column of Estimate Breeding Values stored in RPG_Trait_Data

estimate_BV <- function(RPG,allele_subst_effects,trait,max_trait_value){ ### How to nest only columns relavant to trait specified
  trait_numbers <- c(1:max_trait_value);other_trait_numbers <- trait_numbers[!trait_numbers%in%trait]
  EBV <- RPG %>%
    ungroup() %>%
    select(qtl_genotypes) %>%
    unnest(qtl_genotypes) %>%
    apply(MARGIN = 1,function(x) x%*%allele_subst_effects)
  RPG <- RPG %>%
    ungroup() %>%
    unnest(paste0("trait",trait)) %>%
    mutate(estimated_breeding_values = EBV) %>%
    nest(-c(id,qtl_genotypes,paste0("trait",c(other_trait_numbers)))) %>%
    rename_with(.fn = ~paste0("trait", trait), .cols = "data")
  return(RPG)
}

# visualisation function
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


### Estimate Per Locus Respone to Selection (R hat)
# INPUT:
# geno = matrix of genotypes for reference population of genotypes (RPG)
# pheno = vector of phenotypes for RPG
# all_freq_change = allele frequency change from RPG at time 0 -> RPG at time t
# true_effects = vector of pre-calculated allele substitution effects
# spread = option to output tibble with 1 row per locus, with a column for each statistic
# plot = option to plot estimated R hat & Allele Substitution Effects
# OUTPUT:
# Tibble with allele substiution effects (ASE), standard error of ASE and response to selection (R hat) for each locus
Rhat_func <- function(geno,pheno,all_freq_change,true_effects=NULL,spread=TRUE,plot=TRUE){
  library(latex2exp)
  ### Estimate Allele Substitution Effects for each locus
  blupModel <- mixed.solve(pheno, Z=geno,K=NULL,SE=T,return.Hinv=FALSE,method="ML")
  effects <- blupModel$u #store effect estimates
  effects_se <- blupModel$u.S #store standard error of effect estimates

  if(!is.null(true_effects)==T){ #Calculate per locus response to selection with known/true values
    per_locus_selection_response <- 2*(true_effects*all_freq_change)
  }

  if(is.null(true_effects)==T){ #Calculate per locus response to selection with estimated values
    per_locus_selection_response<- 2*(effects*all_freq_change)
  }

  Rhat_data <- tibble(QTL=c(1:(n_qtl*n_chr)),ASE = effects, SE_ASE = effects_se, Rhat_per_locus=per_locus_selection_response)

  if(spread==F){
    Rhat_data <- Rhat_data %>%
      gather(key = "Measure",value = "Estimate",-QTL)
  } else {}

  if(plot==T){
    #   #Plot R hat
    #   tmp <- qtl_effect_data %>% filter(Measure=="Rhat_per_locus")
    #   plot.Rhat <- ggplot(tmp,aes(x=QTL,y=Estimate)) +
    #     geom_point(size=3) +
    #     guides(alpha = FALSE) +
    #     theme(
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       panel.background = element_blank(),
    #       axis.line = element_line(colour = "black"),
    #       axis.title = element_text(size = rel(1.5)),
    #       axis.text = element_text(size = rel(1.25)),
    #     ) +
    #     scale_x_continuous(breaks = seq(1,(n_qtl*n_chr),1)) +
    #     scale_y_continuous(TeX("$\\hat{\\R}$")) +
    #     ggtitle("Per Locus Response to Selection") +
    #     theme_Publication()
    #
    #   #Plot ASE
    #   tmp <- qtl_effect_data %>% filter(Measure=="ASE") %>%
    #     mutate(SE = effects_se)
    #   plot.ASE <- ggplot(tmp,aes(x=QTL,y=Estimate)) +
    #     geom_point(size=3) +
    #     geom_errorbar(aes(ymin = Estimate - SE,
    #                       ymax = Estimate + SE)) +
    #     guides(alpha = FALSE) +
    #     theme(
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       panel.background = element_blank(),
    #       axis.line = element_line(colour = "black"),
    #       axis.title = element_text(size = rel(1.5)),
    #       axis.text = element_text(size = rel(1.25)),
    #     ) +
    #     scale_x_continuous(breaks = seq(1,(n_qtl*n_chr),1)) +
    #     scale_y_continuous(TeX("$\\hat{\\ASE}$")) +
    #     ggtitle("Allele Substitution Effects") +
    #     theme_Publication()
    #   rm(tmp)
    #   plot.ASE
    #   plot.Rhat
  }

  return(Rhat_data)
}

### Estimate selection coefficients (s hat)
# INPUT:
# geno = matrix of genotypes for selected parents of next generation
# spread = option to output tibble with 1 row per locus, with a column for each statistic
# OUTPUT:
# Tibble with frequencies, fitnesses & selection coefficients for genotypes/alleles at each locus
s_hat_function <- function(geno,n_qtl,n_chr,spread=TRUE){

  per_locus_genotype_counts <- apply(geno,2,function(x) list(table(x)))

  unname(unlist(per_locus_genotype_counts[[1]]))

  P <- unname(unlist(map(per_locus_genotype_counts,function(x) ((x[[1]]["2"])/sum(unlist(x)))))); P[is.na(P)] <- 0
  H <- unname(unlist(map(per_locus_genotype_counts,function(x) ((x[[1]]["1"])/sum(unlist(x)))))); H[is.na(H)] <- 0
  Q <- unname(unlist(map(per_locus_genotype_counts,function(x) ((x[[1]]["0"])/sum(unlist(x)))))); Q[is.na(Q)] <- 0

  p <- P + (H/2);
  q <- Q + (H/2);

  #Genotype Measures
  genotype_fitness <- Q/P
  rev_ref_genotype<-which(genotype_fitness>1)
  genotype_fitness[rev_ref_genotype] <- P[rev_ref_genotype]/Q[rev_ref_genotype]
  genotype_sel_coefficient <- 1 - genotype_fitness

  gametic_contributions_P <- p^2
  gametic_contributions_H <- 2*p*q
  gametic_contributions_Q <- (q^2)*genotype_fitness

  #Allele Measures
  allele_fitness <- q/p
  rev_ref_allele<-which(allele_fitness>1)
  allele_fitness[rev_ref_allele] <- p[rev_ref_allele]/q[rev_ref_allele]
  allele_sel_coefficient <- 1 - allele_fitness

  #assign selection coefficient of 1 when opposzing homozygotes dont exist
  complete_sel <- which(genotype_sel_coefficient==0)
  genotype_sel_coefficient[complete_sel] <- 1 ;# genotype_fitness[complete_sel] <- 0
  allele_sel_coefficient[complete_sel] <- 1 ;# allele_fitness[complete_sel] <- 0

  s_hat_data <- tibble(QTL=c(1:(n_qtl*n_chr)),P, H, Q, genotype_fitness, genotype_sel_coefficient, gametic_contributions_P, gametic_contributions_H,gametic_contributions_Q, p, q, allele_fitness,allele_sel_coefficient)

  if(spread==F){
    s_hat_data <- s_hat_data %>%
      gather(key = "Measure",value = "Estimate",-QTL)
  } else {}
  return(s_hat_data)
}
### Custom ggplot theme ----
#Paper Scaling
theme_Publication <- function(base_size = 20,
                              base_family = "serif") {
  library(grid)
  library(ggthemes)
  (
    theme_foundation(base_size = base_size, base_family = base_family)
    + theme(plot.title = element_text(face = "bold",size = rel(1),margin = margin(b = 20,t = 20,r = 20, l = 20),hjust=0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(0.75)),
            axis.title.y = element_text(angle = 90),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = rel(0.75)),
            axis.line = element_line(colour = NA, size = rel(0.75)),
            axis.line.x = element_line(colour = 'black'),
            axis.line.y = element_line(colour = 'black'),
            axis.ticks = element_line(),
            panel.grid.major.y = element_line(colour = "#f0f0f0"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.text = element_text(size = rel(0.75)),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(1, "cm"),
            #legend.margin = unit(margin(c(0, 0, 0, 0)), "mm"),
            legend.title = element_text(face = "bold", size = rel(0.75)),
            plot.margin = unit(c(1, 2, 3, 1), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold", size = rel(0.75)),
            legend.background=element_blank(),
            #legend.spacing.y = unit(0.5, "mm")
    )
  )
}
