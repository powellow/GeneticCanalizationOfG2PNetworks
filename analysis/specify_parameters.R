
### Specify Population Parameters

n_founders=2 #Change this value to show effect of sampling on covariance/LD between QTL
n_traits = 4
n_qtn = 10 #Number of QTN per Chromosome
n_snp = 0
n_chr= 4 #Number of Chromosomes
ploidy = 2
n_qtn_per_node = (n_qtn*n_chr)/n_traits
n_seg_sites = (n_qtn+n_snp)*n_chr
start_allele_freq = 0.5 #Starting (Alternate) Allele Frequency (ie. 0 Allele )

n.off.f1   = 1
n.off.f2   = 1000
n.selected = 100


### Selection Simulation
scenario = paste0("ConstrainedSampledEffects_",n_qtn_per_node,"QtlPerNode")
n.cycles = 30
nreps = 100
cycles.of.interest = seq(1,n.cycles)

## Trait Parameters
target.trait = "Time to Bud Outgrowth"
genetic_traits = c("Auxin.AdditiveGeneticValue","Cytokinin.AdditiveGeneticValue",
                   "Sucrose.AdditiveGeneticValue","Strigolactone.AdditiveGeneticValue",
                   "Normalised.Auxin.TotalGeneticValue","Normalised.Cytokinin.TotalGeneticValue",
                   "Normalised.Sucrose.TotalGeneticValue","Normalised.Strigolactone.TotalGeneticValue",
                   "Normalised.Auxin.AdditiveGeneticValue","Normalised.Cytokinin.AdditiveGeneticValue",
                   "Normalised.Sucrose.AdditiveGeneticValue","Normalised.Strigolactone.AdditiveGeneticValue",
                   "Normalised.IntegratorSignal.TotalGeneticValue",
                   "TimeToBudOutgrowth.TotalGeneticValue","Normalised.TimeToBudOutgrowth.TotalGeneticValue",target.trait)
phenotypic_traits = c("Integrator Signal","Bud Outgrowth","Alpha","Sucrose")
selection.or.drift = "Selection"
selection.direction = "Min" # "Min" OR "Max"
threshold = "No"
#heritabilities = c(seq(0.5,1.0,0.1))
#heritabilities = c(seq(0.1,0.9,0.1),seq(0.91,1.0,0.01))
#heritabilities = c(seq(0.1,1.0,0.1))
#heritabilities = seq(0.91,0.99,0.01)
heritabilities = c(0.3,1)
heritability_constant = FALSE

## Environment Parameters
Nutrients = "No"
