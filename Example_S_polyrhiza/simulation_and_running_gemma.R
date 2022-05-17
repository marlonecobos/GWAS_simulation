# ------------------------------------------------------------------------------
# Project: Genome wide association of thermal upper limit in Spirodela polyrhiza  
# Authors: Marlon E. Cobos, John Kelly
# Date: 02-04-2022 (dd/mm/yyyy)
# Script title: Simulation of data and running GEMMA
# ------------------------------------------------------------------------------

# Description
# This script simulates phenotypic and genomic data to run GWAs analyses in 
# GEMMA.
# - A total of 50 individuals are considered in this study
# - Phenotypes and genomes are simulated randomly first
# - Phenotypes derive from a normal distribution with mean = 0 and SD = 1
# - Genome is simulated to get 1000 SNPs and each SNP is simulated to have
#   values of 0, 1, or 2. The SNPs are then evenly distributed in 20 chromosomes
# - The first scenario is one of completely random association (no association) 
# - Then phenotypes are changed so one SNP from the genome influence the trend
#   in the phenotypic values, this is done using the data from the initial 
#   simulation. This is the second scenario
# - After that, using the initial simulated data, phenotypic values are changed
#   so five different SNPs influence additively the trend. This is the third 
#   scenario
# - The fourth and fifth scenarios simulate associations for ten and twenty SNPs
#   respectively
# - For each scenario, 100 replicates are produced 
#
# After all the data are prepared, GWAs analyses are run using GEMMA. GWAs are
# run using Linear Mixed Models (likelihood ratio test). The relatedness matrix 
# used to run the GWAs is a centered matrix.
#
# Everything will work from R using the function "gwas_simulation" found in the
# file "functions.R". Internally we are using "system" to interact with
# the terminal and GEMMA directly. Another way to do this is to prepare a 
# batch running and use bash; however, the current way prevent us from writing 
# all the data needed for all simulations (saves some space).


# load function to run analyses ------------------------------------------------
source("Scripts/functions.R")
# ------------------------------------------------------------------------------


# initial parameters -----------------------------------------------------------
# number of simulations
nsim <- 100

# number of individuals
nind <- 50

# number of SNPs
nsnp <- 1000

# genome simulation (totally random association)
## genotype options
genos <- c(0, 1, 2)
# ------------------------------------------------------------------------------


# phenotype and genome simulation ----------------------------------------------
## random association
genos_list <- pheno_geno_simulation(n_individuals = nind, n_snps = nsnp, 
                                   n_simulations = nsim, pheno_mean = 0, 
                                   pheno_sd = 1, geno_options = genos, 
                                   structure = FALSE)

head(genos_list$simulation_1)[, 1:10]



## Phenotype association with one SNP 
### SNP to be associated
snp_one <- 900 + 1 # + 1, because phenotype is in the first column

### changing phenotype according to SNP 900
genos_list_one <- pheno_geno_asociation(pheno_geno_list = genos_list, 
                                        snp_reference = snp_one)

head(genos_list_one$simulation_1)[, 1:10]



## Phenotype association with five SNPs 
### SNPs to be associated
snp_five <- c(9, 90, 199, 499, 900)  + 1 # + 1, as phenotype is in first column

### changing phenotype according to SNPs
genos_list_five <- pheno_geno_asociation(pheno_geno_list = genos_list, 
                                         snp_reference = snp_five)

head(genos_list_five$simulation_1)[, 1:10]


## Phenotype association with ten SNPs 
### SNP to be associated
snp_ten <- c(9, 90, 199, 299, 399, 499, 600, 799, 900, 999)  + 1 # same

### changing phenotype according to SNPs
genos_list_ten <- pheno_geno_asociation(pheno_geno_list = genos_list, 
                                        snp_reference = snp_ten)

head(genos_list_ten$simulation_1)[, 1:10]


## Phenotype association with twenty SNPs 
### SNP to be associated
snp_twenty <- round(seq(9, 995, length.out = 20)) + 1 # same

### changing phenotype according to SNPs
genos_list_twenty <- pheno_geno_asociation(pheno_geno_list = genos_list, 
                                           snp_reference = snp_twenty)

head(genos_list_twenty$simulation_1)[, 1:10]
# ------------------------------------------------------------------------------



# SNP ids (e.g., snp-1-5151352; snp - chr number - SNP position) ---------------
## chromosomes
nchr <- 20

## let's assign an equal number of SNPs to each chromosome
npc <- nsnp / nchr

chr_rep <- rep(1:nchr, each = npc)

snp_id <- paste0("snp-", chr_rep, "-", 1:nsnp)


# bases for each SNIP (GCTA; e.g., G U)
bases <- sapply(1:nsnp, function(x) {
  paste0(sample(c("A", "T", "G", "C"), 2), collapse = " ")
})


# SNP id and bases
snp_base <- paste(snp_id, bases)


# geno map file to run gemma
map_table <- data.frame(snp_id, 1:nsnp, chr_rep)
# ------------------------------------------------------------------------------



# variance from phenotype explained by selected SNPs ---------------------------
## getting variances from E (original phenotype), P (adjusted phenotypes) and
## G (SNPs contributing additively to the variance in P)
## calculation: h2 = VGi / VP; VP = VE + sum(VG)

## proportion of variance explained by SNPs in the examples created
### E = original phenotype in simulations
o_phenotypes <- lapply(genos_list, function(x) {x[, 1]})

### calculations
h2_1 <- calc_h2(orig_phenotype_E = o_phenotypes, 
                pheno_geno_list = genos_list_one, snp_reference = snp_one)

h2_5 <- calc_h2(orig_phenotype_E = o_phenotypes, 
                pheno_geno_list = genos_list_five, snp_reference = snp_five)

h2_10 <- calc_h2(orig_phenotype_E = o_phenotypes, 
                 pheno_geno_list = genos_list_ten, snp_reference = snp_ten)

h2_20 <- calc_h2(orig_phenotype_E = o_phenotypes, 
                 pheno_geno_list = genos_list_twenty, snp_reference = snp_twenty)
# ------------------------------------------------------------------------------


# saving data for later --------------------------------------------------------
save(genos_list, genos_list_one, genos_list_five, genos_list_ten, 
     genos_list_twenty, snp_base, map_table, 
     file = "Data/Simulated.RData")

save(snp_one, snp_five, snp_ten, snp_twenty, 
     file = "Data/SNPs_simulated.RData")

save(h2_1, h2_5, h2_10, h2_20, file = "Data/Prop_variance.RData")
# ------------------------------------------------------------------------------



# running things in GEMMA ------------------------------------------------------
## folders to store data and results from GEMMA
odir <- "gemma_random"
odir_one <- "gemma_one"
odir_five <- "gemma_five"
odir_ten <- "gemma_ten"
odir_twenty <- "gemma_twenty"

dir.create(odir)
dir.create(odir_one)
dir.create(odir_five)
dir.create(odir_ten)
dir.create(odir_twenty)

### complete path for output directories
odir <- normalizePath(odir)
odir_one <- normalizePath(odir_one)
odir_five <- normalizePath(odir_five)
odir_ten <- normalizePath(odir_ten)
odir_twenty <- normalizePath(odir_twenty)

## write map table (same for all scenarios; only written in random directory)
map_file <- paste0(odir, "/map.txt")

write.table(map_table, file = map_file, sep = " ", col.names = FALSE,
            row.names = FALSE, quote = FALSE)

## phenotype files to be written
pheno_file <- "pheno.txt"

## genotype files to be written
geno_file <- "geno.txt"


## running analyses in loops
### base name for relatedness matrix
out_rm <- "dw_rm"

### base name for GWAs result table
out_gwas <- "dw_gwas"


### simulation with random associations
gwas_simulation(pheno_geno_list = genos_list, map_file = map_file, 
                snp_reference = snp_base, phenotype_file = pheno_file, 
                genotype_file = geno_file, rel_matrix_name = out_rm, 
                gwas_results_name = out_gwas, output_directory = odir)


### simulation for associations with one SNPs
gwas_simulation(pheno_geno_list = genos_list_one, map_file = map_file, 
                snp_reference = snp_base, phenotype_file = pheno_file, 
                genotype_file = geno_file, rel_matrix_name = out_rm, 
                gwas_results_name = out_gwas, output_directory = odir_one)


### simulation for associations with five SNPs
gwas_simulation(pheno_geno_list = genos_list_five, map_file = map_file, 
                snp_reference = snp_base, phenotype_file = pheno_file, 
                genotype_file = geno_file, rel_matrix_name = out_rm, 
                gwas_results_name = out_gwas, output_directory = odir_five)


### simulation for associations with ten SNPs
gwas_simulation(pheno_geno_list = genos_list_ten, map_file = map_file, 
                snp_reference = snp_base, phenotype_file = pheno_file, 
                genotype_file = geno_file, rel_matrix_name = out_rm, 
                gwas_results_name = out_gwas, output_directory = odir_ten)


### simulation for associations with twenty SNPs
gwas_simulation(pheno_geno_list = genos_list_twenty, map_file = map_file, 
                snp_reference = snp_base, phenotype_file = pheno_file, 
                genotype_file = geno_file, rel_matrix_name = out_rm, 
                gwas_results_name = out_gwas, output_directory = odir_twenty)
