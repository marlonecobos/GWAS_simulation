pheno_geno_simulation <- function(n_individuals = 100, n_snps = 1000, 
                                  n_simulations = 100, pheno_mean = 0, 
                                  pheno_sd = 1, geno_options = c(0, 1, 2), 
                                  structure = FALSE, 
                                  lowhomo_geno_prob = c(0.6, 0.2, 0.2), 
                                  hete_geno_prob = c(0.2, 0.6, 0.2),
                                  highhomo_geno_prob = c(0.2, 0.2, 0.6)) {
  
  ## number of individuals with low, intermediate, and high traits
  lind <- round(n_individuals / 3)
  iind <- lind
  hind <- n_individuals - lind - iind
  
  ## simulation
  genos_list <- lapply(1:n_simulations, function(x) {
    ### phenotype vector
    set.seed((1 + x))
    pheno <- round(rnorm(n_individuals, mean = pheno_mean, sd = pheno_sd), 2) 
    
    ### genotype table (n_individuals x n_snps)
    if (structure == FALSE) {
      ### genotype table (nind x nsnp)
      set.seed((n_simulations + x))
      genos_table <- matrix(sample(geno_options, size = (n_individuals * n_snps), 
                                   replace = TRUE), ncol = n_snps)
      
    } else {
      genos_table <- sapply(1:n_snps, function(y) {
        set.seed((n_simulations + n_snps + x + y))
        
        snp1 <- sample(geno_options, size = 1)
        
        if (snp1 == 0) {
          snps <- c(sample(geno_options, size = lind, prob = lowhomo_geno_prob, 
                           replace = TRUE),
                    sample(geno_options, size = iind, prob = hete_geno_prob, 
                           replace = TRUE),
                    sample(geno_options, size = hind, prob = highhomo_geno_prob, 
                           replace = TRUE))
        } else {
          snps <- c(sample(geno_options, size = hind, prob = highhomo_geno_prob, 
                           replace = TRUE),
                    sample(geno_options, size = iind, prob = hete_geno_prob, 
                           replace = TRUE),
                    sample(geno_options, size = lind, prob = lowhomo_geno_prob, 
                           replace = TRUE))
        }
      })
    }
    
    ### combine phenotype and genotype
    cbind(Phenotype = pheno, genos_table)
  })
  
  ## naming simulation list
  names(genos_list) <- paste0("simulation_", 1:n_simulations)
  
  return(genos_list)
}



# generating association of phenotypes with snps
pheno_geno_asociation <- function(pheno_geno_list, snp_reference) {
  ## generating association additively
  genos_list_asso <- lapply(1:length(pheno_geno_list), function(x) {
    tab <- pheno_geno_list[[x]]
    tab[, 1] <- apply(tab[, c(1, snp_reference)], 1, sum)
    tab
  })
  
  ## naming simulation list
  names(genos_list_asso) <- names(pheno_geno_list)
  
  return(genos_list_asso)
}





# calculation of proportion of variances explained by SNPs in different cases
calc_h2 <- function(orig_phenotype_E, pheno_geno_list, snp_reference) {
  
  h2 <- sapply(1:length(pheno_geno_list), function(x) {
    ## original phenotype to calculate VE
    ve <- var(orig_phenotype_E[[x]])
    
    ## VP and VG calculation
    vp <- var(pheno_geno_list[[x]][, 1])
    vg <- apply(as.matrix(pheno_geno_list[[x]][, snp_reference]), 2, var)
    
    ## proportion of variance explained by SNPs involved
    vg / (ve + sum(vg))
  })
  
  ## organizing results
  if (class(h2)[1] == "matrix") {
    h2 <- t(h2)
  } else {
    h2 <- as.matrix(h2)
  }
  
  colnames(h2) <- snp_reference - 1
  
  return(h2)
}




# exploring relatedness among groups of individuals
relatedness <- function(rel_matrix, output_file) {
  ## read matrix
  rmat <- read.table(rel_matrix, sep = "\t")
  
  ## excluding diagonal values
  rmm <- as.matrix(rmat)
  diag(rmm) <- NA
  
  ## preparing groups
  cuts <- floor(seq(1, ncol(rmat), length.out = 4))
  
  g1 <- cuts[1]:cuts[2]
  g2 <- (cuts[2]+1):cuts[3]
  g3 <- (cuts[3]+1):cuts[4]
  
  ## preparing results
  df <- data.frame(
    G1_G1 = mean(unlist(rmat[g1, g1]), na.rm = TRUE),
    G1_G2 = mean(unlist(rmat[g1, g2]), na.rm = TRUE),
    G1_G3 = mean(unlist(rmat[g1, g3]), na.rm = TRUE),
    G2_G2 = mean(unlist(rmat[g2, g2]), na.rm = TRUE),
    G2_G3 = mean(unlist(rmat[g2, g3]), na.rm = TRUE),
    G3_G3 = mean(unlist(rmat[g3, g3]), na.rm = TRUE)
  )
  
  ## writing results
  if (!file.exists(output_file)) {
    write.table(df, output_file, row.names = FALSE, sep = "\t")  
  } else {
    write.table(unname(df), output_file, row.names = FALSE, append = TRUE, 
                sep = "\t")
  }
}




# preparing files and running analyses in GEMMA
gwas_simulation <- function(pheno_geno_list, map_file, snp_reference, 
                            phenotype_file, genotype_file, rel_matrix_name, 
                            gwas_results_name, output_directory) {
  
  pheno <- paste0(output_directory, "/", phenotype_file)
  geno <- paste0(output_directory, "/", genotype_file)
  
  nsim <- length(pheno_geno_list)
  
  for (i in 1:nsim) {
    #### writing phenotype file
    write.table(pheno_geno_list[[i]][, 1], file = pheno, sep = " ", 
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    #### writing genotype file
    gf <- data.frame(snp_reference, t(pheno_geno_list[[i]][, -1]))
    
    write.table(gf, file = geno, sep = " ", col.names = FALSE, 
                row.names = FALSE, quote = FALSE)
    
    #### compute realized relatedness matrix from the genotype files
    Sys.sleep(0.1)
    gemma_rmat <- paste("gemma -p", pheno, "-g", geno, "-a", map_file,
                        "-gk 1 -outdir", output_directory, "-o", 
                        rel_matrix_name)
    
    system(gemma_rmat)
    
    #### calculating relatedness stats
    rel_mat <- paste0(output_directory, "/", rel_matrix_name, ".cXX.txt")
    
    relatedness(rel_matrix = rel_mat, 
                output_file = paste0(output_directory, "/relatedness.txt"))
    
    #### LMM-based association statistics for each SNP
    Sys.sleep(0.1)
    out_gwas_r <- paste0(out_gwas, "_r", i)
    
    gemma_gwas <- paste("gemma -p", pheno, "-a", map_file, "-g", geno, 
                        "-notsnp", "-k", rel_mat, "-lmm 2 -outdir", 
                        output_directory, "-o", out_gwas_r)
    
    system(gemma_gwas)
    
    message("\nSimulation step ", i, " of ", nsim, "\n\n")
  }
}



# reading and post-processing results from GWAS analyses
gwas_result_list <- function(gwas_output_files, p_column = 10) {
  # number of simulations considered
  nsim <- length(gwas_output_files)
  
  ### reading the first table of results
  gwsall <- read.table(gwas_output_files[1], header = TRUE, sep = "\t")
  
  for (i in 2:nsim) {
    ### reading other results
    gwscan <- read.table(gwas_output_files[i], header = TRUE, sep = "\t")
    
    ### adding further results to initial table
    gwsall <- cbind(gwsall, gwscan[, p_column])
  }
  
  ## naming columns more appropriately
  colnames(gwsall)[p_column:(p_column + nsim - 1)] <- paste0("Sim_", 1:nsim)
  
  return(gwsall)
}



# identification of significant results
significant_gwas_results <- function(gwas_result_list, p_column = 10) {
  # number of simulations considered
  nsim <- length(gwas_result_list) - p_column - 1
  
  ## detection after FDR correction
  for (i in 1:nsim) {
    ### FDR correction of p-values
    pfdr <- p.adjust(gwas_result_list[, (p_column - 1 + i)], method = "fdr") 
    
    ### checking significant ones
    gwas_result_list[, (p_column - 1 + i)] <- (pfdr <= 0.05) * 1
  }
  
  return(gwas_result_list)
}



# error whiskers
error_whiskers <- function(x, upper, lower, length = 0.1,...){
  if(length(x) !=length(lower) | length(lower) != length(upper))
    stop("vectors must have the same length")
  arrows(x, upper, x, lower, angle = 90, code = 3, length = length, ...)
}
