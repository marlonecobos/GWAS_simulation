# ------------------------------------------------------------------------------
# Project: Genome wide association of thermal upper limit in Spirodela polyrhiza  
# Authors: Marlon E. Cobos, John Kelly
# Date: 23-04-2022 (dd/mm/yyyy)
# Script title: Summary of results from GEMMA on simulated data (structure)
# ------------------------------------------------------------------------------

# Description ------------------------------------------------------------------
# This script helps to prepare a summary from GWAs simulations run under five
# scenarios: 1) random association; 2-5) 1, 5, 10, and 20 SNPs associated with 
# the phenotype (additively).
# ------------------------------------------------------------------------------


# load function to run analyses ------------------------------------------------
source("Scripts/functions.R")
# ------------------------------------------------------------------------------


# loading data -----------------------------------------------------------------
load("Data/SNPs_simulated_structure.RData")
load("Data/Prop_variance_structure.RData")
# ------------------------------------------------------------------------------


# summarizing proportion of variances ------------------------------------------
## vector with means
h2_means <- c(One = unname(apply(h2_1, 2, mean)),
              Five = mean(apply(h2_5, 2, mean)),
              Ten = mean(apply(h2_10, 2, mean)),
              Twenty = mean(apply(h2_20, 2, mean)))

## table with confidence intervals
h2_ci <- as.matrix(data.frame(
  One = unname(apply(h2_1, 2, quantile, prob = c(0.025, 0.975))),
  Five = apply(apply(h2_5, 2, quantile, prob = c(0.025, 0.975)), 1, mean),
  Ten = apply(apply(h2_10, 2, quantile, prob = c(0.025, 0.975)), 1, mean),
  Twenty = apply(apply(h2_20, 2, quantile, prob = c(0.025, 0.975)), 1, mean)
))
# ------------------------------------------------------------------------------


# summarizing genetic relatedness ----------------------------------------------
## reading matrix
relatedness <- read.table(file = "gemma_random_str/relatedness.txt", 
                          header = TRUE, sep = "\t")

colnames(relatedness) <- gsub("_", " & ", colnames(relatedness))
colnames(relatedness)[c(1, 4, 6)] <- gsub(" & .*", "", 
                                          colnames(relatedness)[c(1, 4, 6)])


## getting the absolute values
relatedness[] <- abs(relatedness[])

## mean values for groups
relatedness_mean <- apply(relatedness, 2, mean)

## confidence intervals for groups
relatedness_ci <- apply(relatedness, 2, quantile, prob = c(0.025, 0.975))
# ------------------------------------------------------------------------------


# post-processing results from gwas analyses -----------------------------------
# random associations
## listing results from gwas analyses 
gwscan_files <- list.files(path = "gemma_random_str", pattern = "assoc.txt$", 
                           full.names = TRUE)

## reading results from gwas analyses 
gwsall <- gwas_result_list(gwas_output_files = gwscan_files, p_column = 10)

## detecting significant p-values
gwsall_01 <- significant_gwas_results(gwas_result_list = gwsall, p_column = 10)


# results from associations with one SNP
gwscan_files <- list.files(path = "gemma_one_str", pattern = "assoc.txt$", 
                           full.names = TRUE)

## reading results from gwas analyses 
gwsall_one <- gwas_result_list(gwas_output_files = gwscan_files, p_column = 10)

## detecting significant p-values
gwsall_one_01 <- significant_gwas_results(gwas_result_list = gwsall_one, 
                                          p_column = 10)


# results from associations with five SNPs
gwscan_files <- list.files(path = "gemma_five_str", pattern = "assoc.txt$", 
                           full.names = TRUE)

## reading results from gwas analyses 
gwsall_five <- gwas_result_list(gwas_output_files = gwscan_files, p_column = 10)

## detecting significant p-values
gwsall_five_01 <- significant_gwas_results(gwas_result_list = gwsall_five, 
                                           p_column = 10)


# results from associations with ten SNPs
gwscan_files <- list.files(path = "gemma_ten_str", pattern = "assoc.txt$", 
                           full.names = TRUE)

## reading results from gwas analyses 
gwsall_ten <- gwas_result_list(gwas_output_files = gwscan_files, p_column = 10)

## detecting significant p-values
gwsall_ten_01 <- significant_gwas_results(gwas_result_list = gwsall_ten, 
                                          p_column = 10)


# results from associations with twenty SNPs
gwscan_files <- list.files(path = "gemma_twenty_str", pattern = "assoc.txt$", 
                           full.names = TRUE)

## reading results from gwas analyses 
gwsall_twenty <- gwas_result_list(gwas_output_files = gwscan_files, p_column = 10)

## detecting significant p-values
gwsall_twenty_01 <- significant_gwas_results(gwas_result_list = gwsall_twenty, 
                                             p_column = 10)
# ------------------------------------------------------------------------------



# preparing summary of gwas results --------------------------------------------
## frequency of significance
stats <- apply(gwsall_01[, 10:109], 1, sum)
stats_one <- apply(gwsall_one_01[, 10:109], 1, sum)
stats_five <- apply(gwsall_five_01[, 10:109], 1, sum)
stats_ten <- apply(gwsall_ten_01[, 10:109], 1, sum)
stats_twenty <- apply(gwsall_twenty_01[, 10:109], 1, sum)

## table summarizing stats
summary_gwas <- data.frame(Reference = gwsall_01$rs, 
                           Frequency_random = stats, 
                           Frequency_one = stats_one, 
                           Frequency_five = stats_five, 
                           Frequency_ten = stats_ten,
                           Frequency_twenty = stats_twenty)

## getting only SNPs simulated to be associated and the ones detected 
initial <- unique(c(snp_one, snp_five, snp_ten, snp_twenty)) - 1

significant <- which(summary_gwas$Frequency_random > 0 | 
                       summary_gwas$Frequency_one > 0 |
                       summary_gwas$Frequency_five > 0 |
                       summary_gwas$Frequency_ten > 0 |
                       summary_gwas$Frequency_twenty > 0)

of_interest <- sort(unique(c(initial, significant)))

summary_gwas <- summary_gwas[of_interest, ]

## writing table
write.csv(summary_gwas, "Summary_gwas_results_structure.csv", row.names = FALSE)
# ------------------------------------------------------------------------------



# plotting results -------------------------------------------------------------
## proportion of variance explained by snps
maxh2 <- max(h2_ci[2, ])
ylim <- c(0, (maxh2 + (maxh2 * 0.1)))

png("Figures/Summary_h2_structure.png", width = 80, height = 60, units = "mm", 
    res = 600)

par(mar = c(4.5, 4.5, 1.5, 1), cex = 0.5)

bar_h2 <- barplot(h2_means, space = 0.7, ylim = ylim, xlab = "Number of SNPs", 
                  ylab = expression("VP explained by VG"[i]*" (h"^2*")"), 
                  las = 1)

error_whiskers(x = bar_h2, upper = h2_ci[2, ], lower = h2_ci[1, ])

legend("topright", legend = "Mean and 95% CI", bty = "n")

box(bty = "l", col = "black", lwd = 1.5)

dev.off()



## relatedness among groups
maxrel <- max(relatedness_ci[2, ])
ylim1 <- c(0, (maxrel + (maxrel * 0.1)))

png("Figures/Summary_relatedness_structure.png", width = 80, height = 60, units = "mm",
    res = 600)

par(mar = c(4.5, 4.7, 1.5, 1), cex = 0.5)

bar_rel <- barplot(relatedness_mean, space = 0.7, ylim = ylim1, xlab = "Groups", 
                   ylab = "", las = 1)
title(ylab = "Relatedness", line = 3.5)

error_whiskers(x = bar_rel, upper = relatedness_ci[2, ], 
               lower = relatedness_ci[1, ])

legend("top", legend = "Mean and 95% CI", bty = "n", inset = 0)

box(bty = "l", col = "black", lwd = 2)

dev.off()



## gwas results
## relative position of SNPs
place <- gwsall_01$ps# / 1000

## color corresponding to chromosomes
colf <- as.factor(gwsall$chr %% 2)
cols <- c("#366B9A", "#71B2EB")

## ylim
ylimgwas <- c(0, max(summary_gwas[, -1]))

## point type
pt <- 16

## exporting the figures
### random
png("Figures/Summary_gwas_random_structure.png", width = 166, height = 70, 
    units = "mm", res = 600)
par(mar = c(4, 4.5, 0.5, 0.5),  cex = 0.7)

plot(place, stats, col = cols[colf], pch = pt, las = 1, ylim = ylimgwas,
     xlab = "SNP position", ylab = "Times detected from 100")

dev.off()

### one SNP
png("Figures/Summary_gwas_one_structure.png", width = 166, height = 70, 
    units = "mm", res = 600)
par(mar = c(4, 4.5, 0.5, 0.5),  cex = 0.7)

plot(place, stats_one, col = cols[colf], las = 1, type = "n", ylim = ylimgwas, 
     xlab = "SNP position", ylab = "Times detected from 100")
abline(v = snp_one, col = "gray75")
points(place, stats_one, col = cols[colf], pch = pt)

dev.off()

### five SNPs
png("Figures/Summary_gwas_five_structure.png", width = 166, height = 70, 
    units = "mm", res = 600)
par(mar = c(4, 4.5, 0.5, 0.5),  cex = 0.7)

plot(place, stats_five, col = cols[colf], las = 1, type = "n", ylim = ylimgwas, 
     xlab = "SNP position", ylab = "Times detected from 100")
abline(v = snp_five, col = "gray75")
points(place, stats_five, col = cols[colf], pch = pt)

dev.off()

### ten SNPs
png("Figures/Summary_gwas_ten_structure.png", width = 166, height = 70, 
    units = "mm", res = 600)
par(mar = c(4, 4.5, 0.5, 0.5),  cex = 0.7)

plot(place, stats_ten, col = cols[colf], las = 1, type = "n", ylim = ylimgwas, 
     xlab = "SNP position", ylab = "Times detected from 100")
abline(v = snp_ten, col = "gray75")
points(place, stats_ten, col = cols[colf], pch = pt)

dev.off()


### twenty SNPs
png("Figures/Summary_gwas_twenty_structure.png", width = 166, height = 70, 
    units = "mm", res = 600)
par(mar = c(4, 4.5, 0.5, 0.5),  cex = 0.7)

plot(place, stats_twenty, col = cols[colf], las = 1, type = "n", ylim = ylimgwas, 
     xlab = "SNP position", ylab = "Times detected from 100")
abline(v = snp_twenty, col = "gray75")
points(place, stats_twenty, col = cols[colf], pch = pt)

dev.off()
# ------------------------------------------------------------------------------



# saving results for later -----------------------------------------------------
save(h2_means, h2_ci, relatedness_mean, relatedness_ci,
     gwsall, gwsall_01, gwsall_one, gwsall_one_01, gwsall_five, gwsall_five_01,
     gwsall_ten, gwsall_ten_01, gwsall_twenty, gwsall_twenty_01,
     stats, stats_one, stats_five, stats_ten, stats_twenty, summary_gwas, 
     file = "Summary_all_results_structure.RData")
# ------------------------------------------------------------------------------
