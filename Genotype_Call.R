# Title: Genotype_Call
# Authors: Adam Staadig et al.
# License: MIT License (OSI-approved, see LICENSE file in this repository)
# GitHub repository: https://github.com/astaadig/Genotype_Call

# DISCLAIMER
#This is an in-house R script for genotype calling. 
#This code was used to analyze the data described in Comparative assessment of SNP genotyping assays for challenging forensic samples utilizing ancient DNA methods (Staadig et. al).
# -------------------------------------------------------------------------  

Genotype_Call <- function(het_bal1_cut, het_bal2_cut, het_bal3_cut,
                                  hom_bal_cut, Q_cut,
                                  cov_hom_cut, cov_het_cut, cov_max,
                                  F_R_cut,
                                  indata, outfile) {
  # Read data
  raw_data <- read.table(indata, fill = TRUE, dec = ".", header = TRUE, sep = "\t")
  raw_data[is.na(raw_data)] <- 0

  # Number of SNPs (4 rows per SNP)
  storlek <- dim(raw_data) / 4

  # Output table: Chromosome, Region, rsID, Genotype
  tabell1 <- matrix(0, storlek[1], 4)

  # Loop over SNPs
  for (u in 1:storlek[1]) {

    # Extract 4 rows x 7 columns for this SNP
    temp_data <- data.frame(nrow = 4, ncol = 7)
    for (p in 1:4) {
      for (h in 1:7) {
        temp_data[p, h] <- raw_data[4 * (u - 1) + p, h]
      }
    }

    # Sort alleles by coverage (column 4)
    temp_data <- temp_data[order(-temp_data[, 4]), ]
    status_temp <- 0  # 0=no call, 1=homo, 2=hetero

    # Allele balance (top vs top+second)
    if (temp_data[1, 4] == 0) {
      temp_bal <- 0
    } else {
      temp_bal <- temp_data[1, 4] / (temp_data[1, 4] + temp_data[2, 4])
    }

    # Third-allele check (background)
    if (temp_data[1, 4] == 0) {
      temp_third <- 0
    } else {
      temp_third <- temp_data[3, 4] / temp_data[1, 4]
    }

    # Q-scores, F/R ratios, coverages for top 2 alleles
    tempQ1  <- as.numeric(temp_data[1, 6])
    tempQ2  <- as.numeric(temp_data[2, 6])
    tempF_R1 <- as.numeric(temp_data[1, 5])
    tempF_R2 <- as.numeric(temp_data[2, 5])
    tempC1  <- as.numeric(temp_data[1, 4])
    tempC2  <- as.numeric(temp_data[2, 4])

    # Decide heterozygote vs homozygote from balance
    if (temp_bal >= het_bal1_cut && temp_bal <= het_bal2_cut) {
      status_temp <- 2   # hetero
    } else if (temp_bal >= hom_bal_cut) {
      status_temp <- 1   # homo
    }

    # QC for heterozygotes
    if (status_temp == 2) {
      if (temp_third > het_bal3_cut ||
          tempQ1 < Q_cut || tempQ2 < Q_cut ||
          tempC1 < cov_het_cut || tempC2 < cov_het_cut ||
          tempF_R1 < F_R_cut || tempF_R2 < F_R_cut ||
          tempC1 > cov_max || tempC2 > cov_max) {
        status_temp <- 0
      } else {
        genotyp <- paste(as.character(temp_data[1, 3]),
                         as.character(temp_data[2, 3]),
                         sep = "")
      }
    }

    # QC for homozygotes
    if (status_temp == 1) {
      if (tempQ1 < Q_cut ||
          tempC1 < cov_hom_cut ||
          tempF_R1 < F_R_cut ||
          tempC1 > cov_max) {
        status_temp <- 0
      } else {
        genotyp <- paste(as.character(temp_data[1, 3]),
                         as.character(temp_data[1, 3]),
                         sep = "")
      }
    }

    # Fill output table: chromosome, region, rsID
    tabell1[u, 1] <- temp_data[1, 1]
    tabell1[u, 2] <- temp_data[1, 2]
    tabell1[u, 3] <- as.character(temp_data[1, 7])

    # Genotype or no-call
    if (status_temp == 0) {
      tabell1[u, 4] <- "--"
    } else {
      tabell1[u, 4] <- genotyp
    }
  }


  write.table(
    tabell1,
    paste(outfile, "all_SNPs.txt", sep = "_"),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
}

