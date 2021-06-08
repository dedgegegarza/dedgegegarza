# raw data
myRaw <- read.table(text = "
rs987435        C       G       1       1       1       0       2
rs345783        C       G       0       0       1       0       0
rs955894        G       T       1       1       2       2       1
rs6088791       A       G       1       2       0       0       1
rs11180435      C       T       1       0       1       1       1
rs17571465      A       T       1       2       2       2       2
rs17011450      C       T       2       2       2       2       2
rs6919430       A       C       2       1       2       2       2
rs2342723       C       T       0       2       0       0       0
rs11992567      C       T       2       2       2       2       2")

nIndividuals <- ncol(myRaw) - 3
nSNPs <- nrow(myRaw)

# make map, easy
MAP <- data.frame(
  CHR = 1,
  SNP = myRaw$V1,
  CM = 0,
  BP = seq(nSNPs))

# get first 6 columns of PED, easy
PED6 <- data.frame(
  FID = seq(nIndividuals),
  IID = seq(nIndividuals),
  FatherID = 0,
  MotherID = 0,
  Sex = 1,
  Phenotype = 1)

# convert 0,1,2 to genotypes, a bit tricky
# make helper dataframe for matching alleles
myAlleles <- data.frame(
  AA = paste(myRaw$V2, myRaw$V2),
  AB = paste(myRaw$V2, myRaw$V3),
  BB = paste(myRaw$V3, myRaw$V3))

# make index to match with alleles
PEDsnps <- myRaw[, 4:ncol(myRaw)] + 1

# convert
PEDsnpsAB <- 
  sapply(seq(nSNPs), function(snp)
    sapply(PEDsnps[snp, ], function(ind) myAlleles[snp, ind]))

# column bind first 6 cols with genotypes
PED <- cbind(PED6, PEDsnpsAB)

#output PED and MAP
write.table(PED, "gwas.ped", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
write.table(MAP, "gwas.map", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

# test plink
# plink --file gwas
# PLINK v1.90b3c 64-bit (2 Feb 2015)         https://www.cog-genomics.org/plink2
# (C) 2005-2015 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to plink.log.
# 258273 MB RAM detected; reserving 129136 MB for main workspace.
# .ped scan complete (for binary autoconversion).
# Performing single-pass .bed write (10 variants, 5 people).
# --file: plink.bed + plink.bim + plink.fam written.