# Start off by loading libraries
library(MASS)
library(LDlinkR)
library(rsnps)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(fastman)
library(ggcorrplot)



# Load data into R
GenotypeData <- read.csv("genotypes.csv", header = T, row.names = 1)
PhenotypeData <- read.csv("phenotypes.csv", header = T, row.names = 1)
CovariateData <- read.csv("covars.csv", header = T, row.names = 1)
GenoInfo <- read.csv("gene_info.csv", header = T)
SNPInfo <- read.csv("SNP_info.csv", header = T)



# Edit data
# Re-code: Male - 1 and Female - 0
CovariateData$Sex[CovariateData$Sex == "MALE"] <- 1
CovariateData$Sex[CovariateData$Sex == "FEMALE"] <- 0
CovariateData$Sex <- as.numeric(CovariateData$Sex)

# Re-code population data into indicator random variables with TSI being the 
# reference/base case 
CovariateData$PopCEU <- ifelse(CovariateData$Population == "CEU", 1, 0)
CovariateData$PopFIN <- ifelse(CovariateData$Population == "FIN", 1, 0)
CovariateData$PopGBR <- ifelse(CovariateData$Population == "GBR", 1, 0)
OriginalCovData <- CovariateData
CovariateData <- CovariateData[,2:5]

# Edit eQTL names to symbolic names for ease of use
colnames(PhenotypeData) <- c("ERAP2", "PEX6", "FAHD1", "GFM1", "MARCH7")

# For ploting Manhattan Plot
SNPInfo$color <- ifelse(SNPInfo$chromosome%%2==0, "#276FBF", "#183059")



# Plot Phenotype Data - If we are going to use a normal linear regression model 
# our phenotypes should be ~ normally distributed
perap2 <- ggplot(PhenotypeData, aes(x=ERAP2)) + 
  geom_histogram() + 
  labs(x="ERAP2 Phenotype Values", y="Frequency", 
       title = "ERAP2 Phenotype Histogram")
ppex6 <- ggplot(PhenotypeData, aes(x=PEX6)) + 
  geom_histogram() + 
  labs(x="PEX6 Phenotype Values", y="Frequency", 
       title = "PEX6 Phenotype Histogram")
pfahd1 <- ggplot(PhenotypeData, aes(x=FAHD1)) + 
  geom_histogram() + labs(x="FAHD1 Phenotype Values", y="Frequency", 
                          title = "FAHD1 Phenotype Histogram")
pgfm1 <- ggplot(PhenotypeData, aes(x=GFM1)) + 
  geom_histogram() + labs(x="GFM1 Phenotype Values", y="Frequency", 
                          title = "GFM1 Phenotype Histogram")
pMARCH7 <- ggplot(PhenotypeData, aes(x=MARCH7)) + 
  geom_histogram() + labs(x="MARCH7 Phenotype Values", y="Frequency", 
                          title = "MARCH7 Phenotype Histogram")
plot_grid(plotlist = list(perap2,ppex6,pfahd1,pgfm1,pMARCH7), nrow=2, ncol=3)



# Calculate MAF 
calcMAFs <- function(genotypeDF){
  dimsGeno <- dim(genotypeDF)
  iters <- dimsGeno[2]
  total <- 2*dimsGeno[1]
  mafs <- c()
  for(i in 1:iters){
    maf <- sum(genotypeDF[,i])/total
    if(maf > 0.5){
      maf <- 1 - maf
    }
    mafs <- c(mafs, maf)
  }
  return(mafs)
}

# Plot MAF
MAFs <- calcMAFs(GenotypeData)
mafDF <- data.frame(MAF = MAFs)
mafhist <- ggplot(mafDF, aes(x=MAF)) + geom_histogram() +
  labs(x="MAFs Values", y="Frequency", title = "MAFs Histogram")
mafhist


# Edit GenotypeData

editGenotypeData <- function(genotypeDF){
  dimsGeno <- dim(genotypeDF)
  n <- dimsGeno[1] # Number of individuals
  N <- dimsGeno[2] # Number of genotype/
  
  indivDrops <- c()
  genotypeDrops <- c()
  mafDrops <- c()
  
  for(i in 1:n){
    indiv <- genotypeDF[i,]
    # Drop individuals with > 10% missing data
    percent.missing <- sum(is.na(indiv))/N
    if(percent.missing > 0.10){
      indivDrops <- c(indivDrops, i)
    }
  }
  
  for(i in 1:N){
    locus <- genotypeDF[,i]
    # Drop locus with > 5% missing data
    percent.missing <- sum(is.na(locus))/n
    if(percent.missing > 0.05){
      genotypeDrops <- c(genotypeDrops, i)
    }
    # Drop locus with <= 5% MAF
    maf <- sum(locus)/(2*n)
    if(maf <= 0.05){
      mafDrops <- c(mafDrops, i)
    }
  }
  
  
  drop.list <- list(indivDrops, genotypeDrops, mafDrops)
  return(drop.list)
}

DropsList <- editGenotypeData(GenotypeData)
DropsList[[1]]
DropsList[[2]]
DropsList[[3]]



# PCA of GenotypeData

# Thin GenotypeData -> take every 10th snp
pca.index <- seq(from=1, to=50001, by=10)
iters <- length(pca.index)
pca.index[iters] <- 50000

# Make PCA dataframe
pcaDF <- GenotypeData[pca.index]

# Normalize each SNP
for(i in 1:iters){
  snps <- pcaDF[,i]
  avg.snps <- mean(snps)
  var.snps <- var(snps)
  norm.snps <- (snps - avg.snps)/var.snps
  pcaDF[,i] <- norm.snps
}

# Perform & plot PCA
pca.results <- prcomp(pcaDF)

pca.plotDF <- data.frame(PC1=pca.results$x[,1], PC2=pca.results$x[,2],
                         Population=OriginalCovData$Population)
PCA <- ggplot(pca.plotDF, aes(x=PC1,y=PC2)) + geom_point(aes(color=Population), 
                                                         alpha = 0.75) + 
  labs(title="Population Structure via PCA of Genotype Data")
PCA



# Creating Xa, Xd, Xz matrices

# Create Xa matrix
xa.mat <- as.matrix(GenotypeData)
xa.mat <- xa.mat - 1

# Create Xd matrix
xd.mat <- 1 - 2*abs(xa.mat)

# Create Xz matrix
xz.mat <- as.matrix(CovariateData)




# Create functions needed for GWAS analysis


# MLE calculator
calcMLE <- function(x, y){
  x <- as.matrix(x)
  y <- as.matrix(y)
  betas <- ginv(t(x) %*% x) %*% t(x) %*% y
  return(betas)
}

# Make Beta MLE DF parameters
calcBetas <- function(xa = NULL, xd = NULL, xz = NULL, y, snps=50000){
  indivs <- length(y)
  one.vec <- rep(1, times=indivs)
  if(is.null(xa) & is.null(xd) & is.null(xz)){
    b.mat <- matrix(NA, nrow=snps, ncol=1)
    betas <- calcMLE(one.vec, y)[1,1]
    betas <- rep(betas, times=snps)
    b.mat[,1] <- betas
    BetaDF <- as.data.frame(b.mat)
    colnames(BetaDF) <- c("Bu")
  }
  else if(!is.null(xa) & !is.null(xd) & is.null(xz)){
    b.mat <- matrix(NA, nrow=snps, ncol=3)
    for(i in 1:snps){
      x <- cbind(one.vec, xa[,i], xd[,i])
      betas <- calcMLE(x, y)
      b.mat[i,] <- betas
    }
    BetaDF <- as.data.frame(b.mat, row.names=colnames(xa))
    colnames(BetaDF) <- c("Bu", "Ba", "Bd")
  }
  else if(is.null(xa) & is.null(xd) & !is.null(xz)){
    b.mat <- matrix(NA, nrow=snps, ncol=5)
    x <- cbind(one.vec, xz)
    betas <- calcMLE(x, y)
    for(i in 1:snps){
      b.mat[i,] <- betas
    }
    BetaDF <- as.data.frame(b.mat)
    colnames(BetaDF) <- c("Bu", "Bsex", "Bcue", "Bfin", "Bgbr")
  }
  else if(!is.null(xa) & !is.null(xd) & !is.null(xz)){
    b.mat <- matrix(NA, nrow=snps, ncol=7)
    for(i in 1:snps){
      x <- cbind(one.vec, xa[,i], xd[,i], xz)
      betas <- calcMLE(x,y)
      b.mat[i,] <- betas
    }
    BetaDF <- as.data.frame(b.mat, row.names=colnames(xa))
    colnames(BetaDF) <- c("Bu", "Ba", "Bd", "Bsex", "Bcue", "Bfin", "Bgbr")
  }
  return(BetaDF)
}

# Calculate F-statistic
calcFstat <- function(y, yhat0, yhat1, df0, df1){
  SSE0 <- sum((y-yhat0)^2)
  SSE1 <- sum((y-yhat1)^2)
  num <- (SSE0 - SSE1)/(df0 - df1)
  den <- SSE1/df1
  f <- num/den
  return(f)
}

# Calculating Betas
ERAP2B0 <- calcBetas(y=PhenotypeData$ERAP2)
ERAP2B1 <- calcBetas(xa=xa.mat, xd=xd.mat, y=PhenotypeData$ERAP2)
ERAP2B0Cov <- calcBetas(xz=xz.mat,y=PhenotypeData$ERAP2)
ERAP2B1Cov <- calcBetas(xa=xa.mat, xd=xd.mat, xz=xz.mat, y=PhenotypeData$ERAP2)

PEX6B0 <- calcBetas(y=PhenotypeData$PEX6)
PEX6B1 <- calcBetas(xa=xa.mat, xd=xd.mat, y=PhenotypeData$PEX6)
PEX6B0Cov <- calcBetas(xz=xz.mat,y=PhenotypeData$PEX6)
PEX6B1Cov <- calcBetas(xa=xa.mat, xd=xd.mat, xz=xz.mat, y=PhenotypeData$PEX6)

FAHD1B0 <- calcBetas(y=PhenotypeData$FAHD1)
FAHD1B1 <- calcBetas(xa=xa.mat, xd=xd.mat, y=PhenotypeData$FAHD1)
FAHD1B0Cov <- calcBetas(xz=xz.mat,y=PhenotypeData$FAHD1)
FAHD1B1Cov <- calcBetas(xa=xa.mat, xd=xd.mat, xz=xz.mat, y=PhenotypeData$FAHD1)

GFM1B0 <- calcBetas(y=PhenotypeData$GFM1)
GFM1B1 <- calcBetas(xa=xa.mat, xd=xd.mat, y=PhenotypeData$GFM1)
GFM1B0Cov <- calcBetas(xz=xz.mat,y=PhenotypeData$GFM1)
GFM1B1Cov <- calcBetas(xa=xa.mat, xd=xd.mat, xz=xz.mat, y=PhenotypeData$GFM1)

MARCH7B0 <- calcBetas(y=PhenotypeData$MARCH7)
MARCH7B1 <- calcBetas(xa=xa.mat, xd=xd.mat, y=PhenotypeData$MARCH7)
MARCH7B0Cov <- calcBetas(xz=xz.mat,y=PhenotypeData$MARCH7)
MARCH7B1Cov <- calcBetas(xa=xa.mat, xd=xd.mat, xz=xz.mat, 
                         y=PhenotypeData$MARCH7)

B0 <- list(ERAP2B0, PEX6B0, FAHD1B0, GFM1B0, MARCH7B0)
B1 <- list(ERAP2B1, PEX6B1, FAHD1B1, GFM1B1, MARCH7B1)
B0Cov <- list(ERAP2B0Cov, PEX6B0Cov, FAHD1B0Cov, GFM1B0Cov, MARCH7B0Cov)
B1Cov <- list(ERAP2B1Cov, PEX6B1Cov, FAHD1B1Cov, GFM1B1Cov, MARCH7B1Cov)

# GWAS function
perform.GWAS <- function(xa, xd, xz = NULL, yDF, b0.list, b1.list){
  # Get Beta values
  ERAP2B0 <- b0.list[[1]]
  ERAP2B1 <- b1.list[[1]]
  PEX6B0 <- b0.list[[2]]
  PEX6B1 <- b1.list[[2]]
  FAHD1B0 <- b0.list[[3]]
  FAHD1B1 <- b1.list[[3]]
  GFM1B0 <- b0.list[[4]]
  GFM1B1 <- b1.list[[4]]
  MARCH7B0 <- b0.list[[5]]
  MARCH7B1 <- b1.list[[5]]
  
  # Get Phenotype Values
  ERAP2 <- yDF[,1]
  PEX6 <- yDF[,2]
  FAHD1 <- yDF[,3]
  GFM1 <- yDF[,4]
  MARCH7 <- yDF[,5]
  
  # Get Basic values
  N <- ncol(xa)
  n <- nrow(xa)
  one.vec <- rep(1, times = n)
  gwas.results <- matrix(NA, ncol=5, nrow=N)
  
  # Without covariates
  if(is.null(xz)){
    for(i in 1:N){
      x <- as.matrix(cbind(one.vec, xa[,i], xd[,i]))
      df0 <- n - 1
      df1 <- n - 3
      
      erapB0 <- as.numeric(ERAP2B0[i,])
      erapB1 <- as.numeric(ERAP2B1[i,])
      e.y.h0 <- one.vec %*% as.matrix(erapB0)
      e.y.h1 <- x %*% as.matrix(erapB1)
      e.f.stat <- calcFstat(ERAP2, e.y.h0, e.y.h1, df0, df1)
      e.p.val <- pf(e.f.stat, (df0 - df1), df1, lower.tail=F)
      
      pexB0 <- as.numeric(PEX6B0[i,])
      pexB1 <- as.numeric(PEX6B1[i,])
      p.y.h0 <- one.vec %*% as.matrix(pexB0)
      p.y.h1 <- x %*% as.matrix(pexB1)
      p.f.stat <- calcFstat(PEX6, p.y.h0, p.y.h1, df0, df1)
      p.p.val <- pf(p.f.stat, (df0 - df1), df1, lower.tail=F)
      
      
      fahdB0 <- as.numeric(FAHD1B0[i,])
      fahdB1 <- as.numeric(FAHD1B1[i,])
      f.y.h0 <- one.vec %*% as.matrix(fahdB0)
      f.y.h1 <- x %*% as.matrix(fahdB1)
      f.f.stat <- calcFstat(FAHD1, f.y.h0, f.y.h1, df0, df1)
      f.p.val <- pf(f.f.stat, (df0 - df1), df1, lower.tail=F)
      
      gfmB0 <- as.numeric(GFM1B0[i,])
      gfmB1 <- as.numeric(GFM1B1[i,])
      g.y.h0 <- one.vec %*% as.matrix(gfmB0)
      g.y.h1 <- x %*% as.matrix(gfmB1)
      g.f.stat <- calcFstat(GFM1, g.y.h0, g.y.h1, df0, df1)
      g.p.val <- pf(g.f.stat, (df0 - df1), df1, lower.tail=F)
      
      marchB0 <- as.numeric(MARCH7B0[i,])
      marchB1 <- as.numeric(MARCH7B1[i,])
      m.y.h0 <- one.vec %*% as.matrix(marchB0)
      m.y.h1 <- x %*% as.matrix(marchB1)
      m.f.stat <- calcFstat(MARCH7, m.y.h0, m.y.h1, df0, df1)
      m.p.val <- pf(m.f.stat, (df0 - df1), df1, lower.tail=F)
      
      gwas.results[i,] <- c(e.p.val, p.p.val, f.p.val, g.p.val, m.p.val) 
    } 
  }
  # With covariates
  else{
    for(i in 1:N){
      x0 <- as.matrix(cbind(one.vec, xz))
      x1 <- as.matrix(cbind(one.vec, xa[,i], xd[,i], xz))
      df0 <- n - 5
      df1 <- n - 7
      
      erapB0 <- as.numeric(ERAP2B0[i,])
      erapB1 <- as.numeric(ERAP2B1[i,])
      e.y.h0 <- x0 %*% as.matrix(erapB0)
      e.y.h1 <- x1 %*% as.matrix(erapB1)
      e.f.stat <- calcFstat(ERAP2, e.y.h0, e.y.h1, df0, df1)
      e.p.val <- pf(e.f.stat, (df0 - df1), df1, lower.tail=F)
      
      pexB0 <- as.numeric(PEX6B0[i,])
      pexB1 <- as.numeric(PEX6B1[i,])
      p.y.h0 <- x0 %*% as.matrix(pexB0)
      p.y.h1 <- x1 %*% as.matrix(pexB1)
      p.f.stat <- calcFstat(PEX6, p.y.h0, p.y.h1, df0, df1)
      p.p.val <- pf(p.f.stat, (df0 - df1), df1, lower.tail=F)
      
      
      fahdB0 <- as.numeric(FAHD1B0[i,])
      fahdB1 <- as.numeric(FAHD1B1[i,])
      f.y.h0 <- x0 %*% as.matrix(fahdB0)
      f.y.h1 <- x1 %*% as.matrix(fahdB1)
      f.f.stat <- calcFstat(FAHD1, f.y.h0, f.y.h1, df0, df1)
      f.p.val <- pf(f.f.stat, (df0 - df1), df1, lower.tail=F)
      
      gfmB0 <- as.numeric(GFM1B0[i,])
      gfmB1 <- as.numeric(GFM1B1[i,])
      g.y.h0 <- x0 %*% as.matrix(gfmB0)
      g.y.h1 <- x1 %*% as.matrix(gfmB1)
      g.f.stat <- calcFstat(GFM1, g.y.h0, g.y.h1, df0, df1)
      g.p.val <- pf(g.f.stat, (df0 - df1), df1, lower.tail=F)
      
      marchB0 <- as.numeric(MARCH7B0[i,])
      marchB1 <- as.numeric(MARCH7B1[i,])
      m.y.h0 <- x0 %*% as.matrix(marchB0)
      m.y.h1 <- x1 %*% as.matrix(marchB1)
      m.f.stat <- calcFstat(MARCH7, m.y.h0, m.y.h1, df0, df1)
      m.p.val <- pf(m.f.stat, (df0 - df1), df1, lower.tail=F)
      
      gwas.results[i,] <- c(e.p.val, p.p.val, f.p.val, g.p.val, m.p.val)
    }     
  }
  
  GWASDF <- as.data.frame(gwas.results, row.names=colnames(xa))
  colnames(GWASDF) <- c("ERAP2.p.val", "PEX6.p.val", "FAHD1.p.val", 
                        "GFM1.p.val", "MARCH7.p.val")
  
  return(GWASDF)
} 

# Get GWAS results for Models with & with out covariates
GWASResults <- perform.GWAS(xa=xa.mat, xd=xd.mat, 
                            yDF=PhenotypeData, b0.list=B0, b1.list=B1)
GWASResultsCov <- perform.GWAS(xa=xa.mat, xd=xd.mat, 
                               xz=xz.mat, yDF=PhenotypeData, 
                               b0.list=B0Cov, b1.list=B1Cov)


# Generate Diagnostic Plots

# QQ Plots & Manhattan Plots for each Phenotype with and without covariates

# Basic Calculations
expected_pvals <- qunif(seq(0, 1, length.out = ncol(xa.mat) + 2), 
                        min = 0, max = 1)
expected_pvals <- expected_pvals[expected_pvals != 0 & expected_pvals != 1]
expected_pvals <- -log10(expected_pvals)
expected_pvals <- expected_pvals[order(expected_pvals)]

fwer <- 0.05
bonferonni <- fwer/ncol(xa.mat)

# Plot for ERAP2

erap.obs <- -log10(GWASResults$ERAP2.p.val)
erap.obs <- erap.obs[order(erap.obs)]
erap.df.qqplot <- data.frame(obs = erap.obs, exp = expected_pvals)
erap.qqplot <- ggplot(erap.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'ERAP2 QQ plot')

erap.man.df <- data.frame(index=SNPInfo$chromosome, 
                          p.vals = GWASResults$ERAP2.p.val)
erap.manplot <- ggplot(erap.man.df, aes(x=index, y=-log10(p.vals))) +
  geom_point(alpha=0.5, color=SNPInfo$color) + labs(x="Chromosome", 
                                                    y="-log10(p.vals)", 
                                                    title="ERAP2 Manhattan Plot") 
+ geom_hline(yintercept=-log10(bonferonni), color="red")

erapcov.obs <- -log10(GWASResultsCov$ERAP2.p.val)
erapcov.obs <- erapcov.obs[order(erapcov.obs)]
erapcov.df.qqplot <- data.frame(obs = erapcov.obs, exp = expected_pvals)
erapcov.qqplot <- ggplot(erapcov.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'ERAP2 with Covariates QQ plot')

erapcov.man.df <- data.frame(index=SNPInfo$chromosome, 
                             p.vals = GWASResultsCov$ERAP2.p.val)
erapcov.manplot <- ggplot(erapcov.man.df, aes(x=index, y=-log10(p.vals))) 
+ geom_point(alpha=0.5, color=SNPInfo$color)+ 
  labs(x="Chromosome", y="-log10(p.vals)", 
       title ="ERAP2 with Covariates Manhattan Plot") +
  geom_hline(yintercept=-log10(bonferonni), color="red")

plot_grid(plotlist = list(erap.qqplot,erap.manplot,erapcov.qqplot,
                          erapcov.manplot), ncol=2, nrow=2)
#ggsave("~/Desktop/erap2.png", width = 10, height=5)

# Plot for PEX6

pex6.obs <- -log10(GWASResults$PEX6.p.val)
pex6.obs <- pex6.obs[order(pex6.obs)]
pex6.df.qqplot <- data.frame(obs = pex6.obs, exp = expected_pvals)
pex6.qqplot <- ggplot(pex6.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'PEX6 QQ plot')

pex6.man.df <- data.frame(index=SNPInfo$chromosome, 
                          p.vals = GWASResults$PEX6.p.val)
pex6.manplot <- ggplot(pex6.man.df, aes(x=index, y=-log10(p.vals))) 
+ geom_point(alpha=0.5, color=SNPInfo$color) 
+ labs(x="Chromosome", y="-log10(p.vals)", 
       title="PEX6 Manhattan Plot") 
+ geom_hline(yintercept=-log10(bonferonni), color="red")

pex6cov.obs <- -log10(GWASResultsCov$PEX6.p.val)
pex6cov.obs <- pex6cov.obs[order(pex6cov.obs)]
pex6cov.df.qqplot <- data.frame(obs = pex6cov.obs, exp = expected_pvals)
pex6cov.qqplot <- ggplot(pex6cov.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'PEX6 with Covariates QQ plot')

pex6cov.man.df <- data.frame(index=SNPInfo$chromosome, 
                             p.vals = GWASResultsCov$PEX6.p.val)
pex6cov.manplot <- ggplot(pex6cov.man.df, aes(x=index, y=-log10(p.vals))) 
+ geom_point(alpha=0.5, color = SNPInfo$color) 
+ labs(x="Chromosome", y="-log10(p.vals)"
       ,title="PEX6 with Covariates Manhattan Plot") 
+ geom_hline(yintercept=-log10(bonferonni), color="red")

plot_grid(plotlist = list(pex6.qqplot, pex6.manplot, 
                          pex6cov.qqplot, pex6cov.manplot), ncol=2, nrow=2)
#ggsave("~/Desktop/pex6.png", width = 10, height=5)


# Plot for FAHD1

fahd1.obs <- -log10(GWASResults$FAHD1.p.val)
fahd1.obs <- fahd1.obs[order(fahd1.obs)]
fahd1.df.qqplot <- data.frame(obs = fahd1.obs, exp = expected_pvals)
fahd1.qqplot <- ggplot(fahd1.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'FAHD1 QQ plot')

fahd1.man.df <- data.frame(index=SNPInfo$chromosome, 
                           p.vals = GWASResults$FAHD1.p.val)
fahd1.manplot <- ggplot(fahd1.man.df, aes(x=index, y=-log10(p.vals)))
+ geom_point(alpha=0.5, color = SNPInfo$color) 
+ labs(x="Chromosome", y="-log10(p.vals)", title="FAHD1 Manhattan Plot") 
+ geom_hline(yintercept=-log10(bonferonni), color="red")

fahd1cov.obs <- -log10(GWASResultsCov$FAHD1.p.val)
fahd1cov.obs <- fahd1cov.obs[order(fahd1cov.obs)]
fahd1cov.df.qqplot <- data.frame(obs = fahd1cov.obs, exp = expected_pvals)
fahd1cov.qqplot <- ggplot(fahd1cov.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'FAHD1 with Covariates QQ plot')

fahd1cov.man.df <- data.frame(index=SNPInfo$chromosome, 
                              p.vals = GWASResultsCov$FAHD1.p.val)
fahd1cov.manplot <- ggplot(fahd1cov.man.df, aes(x=index, y=-log10(p.vals)))
+ geom_point(alpha=0.5, color = SNPInfo$color) 
+ labs(x="Chromosome", y="-log10(p.vals)", 
       title="FAHD1 with Covariates Manhattan Plot") 
+ geom_hline(yintercept=-log10(bonferonni), color="red")

plot_grid(plotlist = list(fahd1.qqplot, fahd1.manplot, fahd1cov.qqplot, 
                          fahd1cov.manplot), ncol=2, nrow=2)
#ggsave("~/Desktop/fahd1.png", width = 10, height=5)

# Plot for GFM1

gfm1.obs <- -log10(GWASResults$GFM1.p.val)
gfm1.obs <- gfm1.obs[order(gfm1.obs)]
gfm1.df.qqplot <- data.frame(obs = gfm1.obs, exp = expected_pvals)
gfm1.qqplot <- ggplot(gfm1.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'GFM1 QQ plot')

gfm1.man.df <- data.frame(index=SNPInfo$chromosome, 
                          p.vals = GWASResults$GFM1.p.val)
gfm1.manplot <- ggplot(gfm1.man.df, aes(x=index, y=-log10(p.vals)))
+ geom_point(alpha=0.5, color = SNPInfo$color) + 
  labs(x="Chromosome", y="-log10(p.vals)", title="GFM1 Manhattan Plot") 
+ geom_hline(yintercept=-log10(bonferonni), color="red")

gfm1cov.obs <- -log10(GWASResultsCov$GFM1.p.val)
gfm1cov.obs <- gfm1cov.obs[order(gfm1cov.obs)]
gfm1cov.df.qqplot <- data.frame(obs = gfm1cov.obs, exp = expected_pvals)
gfm1cov.qqplot <- ggplot(gfm1cov.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'GFM1 with Covariates QQ plot')

gfm1cov.man.df <- data.frame(index=SNPInfo$chromosome, 
                             p.vals = GWASResultsCov$GFM1.p.val)
gfm1cov.manplot <- ggplot(gfm1cov.man.df, aes(x=index, y=-log10(p.vals)))
+ geom_point(alpha=0.5, color = SNPInfo$color) +
  labs(x="Chromosome", y="-log10(p.vals)", 
       title="GFM1 with Covariates Manhattan Plot") 
+ geom_hline(yintercept=-log10(bonferonni), color="red")

plot_grid(plotlist = list(gfm1.qqplot, gfm1.manplot,
                          gfm1cov.qqplot, gfm1cov.manplot), ncol=2, nrow=2)
#ggsave("~/Desktop/gfm1testttt.png", width = 10, height=5)

# Plot for MARCH7

march.obs <- -log10(GWASResults$MARCH7.p.val)
march.obs <- march.obs[order(march.obs)]
march.df.qqplot <- data.frame(obs = march.obs, exp = expected_pvals)
march.qqplot <- ggplot(march.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'MARCH7 QQ plot')

march.man.df <- data.frame(index=SNPInfo$chromosome,
                           p.vals = GWASResults$MARCH7.p.val)
march.manplot <- ggplot(march.man.df, aes(x=index, y=-log10(p.vals))) 
+ geom_point(alpha=0.5, color = SNPInfo$color) + 
  labs(x="Chromosome", y="-log10(p.vals)", title="MARCH7 Manhattan Plot") +
  geom_hline(yintercept=-log10(bonferonni), color="red")

marchcov.obs <- -log10(GWASResultsCov$MARCH7.p.val)
marchcov.obs <- marchcov.obs[order(marchcov.obs)]
marchcov.df.qqplot <- data.frame(obs = marchcov.obs, exp = expected_pvals)
marchcov.qqplot <- ggplot(marchcov.df.qqplot, aes(x = exp, y = obs)) +
  geom_point(alpha=0.35) +
  geom_abline(intercept = 0, slope = 1, col="purple") +
  labs(x = '-log10 Expected p-val',
       y = '-log10 Observed p-val',
       title = 'MARCH7 with Covariates QQ plot')

marchcov.man.df <- data.frame(index=SNPInfo$chromosome,
                              p.vals = GWASResultsCov$MARCH7.p.val)
marchcov.manplot <- ggplot(marchcov.man.df, aes(x=index, y=-log10(p.vals))) 
+ geom_point(alpha=0.5, color = SNPInfo$color) 
+ labs(x="Chromosome", y="-log10(p.vals)",
       title="MARCH7 with Covariates Manhattan Plot") 
+ geom_hline(yintercept=-log10(bonferonni), color="red")

plot_grid(plotlist = list(marchcov.qqplot, marchcov.manplot), ncol=2, nrow=1)
ggsave("~/Desktop/MARCH7testttttt.png", width = 10, height=5)




## Final Manhattan Plots

ERAP2cov <- data.frame(CHR=SNPInfo$chromosome, pos=SNPInfo$position, 
                       p.val=GWASResultsCov$ERAP2.p.val, rsID = SNPInfo$id)
fastman(ERAP2cov, chr = "CHR", bp="pos", 
        p="p.val", suggestiveline = F, 
        genomewideline = -log10(bonferonni), col=c("#276FBF", "#183059"), 
        main="ERAP2 Covariates Manhattan Plot", maxP = 100, ylim=c(0,100),
        xlim=c(0,23), chrlabs = names(table(ERAP2cov$CHR)))

PEX6cov <- data.frame(CHR=SNPInfo$chromosome, pos=SNPInfo$position, 
                      p.val=GWASResultsCov$PEX6.p.val, rsID = SNPInfo$id)
fastman(PEX6cov, chr = "CHR", bp="pos", p="p.val", 
        suggestiveline = F, genomewideline = -log10(bonferonni), 
        col=c("#276FBF", "#183059"), main="PEX6 Covariates Manhattan Plot", 
        maxP = 100, ylim=c(0,100), xlim=c(0,23), 
        chrlabs = names(table(ERAP2cov$CHR)))

FAHD1cov <- data.frame(CHR=SNPInfo$chromosome, pos=SNPInfo$position, 
                       p.val=GWASResultsCov$FAHD1.p.val, rsID = SNPInfo$id)
fastman(FAHD1cov, chr = "CHR", bp="pos", p="p.val", 
        suggestiveline = F, genomewideline = -log10(bonferonni), 
        col=c("#276FBF", "#183059"), main="FAHD1 Covariates Manhattan Plot", 
        maxP = 100, ylim=c(0,100), xlim=c(0,23), 
        chrlabs = names(table(ERAP2cov$CHR)))

GFM1cov <- data.frame(CHR=SNPInfo$chromosome, pos=SNPInfo$position, 
                      p.val=GWASResultsCov$GFM1.p.val, rsID = SNPInfo$id)
fastman(GFM1cov, chr = "CHR", bp="pos", p="p.val", suggestiveline = F, 
        genomewideline = -log10(bonferonni), col=c("#276FBF", "#183059"), 
        main="GFM1 Covariates Manhattan Plot", maxP = 100, ylim=c(0,100), 
        xlim=c(0,23), chrlabs = names(table(ERAP2cov$CHR)))

M17cov <- data.frame(CHR=SNPInfo$chromosome, pos=SNPInfo$position, 
                     p.val=GWASResultsCov$MARCH7.p.val, rsID = SNPInfo$id)
fastman(M17cov, chr = "CHR", bp="pos", p="p.val", suggestiveline = F, 
        genomewideline = -log10(bonferonni), col=c("#276FBF", "#183059"), 
      main="MARCH7 Covariates Manhattan Plot", maxP = 100, ylim=c(0,100),
      xlim=c(0,23))




# Lets zoom into our hits and get a better idea of what we are looking at -> 
# trying to identify the resolution of our GWAS
# NOTE: GFM1 & MARCH7 did not have any hits, hence, no plots will be made!

# ERAP2 best hits
erap2hits <- GWASResultsCov$ERAP2.p.val < bonferonni
erap2Hit.list <- list(GWASResultsCov$ERAP2.p.val[erap2hits], 
                      which(erap2hits), rownames(GWASResultsCov[erap2hits,]))
# GET info about our SNPS & the genes they come from via NCBI
specialHitsERAP2df <- ncbi_snp_query(snps=erap2Hit.list[[3]])
# Calculate the Correlation Matrix
erapcorrdf <- LDmatrix(snps = erap2Hit.list[[3]], pop = c("CEU", "GBR", "FIN", 
                                                          "TSI"), 
                       token = "65447c247e34", genome_build = "grch38")
# make df & plot
erap2hitsDF <- data.frame(y=-log10(erap2Hit.list[[1]]),
                          x=specialHitsERAP2df$bp, name=erap2Hit.list[[3]])
ezoom <- ggplot(erap2hitsDF, aes(x=x, y=y)) + 
  geom_point(alpha=0.85, color="purple") + 
  geom_label_repel(aes(label = name),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'blue',
                   max.overlaps = 20) +
  theme_classic() + labs(x="Chromosome 5 (bp)", 
                         y="-log10(p-values)", title = "Zoomed in ERAP2 Hits")
e2corrmat <- as.matrix(erapcorrdf[,-1])
e2corr <- ggcorrplot(e2corrmat, method="square",
                     type = "lower") + labs(title="ERAP2 LD Plot")
ezoom
#ggsave("~/Desktop/ezoom.png", width=15, heigh=8)
e2corr
#ggsave("~/Desktop/e2corr.png", width=15, heigh=8)

# PEX6 best hits
pex6hits <- GWASResultsCov$PEX6.p.val < bonferonni
pex6Hit.list <- list(GWASResultsCov$PEX6.p.val[pex6hits], 
                     which(pex6hits), rownames(GWASResultsCov[pex6hits,]))
# GET info about our SNPS & the genes they come from via NCBI
specialHitsPEX6df <- ncbi_snp_query(snps=pex6Hit.list[[3]])
# Calculate the Correlation Matrix
pex6corrdf <- LDmatrix(snps = pex6Hit.list[[3]][-1], 
                       pop = c("CEU", "GBR", "FIN", "TSI"), 
                       token = "65447c247e34", genome_build = "grch38")
# make df & plot
pex6hitsDF <- data.frame(y=-log10(pex6Hit.list[[1]][-1]), 
                         x=specialHitsPEX6df$bp[-1], name=pex6Hit.list[[3]][-1])
pzoom <- ggplot(pex6hitsDF, aes(x=x, y=y)) + 
  geom_point(alpha=0.85, color="purple") + 
  geom_label_repel(aes(label = name),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'blue',
                   max.overlaps = 25) +
  theme_classic() + labs(x="Chromosome 6 (bp)", 
                         y="-log10(p-values)", title = "Zoomed in PEX6 Hits")
p6corrmat <- as.matrix(pex6corrdf[,-1])
p6corr <- ggcorrplot(p6corrmat, method="square", 
                     type = "lower") + labs(title="PEX6 LD Plot")
pzoom
ggsave("~/Desktop/pzoom.png", width=15, heigh=8)
p6corr
#ggsave("~/Desktop/p6corr.png", width=15, heigh=8)


# FAHD1 best hits
fahd1hits <- GWASResultsCov$FAHD1.p.val < bonferonni
fahd1Hit.list <- list(GWASResultsCov$FAHD1.p.val[fahd1hits], 
                      which(fahd1hits), rownames(GWASResultsCov[fahd1hits,]))
# GET info about our SNPS & the genes they come from via NCBI
specialHitsFAHD1df <- ncbi_snp_query(snps=fahd1Hit.list[[3]])
# Calculate the Correlation Matrix
fahd1corrdf <- LDmatrix(snps = fahd1Hit.list[[3]], pop = c("CEU", "GBR", "FIN", 
                                                           "TSI"), 
                        token = "65447c247e34", genome_build = "grch38")
# make df & plot
fahd1hitsDF <- data.frame(y=-log10(fahd1Hit.list[[1]]), 
                          x=specialHitsFAHD1df$bp, name=fahd1Hit.list[[3]])
fzoom <- ggplot(fahd1hitsDF, aes(x=x, y=y)) + 
  geom_point(alpha=0.85, color="purple") + 
  geom_label_repel(aes(label = name),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'blue',
                   max.overlaps = 30) +
  theme_classic() + labs(x="Chromosome 16 (bp)", 
                         y="-log10(p-values)", title = "Zoomed in FAHD1 Hits")
f1corrmat <- as.matrix(fahd1corrdf[,-1])
f1corr <- ggcorrplot(f1corrmat, method="square", 
                     type = "lower") + labs(title="FAHD1 LD Plot")
fzoom
ggsave("~/Desktop/f1zoom.png", width=15, heigh=8)
f1corr