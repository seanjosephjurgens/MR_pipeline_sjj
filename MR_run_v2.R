#!/usr/bin/env Rscript

library(optparse)
option_list <- list(
  make_option("--sumstats_path_pheno1", type="character",default=""),
  make_option("--sumstats_path_pheno2", type="character",default=""),
  make_option("--name_pheno1", type="character",default=""),
  make_option("--name_pheno2", type="character",default=""),
  make_option("--run_working_directory", type="character",default=""),
  make_option("--CHR_col_pheno1", type="character",default="CHR"),
  make_option("--CHR_col_pheno2", type="character",default="CHR"),
  make_option("--BP_col_pheno1", type="character",default="BP"),
  make_option("--BP_col_pheno2", type="character",default="BP"),
  make_option("--A1_col_pheno1", type="character",default="ALLELE1"),
  make_option("--A1_col_pheno2", type="character",default="Allele1"),
  make_option("--A2_col_pheno1", type="character",default="ALLELE0"),
  make_option("--A2_col_pheno2", type="character",default="Allele2"),
  make_option("--Freq_col_pheno1", type="character",default="A1FREQ"),
  make_option("--Freq_col_pheno2", type="character",default="Freq1"),
  make_option("--INFO_col_pheno1", type="character",default="INFO"),
  make_option("--INFO_col_pheno2", type="character",default="INFO"),
  make_option("--BETA_col_pheno1", type="character",default="BETA"),
  make_option("--BETA_col_pheno2", type="character",default="Effect"),
  make_option("--SE_col_pheno1", type="character",default="SE"),
  make_option("--SE_col_pheno2", type="character",default="StdErr"),
  make_option("--Pval_col_pheno1", type="character",default="P_BOLT_LMM"),
  make_option("--Pval_col_pheno2", type="character",default="P-value"),
  make_option("--distance_cutoff", type="numeric",default=10000),
  make_option("--r2_cutoff", type="numeric",default=0.0005),
  make_option("--LD_reference_prefix", type="character", default="/medpop/afib/sjurgens/UKBB_ldref/merged/v2/UKBB_ldref_chr"),
  make_option("--LD_reference_suffix", type="character", default="_v2")
)
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)
pheno1_file <- opt$sumstats_path_pheno1
pheno2_file <- opt$sumstats_path_pheno2
pheno1 <- opt$name_pheno1
pheno2 <- opt$name_pheno2
wd0 <- opt$run_working_directory
CHR_col_pheno1 <- opt$CHR_col_pheno1
CHR_col_pheno2 <- opt$CHR_col_pheno2
BP_col_pheno1 <- opt$BP_col_pheno1
BP_col_pheno2 <- opt$BP_col_pheno2
A1_col_pheno1 <- opt$A1_col_pheno1
A1_col_pheno2 <- opt$A1_col_pheno2
A2_col_pheno1 <- opt$A2_col_pheno1
A2_col_pheno2 <- opt$A2_col_pheno2
Freq_col_pheno1 <- opt$Freq_col_pheno1
Freq_col_pheno2 <- opt$Freq_col_pheno2
INFO_col_pheno1 <- opt$INFO_col_pheno1
INFO_col_pheno2 <- opt$INFO_col_pheno2
BETA_col_pheno1 <- opt$BETA_col_pheno1
BETA_col_pheno2 <- opt$BETA_col_pheno2
SE_col_pheno1 <- opt$SE_col_pheno1
SE_col_pheno2 <- opt$SE_col_pheno2
Pval_col_pheno1 <- opt$Pval_col_pheno1
Pval_col_pheno2 <- opt$Pval_col_pheno2
distance_cutoff <- opt$distance_cutoff
r2_cutoff <- opt$r2_cutoff
LD_reference_prefix <- opt$LD_reference_prefix
LD_reference_suffix <- opt$LD_reference_suffix

#### CAUSE MR analysis, setup ####
#devtools::install_github("jean997/cause@v1.2.0", lib='~/R/x86_64-pc-linux-gnu-library/4.0')
#devtools::install_github("explodecomputer/genetics.binaRies", lib='~/R/x86_64-pc-linux-gnu-library/4.0')
library(cause)
library(readr)
library(dplyr)
library(data.table)

######library(TwoSampleMR) ### TwoSampleMR is needed, but is loaded when needed to prevent issues with cause
######library(tidyr) ### tidyr is loaded when needed

cat('\n\nStarting analayses: 2Sample-MR + CAUSE MR analysis for testing of causal association between', pheno1, 'and', pheno2, '.\n\n')


################################
# 1. Preprocessing files/data
################################

cat('\n\n\n-------------------------------------------------------\n')
cat('STEP1: Processing input data\n')
cat('-------------------------------------------------------\n\n\n\n')

## Make some directories
try(system(paste0("mkdir ", wd0, 'running/')), silent=TRUE)
try(system(paste0("mkdir ", wd0, 'results/')), silent=TRUE)
try(system(paste0("mkdir ", wd0, 'results/TwoSampleMR/')), silent=TRUE)
try(system(paste0("mkdir ", wd0, 'results/TwoSampleMR/figures/')), silent=TRUE)
try(system(paste0("mkdir ", wd0, 'results/CAUSE/')), silent=TRUE)
try(system(paste0("mkdir ", wd0, 'results/CAUSE/figures/')), silent=TRUE)


cat('Reading in UKBB 5k LD reference curated by Sean...\n')
# Read in LD reference file
ldfile <- NULL
for(CHR in c(1:22)){
        LDreference_prefix_bim <- paste0(LD_reference_prefix, CHR, LD_reference_suffix, '.bim')
	ldfile <- rbind(ldfile, fread(LDreference_prefix_bim, stringsAsFactors=F, data.table=F, header=F))
}
ldfile <- as.data.frame(ldfile$V2)
colnames(ldfile) <- "varid"

setwd(wd0)
# Read in study one sum stats
cat('\nReading in summary statistics for trait 1,', pheno1, '...\n')
X1 <- fread(pheno1_file, stringsAsFactors=F, data.table=F)
rm <- which(is.na(X1[,BETA_col_pheno1]) | is.na(X1[,Pval_col_pheno1]))
if(length(rm)>0){X1 <- X1[-rm,]}
if(INFO_col_pheno1 %in% colnames(X1)){
	cat('\tINFO column found, will subset to INFO>0.3...\n')
	X1$INFO <- X1[,INFO_col_pheno1]
        X1[is.na(X1$INFO), 'INFO'] <- 0
	X1 <- X1[X1$INFO>=0.3,]
}else{
	cat('\tno INFO column found.\n')
}
if(Freq_col_pheno1 %in% colnames(X1)){
	cat('\tFrequency column found, will subset to MAF>=0.01...\n')
	X1$MAF <- X1[,Freq_col_pheno1]
	X1[is.na(X1$MAF), 'MAF'] <- 0
	X1[X1$MAF>0.5,'MAF'] <- 1-X1[X1$MAF>0.5,'MAF']
	X1 <- X1[X1$MAF>=0.01,]
}else{
	cat('\tno Frequency column found.\n')
}
X1 <- X1[,c(CHR_col_pheno1, BP_col_pheno1, A1_col_pheno1, A2_col_pheno1, BETA_col_pheno1, SE_col_pheno1, Pval_col_pheno1)]
colnames(X1) <- c("CHR", "BP", "A1", "A2", "BETA", "SE", "P")
X1$A1 <- toupper(X1$A1)
X1$A2 <- toupper(X1$A2)
X1$CHR <- gsub("chr",	"", X1$CHR)
X1_1 <- X1
X1_1$varid <- paste0(X1_1$CHR, ":", X1_1$BP, "_", X1_1$A1, "_", X1_1$A2)
X1_2 <- X1
X1_2$varid <- paste0(X1_2$CHR, ":", X1_2$BP, "_", X1_2$A2, "_", X1_2$A1)
X1 <- rbind(X1_1, X1_2)
X1 <- merge(X1, ldfile, by="varid", all=F)

# Read in study two 
cat('\nReading in summary statistics for trait 2,', pheno2, '...\n')
X2 <- fread(pheno2_file, stringsAsFactors=F, data.table=F)
rm <- which(is.na(X2[,BETA_col_pheno2]) | is.na(X2[,Pval_col_pheno2]))
if(length(rm)>0){X2 <- X2[-rm,]}
if(INFO_col_pheno2 %in% colnames(X2)){
        cat('\tINFO column found, will subset to INFO>0.3...\n')
        X2$INFO <- X2[,INFO_col_pheno2]
        X2[is.na(X2$INFO), 'INFO'] <- 0
        X2 <- X2[X2$INFO>=0.3,]
}else{
      	cat('\tno INFO column found.\n')
}
if(Freq_col_pheno2 %in% colnames(X2)){
        cat('\tFrequency column found, will subset to MAF>=0.01...\n')
        X2$MAF <- X2[,Freq_col_pheno2]
        X2[is.na(X2$MAF), 'MAF'] <- 0
        X2[X2$MAF>0.5,'MAF'] <- 1-X2[X2$MAF>0.5,'MAF']
        X2 <- X2[X2$MAF>=0.01,]
}else{
      	cat('\tno Frequency column found.\n')
}
X2 <- X2[,c(CHR_col_pheno2, BP_col_pheno2, A1_col_pheno2, A2_col_pheno2, BETA_col_pheno2, SE_col_pheno2, Pval_col_pheno2)]
colnames(X2) <-	c("CHR", "BP", "A1", "A2", "BETA", "SE", "P")
X2$A1 <- toupper(X2$A1)
X2$A2 <- toupper(X2$A2)
X2$CHR <- gsub("chr",   "", X2$CHR)
X2_1 <- X2
X2_1$varid <- paste0(X2_1$CHR, ":", X2_1$BP, "_", X2_1$A1, "_", X2_1$A2)
X2_2 <- X2
X2_2$varid <- paste0(X2_2$CHR, ":", X2_2$BP, "_", X2_2$A2, "_", X2_2$A1)
X2 <- rbind(X2_1, X2_2)
X2 <- merge(X2, ldfile, by="varid", all=F)


# Merge data
cat('\nMerging summary statistics...\n')
X <- gwas_merge(X1, X2, snp_name_cols = c("varid", "varid"), 
                       beta_hat_cols = c("BETA", "BETA"), 
                       se_cols = c("SE", "SE"), 
                       A1_cols = c("A1", "A1"), 
                       A2_cols = c("A2", "A2") 
                       #,pval_cols = c("P", "P")
)

################################
# 2. Clumping for 2-Sample MR
################################

cat('\n\n\n-------------------------------------------------------\n')
cat('STEP2: LD-clumping for 2-Sample MR\n')
cat('-------------------------------------------------------\n\n\n\n')

##### LD pruning: Will use integrated LD reference for now, can use UKBB reference in more advanced runs later
# Save clump readible file
cat('\nRunning LD clumping using UKBB EUR LD reference, for running 2Sample-MR...\n')
wd <- paste0(wd0, 'running/')

r2_thresh = r2_cutoff
pval_thresh = 5e-8
distance = distance_cutoff
clumped_rez <- NULL
cat('\t\t using r^2<', r2_thresh, ', pval<', pval_thresh, ' and distance>', distance/1000, 'Mb\n')

phenotype <- paste0(pheno1, "__", pheno2)
X_t <- merge(X, X1[,c("varid", "P")], by.x="snp", by.y="varid", all=F)
Clumpfile <- X_t[,c("snp", "A1", "beta_hat_1", "P")]
colnames(Clumpfile) <- c("SNP", "A1", "BETA", "P")
write.table(Clumpfile, file=paste0(wd,phenotype,'_Clumpfile.txt'), col.names=T, row.names=F, quote=F)

for(CHR in c(1:22)){
        cat('\t\tBusy with chr', CHR, '...\n')
        #LDreference_prefix_bed <- paste0('/medpop/afib/sjurgens/UKBB_ldref/merged/v2/UKBB_ldref_chr',CHR, '_v2.bed')
        #LDreference_prefix_fam <- paste0('/medpop/afib/sjurgens/UKBB_ldref/merged/v2/UKBB_ldref_chr',CHR, '_v2.fam')
        #LDreference_prefix_bim <- paste0('/medpop/afib/sjurgens/UKBB_ldref/merged/v2/UKBB_ldref_chr',CHR, '_v2.bim')
	LDreference_prefix_bed <- paste0(LD_reference_prefix, CHR, LD_reference_suffix, '.bed')
	LDreference_prefix_fam <- paste0(LD_reference_prefix, CHR, LD_reference_suffix, '.fam')
	LDreference_prefix_bim <- paste0(LD_reference_prefix, CHR, LD_reference_suffix, '.bim')
        cat('   ... p==',pval_thresh,'and r2==',r2_thresh,'...\n')
        system(paste0('/medpop/afib/software/plink1.9/Aug16_2016/plink --bed ',LDreference_prefix_bed, '  --fam ', LDreference_prefix_fam, '  --bim ', LDreference_prefix_bim,
              ' --clump ',wd,phenotype,'_Clumpfile.txt --clump-field P --clump-p1 ', pval_thresh, ' --clump-p2 ', pval_thresh,
              ' --clump-r2  ', r2_thresh, '  --clump-kb ', distance, ' --out ',wd,phenotype,'_chr', CHR, '_Clumped_p',pval_thresh,'_r2', r2_thresh, '.txt'))
        if(file.exists(paste0(wd,phenotype,'_chr', CHR, '_Clumped_p',pval_thresh,'_r2', r2_thresh, '.txt.clumped'))){
                interres <- fread(paste0(wd,phenotype,'_chr', CHR, '_Clumped_p', pval_thresh, '_r2', r2_thresh, '.txt.clumped'), stringsAsFactors=F)
                interres <- as.data.frame(interres$SNP)
                colnames(interres)<-"SNP"
                Clumpedfile <- merge(interres, Clumpfile, by="SNP", all=F)
                system(paste0('rm ', wd,phenotype,'_chr', CHR, '_Clumped_p',pval_thresh,'_r2', r2_thresh, '.txt.clumped'))
		clumped_rez <- rbind(clumped_rez, Clumpedfile)
        }
        system(paste0('rm ', wd,phenotype,'_chr', CHR, '_Clumped_p',pval_thresh,'_r2', r2_thresh, '.txt.log'))
}
system(paste0('rm ',wd,phenotype,'_Clumpfile.txt'))
top_vars <- clumped_rez$SNP

cat("\t\tfound ", length(top_vars), "lead variants for 2Sample-MR.\n")


################################
# 3. Two-sample MR running
################################

cat('\n\n\n-------------------------------------------------------\n')
cat('STEP3: Running 2-Sample MR\n')
cat('-------------------------------------------------------\n\n\n\n')

##### Running 2-Sample MR
cat("Start running 2-Sample MR.\n\n")
wd <- paste0(wd0, 'TwoSampleMR/')
X_2SMR <- X[X$snp %in% top_vars, ]
X_2SMR$snp2 <- X_2SMR$snp
X_2SMR <- tidyr::separate(X_2SMR, col=snp2, into=c("chr", "pos_a1_a2"), sep=":")
X_2SMR$id1 <- pheno1
X_2SMR$id2 <- pheno2
X_2SMR$pval1 <- pnorm(abs(X_2SMR$beta_hat_1 / X_2SMR$seb1), lower.tail=FALSE)
X_2SMR$pval2 <-	pnorm(abs(X_2SMR$beta_hat_2 / X_2SMR$seb2), lower.tail=FALSE) 
exposures_mr <- TwoSampleMR::format_data(X_2SMR,
                type = "exposure",
                snp_col = "snp",
                beta_col = "beta_hat_1",
                se_col = "seb1",
                pval_col = "pval1",
                log_pval=FALSE,
                effect_allele_col = "A1",
                other_allele_col = "A2",
                phenotype_col = "id1")
outcome_mr <- TwoSampleMR::format_data(X_2SMR,
                type = "outcome",
                snp_col = "snp",
                beta_col = "beta_hat_2",
                se_col = "seb2",
                effect_allele_col = "A1",
                other_allele_col = "A2",
                pval_col = "pval2",
		phenotype_col = "id2"
)
TwoSMR <- TwoSampleMR::harmonise_data(
            exposure_dat = exposures_mr,
            outcome_dat = outcome_mr,
            action = 2
)
res_TwoSMR <- TwoSampleMR::mr(TwoSMR)
res_TwoSMR <-  TwoSampleMR::generate_odds_ratios(res_TwoSMR)
res_TwoSMR_eggerintercept <- TwoSampleMR::mr_pleiotropy_test(TwoSMR)
res_TwoSMR_het <- TwoSampleMR::mr_heterogeneity(TwoSMR, method_list = c("mr_ivw"))
save(res_TwoSMR, res_TwoSMR_eggerintercept, res_TwoSMR_het, file=paste0(wd0, 'results/TwoSampleMR/MR_TwoSampleMR_results_', phenotype, '.RData'))

res_TwoSMR_betaIVW <- res_TwoSMR[res_TwoSMR$method=="Inverse variance weighted", 'b']
res_TwoSMR_pvalIVW <- res_TwoSMR[res_TwoSMR$method=="Inverse variance weighted", 'pval']
res_TwoSMR_betaWM <- res_TwoSMR[res_TwoSMR$method=="Weighted median", 'b']
res_TwoSMR_pvalWM <- res_TwoSMR[res_TwoSMR$method=="Weighted median", 'pval']
cat("\nInverse Variance Weighted: beta=", res_TwoSMR_betaIVW, ", P=", res_TwoSMR_pvalIVW, "\n")
cat("Weighted Median: beta=", res_TwoSMR_betaWM, ", P=", res_TwoSMR_pvalWM, "\n")
res_TwoSMR_betaEgger <- res_TwoSMR[res_TwoSMR$method=="MR Egger", 'b']
res_TwoSMR_pvalEgger <- res_TwoSMR[res_TwoSMR$method=="MR Egger", 'pval']
res_TwoSMR_egger_int <- res_TwoSMR_eggerintercept[, 'egger_intercept']
res_TwoSMR_egger_pval <- res_TwoSMR_eggerintercept[, 'pval']
cat("\nMR Egger: beta=", res_TwoSMR_betaEgger, ", P=", res_TwoSMR_pvalEgger, "\n")
cat("MR Egger: intercept=", res_TwoSMR_egger_int, ", P=", res_TwoSMR_egger_pval, "\n")

### Plotting
cat("Plotting 2-Sample MR results.\n\n")
p1 <- TwoSampleMR::mr_scatter_plot(res_TwoSMR, TwoSMR)
res_TwoSMR_single <- TwoSampleMR::mr_singlesnp(TwoSMR, all_method=c("mr_weighted_median", "mr_egger_regression"))
p2 <- TwoSampleMR::mr_forest_plot(res_TwoSMR_single)
res_TwoSMR_loo <- TwoSampleMR::mr_leaveoneout(TwoSMR)
p3 <- TwoSampleMR::mr_leaveoneout_plot(res_TwoSMR_loo)
sin_TwoSMR <- TwoSampleMR::mr_singlesnp(TwoSMR, all_method=c("mr_weighted_median", "mr_egger_regression"))
p4 <- TwoSampleMR::mr_funnel_plot(sin_TwoSMR)
pdf(paste0(wd0, 'results/TwoSampleMR/figures/MR_TwoSampleMR_results_', phenotype, '_scatterplots.pdf'))
pdf.options(width = 7, height = 5)
for (i in 1:length(p1)){
	print(p1[[i]])
}
dev.off()
pdf(paste0(wd0, 'results/TwoSampleMR/figures/MR_TwoSampleMR_results_', phenotype, '_forestplots.pdf'))
pdf.options(width = 9, height = 10)
for (i in 1:length(p2)){
	print(p2[[i]])
}
dev.off()
pdf(paste0(wd0, 'results/TwoSampleMR/figures/MR_TwoSampleMR_results_', phenotype, '_LOVOplots.pdf'))
pdf.options(width = 9, height = 10)
for (i in 1:length(p3)){
	print(p3[[i]])
}
dev.off()
pdf(paste0(wd0, 'results/TwoSampleMR/figures/MR_TwoSampleMR_results_', phenotype, '_funnelplots.pdf'))
pdf.options(width = 9, height = 7)
for (i in 1:length(p4)){
	print(p4[[i]])
}
dev.off()



################################
# 4. CAUSE: step 1, nuissance param  
################################

cat('\n\n\n-------------------------------------------------------\n')
cat('STEP4: CAUSE nuissance paramaters\n')
cat('-------------------------------------------------------\n\n\n\n')

###### Step2: Calculate nuisance parameters, using 1M variants
cat('\nCompute CAUSE nuisance paramaters using 1M random variants...\n')
setwd(wd0)
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)

class(params)
params$rho
head(params$mix_grid)


################################
# 5. CAUSE: step 2, clumping
################################

cat('\n\n\n-------------------------------------------------------\n')
cat('STEP5: LD-clumping for CAUSE\n')
cat('-------------------------------------------------------\n\n\n\n')

##### LD pruning: Will use integrated LD reference for now, can use UKBB reference in more advanced runs later
# Save clump readible file
cat('\nRunning LD clumping using UKBB EUR LD reference, for CAUSE...\n')
wd <- paste0(wd0, 'running/')

r2_thresh = r2_cutoff
pval_thresh = 1e-3
distance = distance_cutoff
clumped_rez <- NULL
cat('\t\t using r^2<', r2_thresh, ', pval<', pval_thresh, ' and distance>', distance/1000, 'Mb\n')

phenotype <- paste0(pheno1, "__", pheno2)
X_t <- merge(X, X1[,c("varid", "P")], by.x="snp", by.y="varid", all=F)
Clumpfile <- X_t[,c("snp", "A1", "beta_hat_1", "P")]
colnames(Clumpfile) <- c("SNP", "A1", "BETA", "P")
write.table(Clumpfile, file=paste0(wd,phenotype,'_Clumpfile.txt'), col.names=T, row.names=F, quote=F)

for(CHR in c(1:22)){
	cat('Busy with chr', CHR, '...\n')
        #LDreference_prefix_bed <- paste0('/medpop/afib/sjurgens/UKBB_ldref/merged/v2/UKBB_ldref_chr',CHR, '_v2.bed')
        #LDreference_prefix_fam <- paste0('/medpop/afib/sjurgens/UKBB_ldref/merged/v2/UKBB_ldref_chr',CHR, '_v2.fam')
        #LDreference_prefix_bim <- paste0('/medpop/afib/sjurgens/UKBB_ldref/merged/v2/UKBB_ldref_chr',CHR, '_v2.bim')
	LDreference_prefix_bed <- paste0(LD_reference_prefix, CHR, LD_reference_suffix, '.bed')
	LDreference_prefix_fam <- paste0(LD_reference_prefix, CHR, LD_reference_suffix, '.fam')
	LDreference_prefix_bim <- paste0(LD_reference_prefix, CHR, LD_reference_suffix, '.bim')
	cat('   ... p==',pval_thresh,'and r2==',r2_thresh,'...\n')
	system(paste0('/medpop/afib/software/plink1.9/Aug16_2016/plink --bed ',LDreference_prefix_bed, '  --fam ', LDreference_prefix_fam, '  --bim ', LDreference_prefix_bim,
              ' --clump ',wd,phenotype,'_Clumpfile.txt --clump-field P --clump-p1 ', pval_thresh, ' --clump-p2 ', pval_thresh,
              ' --clump-r2  ', r2_thresh, '  --clump-kb ', distance, ' --out ',wd,phenotype,'_chr', CHR, '_Clumped_p',pval_thresh,'_r2', r2_thresh, '.txt'))
	if(file.exists(paste0(wd,phenotype,'_chr', CHR, '_Clumped_p',pval_thresh,'_r2', r2_thresh, '.txt.clumped'))){
        	interres <- fread(paste0(wd,phenotype,'_chr', CHR, '_Clumped_p', pval_thresh, '_r2', r2_thresh, '.txt.clumped'), stringsAsFactors=F)
        	interres <- as.data.frame(interres$SNP)
        	colnames(interres)<-"SNP"
        	Clumpedfile <- merge(interres, Clumpfile, by="SNP", all=F)
        	system(paste0('rm ', wd,phenotype,'_chr', CHR, '_Clumped_p',pval_thresh,'_r2', r2_thresh, '.txt.clumped'))
	}
	clumped_rez <- rbind(clumped_rez, Clumpedfile)
        system(paste0('rm ', wd,phenotype,'_chr', CHR, '_Clumped_p',pval_thresh,'_r2', r2_thresh, '.txt.log'))
}
system(paste0('rm ',wd,phenotype,'_Clumpfile.txt'))
top_vars <- clumped_rez$SNP


################################
# 6. CAUSE: step 3, actual algo       
################################

cat('\n\n\n-------------------------------------------------------\n')
cat('STEP6: Running CAUSE\n')
cat('-------------------------------------------------------\n\n\n\n')

# Run CAUSE
cat('\nRunning CAUSE...\n')
cat('\tusing ', length(top_vars), 'variants...\n')
res <- cause(X=X, variants = top_vars, param_ests = params)

# Diagnostics
cat('\nRunning diagnostics for model...\n')
res$loos[[2]]
loo::pareto_k_table(res$loos[[2]])
res$loos[[3]]
loo::pareto_k_table(res$loos[[3]])

# Results
cat('\nSummary results:\n')
summary(res, ci_size=0.95)

cat('\nSaving...\n')
save(res, file=paste0(wd0, 'results/CAUSE/MR_CAUSE_results_', phenotype, '.RData'))

# Making plots
cat('\nPlotting results...\n')
probability_density_per_model_plots <- plot(res)
variant_probability_per_model_plots <- plot(res, type="data")
save(probability_density_per_model_plots, variant_probability_per_model_plots, file=paste0(wd0, 'results/CAUSE/figures/MR_CAUSE_results_', phenotype, '_figures.RData'))

cat('\nDone!')
