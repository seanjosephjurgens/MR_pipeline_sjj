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
pheno1_file_vec <- opt$sumstats_path_pheno1
pheno2_file <- opt$sumstats_path_pheno2
pheno1_vec <- opt$name_pheno1
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

cat('\n\nStarting analayses: 2Sample-MR analysis for testing of multivariable causal associations for', pheno2, ' and several outcomes (', pheno1_vec, ').\n\n')


################################
# 1. Preprocessing files/data
################################

cat('\n\n\n-----------------------------------------------------------------\n')
cat('STEP1: Processing input data\n')
cat('-----------------------------------------------------------------\n\n\n\n')

## Make some directories
try(system(paste0("mkdir ", wd0, 'running/')), silent=TRUE)
try(system(paste0("mkdir ", wd0, 'results/')), silent=TRUE)
try(system(paste0("mkdir ", wd0, 'results/TwoSampleMR/')), silent=TRUE)
try(system(paste0("mkdir ", wd0, 'results/TwoSampleMR/figures/')), silent=TRUE)
try(system(paste0("mkdir ", wd0, 'results/TwoSampleMR/MVMR/')), silent=TRUE)


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
# Read in exposure study sum stats
pheno1_file_vec <- as.vector(base::strsplit(pheno1_file_vec, split=";")[[1]])
pheno1_vec	<- as.vector(base::strsplit(pheno1_vec, split=";")[[1]])
X1_list <- list()
for(pheno1_num in 1:length(pheno1_file_vec)){
	pheno1_file <- pheno1_file_vec[pheno1_num]
	pheno1<- pheno1_vec[pheno1_num]
	X1 <- NULL

	cat('\nReading in summary statistics for exposure trait,', pheno1, '...\n')
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
	X1_list[[pheno1]] <- X1
}

# Read in study two (outcome)
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


#######################################
# 2. Merging+Clumping for each exposure
#######################################

cat('\n\n\n-----------------------------------------------------------------\n')
cat('STEP2: Merging and LD-clumping for 2-Sample multivariable MR\n')
cat('-----------------------------------------------------------------\n\n\n\n')

X_list <- list()
#top_vars_tot <- NULL
overlapping_variants <- NULL
for(pheno1_num in 1:length(pheno1_file_vec)){
	#top_vars <- NULL
	X <- NULL
	
	pheno1<- pheno1_vec[pheno1_num]
	cat(paste0('\n\n\nMerging exposure-and-outcome data for ', pheno1, '...\n\n\n')) 
        cat('\nMerging summary statistics...\n')
	X1 <- X1_list[[pheno1]]
	X <- gwas_merge(X1, X2, snp_name_cols = c("varid", "varid"), 
        	                beta_hat_cols = c("BETA", "BETA"), 
                	        se_cols = c("SE", "SE"), 
                       	 	A1_cols = c("A1", "A1"), 
                       		A2_cols = c("A2", "A2") 
                       		#,pval_cols = c("P", "P")
	)
	if(pheno1_num==1){
		overlapping_variants <- X$snp
	}else{
		overlapping_variants <- overlapping_variants[overlapping_variants %in% X$snp]
	}
        cat(paste0('\nnumber of overlapping variants with previous combination data: ', length(overlapping_variants), '\n\n'))
	X_list[[pheno1]] <- X

}

cat(paste0('\n\n\nFiltering to ', length(overlapping_variants), ' overlapping variants...\n\n\n'))
for(pheno1_num in 1:length(pheno1_file_vec)){
        #top_vars <- NULL
	
	pheno1<- pheno1_vec[pheno1_num]
	X_list[[pheno1]] <- X_list[[pheno1]][ X_list[[pheno1]]$snp %in% overlapping_variants ,  ]
	X1_list[[pheno1]] <- X1_list[[pheno1]][ X1_list[[pheno1]]$varid %in% overlapping_variants ,  ]
	
}


# Save clump readible file
cat('\nMaking a combined clump file ordered by P-values in any of the exposures...\n')
X_t <- NULL
for(pheno1_num in 1:length(pheno1_file_vec)){
        pheno1<- pheno1_vec[pheno1_num]
	X <- X_list[[pheno1]]
	#X1 <- X1_list[[pheno1]]
	#X_t <- rbind(X_t, merge(X, X1[,c("varid", "P")], by.x="snp", by.y="varid", all=F))
	X$P <- pnorm(abs(X$beta_hat_1 / X$seb1), lower.tail=FALSE)
	X_t <- rbind(X_t, X)
}
X_t <- X_t[order(X_t$P), ]
X_t <- X_t[-which(duplicated(X_t$snp)), ]

cat('\nRunning LD clumping using UKBB EUR LD reference, for combined exposure file...\n')
r2_thresh = r2_cutoff
pval_thresh = 5e-8
distance = distance_cutoff
clumped_rez <- NULL
cat('\t\t using r^2<', r2_thresh, ', pval<', pval_thresh, ' and distance>', distance/1000, 'Mb\n')

clumped_rez <- NULL
wd <- paste0(wd0, 'running/')
phenotype <- NULL
for(pheno1_num in 1:length(pheno1_file_vec)){
	pheno1<- pheno1_vec[pheno1_num]
	phenotype <- paste0(phenotype, "_", pheno1)
}	
phenotype <- sub("_", "", paste0(phenotype, "__", pheno2))
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

cat("\t\tfound ", length(top_vars), "lead variants across exposure phenotypes to be used in multivariable 2Sample-MR.\n")


################################
# 3. Multi-variable MR running
################################

cat('\n\n\n-----------------------------------------------------------------\n')
cat('STEP3: Multivariable MR running\n')
cat('-----------------------------------------------------------------\n\n\n\n')

##### Running 2-Sample MR
cat("Collect all exposure variants for multivariable MR")
wd <- paste0(wd0, 'results/TwoSampleMR/MVMR/')
try(setwd(wd))

X_tot <- NULL
for(pheno1_num in 1:length(pheno1_file_vec)){
        pheno1<- pheno1_vec[pheno1_num]
	
	X <- X_list[[pheno1]]
	X_2SMR <- X[X$snp %in% top_vars, ]
 	X_2SMR$snp2 <- X_2SMR$snp
	X_2SMR <- tidyr::separate(X_2SMR, col=snp2, into=c("chr", "pos_a1_a2"), sep=":")
	X_2SMR$id1 <- pheno1
	X_2SMR$id2 <- pheno2
	X_2SMR$pval1 <- pnorm(abs(X_2SMR$beta_hat_1 / X_2SMR$seb1), lower.tail=FALSE)
	X_2SMR$pval2 <- pnorm(abs(X_2SMR$beta_hat_2 / X_2SMR$seb2), lower.tail=FALSE)
	X_tot <- rbind(X_tot, X_2SMR)

}
X_2SMR <- X_tot

cat("Start running 2-Sample multivariable MR.\n\n")
library(TwoSampleMR)

X_2SMR$eaf.exposure <- NA
X_2SMR <- X_2SMR[,c("snp", "id1", "id1", "A1", "A2", "eaf.exposure", "beta_hat_1", "seb1", "pval1")]
colnames(X_2SMR) <- c("SNP", "exposure", "id.exposure",    
		      "effect_allele.exposure", "other_allele.exposure", "eaf.exposure",
		      "beta.exposure", "se.exposure", "pval.exposure"       
)

X_tot$pos <- gsub("_.*", "", X_tot$pos)
X_tot$samplesize.outcome <- NA
X_tot$eaf.outcome <- NA
X_tot$mr_keep.outcome <- TRUE
X_tot$data_source.outcome <- "Sean_Pipeline"
X_tot <- X_tot[,c("snp", "chr", "pos", "beta_hat_2", "seb2", "samplesize.outcome", "pval2", 
		  "eaf.outcome", "A1", "A2", "id2", "id2", "id2", "id2",
		  "mr_keep.outcome", "data_source.outcome")
]
colnames(X_tot) <- c("SNP", "chr", "pos", "beta.outcome", "se.outcome", "samplesize.outcome",
		     "pval.outcome", "eaf.outcome", "effect_allele.outcome", "other_allele.outcome",
		     "outcome", "id.outcome", "originalname.outcome", "outcome.deprecated", 
		     "mr_keep.outcome", "data_source.outcome"
)

mvdat <- TwoSampleMR::mv_harmonise_data(X_2SMR, X_tot)

res <- NULL
res_wintercept <- NULL
res_residualized <- NULL

try(res <- TwoSampleMR::mv_multiple(mvdat))
try(res_wintercept <- TwoSampleMR::mv_multiple(mvdat, intercept=T))
try(res_residualized <- TwoSampleMR::mv_residual(mvdat))
save(res, res_wintercept, res_residualized, file=paste0('MVMR_results__', phenotype, '.RData'))

cat("\n\n\nMultivariable MR results:\n\n")
try(print(res))
cat("\n\n\nMultivariable MR results, including non-zero intercept:\n\n")
try(print(res_wintercept))
cat("\n\n\nMultivariable MR results, using original residualization method:\n\n")
try(print(res_residualized$result))


cat('\nDone!')
