require("reticulate")

#Perform two degree of freedom test on the maternal and fetal data
 
rm(list=ls())

 args = commandArgs(TRUE)

#Read in GWAS of birthweight gzipped files downloaded from EGG consortium website
source_python("/home/fulminatingmoat/dev/projects/uq_imb/scripts/pickle_reader.py")
maternal <- read_pickle_file(args[1])
fetal <- read_pickle_file(args[2])
int <- args[3]
result_file <- args[4]
print("a")
#Column names from file
colnames(fetal) <- c("SNP", "chr","pos", "ea_fetal", "nea_fetal", "eaf_fetal", "Beta_fetal", "SE_fetal","p_fetal","het_p_fetal","n_fetal","rsid")
colnames(maternal) <- c("rsid","chr","pos","ea_maternal","nea_maternal","eaf_maternal","Beta_maternal","SE_maternal","p_maternal","het_p_maternal","n_maternal")
 
#Merge datasets
#NB! There appears to be many non-overlapping SNPs - not sure why this is the case? Looks like there are issues with just merging on rs number and there would be better ways to do this
data <- merge(fetal,maternal, by=c("rsid","chr","pos"))
print("b")
#Need to ensure that effect alleles are the same
#If alleles do not match, switch maternal alleles and change the sign of the beta coefficient
data$Beta_maternal_new <- ifelse(data$ea_fetal == data$nea_maternal, data$Beta_maternal*(-1), data$Beta_maternal)
data$ea_maternal_new <- ifelse(data$ea_fetal == data$nea_maternal, data$nea_maternal , data$ea_maternal)
data$nea_maternal_new <- ifelse(data$ea_fetal == data$nea_maternal, data$ea_maternal , data$nea_maternal)
data$eaf_maternal_new <- ifelse(data$ea_fetal == data$nea_maternal, 1-data$eaf_maternal , data$eaf_maternal)
###
 
data$Beta_maternal <- data$Beta_maternal_new
data$ea_maternal <- data$ea_maternal_new
data$nea_maternal <- data$nea_maternal_new
data$eaf_maternal <- data$eaf_maternal_new
 
#NB Insert intercept from LD score regression
int <- 0.1352

#Calculate conditional fetal and maternal regression coefficients
data$fetal_beta_adjusted <- ((4/3)*data$Beta_fetal) - ((2/3)*data$Beta_maternal)
data$maternal_beta_adjusted <- ((4/3)*data$Beta_maternal) - ((2/3)*data$Beta_fetal)
 
#Calculate variance of adjusted fetal and maternal regression coefficients
data$fetal_var_adjusted <- ((16/9)*(data$SE_fetal)^2) + ((4/9)*(data$SE_maternal)^2) - ((16/9)*int*data$SE_fetal*data$SE_maternal)
data$maternal_var_adjusted <- ((16/9)*(data$SE_maternal)^2) + ((4/9)*(data$SE_fetal)^2) - ((16/9)*int*data$SE_fetal*data$SE_maternal)
 
#Calculate covariance of adjusted fetal and maternal regression coefficients
data$covar <- ((20/9)*int*data$SE_fetal*data$SE_maternal) - ((8/9)*(data$SE_fetal)^2)-((8/9)*(data$SE_maternal)^2)

#Perform 2 df test 
effects <- matrix(nrow = 1, ncol = 2)
sigma <- matrix(nrow = 2, ncol = 2)
print("c")
chisq_2df <- vector(length= dim(data)[1])
 
for (i in 1:dim(data)[1]) {
effects <- matrix(c(data$fetal_beta_adjusted[i], data$maternal_beta_adjusted[i]), nrow=1, ncol=2, byrow=TRUE)
sigma <- matrix(c(data$fetal_var_adjusted[i], data$covar[i], data$covar[i], data$maternal_var_adjusted[i]), nrow=2, ncol=2, byrow=TRUE)
chisq_2df[i] <- effects%*%solve(sigma) %*%t(effects)
}

print("d")

pval_2df <- pchisq(chisq_2df, df = 2, ncp = 0, lower.tail = FALSE, log.p = FALSE)

#Perform conditional tests for maternal and fetal genome
chisq_fetal_adj <- (data$fetal_beta_adjusted^2)/(data$fetal_var_adjusted)
chisq_maternal_adj <- (data$maternal_beta_adjusted^2)/(data$maternal_var_adjusted)
pval_fetal_adj <- pchisq(chisq_fetal_adj, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
pval_maternal_adj <- pchisq(chisq_maternal_adj, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)

print("e")

#Output complete results in a nice format
complete <- subset(data, select = c(chr,rsid,pos,ea_fetal,nea_fetal,eaf_fetal,Beta_fetal,SE_fetal,p_fetal,Beta_maternal, SE_maternal, p_maternal, fetal_beta_adjusted, fetal_se_adjusted, pval_fetal_adj, maternal_beta_adjusted, maternal_se_adjusted, pval_maternal_adj, chi_2df, pval_2df))
colnames(complete) <- c("CHR","SNP","BP","A1","A2","FREQ","BETA_F","SE_F","PVAL_F","BETA_M","SE_M","PVAL_M","BETA_F_ADJ","SE_F_ADJ","PVAL_F_ADJ","BETA_M_ADJ","SE_M_ADJ","PVAL_M_ADJ","CHISQ_2DF","PVAL_2DF")

print("f")
write.table(complete, file=result_file, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


