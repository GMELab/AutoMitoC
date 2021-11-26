args = commandArgs(trailingOnly=TRUE)
start_time<-Sys.time()
print(paste0("AutoMitoC Script v1.0 Initialized @ ",start_time))
# Load R packages # 
if (!require(data.table)) {install.packages('data.table')}
if (!require(ggplot2)) {install.packages('ggplot2')}
if (!require(parallel)) {install.packages('parallel')}
library(data.table)
library(ggplot2)
library(parallel)

print("Reading in Input Arguments...")
AUTO_LRR_MATRIX<-as.data.frame(fread(args[1],sep=",")) # AUTO_LRR.csv
MITO_LRR_MATRIX<-as.data.frame(fread(args[2],sep=",")) # MT_LRR.csv
SAMPLES_FILE<-as.data.frame(fread(args[3],sep=",")) # Column 1: Sample IDs ; Column 2: Gender (Male=1; Female=2) : AGE - numeric # Must be in the same order as Auto and Mito LRR Matrices # SAMPLE_FILE.csv
names(SAMPLES_FILE)[c(1:3)]<-c("SAMPLE","SEX","AGE") # Samples File can optionally include secondary mtDNA-CN estimates to benchmark in the fourth column # 
if (ncol(SAMPLES_FILE)==4){names(SAMPLES_FILE)[4]<-"MITO_BENCHMARK"}
NUM_CORES<-paste0(args[4]) # for mclapply 
OUTPUT_PREFIX<-paste0(args[5])

################### STEP1: Initial PCA Adjustment ####################
# i) Calculate initial autosomal PCs 
print("Step 1. Preliminary PCA Correction of Autosomal Probe Intensities ... ")
AUTO_PCA<-prcomp(AUTO_LRR_MATRIX,scale=TRUE,rank.=500) # 10 minutes # 
VAR_EXP<-as.data.frame(cbind(1:length(AUTO_PCA$sdev),AUTO_PCA$sdev^2/sum(AUTO_PCA$sdev^2))) # variance explained 
VAR_EXP[,3]<-cumsum(VAR_EXP[,2])
names(VAR_EXP)<-c("PC","IND_VAR_EXP","CUM_VAR_EXP")
# ii) Find K PCs corresponding to 70% of total variance explained 
K70<-which(abs(VAR_EXP$CUM_VAR_EXP-0.70)==min(abs(VAR_EXP$CUM_VAR_EXP-0.70)))
print(paste0("The Top ", K70," Autosomal Principal Components Explain ~70% variance in Probe Intensities and will be corrected for "))
#Plot PCA
pdf(paste0(OUTPUT_PREFIX,".STEP1.PCA.CUM_VAR.pdf"))
ggplot(VAR_EXP[1:min(500,dim(VAR_EXP)[1]),], aes(x = PC, y = CUM_VAR_EXP)) + geom_line() + geom_point() +ggtitle(paste0("PCA1: Cumulative Variance in Autosomal Probe Intensities Explained by Top PCs \n 70% @ ",K70, "PCs")) + theme(axis.text.x = element_text(angle=90, hjust = 1)) + xlab("# Top PCs") + ylab(" Cumulative Variance Explained") + geom_vline(xintercept = K70, linetype="dashed",color = "red")
dev.off()
# iii) Find sex-associated PCs and don't residualize out #### 
SEX_ASSOC_PC<-NULL
for (i in 1:K70) {if (summary(lm(SAMPLES_FILE$SEX~AUTO_PCA$x[,i]))$r.squared > 0.01) {SEX_ASSOC_PC<-c(SEX_ASSOC_PC,i)} } # R>0.1
# iv) Residualize out K PCs 
AUTO_LRR_LIST<-as.list(AUTO_LRR_MATRIX)
MITO_LRR_LIST<-as.list(MITO_LRR_MATRIX)
print(paste0("... Running Probe Correction with ",NUM_CORES," core(s)"))
if (class(SEX_ASSOC_PC)=="numeric" | class(SEX_ASSOC_PC)=="integer"){
RESID_AUTO<-as.data.frame(mclapply(AUTO_LRR_LIST,function(x){resid(lm(x~AUTO_PCA$x[,c(1:K70)[-SEX_ASSOC_PC]]))},mc.cores=NUM_CORES)) # 1 core = 6 hours; 12 cores ~ 50 minutes # 30 cores ~ 200 GB RAM ~ 
RESID_MITO<-as.data.frame(mclapply(MITO_LRR_LIST,function(x){resid(lm(x~AUTO_PCA$x[,c(1:K70)[-SEX_ASSOC_PC]]))},mc.cores=NUM_CORES)) 
} else {
RESID_AUTO<-as.data.frame(mclapply(AUTO_LRR_LIST,function(x){resid(lm(x~AUTO_PCA$x[,c(1:K70)]))},mc.cores=NUM_CORES)) # 1 core = 6 hours; 12 cores ~ 50 minutes # 30 cores ~ 200 GB RAM ~ 
RESID_MITO<-as.data.frame(mclapply(MITO_LRR_LIST,function(x){resid(lm(x~AUTO_PCA$x[,c(1:K70)]))},mc.cores=NUM_CORES)) 
}
######################################################################

#################### STEP2: Identification of Cross-Hybridizing Probes ####################
# i) Obtain median of adjusted mito and auto intensities # 
print("Step 2. Probe Cross-hybridization Check ... ")

MEDIAN_VALUES<-as.data.frame(cbind(SAMPLES_FILE[,c(1:3)],apply(RESID_MITO,1,median,na.rm=T),apply(RESID_AUTO,1,median,na.rm=T)))
names(MEDIAN_VALUES)<-c("SAMPLE","SEX","AGE","MEDIAN_MITO","MEDIAN_AUTO")
# ii) Association of autosomal probes with sex # 5 minutes # 

write.table(cbind("VARIANT_INDEX","R","PVALUE"),paste0(OUTPUT_PREFIX,".AUTO_SNP.SEX_ASSOC.csv"),sep=",",col.names=F,row.names=F,quote=F,append=F)
for (i in 1:ncol(RESID_AUTO)){
write.table(cbind(i,sqrt(summary(lm(MEDIAN_VALUES$SEX~RESID_AUTO[,i]))$r.squared),signif(summary(lm(MEDIAN_VALUES$SEX~RESID_AUTO[,i]))$coefficients[2,4],3)),paste0(OUTPUT_PREFIX,".AUTO_SNP.SEX_ASSOC.csv"),sep=",",col.names=F,row.names=F,quote=F,append=T)
}

# iii) association of autosomal probes with median mito # 5 minutes # 
write.table(cbind("VARIANT_INDEX","R","PVALUE"),paste0(OUTPUT_PREFIX,".AUTO_SNP.MED_MITO_ASSOC.csv"),sep=",",col.names=F,row.names=F,quote=F,append=F)
for (i in 1:ncol(RESID_AUTO)){
write.table(cbind(i,sqrt(summary(lm(MEDIAN_VALUES$MEDIAN_MITO~RESID_AUTO[,i]))$r.squared),signif(summary(lm(MEDIAN_VALUES$MEDIAN_MITO~RESID_AUTO[,i]))$coefficients[2,4],3)),paste0(OUTPUT_PREFIX,".AUTO_SNP.MED_MITO_ASSOC.csv"),sep=",",col.names=F,row.names=F,quote=F,append=T)
}
# iv) association of mito probes with median autosome # 
write.table(cbind("VARIANT_INDEX","R","PVALUE"),paste0(OUTPUT_PREFIX,".MITO_SNP.MED_AUTO_ASSOC.csv"),sep=",",col.names=F,row.names=F,quote=F,append=F)
for (i in 1:ncol(RESID_MITO)){
write.table(cbind(i,sqrt(summary(lm(MEDIAN_VALUES$MEDIAN_AUTO~RESID_MITO[,i]))$r.squared),signif(summary(lm(MEDIAN_VALUES$MEDIAN_AUTO~RESID_MITO[,i]))$coefficients[2,4],3)),paste0(OUTPUT_PREFIX,".MITO_SNP.MED_AUTO_ASSOC.csv"),sep=",",col.names=F,row.names=F,quote=F,append=T)
}
# v) association of mito probes with sex # since we expect an association with sex, increase the r correlation coefficient threshold to be more stringent  
write.table(cbind("VARIANT_INDEX","R","PVALUE"),paste0(OUTPUT_PREFIX,".MITO_SNP.SEX_ASSOC.csv"),sep=",",col.names=F,row.names=F,quote=F,append=F)
for (i in 1:ncol(RESID_MITO)){
write.table(cbind(i,sqrt(summary(lm(MEDIAN_VALUES$SEX~RESID_MITO[,i]))$r.squared),signif(summary(lm(MEDIAN_VALUES$SEX~RESID_MITO[,i]))$coefficients[2,4],3)),paste0(OUTPUT_PREFIX,".MITO_SNP.SEX_ASSOC.csv"),sep=",",col.names=F,row.names=F,quote=F,append=T)
}

AUTO_SNP_SEX_ASSOC<-read.csv(paste0(OUTPUT_PREFIX,".AUTO_SNP.SEX_ASSOC.csv"),header=T)
AUTO_SNP_MITO_ASSOC<-read.csv(paste0(OUTPUT_PREFIX,".AUTO_SNP.MED_MITO_ASSOC.csv"),header=T)
MITO_SNP_AUTO_ASSOC<-read.csv(paste0(OUTPUT_PREFIX,".MITO_SNP.MED_AUTO_ASSOC.csv"),header=T)
MITO_SNP_SEX_ASSOC<-read.csv(paste0(OUTPUT_PREFIX,".MITO_SNP.SEX_ASSOC.csv"),header=T)

# vi) Exclude probes exhibiting at least moderate correlation (r>0.05) # 
AUTO_SNP_RM <- which(AUTO_SNP_SEX_ASSOC$R > 0.05 | AUTO_SNP_MITO_ASSOC$R > 0.05) # 1605 SNPs total (1558 Sex; 48 Mito)
MITO_SNP_RM <- which(MITO_SNP_AUTO_ASSOC$R > 0.05 | MITO_SNP_SEX_ASSOC$R > 0.20) # 1 SNP

if (sum(AUTO_SNP_RM > 0)) {print(paste0("Removing ",length(AUTO_SNP_RM)," Autosomal probes with potential cross-hybridization to the sex (R > 0.05) or mitochondrial (R > 0.05) genomes "))}
if (sum(MITO_SNP_RM > 0)) {print(paste0("Removing ",length(MITO_SNP_RM)," Mitochondrial probes with potential cross-hybridization to the sex (R > 0.20) or autosomal (R > 0.05) genomes "))
print("*** Note that due to the robust epidemiological association between mtDNA-CN and sex, there is an expectation that mitochondrial probes will be associated with sex and therefore the correlation coefficient threshold for removing mitochondrial probes with evidence of cross-hybridization to sex chromosomes is more stringent")
}

######################################################################

#################### STEP3: Final PCA Adjustment ##################### 

# i) Redo PCA 
print("Step 3. Final PCA of Clean Autosomal Probe Intensities ...")
if (sum(AUTO_SNP_RM > 0)) {
AUTO_PCA2<-prcomp(AUTO_LRR_MATRIX[,-AUTO_SNP_RM],scale=TRUE,rank.=500)
VAR_EXP2<-as.data.frame(cbind(1:length(AUTO_PCA2$sdev),AUTO_PCA2$sdev^2/sum(AUTO_PCA2$sdev^2))) # variance explained 
VAR_EXP2[,3]<-cumsum(VAR_EXP2[,2])
names(VAR_EXP2)<-c("PC","IND_VAR_EXP","CUM_VAR_EXP")
# ii) Find K PCs corresponding to 70% of total variance explained 
K70_2<-which(abs(VAR_EXP2$CUM_VAR_EXP-0.70)==min(abs(VAR_EXP2$CUM_VAR_EXP-0.70)))
print(paste0("The Top ", K70_2," Autosomal Principal Components Explain ~70% variance in Probe Intensities and will be corrected for "))

#Plot PCA
pdf(paste0(OUTPUT_PREFIX,".STEP2.PCA.CUM_VAR.pdf"))
print(ggplot(VAR_EXP2[1:min(500,dim(VAR_EXP2)[1]),], aes(x = PC, y = CUM_VAR_EXP)) + geom_line() + geom_point() +ggtitle(paste0("PCA2: Autosomal LRR Variance Explained by PCs \n 70% @ ",K70_2, "PCs")) + theme(axis.text.x = element_text(angle=90, hjust = 1)) + xlab("# Top PCs") + ylab(" Cumulative Variance Explained") + geom_vline(xintercept = K70_2, linetype="dashed",color = "red"))
dev.off()
} else {
AUTO_PCA2 <- AUTO_PCA
}
# iii) Find strong sex-associated PCs and make sure to not residualize out #### 
SEX_ASSOC_PC2<-NULL
for (i in 1:K70_2) {if ((summary(lm(SAMPLES_FILE$SEX~AUTO_PCA2$x[,i]))$r.squared) > 0.01) {SEX_ASSOC_PC2<-c(SEX_ASSOC_PC2,i)} } # R>0.1
# iv) Residualize out K70 PCs from Mitochondrial SNPs # 
if (class(SEX_ASSOC_PC2)=="numeric"| class(SEX_ASSOC_PC2)=="integer"){
if (sum(MITO_SNP_RM > 0)) {
RESID_MITO2<-apply(MITO_LRR_MATRIX[,-MITO_SNP_RM],2,function(x){resid(lm(x~AUTO_PCA2$x[,c(1:K70_2)[-SEX_ASSOC_PC2]]))})
} else {
RESID_MITO2<-apply(MITO_LRR_MATRIX,2,function(x){resid(lm(x~AUTO_PCA2$x[,c(1:K70_2)[-SEX_ASSOC_PC2]]))})
}
}

if (class(SEX_ASSOC_PC2) == "NULL"){
if (sum(MITO_SNP_RM > 0)) {
RESID_MITO2<-apply(MITO_LRR_MATRIX[,-MITO_SNP_RM],2,function(x){resid(lm(x~AUTO_PCA2$x[,c(1:K70_2)]))}) 
} else {
RESID_MITO2<-apply(MITO_LRR_MATRIX,2,function(x){resid(lm(x~AUTO_PCA2$x[,c(1:K70_2)]))})
}
}

# v) Calculate and output final mtDNA-CN # 
RESULTS<-as.data.frame(cbind(SAMPLES_FILE[,c(1:3)],scale(apply(RESID_MITO2,1,median)),scale(prcomp(RESID_MITO2,scale=TRUE)$x[,1])))

names(RESULTS)<-c("SAMPLE","SEX","AGE","MEDIAN_MTDNACN","PC_MTDNACN")
# vi) Calibrate sign based on mito benchmark or if not available association with age (expect inverse association)
if (ncol(SAMPLES_FILE)==4){
ASSOC_SIGN_MED<-sign(summary(lm(RESULTS$MEDIAN_MTDNACN~SAMPLES_FILE$MITO_BENCHMARK))$coefficients[2,1])
ASSOC_SIGN_PCA<-sign(summary(lm(RESULTS$PC_MTDNACN~SAMPLES_FILE$MITO_BENCHMARK))$coefficients[2,1])
} else {
ASSOC_SIGN_MED<-sign(summary(lm(RESULTS$MEDIAN_MTDNACN~RESULTS$AGE))$coefficients[2,1])*-1
ASSOC_SIGN_PCA<-sign(summary(lm(RESULTS$PC_MTDNACN~RESULTS$AGE))$coefficients[2,1])*-1
}
RESULTS$MEDIAN_MTDNACN<-RESULTS$MEDIAN_MTDNACN*ASSOC_SIGN_MED
RESULTS$PC_MTDNACN<-RESULTS$PC_MTDNACN*ASSOC_SIGN_PCA
write.table(RESULTS,paste0(OUTPUT_PREFIX,".MTDNA_CN_RESULTS.csv"),sep=",",quote=F,row.names=F,col.names=T)
print("... AutoMitoC mtDNA-CN estimates are finished!")

######################################################################

#################### OPTIONAL STEP4: Benchmark mtDNA-CN against known estimates (OPTIONAL) ####################
if (ncol(SAMPLES_FILE)==4){
RESULTS2<-merge(RESULTS[,c("SAMPLE","MEDIAN_MTDNACN","PC_MTDNACN")],SAMPLES_FILE,by="SAMPLE")
benchmark_r<-round(sqrt(summary(lm(PC_MTDNACN~MITO_BENCHMARK,data=RESULTS2))$r.squared),3)
benchmark_p<-signif(summary(lm(PC_MTDNACN~MITO_BENCHMARK,data=RESULTS2))$coefficients[2,4],3)
if (benchmark_p ==0) {benchmark_p <- "< 5e-324"} 
pdf(paste0(OUTPUT_PREFIX,".MTDNA_CN_BENCHMARK.pdf"))
ggplot(RESULTS2,aes(y=scale(MITO_BENCHMARK),x=PC_MTDNACN))+geom_point()+geom_smooth(method="lm")+xlab("AutoMitoC Estimate (Standardized to Mean=0;SD=1)")+ylab("Reference Estimate (Standardized to Mean=0;SD=1)")+ggtitle(paste0("AutoMitoC Benchmarking Comparison \n R = ",benchmark_r))+xlim(-3.5,3.5)+ylim(-3.5,3.5)
dev.off()
print("... Benchmarking Complete!")
print(paste0("AutoMitoC vs. user-supplied mtDNA-CN measurement performance: " , "R = ",benchmark_r,"; Association P-value = ",benchmark_p ))
}

end_time<-Sys.time()
duration<-as.numeric(end_time)-as.numeric(start_time) # in seconds 
hours<- floor(duration / 3600)
minutes <- floor((duration-hours*3600) / 60)
seconds <- round(duration-hours*3600-minutes*60,2)
print(paste0("AutoMitoC Script Finished @ ",end_time))
print(paste0("Script Duration: ",hours ," hrs, ", minutes ," mins, ",seconds," secs ..."))

######################################################################
