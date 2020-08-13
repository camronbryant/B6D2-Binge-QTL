library(qtl)
library(qtlcharts)
library(XLConnect)
library(bestNormalize)

QT <- function(x) orderNorm(x)$x.t

#-----------------------------------------------------------------------------
B6D2_quantile <- phenotype_data
B6D2_quantile <- apply(B6D2_quantile, 2, QT)
writeWorksheetToFile(B6D2_quantile, file="quantile data.xlsx", sheet="Sheet1")
#-----------------------------------------------------------------------------

B6D2_quantile<-read.cross(format=c("csv"), file="B6D2_Quantile.csv", na.strings=c("-"),
                          genotypes=c("B","H","D"),alleles=c("B","D"),estimate.map=TRUE,
                          convertXdata=TRUE, error.prob=0.0001)

B6D2_quantile<-jittermap(B6D2_quantile)

#------------------------------------------
#DELETING MARKERS W CALL RATE < 95%
#total # individuals=131. 5% of 131=6.55. 1% of 131=1.31
na<-nmissing(B6D2_quantile, what="mar")
na<-as.data.frame(na)
na<-data.frame(names = row.names(na), na)
row.names(na)=NULL
na$names<-as.character(na$names)
i = 1
for (i in 1:nrow(na)){
  if(na$na[i] > 6.55){
    B6D2_quantile<-drop.markers(B6D2_quantile, na$names[i])
  } else{NULL}
}

#Counting number of crossover events. identifying and removing outliers
nxo <- countXO(B6D2_quantile)
plot(nxo, ylab="# crossovers")
nxo[nxo>500]
B6D2_quantile<-B6D2_quantile[ , -c(81,85,122)]
mean(nxo)

#Identify genotyping errors (double crossover)
B6D2_quantile<-calc.errorlod(B6D2_quantile, error.prob = 0.05)
top<-top.errorlod(B6D2_quantile, cutoff = 5)
top$id<-as.character(top$id)
plotGeno(B6D2_quantile, chr=14, ind=top$id[top$chr=="14"], cutoff=5)
top <- top[!duplicated(top$marker),]
B6D2_quantile<-drop.markers(B6D2_quantile, top$marker)

howmany<-markernames(B6D2_quantile)
length(howmany)
#------------------------------------------

qmap<-pull.map(B6D2_quantile)
par(font.lab=2, lwd=3, cex.axis=1)
plot(qmap,alternate.chrid = TRUE)
axis(side=1,lwd=3,labels=FALSE, at=c(0:20))

#effect plot
par(cex=1.25,mar=c(5,4,4,8),ann=TRUE,tck=NA,xaxt="s",yaxt="s")
effectplot(B6D2_quantile, mname1="13@30.792", pheno.col=5, ylab="",xlab="", col="red", ylim=c(-0.4,0.5),
           main= "ch13 30.792 cM")
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(B6D2_quantile, mname1="13@30.792", pheno.col=6, ylim=c(-0.4,0.5), col="orange")
effectplot(B6D2_quantile, mname1="13@30.792", pheno.col=7, ylim=c(-0.4,0.5), col="yellow")
effectplot(B6D2_quantile, mname1="13@30.792", pheno.col=8, ylim=c(-0.4,0.5), col="green")
effectplot(B6D2_quantile, mname1="13@30.792", pheno.col=9, ylim=c(-0.4,0.5), col="cyan")
effectplot(B6D2_quantile, mname1="13@30.792", pheno.col=10, ylim=c(-0.4,0.5), col="blue")
effectplot(B6D2_quantile, mname1="13@30.792", pheno.col=12, ylim=c(-0.4,0.5), col="purple")
par(cex=1)
legend("topright", legend=c("D2","D4","D9","D11","D16","D18","D23"), lwd=2, 
       col=c("red","orange","yellow","green","cyan","blue","purple"),bty="n",xpd=TRUE,inset=c(-0.50,0))

#Shapiro Wilk Test for normality
shapiro.test(B6D2_quantile$pheno$D2)
hist(B6D2_quantile$pheno$D2)


#calculating probability of genotypes
B6D2_quantile<-calc.genoprob(B6D2_quantile, step=1)
sex<-as.matrix(B6D2_quantile$pheno$Sex)
coh<-as.matrix(B6D2_quantile$pheno$Study)
B6D2_q_coh<-scanone(B6D2_quantile,pheno.col=c(25:27),addcovar=coh,intcovar=sex,method="hk")
B6D2_q_cohperm<-scanone(B6D2_quantile,pheno.col=c(25:27),addcovar=coh,intcovar=sex,n.perm=1000,method="hk")
writeWorksheetToFile(summary(B6D2_q_coh,alpha=0.05,format='allpheno',pvalues=TRUE,perms=B6D2_q_cohperm),file="B6D2_quantile_summary.xlsx",sheet="coh add, sex int")


bayesint(B6D2_q_coh, chr=13, 0.95, lodcolumn=30)
lodint(B6D2_q_coh, chr=13, 1.5, lodcolumn=30)

#plot intake phenotypes genome wide
par(cex=1.25,mar=c(5,4,4,8))
plot(B6D2_q_coh, lodcolumn=c(21:23), ylim = c(0,6), alternate.chrid=TRUE,col=c("red","blue","black"),ylab="LOD score",lwd=1)
plot(B6D2_q_coh, lodcolumn=c(30), ylim = c(0,6), add=TRUE,col=c("purple"),lwd=2,chr=13)
plot(B6D2_q_coh, lodcolumn=c(27:29), ylim = c(0,6), add=TRUE,col=c("green","cyan","blue"),lwd=2,chr=13)

add.threshold(B6D2_q_coh, perms=B6D2_q_cohperm, alpha=0.05)    
add.threshold(B6D2_q_coh, perms=B6D2_q_cohperm, alpha=0.63, lty=2)  
par(cex=1)
legend("topright",
       legend=c("D1 Right Time", "D22 Right Time", "D22-D1 Right Time"), lwd = 2, inset=c(-0.43,0),
       col = c("red","blue","black"),xpd=TRUE,bty="n")



#males and females separate
B6D2_quantile.hk<-scanone(B6D2_quantile, pheno.col=c(5:35), method="hk")
B6D2_quantile.hkperm<-scanone(B6D2_quantile, pheno.col=c(5:35), n.perm=1000, method="hk")

quantout.f<-scanone(subset(B6D2_quantile,ind=sex==0),
                    pheno.col=c(5:13), method="hk")
quantout.m<-scanone(subset(B6D2_quantile,ind = sex==1),
                    pheno.col=c(5:13), method="hk")
par(cex=1.25)
plot(out.m.coh, out.f.coh, lodcolumn=30, ylim=c(0,6),col = c("blue","red"),main='BW D23, chr. 13',ylab='LOD score', alternate.chrid=TRUE,lwd=2, ch=13)
add.threshold(out.f.coh, perms=out.f.cohperm, alpha=0.05)
legend("topright",
       legend=c("Males","Females"), lwd = 2, col = c("blue", "red"),xpd=TRUE,bty="n")

bayesint(B6D2_quantile.hki, chr=6, 0.95, lodcolumn=8)
lodint(B6D2_quantile.hki, chr=13, 1.5, lodcolumn=30)
writeWorksheetToFile(B6D2_quantile.hki,file="temp interval search.xlsx",sheet="ch13, sex combined")



#----------------QUANTILE TRANSFORMING B6D2 MALES--------------------
B6D2_quant_m <- phenotype_data
B6D2_quant_m <- apply(B6D2_quant_m, 2, QT)

writeWorksheetToFile(B6D2_quant_m, file="B6D2_quant_m.xlsx", sheet="Sheet1")


#------------------------------------------
out.m<-subset(B6D2_quantile, ind=sex==1)

m<-read.cross(format=c("csv"), file="B6D2_Quantile_M.csv", na.strings=c("-"),
              genotypes=c("B","H","D"),alleles=c("B","D"),estimate.map=TRUE,
              convertXdata=TRUE, error.prob=0.0001)

out.m$pheno<-m$pheno

#------------------------------------------

#calculating probability of genotypes
out.m<-calc.genoprob(out.m, step=1)


#coh additive
m.coh<-as.matrix(out.m$pheno$Study)
out.m.coh<-scanone(out.m, pheno.col=5:35, addcovar=m.coh, method="hk")
out.m.cohperm<-scanone(out.m, pheno.col=5:35, n.perm=1000, addcovar=m.coh, method="hk")
writeWorksheetToFile(summary(out.m.coh,alpha=0.05,format='allpheno',pvalues=TRUE,perms=out.m.cohperm),file="B6D2_quantile_summary.xlsx",sheet="MALES coh add")

lodint(out.m.coh, chr=6, 1.5, lodcolumn=8)
bayesint(out.m.coh, chr=6, 0.95, lodcolumn=8)
find.marker(out.m, chr=6, pos=52.86)

#effect plot
par(cex=1.25,mar=c(5,4,4,5),ann=TRUE,tck=NA,xaxt="s",yaxt="s",lwd=1)
effectplot(out.m, mname1="6@52.86", pheno.col=9, ylab="",xlab="",ylim=c(-1.1,1.1), col="cyan")
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(out.m, mname1="6@52.86", pheno.col=12, ylim=c(-1.1,1.1), col="black")

par(cex=1)
legend("topleft", legend=c("D16","D23"), lwd=2, col=c("cyan", "black"),bty="n",xpd=TRUE,y.intersp = 0.9)

#effect plot over days
par(cex=1.25,mar=c(5,4,4,7),ann=TRUE,tck=NA,xaxt="s",yaxt="s")
effectplot(out.f, mname1="13@32.98", pheno.col=28, ylab="",xlab="", col="red", ylim=c(-1.2,-0.03),
           main= "ch13 32.98 cM")
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(out.f, mname1="13@32.98", pheno.col=29, ylim=c(-1.2,-0.03), col="orange")
effectplot(out.f, mname1="13@32.98", pheno.col=30, ylim=c(-1.2,-0.03), col="yellow")
effectplot(out.f, mname1="13@32.98", pheno.col=31, ylim=c(-1.2,-0.03), col="green")
effectplot(out.f, mname1="13@32.98", pheno.col=32, ylim=c(-1.2,-0.03), col="cyan")
effectplot(out.f, mname1="13@32.98", pheno.col=33, ylim=c(-1.2,-0.03), col="blue")
effectplot(out.f, mname1="13@32.98", pheno.col=34, ylim=c(-1.2,-0.03), col="purple")
par(cex=1)
legend("topright", legend=c("BW2","BW4","BW9","BW11","BW16","BW18","BW23"), lwd=2, 
       col=c("red","orange","yellow","green","cyan","blue","purple"),bty="n",xpd=TRUE,inset=c(-0.50,0))


male.summary<-summary(out.m.coh, alpha=0.05, format='allpheno', pvalues=TRUE, perms=out.m.cohperm)
writeWorksheetToFile(male.summary,file="B6D2_quantile_summary.xlsx",sheet="males")


#calculating Bayes and 1.5 lod interval
bayesint(out.m.hk, chr=6, 0.95, lodcolumn=5)
lodint(out.f.coh, chr=13, 1.5, lodcolumn=30)


#plot all phenotypes genome wide
par(cex=1.25,mar=c(5,4,4,8))
plot(out.m.coh, lodcolumn=c(21:23), ylim = c(0,6), alternate.chrid=TRUE,col=c("red","blue","black"),lwd=1,ylab="LOD score")
plot(out.m.coh, lodcolumn=c(7:9), ylim = c(0,6), add=TRUE,col=c("purple","black","gray"),lwd=2,ch=6)
plot(out.m.coh, lodcolumn=c(4:6), ylim = c(0,6), add=TRUE,col=c("green","cyan","blue"),lwd=2,ch=6)

add.threshold(out.m.coh, perms=out.m.cohperm, alpha=0.05)
add.threshold(out.m.coh, perms=out.m.cohperm, alpha=0.63, lty=2)
par(cex=1)
legend("topright",
       c("D1 Right Time", "D22 Right Time", "D22-D1 Right Time"), lwd = 2, inset=c(-0.43,0),
       col = c("red", "blue", "black"),xpd=TRUE,bty="n")



#----------------QUANTILE TRANSFORMING B6D2 FEMALES--------------------
B6D2_quant_f <- phenotype_data
B6D2_quant_f <- apply(B6D2_quant_f, 2, QT)

writeWorksheetToFile(B6D2_quant_f, file="B6D2_Quantile_F.xlsx", sheet="Sheet1")
#------------------------------------------
out.f<-subset(B6D2_quantile, ind=sex==0)

f<-read.cross(format=c("csv"), file="B6D2_Quantile_F.csv", na.strings=c("-"),
              genotypes=c("B","H","D"),alleles=c("B","D"),estimate.map=TRUE,
              convertXdata=TRUE, error.prob=0.0001)

out.f$pheno<-f$pheno

#------------------------------------------
par(cex=1.25)
effectplot(out.f, mname1="6@69.24", pheno.col=27, ylab="",xlab="", main="D22-D1 Right Time
           chr.6 69.24 cM")


#multiple phenotype effect plot
par(cex=1.25,mar=c(5,4,4,5),ann=TRUE,tck=NA,xaxt="s",yaxt="s",lwd=1)
effectplot(out.f, mname1="6@10.65", pheno.col=10, ylab="",xlab="",ylim=c(-0.3,1), col="blue")
par(new=TRUE,ann=FALSE,tck=0,xaxt="n",yaxt="n")
effectplot(out.f, mname1="5@10.658", pheno.col=11, ylim=c(-0.3,1), col="purple")

par(cex=1)
legend("topleft", legend=c("D2","D18","Summed"), lwd=2, col=c("red","blue","purple"),bty="n",xpd=TRUE,y.intersp = 0.9)

#calculating probability of genotypes
out.f<-calc.genoprob(out.f, step=1)

#cohort additive
fcoh<-as.matrix(out.f$pheno$Study)
out.f.coh<-scanone(out.f, pheno.col=5:35, addcovar=fcoh, method="hk")
out.f.cohperm<-scanone(out.f, pheno.col=5:35, addcovar=fcoh, n.perm=1000, method='hk')
writeWorksheetToFile(summary(out.f.coh,alpha=0.05,format='allpheno',pvalues=TRUE,perms=out.f.cohperm),file="B6D2_quantile_summary.xlsx",sheet="FEMALES coh add 0.05")

#calculating Bayes and 1.5 lod interval
bayesint(out.f.coh, chr=6, 0.95, lodcolumn=21)
lodint(out.f.coh, chr=6, 1.5, lodcolumn=21)
find.marker(out.f, chr=6, pos=52.86)

#plot all phenotypes genome wide
par(cex=1.25,mar=c(5,4,4,5.5))
plot(out.f.coh, lodcolumn=c(1:3), ylim = c(0,6), alternate.chrid=TRUE,col=c("red","orange","yellow"),lwd=2,ylab="LOD score",ch=18)
plot(out.f.coh, lodcolumn=c(7:9), ylim = c(0,6), add=TRUE,col=c("purple","black","gray"),lwd=2,ch=18)
plot(out.f.coh, lodcolumn=c(4:6), ylim = c(0,6), add=TRUE,col=c("green","cyan","blue"),lwd=2,ch=18)

add.threshold(out.f.coh, perms=out.f.cohperm, alpha=0.05)    
add.threshold(out.f.coh, perms=out.f.cohperm, alpha=0.63, lty=2)  
par(cex=1)
legend("topright",
       colnames(out.f$pheno[c(5:13)]), lwd = 2, inset=c(-0.32,0),
       col = c("red","orange","yellow","green","cyan","blue","purple","black","gray"),xpd=TRUE,bty="n")


#--------------------MAKE QTL AND FIT QTL--------------------
all_5<-makeqtl(B6D2_quantile, 5, 25.824, what='prob')
all_6<-makeqtl(B6D2_quantile, 6, 52.863, what='prob')
all_13<-makeqtl(B6D2_quantile, 13, 30.792, what='prob')
cov<-cbind(coh,sex)
colnames(cov)<-c('coh','sex')

male_5_10<-makeqtl(B6D2_quant_m, 5, 10.658, what='prob')
male_5_24<-makeqtl(B6D2_quant_m, 5, 24.452, what='prob')
male_6_52<-makeqtl(B6D2_quant_m, 6, 52.863, what='prob')
colnames(m.coh)<-'coh'

female_18_23<-makeqtl(B6D2_quant_f, 18, 23.808, what='prob')
female_6_69<- makeqtl(out.f, 6, 69.24, what='prob')
colnames(fcoh)<-'coh'

summary(fitqtl(out.f, pheno.col=25, female_6_69, covar=fcoh, method='hk', formula=y~Q1+coh,
               get.ests=TRUE, run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE))
