#-------#-------#-------#-------#-------#-------#-------#-------#-------#-------#
#2018/07/05
#Raquel Aoki
#-------#-------#-------#-------#-------#-------#-------#-------#-------#-------#

#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#Script to get genes from gistic output
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#---- WORK DIRETORY
rm(list=ls(all=TRUE)) 
#setwd("C:\\Users\\raoki\\Documents\\Project 1")
setwd('C:\\Users\\raque\\Google Drive\\SFU\\CMPT884 - Course Project and presentation')

#---- PACKAGES
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)

#---- GET GENES POSITION 
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
cn = c("1","2","3","4","5","6","7","8","9","10","11","12",
       "13","14","15","16","17","18","19","20","21","22","X","Y")
genes <- getBM(attributes=c('ensembl_gene_id','gene_biotype','hgnc_symbol',
                            'chromosome_name','start_position','end_position'),
               filters = 'chromosome_name', values =cn, mart = ensembl)

genes[genes=='X' | genes=='Y'] = 23
genes = subset(genes, gene_biotype=='protein_coding')
genes = subset(genes, !is.na(hgnc_symbol) & hgnc_symbol!='')
rownames(genes) = NULL
names(genes)[c(5,6)] = c('start','end')
dim(genes)
head(genes)

#---- LOAD GISTIC OUTPUT
sg = read.table("scores_gistic.txt", sep='\t', header=T)
names(sg)=c('type','chromosome_name','start','end','log10qvalue','g_score',
            'average_amplitude','frequency')
sg = subset(sg, type=='Amp')


#---- LOOKING IF IN THE GENE POSITION THERE IS AN AMPLIFICATION
#I dont need all genes, only the ones with amplification areas
#Paper work only with amplifications
genes$g_score = c()
for(i in 1:dim(genes)[1]){
g_score = subset(sg, (sg$start<genes$start[i] & genes$start[i]<sg$end & genes$chromosome_name[i]==sg$chromosome_name) |
           (sg$start<genes$end[i] & genes$end[i]<sg$end & genes$chromosome_name[i]==sg$chromosome_name))#$g_score
if(dim(g_score)[1]!=0){
  genes$g_score[i] = sum(g_score$g_score*g_score$frequency)/sum(g_score$frequency)
}else{
  genes$g_score[i] = NA
}
}


summary(genes$g_score)
dim(subset(genes, is.na(g_score)))
dim(subset(genes, !is.na(g_score)))

write.table(genes, "g_score_by_gene.txt",sep=';',row.names = F)

#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#   join g-score and features 
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#

gsd = read.table('g_score_by_gene.txt',sep=';',header=T)
features = read.table('features_new_sirna.csv',sep=',',header=T)

gsd = subset(gsd, !duplicated(hgnc_symbol))
features = subset(features, !duplicated(gene))

data = merge(gsd,features,by.x='hgnc_symbol',by.y='gene',all.x=T)
dim(gsd)
dim(features)
dim(data)
names(data)[1] =  'gene'

write.table(data, 'data_gscore_features_bygene.csv', row.names = F, sep=';')

#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#   graphic presentation and report
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#
#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#-----#

#Features plot 
features = read.csv('features_new_sirna.csv', sep=',')
#features2 = read.csv('features_normalized.csv',sep=',')
head(features)

par(mfrow=c(3,1))
par(mar= c(3.1, 3.1, 1.1, 1.1))
hist(features$mutsig,col='#2a3990',main='Mutsig',xlab='',ylab='',cex.main=1.5)
hist(features$avg_expr,col='#2a3990',main='Expression',xlab='',ylab='',cex.main=1.5)
hist(features$sirna,col='#2a3990',main='SiRNA',xlab='',ylab='',cex.main=1.5)

#chain
cchain1= read.csv('complete_chain1_T.csv',sep=';', header = F)
cchain2= read.csv('complete_chain2_T.csv',sep=';', header = F)
cchainT= read.csv('complete_matrix_T.csv',sep=';', header = F)
cchainT = cchainT[-c(1:1600),]

schain2= read.csv('simplified_chain2_T.csv',sep=';', header = F)
schainT= read.csv('simplified_matrix_T.csv',sep=';', header = F)

#complete model
par(mfrow=c(1,3))
w_min_max = c(min(cchain1), max(cchain1)*1.5)
names(cchain1) = c('w0','w1','w2','w3')
plot(cchain1$w0,type='l', lwd =2,main='',ylim=w_min_max,ylab='Wi value',
     cex.axis=1.5,cex.main=1.5,cex.lab=1.5,xlab='iterations')
points(cchain1$w1,type='l', lwd =2,main='',ylim=w_min_max,col='darkblue')
points(cchain1$w2,type='l', lwd =2,main='',ylim=w_min_max,col='darkgreen')
points(cchain1$w3,type='l', lwd =2,main='',ylim=w_min_max,col='darkred')
abline(v=1600,lty=2)
legend('topright',lwd = c(3,3,3,3),col=c('black','darkblue','darkgreen','darkred'),
       legend = c('w0','w1','w2','w3'),bty='n', cex=1.5,ncol=2)

mu_min_max = c(min(cchain2[,c(1,2)]), 1.15*max(cchain2[,c(1,2)]))
names(cchain2) = c('mu0','mu1','var0','var1')
plot(cchain2$mu1,type='l', lwd =2,main='',ylim=mu_min_max,xlab='iterations', col = 'darkred',
     ylab=expression(paste(mu,' value',sep='')),cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
points(cchain2$mu0,type='l', lwd =2,main='w',ylim=mu_min_max,col='darkblue')
abline(v=1600,lty=2)
legend('topright',lwd = c(3,3),col=c('darkred','darkblue'), legend = c('Driver','Passenger'),
       bty='n', cex=1.5)

var_min_max = c(min(cchain2[,c(3,4)]), 1.15*max(cchain2[,c(3,4)]))
plot(cchain2$var1,type='l', lwd =2,main='',ylim=var_min_max,xlab='iterations', col = 'darkred',
     ylab=expression(paste(sigma,' value',sep='')),cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
points(cchain2$var0,type='l', lwd =2,main='w',ylim=var_min_max,col='darkblue')
abline(v=1600,lty=2)
legend('topright',lwd = c(3,3),col=c('darkred','darkblue'), legend = c('Driver','Passenger'),
       bty='n', cex=1.5)

#simplified model
par(mfrow=c(2,1))
par(mar=c(4,4.2,1,1))
w_min_max = c(min(schain2[,c(1:4)]), max(schain2[,c(1:4)])*2)
names(schain2) = c('w0','w1','w2','w3','mu0','mu1')
plot(schain2$w0,type='l', lwd =2,main='',ylim=w_min_max,ylab='Wi value',
     cex.axis=1.5,cex.main=1.5,cex.lab=1.5,xlab='iterations')
points(schain2$w1,type='l', lwd =2,main='',ylim=w_min_max,col='darkblue')
points(schain2$w2,type='l', lwd =2,main='',ylim=w_min_max,col='darkgreen')
points(schain2$w3,type='l', lwd =2,main='',ylim=w_min_max,col='darkred')
abline(v=1600,lty=2)
legend('topright',lwd = c(3,3,3,3),col=c('black','darkblue','darkgreen','darkred'),
       legend = c('w0','w1','w2','w3'),bty='n', cex=1.5,ncol=2)
par(mar=c(4,4.1,1,1))
mu_min_max = c(min(schain2[,c(5,6)]), 1.15*max(schain2[,c(5,6)]))
plot(schain2$mu1,type='l', lwd =2,main='',ylim=mu_min_max,xlab='iterations', col = 'darkred',
     ylab=expression(paste(mu,' value',sep='')),cex.axis=1.5,cex.main=1.5,cex.lab=1.5)
points(schain2$mu0,type='l', lwd =2,main='w',ylim=mu_min_max,col='darkblue')
abline(v=1600,lty=2)
legend('topright',lwd = c(3,3),col=c('darkred','darkblue'), legend = c('Driver','Passenger'),
       bty='n', cex=1.5)


#--------- Normal Mix
par(mar=c(1,1,1,1))
plot(  density(rnorm(1000000,0.08,0.015)),main='', col='red',lwd=4,   ylab='',xlab='',axes=F, xlim=c(0,0.15))
points(density(rnorm(1000000,0.04,0.02)) ,main='', col= 'blue', lwd=4,type='l')
text(0.08,-0.5,'G Score')
legend('topright',lwd=c(4,4),col=c('red','blue'),legend=c('Driver Genes','Passenger Genes'),bty='n')

