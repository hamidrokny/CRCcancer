{
  library(fitdistrplus)
  library(logspline)
  library("mclust")
  library(tidyr)
  library(dplyr)
  library('extraDistr')
  library('VGAM')
  library('goftest')
  library('RVAideMemoire')
  library(e1071) 
}

cg<-read.csv('./Desktop/allllll/Colorectal/Colorectal_gene.csv')
cg<-cg[,-1]
cg$gm<-paste0(cg$geneID,'_',cg$motif)
gene<-cg[,c(1,9)]
gene<-cg[,c(1,5)]
gene<-unique(gene)
g<-as.data.frame(table(gene$geneID))
g<-as.data.frame(table(gene$gm))

descdist(g$Freq,boot = 1000, discrete = TRUE)
descdist(g$Freq, boot = 100, discrete = FALSE)

fit.gamma <- fitdist(g$Freq, "gamma")
fit.ln <- fitdist(g$Freq, "lnorm")
fit.w <- fitdist(g$Freq, "weibull")
fit.e <- fitdist(g$Freq, "exp")
fit.c <- fitdist(g$Freq, "pois")
fitdistr (g$Freq,'negative binomial')


hist(g$Freq, pch=2, breaks=seq(-0.5,301+0.5,1), prob=TRUE, main="",xlim = c(1,100),ylim = c(0,0.1),ylab="Probability",xlab="Number of mutation" )
curve(dgamma(x,2.22, 1.23 ), col="red", lwd=2, add=T)
curve(dexp(x,0.555 ), col="blue", lwd=2, add=T)
curve(dpois(x,lambda = 1.8), col="green",lwd=2, add=T)
curve(dnorm(x, 1.8, 2.29),     col="yellow", lwd=2, add=T) 
curve(dlnorm(x,0.348,0.578),     col="brown", lwd=2, add=T)# optim fit
curve(dweibull(x,1.23,1.96  ),     col="pink", lwd=2, add=T)
curve(dnbinom(x,2871.2,1.8  ),     col="pink", lwd=2, add=T)
curve(dbetabinom(x,536,0.03203024,0.02016863 ),     col="pink", lwd=2, add=T)
text(30, 0.80, "Gamma",col='red', cex = 1.1)
text(30, 0.75, "Exponential",col='blue', cex = 1.1)
text(30, 0.70, "Poisson",col='green', cex = 1.1)
text(30, 0.65, "Gaussian",col='yellow', cex = 1.1)
text(30, 0.60, "Lognormal",col='brown', cex = 1.1)
text(30, 0.55, "Weibull",col='pink', cex = 1.1)


bdata <- data.frame(N = 10, mu = 0.5, rho = 0.8)
bdata <- transform(bdata,
                   y = rbetabinom(100, size = N, prob = mu, rho = rho))
fit <- vglm(cbind(y, N-y) ~ 1, betabinomial, data = bdata, trace = TRUE)

ks.test(g$Freq,'ppois',17.1 ,alternative = 'two.sided')

fit <- vglm(cbind(g$Freq, 536 - g$Freq) ~ 1, betabinomial, data = g, trace = TRUE)
Coef(fit)
ks.test(g$Freq,rbetabinom(151,size=151,prob =0.011956986,rho = 0.002715817 ),alternative = 'two.sided')
ks.test(g$Freq,'pbetabinom', 151, mu =0.3016074, rho = 0.1545593 )
ks.test(g$Freq,'pbetabinom',301,0.05769866, 0.03791390 ,alternative = 'two.sided')
ks.test(g$Freq,pbetabinom(270479,151,0.01193757, 0.002715817),alternative = 'two.sided')
ks.test(g$Freq,'pbetabinom',301,0.05769866, 0.03791390 ,alternative = 'two.sided')
cvm.test(g$Freq,'pbetabinom',301,0.05769866, 0.03791390  )
ad.test(g$Freq,'pbetabinom', 301,0.05769866, 0.03791390)

fit <- vglm(cbind(g$Freq,536-g$Freq) ~ 1,  binomialff, data = g, trace = TRUE)
Coef(fit)
ks.test(g$Freq,'pbinom',301,0.05667386  ,alternative = 'two.sided')
cvm.test(g$Freq,'pbinom',301,0.05667386  )
ad.test(g$Freq,'pbinom', 301,0.05667386)
ad.test(g$Freq,'pbinom', 301,0.01193757)

appletree<-as.data.frame(table(g$Freq))
colnames(appletree)<-c('y','w')
appletree$y<-as.numeric(appletree$y)
fit <- vglm(y ~ 1, negbinomial(deviance = TRUE), data = appletree,
            weights = w, crit = "coef") 
Coef(fit)
ks.test(g$Freq,'pnbinom',size= 1.679114 ,mu= 16.923100)
cvm.test(g$Freq,'pnbinom',size= 1.679114 ,mu= 16.923100)
ad.test(g$Freq,'pnbinom', size= 1.679114 ,mu= 16.923100)
ad.test(g$Freq,'pnbinom', size= 3.211010 ,mu= 1.802572)

ks.test(g$Freq,'ppois',17.1  ,alternative = 'two.sided')
cvm.test(g$Freq,'ppois',17.1 )
ad.test(g$Freq,'ppois', 17.1)

g$pval<- 1- pbetabinom(g$Freq,151,0.01193757, 0.002715817) #gene-motif-FINAL
g$pval<- 1- pnbinom(g$Freq,size= 3.211010 ,mu= 1.802572) #gene-motif

g$pval<- 1- pnbinom(g$Freq,size =  1.679114 , mu = 16.923100 ) #gene-FINAL
g$pval<- 1- pnbinom(g$Freq,size =  1.64 , mu = 17.06 ) #gene
g$pval<- 1- pbetabinom(g$Freq,size = 301,prob =0.05769866, rho = 0.03791390 ) #gene
g$pval<- 1- pbinom(g$pval,301,0.05667386)
g$pval<- 1- pgamma(g$Freq,0.0586) #gene
g$pval<- 1- pgamma(g$Freq, 1.4571,0.0854) #gene
g$pval<- 1- plnorm(g$Freq, 2.456,0.934) #gene
g$pval<- 1- pweibull(g$Freq, 1.18,18.15) #gene
sgngm<-subset(g,g$Freq>2)
sgngm<-subset(g,g$Var1 %in% ag$`unique(b$geneID)`)
sgngm<-subset(g,g$pval<0.01)
sgngm<-sgngm[order(-sgngm$pval),]
tt<-sgngm
gt<-subset(cg,cg$gm %in% sgngm$Var1)
gt<-subset(cg,cg$geneID %in% tt$Var1)
gtt<-subset(gt,gt$geneID %in% tt$Var1)
gtt<-subset(gt,gt$gm %in% sgngm$Var1)
gtt<-subset(cg,cg$gm %in% sgngm$Var1 & cg$geneID %in% tt$Var1)
g<-gtt[,c(1,3,4,9)]
g<-unique(g)
g<-g[,c(1,4)]
colnames(g)<-c('Sample','motif')

g<-g[1:100,]

gene_sample<-g  %>% mutate(value=1) %>% complete(motif,Sample,fill=list(value=0)) %>%
  mutate(key=paste0(motif)) %>%
  group_by(Sample,key) %>%
  summarize(value = sum(value)) %>%
  spread(key,value) %>% 
  as.data.frame
write.csv(gene_sample,'./Desktop/table/motif.csv')
#=======histogram===========
library(gridExtra)
vegLengths<-rbind(g$Freq,rbetabinom(19870,536,0.03203024,0.02016863),rbinom(19870,301,0.05667386),rnbinom(19870,size= 1.679114 ,mu= 16.923100))
vegLengths<-rbind(g$Freq,rbetabinom(270479,size=151,prob =0.011956986,rho = 0.002715817),rnbinom(270479,size= 3.211010 ,mu= 1.802572))


real <- data.frame(length=g$Freq)
betabinomial <- data.frame(length=rbetabinom(19870,536,0.03203024,0.02016863))
binomial <- data.frame(length=rbinom(19870,536,0.03182618))
negativebinomial <- data.frame(length=rnbinom(19870,size= 1.679114 ,mu= 16.923100))
poisson <- data.frame(length=rpois(19870,17.1))
weibull <- data.frame(length=rweibull(19870,1.18,18.15))
gamma <- data.frame(length=rgamma(19870,1.4571,0.0854))

betabinomial <- data.frame(length=rbetabinom(270479,size=151,prob =0.011956986,rho = 0.002715817))
negativebinomial <- data.frame(length=rnbinom(270479,size= 3.211010 ,mu= 1.802572))

betabinomial$Distribution <- 'Beta binomial'
real$Distribution <- 'Real'
binomial$Distribution <- 'Binomial'
negativebinomial$Distribution <- 'Negative binomial'
poisson$Distribution <- 'Poisson'
weibull$Distribution <- 'Weibull'
gamma$Distribution <- 'Gamma'

vegLengths <- rbind( real,binomial,negativebinomial,poisson,betabinomial,weibull,gamma)
plot1<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
grid.arrange( plot1, nrow=1, ncol=1)

vegLengths <- rbind( real,betabinomial,binomial, negativebinomial,poisson,gamma,weibull)
plot1<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 301, aes(y = ..density..), position = 'identity')
grid.arrange( plot1, nrow=1, ncol=1)


vegLengths <- rbind( real,betabinomial)
plot1<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 301, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,binomial)
plot2<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 301, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,negativebinomial)
plot3<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 301, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,gamma)
plot4<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 301, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,poisson)
plot5<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 301, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,weibull)
plot6<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 301, aes(y = ..density..), position = 'identity')

grid.arrange( plot1,plot2,plot3,plot4,plot5,plot6, nrow=3, ncol=2)


vegLengths <- rbind( real,betabinomial)
plot1<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,binomial)
plot2<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,negativebinomial)
plot3<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,gamma)
plot4<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,poisson)
plot5<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,weibull)
plot6<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)

grid.arrange( plot1,plot2,plot3,plot4,plot5,plot6, nrow=3, ncol=2)
##########################Gene-motif
real <- data.frame(length=g$Freq)
betabinomial <- data.frame(length=rbetabinom(270479,536,0.0033643966, 0.0007376344))
binomial <- data.frame(length=rlnorm(270479, 0.368,0.609))
negativebinomial <- data.frame(length=rnbinom(270479,size= 1502.3462629 ,mu= 1.9124377))
poisson <- data.frame(length=rpois(270479,1.91))
weibull <- data.frame(length=rweibull(270479,1.15,2.04))
gamma <- data.frame(length=rgamma(270479,1.93,1.01))

betabinomial$Distribution <- 'Beta binomial'
real$Distribution <- 'Real'
binomial$Distribution <- 'Log-normal'
negativebinomial$Distribution <- 'Negative binomial'
poisson$Distribution <- 'Poisson'
weibull$Distribution <- 'Weibull'
gamma$Distribution <- 'Gamma'

vegLengths <- rbind( real,binomial,negativebinomial,poisson,betabinomial,weibull,gamma)
plot1<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
grid.arrange( plot1, nrow=1, ncol=1)

vegLengths <- rbind( real,betabinomial,binomial, negativebinomial,poisson,gamma,weibull)
plot1<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 301, aes(y = ..density..), position = 'identity')
grid.arrange( plot1, nrow=1, ncol=1)


vegLengths <- rbind( real,betabinomial)
plot1<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 151, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,binomial)
plot2<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 151, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,negativebinomial)
plot3<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 151, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,gamma)
plot4<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 151, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,poisson)
plot5<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 151, aes(y = ..density..), position = 'identity')
vegLengths <- rbind( real,weibull)
plot6<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_histogram(alpha = 0.5,bins = 151, aes(y = ..density..), position = 'identity')

grid.arrange( plot1,plot2,plot3,plot4,plot5,plot6, nrow=3, ncol=2)


vegLengths <- rbind( real,betabinomial)
plot1<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,binomial)
plot2<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,negativebinomial)
plot3<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,gamma)
plot4<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,poisson)
plot5<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)
vegLengths <- rbind( real,weibull)
plot6<-ggplot(vegLengths, aes(length, fill = Distribution)) + geom_density(alpha = 0.2)

grid.arrange( plot1,plot2,plot3,plot4,plot5,plot6, nrow=3, ncol=2)

