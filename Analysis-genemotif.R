{
  library('clues')
  library("mclust")
  library("ggplot2")
  library("factoextra")
  library(cluster)
  library(RColorBrewer)
  library(tidyr)
  library(dplyr)
  library(vegan)
  library(data.table)
  library("dbscan")
  library('ARTP2')
  library(rtracklayer)
  library(ggfortify)
  library(KernSmooth)
  library("fpc")
  library("gridExtra")
}
gene_sample<-read.csv('../soudbeh/Desktop/Gene-motif/gnmtf1.csv')
gene_sample<-read.csv('../soudbeh/Desktop/Gene-motif/gngnmt1.csv')
gene_sample<-read.csv('../hamed/Desktop/allllll/table/gngnmt1.csv')
gene_sample<-read.csv('./Desktop/gsgms-sample.csv')
gene_sample<-gene_sample[,-1]
row.names(gene_sample)<-gene_sample$Sample
gene_sample<-gene_sample[,-1]

t<-gene_sample[apply(gene_sample[,-1], 1, function(x) !all(x==0)),]
t<-gene_sample[apply(gene_sample[,-1], 1, function(x) all(x==0)),]
t<-gene_sample[,which(!apply(gene_sample,2,FUN = function(x){all(x == 0)}))]
{
  y1<-as.data.frame(cl1)
  y1<-as.data.frame(cl2)
}
s<-cbind(cl1,cl2,cl3,cl4)
y<-as.data.frame(y1)
t<-subset(gene_sample,row.names(gene_sample) %in% row.names(y))
#t<-t[,which(!apply(t,2,FUN = function(x){all(x == 0)}))]
y<-as.matrix(t)
#y<-as.matrix(gene_sample)
mc12 <- Mclust(y)            # Model-based-clustering
summary(mc12)
class<-mc12$classification
cl1<-subset(class,class==1)
cl2<-subset(class,class==2)

cl21<-subset(class,class==1)
cl22<-subset(class,class==2)

cl21<-subset(class,class==1)
cl21<-subset(class,class==2)


{
  y1<-as.data.frame(t[,1])
  row.names(y1)<-row.names(t)
  y2<-as.data.frame(cl1)
  y3<-as.data.frame(cl2)
  y4<-as.data.frame(cl3)
  y5<-as.data.frame(cl4)
  colnames(y1)<-'1'
  colnames(y2)<-'1'
  colnames(y3)<-'1'
  colnames(y4)<-'1'
  colnames(y5)<-'1'
  y1$`1`<-1
  y2$`1`<-2
  y3$`1`<-3
  y4$`1`<-4
  y5$`1`<-5
  y<-rbind(y1,y2,y3,y4,y5)
}
{
  y1<-as.data.frame(cl11)
  y2<-as.data.frame(cl12)
  y3<-as.data.frame(cl13)
  y4<-as.data.frame(cl14)
  y5<-as.data.frame(cl21)
  y6<-as.data.frame(cl22)
  y7<-as.data.frame(cl23)
  colnames(y1)<-'1'
  colnames(y2)<-'1'
  colnames(y3)<-'1'
  colnames(y4)<-'1'
  colnames(y5)<-'1'
  colnames(y6)<-'1'
  colnames(y7)<-'1'
  y1$`1`<-1
  y2$`1`<-2
  y3$`1`<-3
  y4$`1`<-4
  y5$`1`<-5
  y6$`1`<-6
  y7$`1`<-7
  y<-rbind(y1,y2,y3,y4,y5,y6,y7)
}
s<-as.matrix(y)
s<-s[order(row.names(s)),]
fviz_mclust(mc12, "BIC", palette = "jco")
fviz_mclust(mc12, "classification", geom = "point", 
            pointsize = 1.5, palette = "jco")
fviz_mclust(mc12, "uncertainty", palette = "jco")
plot(mc12, "density")
cs = cluster.stats(gene_sample, mc12$classification)
cs[c("within.cluster.ss","avg.silwidth")]
ss <- silhouette(s, dist(df))
plot(ss,col=y$X1)

pca<-prcomp(t(gene_sample))

pres <- as.data.frame(prcomp(t(gene_sample))[[2]])[,c(1,2)]
plot(pres$PC1, pres$PC2)
library(ggplot2)
identical(as.character(y$X), rownames(pres))
pres<-pres[match(as.character(y$X),rownames(pres)),]
pres$color <- factor(y$X1)


ggplot(pres[which(pres$PC1 < 0.05 & pres$PC2<0.05 & pres$PC1 > -0.05 & pres$PC2>-0.05),]
       , aes(PC1,PC2), fill = color)+
  geom_point(aes(colour = color))


y24<-subset(y,y$X1 ==2  | y$X1 ==4)
gene_sample24<- subset(gene_sample,rownames(gene_sample) %in% y234$X)
pres <- as.data.frame(prcomp(t(gene_sample24))[[2]])[,c(1,2)]
plot(pres$PC1, pres$PC2)

identical(as.character(y24$X), rownames(pres))
pres<-pres[match(as.character(y24$X),rownames(pres)),]
pres$color <- factor(y24$X1)
ggplot(pres, aes(PC1,PC2), fill = color)+
  geom_point(aes(colour = color))
#============dbscan=========
{
  cl <- hdbscan(y, minPts = 50)
  cl
}
{
  cl<-dbscan(y,MinPts = 20,eps = 2)
  cl
}
#===========k-means==========
cl <- kmeans(gene_sample, 9)
fviz_cluster(list(data = gene_sample, cluster = s),geom = ('point'),ellipse.type = 't',show.clust.cent=FALSE, main = "Clustering results in TWO first PCAs vector")

#==========Hirarchical===============
cl <- hclust(dist(gene_sample, method = "binary"), method = 'average')
plot(cl)
cl <- hclust(dist(gene_sample, method = "euclidean"), method = 'average')
plot(cl)
cl <- hclust(dist(gene_sample, method = "euclidean"), method = 'complete')
plot(cl)
cl <- hclust(dist(gene_sample, method = "euclidean"), method = 'ward.D')
plot(cl)
clusterCut <- cutree(cl, 5)
{
  clust <- cutree(cl, k = 6)
  fviz_cluster(list(data = gene_sample, cluster = clust),geom = ('point'),ellipse.type = 't',show.clust.cent=FALSE, main = "Clustering results in TWO first PCAs vector")  ## from ‘factoextra’ package
}

#========Shinking========
cl<-clues(y, n0 = 5, alpha = 0.05, eps = 1.0e-4, itmax = 50,
          K2.vec = NULL, strengthMethod = "sil", strengthIni = -3,
          disMethod ="Euclidean", quiet = TRUE)
summary(cl)
plotClusters(y, cl$mem, plot.dim = c(1, 1, 9), cex.points = 0.01)
#=========================intersect==============================
as<-as.data.frame(y1$X1)
row.names(as)<-y1$X
as<-as.data.frame(mc12$classification)
ah<-as.data.frame(clusterCut)
colnames(ah)<-c('ah')
colnames(as)<-c('as')
bkm<-matrix(0,nrow = max(as$as)+1 ,ncol = max(ah$ah)+1)
for (i in 1:max(as$as)) {
  for (j in 1: max(ah$ah)) {
    k1<-subset(as,as==i)
    k2<-subset(ah,ah==j)
    b<-as.data.frame(subset(k1,row.names(k1) %in% row.names(k2)))
    bkm[i+1,j+1]<-nrow(b)
  }
}
for (i in 1:max(as$as)) {
  bkm[i+1:1]<-nrow(subset(as,as==i))
}
for (j in 1:max(ah$ah)) {
  bkm[1,j+1]<-nrow(subset(ah,ah==j))
}

####Figure S1-mutational load ####

t<-var1[,c(1,6,7)]
t<-unique(t)

df<-as.data.frame(table(t$icgc_samlpe_id))
colnames(df)<-c("Samples", "Frequency")

g <- ggplot()+theme_bw()+
  geom_bar(stat = 'identity', data = df, aes(x=Samples, y=Frequency))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.05*max(df$Freq)))+
  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
plot(g)

####Figure S4 - Samples in PCAs ####
gene_sample<-read.csv('../hamed/Desktop/allllll/table/gngnmt1.csv')
gene_sample<-gene_sample[,-1]
row.names(gene_sample)<-gene_sample$Sample
gene_sample<-gene_sample[,-1]

pca<-prcomp(t(gene_sample))
pres <- as.data.frame(prcomp(t(gene_sample))[[2]])[,c(1,2)]
pres<-pres[match(as.character(y$X),rownames(pres)),]
identical(as.character(y$X), rownames(pres))
pres$subtype <- factor(y$X1)

library(tibble)
library(ggpmisc)

p <-
  ggplot(pres ,  aes(PC1,PC2), fill = color)+
  geom_point(aes(colour = subtype))

df <- tibble(x = 0.25, y = 0.9,
             plot = list(p +
                           coord_cartesian(xlim = c(-0.003, 0.003),
                                           ylim = c(-0.003, 0.005)) +
                           labs(x = NULL, y = NULL) +
                           theme_classic()))


p +
  expand_limits(x = 0, y = 0) +
  geom_plot_npc(data = df, aes(npcx = x, npcy = y, label = plot))


###### Figure S7 - difference in features load #####
cg<-read.csv('~/Desktop/allllll//Colorectal/Colorectal_gene.csv')
cg<-cg[,-1]
cg$gm<-paste0(cg$geneID,'_',cg$motif)
t<-merge(genename,cg, all=TRUE)
t<-subset(t, t$sample_id %in% cg$sample_id)
t$gene_name <- ifelse(is.na(t$gene_name), t$geneID, t$gene_name)
t$cgn<- paste0(t$gene_name,'_',t$motif)
{
  gnmtf<-t[,c(3,10,11)]
  gnmtf$gm<-gsub("-", ".", gnmtf$gm)
  gnmtf<-subset(gnmtf,gnmtf$gm %in% colnames(gene_sample))
  gnmtf<-unique(gnmtf)
  gnmtf<-gnmtf[,c(1,3)]
  colnames(gnmtf)<-c("sample_id","gm")
  
  gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
  gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
  gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
  gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
  gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
  gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
  gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))
  
  gnmtf<-as.data.frame(table(gnmtf$gm))
  gnmtf1<-as.data.frame(table(gnmtf1$gm))
  gnmtf2<-as.data.frame(table(gnmtf2$gm))
  gnmtf3<-as.data.frame(table(gnmtf3$gm))
  gnmtf4<-as.data.frame(table(gnmtf4$gm))
  gnmtf5<-as.data.frame(table(gnmtf5$gm))
  gnmtf6<-as.data.frame(table(gnmtf6$gm))
  gnmtf7<-as.data.frame(table(gnmtf7$gm))
  
  {
    gnmtf$Freq<-gnmtf$Freq/nrow(y)
    gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
    gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
    gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
    gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
    gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
    gnmtf6$Freq<-gnmtf6$Freq/nrow(y6)
    gnmtf7$Freq<-gnmtf7$Freq/nrow(y7)
  }
  {
    t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf1$Var1))
    t$Freq<-0
    colnames(t)<-c('Var1', 'Freq')
    gnmtf1<-rbind(gnmtf1,t)
    t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf2$Var1))
    t$Freq<-0
    colnames(t)<-c('Var1', 'Freq')
    gnmtf2<-rbind(gnmtf2,t)
    t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf3$Var1))
    t$Freq<-0
    colnames(t)<-c('Var1', 'Freq')
    gnmtf3<-rbind(gnmtf3,t)
    t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf4$Var1))
    t$Freq<-0
    colnames(t)<-c('Var1', 'Freq')
    gnmtf4<-rbind(gnmtf4,t)
    t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf5$Var1))
    t$Freq<-0
    colnames(t)<-c('Var1', 'Freq')
    gnmtf5<-rbind(gnmtf5,t)
    t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf6$Var1))
    t$Freq<-0
    colnames(t)<-c('Var1', 'Freq')
    gnmtf6<-rbind(gnmtf6,t)
    t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf7$Var1))
    t$Freq<-0
    colnames(t)<-c('Var1', 'Freq')
    gnmtf7<-rbind(gnmtf7,t)
    
    gnmtf<-gnmtf[order(gnmtf$Var1),]
    gnmtf1<-gnmtf1[order(as.character(gnmtf1$Var1)),]
    gnmtf2<-gnmtf2[order(as.character(gnmtf2$Var1)),]
    gnmtf3<-gnmtf3[order(as.character(gnmtf3$Var1)),]
    gnmtf4<-gnmtf4[order(as.character(gnmtf4$Var1)),]
    gnmtf5<-gnmtf5[order(as.character(gnmtf5$Var1)),]
    gnmtf6<-gnmtf6[order(as.character(gnmtf6$Var1)),]
    gnmtf7<-gnmtf7[order(as.character(gnmtf7$Var1)),]
  }
  
}

  for (i in 1:3109) {
    if(gnmtf[i,2]<0.1){
      gnmtf[i,2]=0
    }
    if(gnmtf1[i,2]<0.1){
      gnmtf1[i,2]=0
    }
    if(gnmtf2[i,2]<0.1){
      gnmtf2[i,2]=0
    }
    if(gnmtf3[i,2]<0.1){
      gnmtf3[i,2]=0
    }
    if(gnmtf4[i,2]<0.1){
      gnmtf4[i,2]=0
    }
    if(gnmtf5[i,2]<0.1){
      gnmtf5[i,2]=0
    }
    if(gnmtf6[i,2]<0.1){
      gnmtf6[i,2]=0
    }
    if(gnmtf7[i,2]<0.1){
      gnmtf7[i,2]=0
    }
  }


{
  df <- gnmtf
  df$x = 1:dim(df)[1]
  df_text <- df[c(2057,190),]
  p1 <- ggplot()+theme_bw()+
    geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
    xlab("Significant gene-motifs")+ ylab("Mutational load")+
    geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  df <- gnmtf1
  df$x = 1:dim(df)[1]
  df_text <- df[order(df$Freq,decreasing = T)[c(1,7)],]
  p2 <- ggplot()+theme_bw()+
    geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
    xlab("Significant gene-motifs")+ ylab("Mutational load")+
    geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  df <- gnmtf2
  df$x = 1:dim(df)[1]
  df_text <- df[order(df$Freq,decreasing = T)[1:2],]
  p3 <- ggplot()+theme_bw()+
    geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
    xlab("Significant gene-motifs")+ ylab("Mutational load")+
    geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=90), size=3, vjust=0.5, hjust = -0.1)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  df <- gnmtf3
  df$x = 1:dim(df)[1]
  df_text <- df[order(df$Freq,decreasing = T)[c(1,4,5)],]
  p4 <- ggplot()+theme_bw()+
    geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
    xlab("Significant gene-motifs")+ ylab("Mutational load")+
    geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  df <- gnmtf4
  df$x = 1:dim(df)[1]
  df_text <- df[order(df$Freq,decreasing = T)[c(1,9,17)],]
  p5 <- ggplot()+theme_bw()+
    geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
    xlab("Significant gene-motifs")+ ylab("Mutational load")+
    geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  df <- gnmtf5
  df$x = 1:dim(df)[1]
  df_text <- df[order(df$Freq,decreasing = T)[c(1,22,23)],]
  p6 <- ggplot()+theme_bw()+
    geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
    xlab("Significant gene-motifs")+ ylab("Mutational load")+
    geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=90), size=3, vjust=0.5, hjust = -0.1)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  df <- gnmtf6
  df$x = 1:dim(df)[1]
  df_text <- df[order(df$Freq,decreasing = T)[c(1,3,5)],]
  p7 <- ggplot()+theme_bw()+
    geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
    xlab("Significant gene-motifs")+ ylab("Mutational load")+
    geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  df <- gnmtf7
  df$x = 1:dim(df)[1]
  df_text <- df[order(df$Freq,decreasing = T)[c(1,6)],]
  p8 <- ggplot()+theme_bw()+
    geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,1.21))+
    xlab("Significant gene-motifs")+ ylab("Mutational load")+
    geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
}
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2)


##### Figure 2-b: motif rate in different gene######
library(networkD3)
library(tidyverse)
library(ggalluvial)
t<- subset(var, var$gene %in% c("ENSG00000134982","ENSG00000141510", "ENSG00000133703","ENSG00000157764",
                                "ENSG00000185008","ENSG00000188580","ENSG00000165323") )
gene_name_selected <- subset(chr1_genes, chr1_genes$ensembl_gene_id %in% c("ENSG00000134982","ENSG00000141510", "ENSG00000133703","ENSG00000157764",
                                                                           "ENSG00000185008","ENSG00000188580","ENSG00000165323"))
gnmtf1$Subtype <- "SC1"
gnmtf2$Subtype <- "SC2"
gnmtf3$Subtype <- "SC3"
gnmtf4$Subtype <- "SC4"
gnmtf5$Subtype <- "SC5"
gnmtf6$Subtype <- "SC6"
gnmtf7$Subtype <- "SC7"

gnfreq<-rbind(gnmtf1,gnmtf2,gnmtf3,gnmtf4,gnmtf5,gnmtf6,gnmtf7)
gnfreqs<- subset(gnfreq, gnfreq$Var1 %in% c("ENSG00000134982","ENSG00000141510", "ENSG00000133703","ENSG00000157764",
                                            "ENSG00000185008","ENSG00000188580","ENSG00000165323"))
colnames(gnfreqs)<-c("ensembl_gene_id","Freq", "Subtype")
gnfreq<- merge(gnfreqs,gene_name_selected)

{
  var<-subset(var1, var1$gene == "ENSG00000157764")
  gnmtf<-var[,c(1,9)]
  colnames(gnmtf)<-c('sample_id','gm')
  
  gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
  gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
  gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
  gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
  gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
  gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
  gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))
  
  gnmtf1<-as.data.frame(table(gnmtf1$gm))
  gnmtf2<-as.data.frame(table(gnmtf2$gm))
  gnmtf3<-as.data.frame(table(gnmtf3$gm))
  gnmtf4<-as.data.frame(table(gnmtf4$gm))
  gnmtf5<-as.data.frame(table(gnmtf5$gm))
  gnmtf6<-as.data.frame(table(gnmtf6$gm))
  gnmtf7<-as.data.frame(table(gnmtf7$gm))
  
  {
    gnmtf1$Freq<-gnmtf1$Freq/sum(gnmtf1$Freq)
    gnmtf2$Freq<-gnmtf2$Freq/sum(gnmtf2$Freq)
    gnmtf3$Freq<-gnmtf3$Freq/sum(gnmtf3$Freq)
    gnmtf4$Freq<-gnmtf4$Freq/sum(gnmtf4$Freq)
    gnmtf5$Freq<-gnmtf5$Freq/sum(gnmtf5$Freq)
    gnmtf6$Freq<-gnmtf6$Freq/sum(gnmtf6$Freq)
    gnmtf7$Freq<-gnmtf7$Freq/sum(gnmtf7$Freq)
  }
  
  gnmtf1$Subtype <- "SC1"
  gnmtf2$Subtype <- "SC2"
  gnmtf3$Subtype <- "SC3"
  gnmtf4$Subtype <- "SC4"
  gnmtf5$Subtype <- "SC5"
  gnmtf6$Subtype <- "SC6"
  gnmtf7$Subtype <- "SC7"
  
  gnmtf<-rbind(gnmtf1,gnmtf2,gnmtf3,gnmtf4,gnmtf5,gnmtf6,gnmtf7)
  gnmtf$Gene<-"BRAF"
}
df<-rbind(df,gnmtf)
df$real<-0
for (i in 1:2688) {
  t<-subset(gnfreq, gnfreq$external_gene_name == df[i,4] & gnfreq$Subtype== df[i,3]) 
  df[i,5] = as.numeric(df[i,2])* t[1,2]*100
}
df[is.na(df)] <- 0
colnames(df)<-c("Motif","Freq","Subtype","Gene","Value")
df1<-subset(df, df$Gene %in% c("KRAS","BRAF","TP53","APC"))
df1<-subset(df1, df1$Motif %in% c("CT-C.G","TG-T.T","TA-G.G","CT-T.G","CT-A.C","CA-T.T","CT-G.G","CT-T.G"))

ggplot(df1,
       aes(y = Value, axis1 = Subtype, axis2 = Motif)) +
  geom_alluvium(aes(fill = Gene), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("Subtype", "Motif"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") 

ggplot(df1,
       aes(y = Value, axis1 = Gene, axis2 = Motif)) +
  geom_alluvium(aes(fill = Subtype), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("Gene", "Motif"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") 


#####Figure S9- hyper mutated in protien coding and lncRNA regions######

gnmtf<-cg[,c(1,5)]
gnmtf<-unique(gnmtf)
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

lnc<-read.csv('./Desktop/allllll/Colorectal/Colorectal_lncRNA.csv')
lnc<-lnc[,-1]
gnmtf<-lnc[,c(1,5)]
gnmtf<-unique(gnmtf)
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

mut<-data.frame(Subtype = c("SC1","SC2","SC3","SC4","SC5","SC6","SC7"),
                Frequency= c(64606,14726,37372,35347,81336,68145,37427))
mut$Type<-"Protein-coding"
mut<-data.frame(Subtype = c("SC1","SC2","SC3","SC4","SC5","SC6","SC7"),
                Frequency= c(7711,1844,5734,4149,13605,12637,24565))
mut$Type<-"lncRNA"
all<-rbind(mut,all)
all$Relative<-0
for (i in 1:14) {
  t<-nrow(subset(y,y$X1==as.numeric(substr(all[i,1],3,3))))
  all[i,4]=all[i,2]/t
}
ggplot(all, aes(fill=Type, y=Relative, x=Subtype)) + 
  geom_bar( stat="identity")+
  ylab("Averge mutation rates")


##### Figure S9+S10 lncRNA and Gene rate ####
{
  gnmtf<-cg[,c(1,5)]
  gnmtf<-unique(gnmtf)
  colnames(gnmtf)<- c("sample_id","ensembl_gene_id")
  t<-merge(chr1_genes,gnmtf,all= TRUE)
  t<-subset(t, t$sample_id %in% cg$sample_id)

  t$external_gene_name <- ifelse(is.na(t$external_gene_name), t$ensembl_gene_id, t$external_gene_name)
  
  gnmtf<-t[,c(3,2)]
  
  colnames(gnmtf)<-c('sample_id','gm')
  
  gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
  gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
  gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
  gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
  gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
  gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
  gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))
  
  gnmtf<-as.data.frame(table(gnmtf$gm))
  gnmtf1<-as.data.frame(table(gnmtf1$gm))
  gnmtf2<-as.data.frame(table(gnmtf2$gm))
  gnmtf3<-as.data.frame(table(gnmtf3$gm))
  gnmtf4<-as.data.frame(table(gnmtf4$gm))
  gnmtf5<-as.data.frame(table(gnmtf5$gm))
  gnmtf6<-as.data.frame(table(gnmtf6$gm))
  gnmtf7<-as.data.frame(table(gnmtf7$gm))
  
  {
    gnmtf$Freq<-gnmtf$Freq/nrow(y)
    gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
    gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
    gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
    gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
    gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
    gnmtf6$Freq<-gnmtf6$Freq/nrow(y6)
    gnmtf7$Freq<-gnmtf7$Freq/nrow(y7)
  }
}

  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf1$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf1<-rbind(gnmtf1,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf2$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf2<-rbind(gnmtf2,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf3$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf3<-rbind(gnmtf3,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf4$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf4<-rbind(gnmtf4,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf5$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf5<-rbind(gnmtf5,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf6$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf6<-rbind(gnmtf6,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf7$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf7<-rbind(gnmtf7,t)
  
  gnmtf<-gnmtf[order(gnmtf$Var1),]
  gnmtf1<-gnmtf1[order(as.character(gnmtf1$Var1)),]
  gnmtf2<-gnmtf2[order(as.character(gnmtf2$Var1)),]
  gnmtf3<-gnmtf3[order(as.character(gnmtf3$Var1)),]
  gnmtf4<-gnmtf4[order(as.character(gnmtf4$Var1)),]
  gnmtf5<-gnmtf5[order(as.character(gnmtf5$Var1)),]
  gnmtf6<-gnmtf6[order(as.character(gnmtf6$Var1)),]
  gnmtf7<-gnmtf7[order(as.character(gnmtf7$Var1)),]


{
  for (i in 1:19805) {
    if(gnmtf[i,2]<0.5){
      gnmtf[i,2]=0
    }
    if(gnmtf1[i,2]<0.5){
      gnmtf1[i,2]=0
    }
    if(gnmtf2[i,2]<0.3){
      gnmtf2[i,2]=0
    }
    if(gnmtf3[i,2]<0.5){
      gnmtf3[i,2]=0
    }
    if(gnmtf4[i,2]<0.5){
      gnmtf4[i,2]=0
    }
    if(gnmtf5[i,2]<0.5){
      gnmtf5[i,2]=0
    }
    if(gnmtf6[i,2]<0.5){
      gnmtf6[i,2]=0
    }
    if(gnmtf7[i,2]<0.5){
      gnmtf7[i,2]=0
    }
  }
  
  {
    df <- gnmtf
    df$x = 1:dim(df)[1]
    df_text <- df[order(df$Freq,decreasing = T)[c(1,3,5)],]
    p1 <- ggplot()+theme_bw()+
      geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq),width = 10)+
      scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
      xlab("all protein-coding genes")+ ylab("Mutational load")+
      ggtitle("all CRC patients")+
      geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    df <- gnmtf1
    df$x = 1:dim(df)[1]
    df_text <- df[order(df$Freq,decreasing = T)[c(1,17,18,22)],]
    p2 <- ggplot()+theme_bw()+
      geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
      scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
      ggtitle("SC1")+
      xlab("all protein-coding genes")+ ylab("Mutational load")+
      geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=90), size=3, vjust=0.5, hjust = -0.1)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    df <- gnmtf2
    df$x = 1:dim(df)[1]
    df_text <- df[order(df$Freq,decreasing = T)[c(1:3)],]
    p3 <- ggplot()+theme_bw()+
      geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq), fill="black", width = 10)+
      scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
      ggtitle("SC2")+
      xlab("all protein-coding genes")+ ylab("Mutational load")+
      geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    df <- gnmtf3
    df$x = 1:dim(df)[1]
    df_text <- df[order(df$Freq,decreasing = T)[c(9,10,5,5,1)],]
    p4 <- ggplot()+theme_bw()+
      geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
      scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
      ggtitle("SC3")+
      xlab("all protein-coding genes")+ ylab("Mutational load")+
      geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    df <- gnmtf4
    df$x = 1:dim(df)[1]
    df_text <- df[order(df$Freq,decreasing = T)[c(1:3)],]
    p5 <- ggplot()+theme_bw()+
      geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq), fill="black", width = 10)+
      scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
      ggtitle("SC4")+
      xlab("all protein-coding genes")+ ylab("Mutational load")+
      geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    df <- gnmtf5
    df$x = 1:dim(df)[1]
    df_text <- df[order(df$Freq,decreasing = T)[c(11,2,5,25,31)],]
    p6 <- ggplot()+theme_bw()+
      geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
      scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
      ggtitle("SC5")+
      xlab("all protein-coding genes")+ ylab("Mutational load")+
      geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=45), size=3, vjust=0.5, hjust = -0.1)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    df <- gnmtf6
    df$x = 1:dim(df)[1]
    df_text <- df[order(df$Freq,decreasing = T)[c(9,1,2)],]
    p7 <- ggplot()+theme_bw()+
      geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq), fill="black", width = 5)+
      scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
      ggtitle("SC6")+
      xlab("all protein-coding genes")+ ylab("Mutational load")+
      geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=90), size=3, vjust=0.5, hjust = -0.1)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    df <- gnmtf7
    df$x = 1:dim(df)[1]
    df_text <- df[order(df$Freq,decreasing = T)[c(30,54,58,13)],]
    p8 <- ggplot()+theme_bw()+
      geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq), fill="black", width = 3)+
      scale_y_continuous(expand = c(0, 0), limits = c(0,1.2))+
      ggtitle("SC7")+
      xlab("all protein-coding genes")+ ylab("Mutational load")+
      geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=90), size=3, vjust=0.5, hjust = -0.1)+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
  }
  grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2)
}

#####Figure 2-c: transcript rate ####  
  var<-read.csv('./Desktop/allllll/Colorectal/var.csv')
  var<-var[,-1]
  var<-unique(var)
  {
    t<-var[,c(1,3,4)]
    t<-unique(t)
    colnames(t)<-c('sample_id',"gene",'gm')
    
    t<- subset(t,t$gene%in% c("ENSG00000124092"))
    
    gnmtf<-t[,c(1,3)]
    
    gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
    gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
    gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
    gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
    gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
    gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
    gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))
    
    gnmtf<-as.data.frame(table(gnmtf$gm))
    gnmtf1<-as.data.frame(table(gnmtf1$gm))
    gnmtf2<-as.data.frame(table(gnmtf2$gm))
    gnmtf3<-as.data.frame(table(gnmtf3$gm))
    gnmtf4<-as.data.frame(table(gnmtf4$gm))
    gnmtf5<-as.data.frame(table(gnmtf5$gm))
    gnmtf6<-as.data.frame(table(gnmtf6$gm))
    gnmtf7<-as.data.frame(table(gnmtf7$gm))
    
    
    {
      gnmtf<-subset(gnmtf,gnmtf$Var1 %in% t$gm)
      gnmtf1<-subset(gnmtf1,gnmtf1$Var1 %in% t$gm)
      gnmtf2<-subset(gnmtf2,gnmtf2$Var1 %in% t$gm)
      gnmtf3<-subset(gnmtf3,gnmtf3$Var1 %in% t$gm)
      gnmtf4<-subset(gnmtf4,gnmtf4$Var1 %in% t$gm)
      gnmtf5<-subset(gnmtf5,gnmtf5$Var1 %in% t$gm)
      gnmtf6<-subset(gnmtf6,gnmtf6$Var1 %in% t$gm)
      gnmtf7<-subset(gnmtf7,gnmtf7$Var1 %in% t$gm)
    }
    
    {
      gnmtf$Freq<-gnmtf$Freq/nrow(y)
      gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
      gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
      gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
      gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
      gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
      gnmtf6$Freq<-gnmtf6$Freq/nrow(y6)
      gnmtf7$Freq<-gnmtf7$Freq/nrow(y7)
      gnmtf1$subtype<-"SC1"
      gnmtf2$subtype<-"SC2"
      gnmtf3$subtype<-"SC3"
      gnmtf4$subtype<-"SC4"
      gnmtf5$subtype<-"SC5"
      gnmtf6$subtype<-"SC6"
      gnmtf7$subtype<-"SC7"
    }
    gnmtf<-rbind(gnmtf1,gnmtf2,gnmtf3,gnmtf4,gnmtf5,gnmtf6,gnmtf7)
    colnames(gnmtf)<-c("individual","value","group")
    gnmtf$value=gnmtf$value*100
    
  }
  {
    # Set a number of 'empty bar' to add at the end of each group
    data<-gnmtf[,c(1,3,2)]
    empty_bar= 5
    to_add = data.frame( matrix(NA, empty_bar*nrow(as.data.frame(unique(data$group))), ncol(data)) )
    colnames(to_add) = colnames(data)
    
    to_add$group=rep(unique(data$group), each=empty_bar)
    data=rbind(data, to_add)
    data=data %>% arrange(group)
    data$id=seq(1, nrow(data))
    
    # Get the name and the y position of each label
    label_data=data
    number_of_bar=nrow(label_data)
    angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust<-ifelse( angle < -90, 1, 0)
    label_data$angle<-ifelse(angle < -90, angle+180, angle)
    
    # prepare a data frame for base lines
    base_data=data %>% 
      group_by(group) %>% 
      summarize(start=min(id), end=max(id) - empty_bar) %>% 
      rowwise() %>% 
      mutate(title=mean(c(start, end)))
    
    # prepare a data frame for grid (scales)
    grid_data = base_data
    grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
    grid_data$start = grid_data$start - 1
    grid_data=grid_data[-1,]
    
    # Make the plot
    p1 = ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
      
      geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
      
      # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
      geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.2 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.2 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.2 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.2 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.2 , inherit.aes = FALSE ) +
      #geom_segment(data=grid_data, aes(x = end, y = 10, xend = start, yend = 10), colour = "grey", alpha=1, size=0.2 , inherit.aes = FALSE ) +
      #geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.2 , inherit.aes = FALSE ) +
      
      # Add text showing the value of each 100/75/50/25 lines
      annotate("text", x = rep(max(data$id),5), y = c(20, 40, 60, 80,100), label = c("20", "40", "60", "80","100") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
      
      geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
      ylim(-90,130) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,7), "cm") 
      ) +
      coord_polar() + 
      geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
      
      # Add base line information
      geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
      geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(0.8,1,1,1,0,0,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)+
      annotate(geom="text", x=63, y=-90, label = "MUC6", size=5, fontface="bold",
               colour = "red") 
  }
  plot(p1)
  grid.arrange(p2, p6 ,p4, p5, nrow = 2)
  grid.arrange( p6 ,p4, nrow = 1)
#========gene-motif=====
cg<-read.csv('~/Desktop/allllll//Colorectal/Colorectal_gene.csv')
cg<-cg[,-1]
cg$gm<-paste0(cg$geneID,'_',cg$motif)

y<-read.csv('./Desktop/allllll/table/cluster7.csv')
row.names(y)<-y$X
y1<-subset(y,y$X1==1)
y2<-subset(y,y$X1==2)
y3<-subset(y,y$X1==3)
y4<-subset(y,y$X1==4)
y5<-subset(y,y$X1==5)
y6<-subset(y,y$X1==6)
y7<-subset(y,y$X1==7)

gnmtf<-cg[,c(1,9)]
gnmtf$gm<-gsub("-", ".", gnmtf$gm)
gnmtf<-subset(gnmtf,gnmtf$gm %in% colnames(gene_sample))
gnmtf<-unique(gnmtf)

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))
gnmtf6<-as.data.frame(table(gnmtf6$gm))
gnmtf7<-as.data.frame(table(gnmtf7$gm))

{
  gnmtf$Freq<-gnmtf$Freq/nrow(y)
  gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
  gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
  gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
  gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
  gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
  gnmtf6$Freq<-gnmtf6$Freq/nrow(y6)
  gnmtf7$Freq<-gnmtf7$Freq/nrow(y7)
}
{
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf1$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf1<-rbind(gnmtf1,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf2$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf2<-rbind(gnmtf2,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf3$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf3<-rbind(gnmtf3,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf4$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf4<-rbind(gnmtf4,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf5$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf5<-rbind(gnmtf5,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf6$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf6<-rbind(gnmtf6,t)
  t<- as.data.frame(setdiff(gnmtf$Var1,gnmtf7$Var1))
  t$Freq<-0
  colnames(t)<-c('Var1', 'Freq')
  gnmtf7<-rbind(gnmtf7,t)
  
  gnmtf<-gnmtf[order(gnmtf$Var1),]
  gnmtf1<-gnmtf1[order(as.character(gnmtf1$Var1)),]
  gnmtf2<-gnmtf2[order(as.character(gnmtf2$Var1)),]
  gnmtf3<-gnmtf3[order(as.character(gnmtf3$Var1)),]
  gnmtf4<-gnmtf4[order(as.character(gnmtf4$Var1)),]
  gnmtf5<-gnmtf5[order(as.character(gnmtf5$Var1)),]
  gnmtf6<-gnmtf6[order(as.character(gnmtf6$Var1)),]
  gnmtf7<-gnmtf7[order(as.character(gnmtf7$Var1)),]
}

par(mfrow=c(2,4))
barplot(gnmtf$Freq, main = 'Whole patients',ylim=c(0, 1))
barplot(gnmtf1$Freq, main = 'SC1',ylim=c(0, 1),col = 'yellow',border = 'yellow')
barplot(gnmtf2$Freq, main = 'SC2',ylim=c(0, 1),col = 'blue',border = 'blue')
barplot(gnmtf3$Freq, main = 'SC3',ylim=c(0, 1),col = 'green',border = 'green')
barplot(gnmtf4$Freq, main = 'SC4',ylim=c(0, 1),col = 'orange',border = 'orange')
barplot(gnmtf5$Freq, main = 'SC5',ylim=c(0, 1),col = 'pink',border = 'pink')
barplot(gnmtf6$Freq, main = 'SC6',ylim=c(0, 1),col = 'brown',border = 'brown')
barplot(gnmtf7$Freq, main = 'SC7',ylim=c(0, 1),col = 'purple',border = 'purple')
#========gene ============
gnmtf<-cg[,c(1,5)]
gnmtf<-unique(gnmtf)
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))
gnmtf6<-as.data.frame(table(gnmtf6$gm))
gnmtf7<-as.data.frame(table(gnmtf7$gm))

{
  gnmtf$Freq<-gnmtf$Freq/nrow(y)
  gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
  gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
  gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
  gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
  gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
  gnmtf6$Freq<-gnmtf6$Freq/nrow(y6)
  gnmtf7$Freq<-gnmtf7$Freq/nrow(y7)
}


df <- gnmtf4

df$x = 1:dim(df)[1]

df_text <- df[order(df$Freq,decreasing = T)[1:3],]

g <- ggplot()+theme_bw()+
  geom_bar(stat = 'identity', data = df, aes(x=x, y=Freq))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.5*max(df$Freq)))+
  
  geom_text(data = df_text, aes(x=x, y=Freq, label=Var1, angle=90), size=4, vjust=0.5, hjust = -0.1)+
  
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
plot(g)




par(mfrow=c(2,4 )) 
barplot(gnmtf$Freq, main = 'Whole patients',ylim=c(0, 1))
barplot(gnmtf1$Freq, main = 'SC1',ylim=c(0, 1),col = 'yellow',border = 'yellow')
barplot(gnmtf2$Freq, main = 'SC2',ylim=c(0, 1),col = 'blue',border = 'blue')
barplot(gnmtf3$Freq, main = 'SC3',ylim=c(0, 1),col = 'green',border = 'green')
barplot(gnmtf4$Freq, main = 'SC4',ylim=c(0, 1),col = 'orange',border = 'orange')
barplot(gnmtf5$Freq, main = 'SC5',ylim=c(0, 1),col = 'pink',border = 'pink')
barplot(gnmtf6$Freq, main = 'SC6',ylim=c(0, 1),col = 'brown',border = 'brown')
barplot(gnmtf7$Freq, main = 'SC7',ylim=c(0, 1),col = 'purple',border = 'purple')

#=========Gene association==========
gnmtf<-cg[,c(1,5)]
gnmtf<-unique(gnmtf)
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))
gnmtf6<-as.data.frame(table(gnmtf6$gm))
gnmtf7<-as.data.frame(table(gnmtf7$gm))

gnmtf1$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq
gnmtf1$fisher<-0
t<-as.matrix(gnmtf1[,2:4])
for (i in 1:19870) {
  tt<-matrix(c(t[i,1], t[i,2], 44-t[i,1], 492-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf1$fisher<-t[,3]

gnmtf2$non<-gnmtf1$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq
gnmtf2$fisher<-0
t<-as.matrix(gnmtf2[,2:4])
for (i in 1:19870) {
  tt<-matrix(c(t[i,1], t[i,2], 137-t[i,1], 396-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf2$fisher<-t[,3]

gnmtf3$non<-gnmtf2$Freq+gnmtf1$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq
gnmtf3$fisher<-0
t<-as.matrix(gnmtf3[,2:4])
for (i in 1:19870) {
  tt<-matrix(c(t[i,1], t[i,2], 64-t[i,1], 472-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf3$fisher<-t[,3]

gnmtf4$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf1$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq
gnmtf4$fisher<-0
t<-as.matrix(gnmtf4[,2:4])
for (i in 1:19870) {
  tt<-matrix(c(t[i,1], t[i,2], 205-t[i,1], 331-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf4$fisher<-t[,3]

gnmtf5$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf1$Freq+gnmtf6$Freq+gnmtf7$Freq
gnmtf5$fisher<-0
t<-as.matrix(gnmtf5[,2:4])
for (i in 1:19870) {
  tt<-matrix(c(t[i,1], t[i,2], 13-t[i,1], 523-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf5$fisher<-t[,3]

gnmtf6$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf1$Freq+gnmtf7$Freq
gnmtf6$fisher<-0
t<-as.matrix(gnmtf6[,2:4])
for (i in 1:19870) {
  tt<-matrix(c(t[i,1], t[i,2], 49-t[i,1], 487-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf6$fisher<-t[,3]

gnmtf7$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf1$Freq
gnmtf7$fisher<-0
t<-as.matrix(gnmtf7[,2:4])
for (i in 1:19870) {
  tt<-matrix(c(t[i,1], t[i,2], 21-t[i,1], 512-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf7$fisher<-t[,3]

{
  gnmtf1<-gnmtf1[order(gnmtf1$fisher),]
  gnmtf2<-gnmtf2[order(gnmtf2$fisher),]
  gnmtf3<-gnmtf3[order(gnmtf3$fisher),]
  gnmtf4<-gnmtf4[order(gnmtf4$fisher),]
  gnmtf5<-gnmtf5[order(gnmtf5$fisher),]
  gnmtf6<-gnmtf6[order(gnmtf6$fisher),]
  gnmtf7<-gnmtf7[order(gnmtf7$fisher),]
}
{
  gnmtf1<-gnmtf1[order(-gnmtf1$Freq),]
  gnmtf2<-gnmtf2[order(-gnmtf2$Freq),]
  gnmtf3<-gnmtf3[order(-gnmtf3$Freq),]
  gnmtf4<-gnmtf4[order(-gnmtf4$Freq),]
  gnmtf5<-gnmtf5[order(-gnmtf5$Freq),]
}
{
  gnmtf1$BH<- p.adjust(gnmtf1$fisher, "BH")
  gnmtf1$bonferroni<- p.adjust(gnmtf1$fisher, "bonferroni")
  gnmtf2$BH<- p.adjust(gnmtf2$fisher, "BH")
  gnmtf2$bonferroni<- p.adjust(gnmtf2$fisher, "bonferroni")
  
  gnmtf3$BH<- p.adjust(gnmtf3$fisher, "BH")
  gnmtf3$bonferroni<- p.adjust(gnmtf3$fisher, "bonferroni")
  gnmtf4$BH<- p.adjust(gnmtf4$fisher, "BH")
  gnmtf4$bonferroni<- p.adjust(gnmtf4$fisher, "bonferroni")
  gnmtf5$BH<- p.adjust(gnmtf5$fisher, "BH")
  gnmtf5$bonferroni<- p.adjust(gnmtf5$fisher, "bonferroni")
  gnmtf6$BH<- p.adjust(gnmtf6$fisher, "BH")
  gnmtf6$bonferroni<- p.adjust(gnmtf6$fisher, "bonferroni")
  gnmtf7$BH<- p.adjust(gnmtf7$fisher, "BH")
  gnmtf7$bonferroni<- p.adjust(gnmtf7$fisher, "bonferroni")
}
{
  write.csv(gnmtf1,'./Desktop/table/gene-sc1.csv')
  write.csv(gnmtf2,'./Desktop/table/gene-sc2.csv')
  write.csv(gnmtf3,'./Desktop/table/gene-sc3.csv')
  write.csv(gnmtf4,'./Desktop/table/gene-sc4.csv')
  write.csv(gnmtf5,'./Desktop/table/gene-sc5.csv')
  write.csv(gnmtf6,'./Desktop/table/gene-sc6.csv')
  write.csv(gnmtf7,'./Desktop/table/gene-sc7.csv')
}
#--------permutation----------
p1<- matrix( ,nrow = 536,ncol = 2)
mx<-nrow(y1)
nry1<-nrow(y1)
np<- 100000
pr1<- matrix(0, nrow = np, ncol = mx)
t<-nrow(y1)+1
p1[1:nrow(y1),1]<-'case'
p1[t:536,1]<- 'control'
for (i in 1:mx) {
  p1[,2]<- runif(536,0,1)
  p1<-p1[order(p1[,2]),]
  for (j in 1:np) {
    r<- sample(1:536, 1)
    nca<- nrow(subset(p1[1:r,],p1[1:r,1]=='case'))
    nco<- nrow(subset(p1[1:r,],p1[1:r,1]=='control'))
    if(r==1){
      if(p1[1,1]=='case'){
        nca=1
        nco=0
      }
      nco=1
      nca=0
    }
    tt<-matrix(c(nca, nco, nry1-nca, 536-nry1-nco), nrow = 2)
    pr1[j,i] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
  }
}

#=======motif =========
var<-read.csv('./Desktop/allllll/Colorectal/var.csv')
var<-var[,-1]
var<-unique(var)

gnmtf<-var[,c(1,9)]
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))
gnmtf6<-as.data.frame(table(gnmtf6$gm))
gnmtf7<-as.data.frame(table(gnmtf7$gm))

{
  gnmtf$Freq<-gnmtf$Freq/sum(gnmtf$Freq)
  gnmtf1$Freq<-gnmtf1$Freq/sum(gnmtf1$Freq)
  gnmtf2$Freq<-gnmtf2$Freq/sum(gnmtf2$Freq)
  gnmtf3$Freq<-gnmtf3$Freq/sum(gnmtf3$Freq)
  gnmtf4$Freq<-gnmtf4$Freq/sum(gnmtf4$Freq)
  gnmtf5$Freq<-gnmtf5$Freq/sum(gnmtf5$Freq)
  gnmtf6$Freq<-gnmtf6$Freq/sum(gnmtf6$Freq)
  gnmtf7$Freq<-gnmtf7$Freq/sum(gnmtf7$Freq)
}

clr<-c(1:96)
clr[1:16]<-1
clr[17:32]<-2
clr[33:48]<-3
clr[49:64]<-4
clr[65:80]<-5
clr[81:96]<-6
par(mfrow=c(2,4 ))
barplot(gnmtf$Freq,col = clr, main = 'Whole patients',ylim=c(0, 0.2))
barplot(gnmtf1$Freq,col = clr, main = 'SC1',ylim=c(0, 0.2))
barplot(gnmtf2$Freq,col = clr, main = 'SC2',ylim=c(0, 0.2))
barplot(gnmtf3$Freq,col = clr, main = 'SC3',ylim=c(0, 0.2))
barplot(gnmtf4$Freq,col = clr, main = 'SC4',ylim=c(0, 0.2))
barplot(gnmtf5$Freq,col = clr, main = 'SC5',ylim=c(0, 0.2))
barplot(gnmtf6$Freq,col = clr, main = 'SC6',ylim=c(0, 0.2))
barplot(gnmtf7$Freq,col = clr, main = 'SC7',ylim=c(0, 0.2))

#---signature similarity------


signature<-read.csv('./Desktop/allllll/table/signatures_probabilities.csv')
signature<-signature[order(signature$Substitution.Type),]
t<-signature[order(signature$Substitution.Type),]
tt<-as.matrix(signature[,c(4:33)])
mysignature<-read.csv("./Desktop/allllll/table/signatures7.csv")
mysignature<-mysignature[,-1]
mysignature<-as.matrix(mysignature)

similarity<-matrix(0,nrow = 7, ncol = 30)
library(lsa)
library(proxy)
library(pheatmap)
#library(pracma)
for (i in 1:7) {
  for (j in 1:27) {
  similarity[i,j]<-simil(list(as.numeric(mysignature[,i]),as.numeric(signatures_probabilities[,j])), 
                         method="cosine")
  # t<-as.numeric(mysignature[,i])-as.numeric(signatures_probabilities[,j])
   #t<-as.matrix(t)
   #similarity[i,j]<- Norm(t,p=2)
  #similarity[i,j]<-cor.test(as.numeric(mysignature[,i]),as.numeric(tt[,j]),method = 'kendall')
  }
}
pheatmap(similarity, display_numbers = T, cluster_rows = F, cluster_cols = F,
         main = "Figure S12. Correlation between our finding signatures and Alexandrov's signatures."
         , labels_row=paste0("signature", 1:7), labels_col = paste0("Alexandrov", 1:30),cellheight = 17)

library(gplots)
heatmap.2(similarity,cellnote=similarity,
        main = "Correlation between our signatures and Alexandrov",
        Colv = NA, Rowv = NA, scale="column")
#=======transcript =========
var<-var1

gnmtf<-var[,c(1,4)]
gnmtf<-unique(gnmtf)
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))
gnmtf6<-as.data.frame(table(gnmtf6$gm))
gnmtf7<-as.data.frame(table(gnmtf7$gm))

{
  gnmtf$Freq<-gnmtf$Freq/nrow(y)
  gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
  gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
  gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
  gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
  gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
  gnmtf6$Freq<-gnmtf6$Freq/nrow(y6)
  gnmtf7$Freq<-gnmtf7$Freq/nrow(y7)
}

par(mfrow=c(2,4 ))
barplot(gnmtf$Freq, main = 'Whole patients',ylim=c(0, 1))
barplot(gnmtf1$Freq, main = 'SC1',ylim=c(0, 1),col = 'yellow',border = 'yellow')
barplot(gnmtf2$Freq, main = 'SC2',ylim=c(0, 1),col = 'blue',border = 'blue')
barplot(gnmtf3$Freq, main = 'SC3',ylim=c(0, 1),col = 'green',border = 'green')
barplot(gnmtf4$Freq, main = 'SC4',ylim=c(0, 1),col = 'orange',border = 'orange')
barplot(gnmtf5$Freq, main = 'SC5',ylim=c(0, 1),col = 'pink',border = 'pink')
barplot(gnmtf6$Freq, main = 'SC6',ylim=c(0, 1),col = 'brown',border = 'brown')
barplot(gnmtf7$Freq, main = 'SC7',ylim=c(0, 1),col = 'purple',border = 'purple')

colnames(gnmtf1)<-c("transcript","Freq")
colnames(gnmtf2)<-c("transcript","Freq")
colnames(gnmtf3)<-c("transcript","Freq")
colnames(gnmtf4)<-c("transcript","Freq")
colnames(gnmtf5)<-c("transcript","Freq")
colnames(gnmtf6)<-c("transcript","Freq")
colnames(gnmtf7)<-c("transcript","Freq")

t<-var[,c(3,4)]
t<-unique(t)
tt<-merge(t,gnmtf)
tt1<-merge(t,gnmtf1)
tt2<-merge(t,gnmtf2)
tt3<-merge(t,gnmtf3)
tt4<-merge(t,gnmtf4)
tt5<-merge(t,gnmtf5)
tt6<-merge(t,gnmtf6)
tt7<-merge(t,gnmtf7)

Transcript<-cbind(tt1,tt2$Freq,tt3$Freq,tt4$Freq,tt5$Freq,tt6$Freq,tt7$Freq)
colnames(Transcript)<-c("Transcript","ensembl_gene_id","SC1","SC2","SC3","SC4","SC5","SC6","SC7")

TN<-merge(Transcript,genename)
TN<-TN[,c(10,1:9)]
write.csv(TN,"~/Desktop/allllll/table/Transcript_mutation_rate.csv")
#=======lncRNA =======
lnc<-read.csv('./Desktop/allllll/Colorectal/Colorectal_lncRNA.csv')
lnc<-lnc[,-1]
gnmtf<-lnc[,c(1,5)]
gnmtf<-unique(gnmtf)
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))
gnmtf6<-as.data.frame(table(gnmtf6$gm))
gnmtf7<-as.data.frame(table(gnmtf7$gm))


{
  gnmtf$Freq<-gnmtf$Freq/nrow(y)
  gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
  gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
  gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
  gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
  gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
  gnmtf6$Freq<-gnmtf6$Freq/nrow(y6)
  gnmtf7$Freq<-gnmtf7$Freq/nrow(y7)
}

par(mfrow=c(2,4 ))
barplot(gnmtf$Freq, main = 'Whole patients',ylim=c(0, 1))
barplot(gnmtf1$Freq, main = 'SC1',ylim=c(0, 1),col = 'yellow',border = 'yellow')
barplot(gnmtf2$Freq, main = 'SC2',ylim=c(0, 1),col = 'blue',border = 'blue')
barplot(gnmtf3$Freq, main = 'SC3',ylim=c(0, 1),col = 'green',border = 'green')
barplot(gnmtf4$Freq, main = 'SC4',ylim=c(0, 1),col = 'orange',border = 'orange')
barplot(gnmtf5$Freq, main = 'SC5',ylim=c(0, 1),col = 'pink',border = 'pink')
barplot(gnmtf6$Freq, main = 'SC6',ylim=c(0, 1),col = 'brown',border = 'brown')
barplot(gnmtf7$Freq, main = 'SC7',ylim=c(0, 1),col = 'purple',border = 'purple')

#========consequnce=======
gnmtf<-var[,c(1,5)]
gnmtf<-unique(gnmtf)
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))
gnmtf6<-as.data.frame(table(gnmtf6$gm))
gnmtf7<-as.data.frame(table(gnmtf7$gm))

{
  gnmtf$Freq<-gnmtf$Freq/sum(gnmtf$Freq)
  gnmtf1$Freq<-gnmtf1$Freq/sum(gnmtf1$Freq)
  gnmtf2$Freq<-gnmtf2$Freq/sum(gnmtf2$Freq)
  gnmtf3$Freq<-gnmtf3$Freq/sum(gnmtf3$Freq)
  gnmtf4$Freq<-gnmtf4$Freq/sum(gnmtf4$Freq)
  gnmtf5$Freq<-gnmtf5$Freq/sum(gnmtf5$Freq)
  gnmtf6$Freq<-gnmtf6$Freq/sum(gnmtf6$Freq)
  gnmtf7$Freq<-gnmtf7$Freq/sum(gnmtf7$Freq)
}
par(mfrow=c(2,4 ))
{
  mids <- barplot(gnmtf$Freq,col = c(1:19), main = 'Whole patients',ylim=c(0, 1))
  axis(1, at=mids, labels=c('3P-UTR','5P-UTR_S','5P-UTR','DG','exon','initiator','intergenic_region','intragenic_variant','intron','missense','SA','SD','SR','start_lost','stop_gained','stop_lost','stop_retained','synonymous','UG'), las=3)
  mids <- barplot(gnmtf1$Freq,col = c(1:19), main = 'SC1',ylim=c(0, 1))
  axis(1, at=mids, labels=c('3P-UTR','5P-UTR_S','5P-UTR','DG','exon','initiator','intergenic_region','intragenic_variant','intron','missense','SA','SD','SR','start_lost','stop_gained','stop_lost','stop_retained','synonymous','UG'), las=3)
  mids <- barplot(gnmtf2$Freq,col = c(1:19), main = 'SC2',ylim=c(0, 1))
  axis(1, at=mids, labels=c('3P-UTR','5P-UTR_S','5P-UTR','DG','exon','initiator','intergenic_region','intragenic_variant','intron','missense','SA','SD','SR','start_lost','stop_gained','stop_lost','stop_retained','synonymous','UG'), las=3)
  mids <- barplot(gnmtf3$Freq,col = c(1:19), main = 'SC3',ylim=c(0, 1))
  axis(1, at=mids, labels=c('3P-UTR','5P-UTR_S','5P-UTR','DG','exon','initiator','intergenic_region','intragenic_variant','intron','missense','SA','SD','SR','start_lost','stop_gained','stop_lost','stop_retained','synonymous','UG'), las=3)
  mids <- barplot(gnmtf4$Freq,col = c(1:19), main = 'SC4',ylim=c(0, 1))
  axis(1, at=mids, labels=c('3P-UTR','5P-UTR_S','5P-UTR','DG','exon','initiator','intergenic_region','intragenic_variant','intron','missense','SA','SD','SR','start_lost','stop_gained','stop_lost','stop_retained','synonymous','UG'), las=3)
  mids <- barplot(gnmtf5$Freq,col = c(1:19), main = 'SC5',ylim=c(0, 1))
  axis(1, at=mids, labels=c('3P-UTR','5P-UTR_S','5P-UTR','DG','exon','initiator','intergenic_region','intragenic_variant','intron','missense','SA','SD','SR','start_lost','stop_gained','stop_lost','stop_retained','synonymous','UG'), las=3)
  mids <- barplot(gnmtf6$Freq,col = c(1:19), main = 'SC6',ylim=c(0, 1))
  axis(1, at=mids, labels=c('3P-UTR','5P-UTR_S','5P-UTR','DG','exon','initiator','intergenic_region','intragenic_variant','intron','missense','SA','SD','SR','start_lost','stop_gained','stop_lost','stop_retained','synonymous','UG'), las=3)
  mids <- barplot(gnmtf7$Freq,col = c(1:19), main = 'SC7',ylim=c(0, 1))
  axis(1, at=mids, labels=c('3P-UTR','5P-UTR_S','5P-UTR','DG','exon','initiator','intergenic_region','intragenic_variant','intron','missense','SA','SD','SR','start_lost','stop_gained','stop_lost','stop_retained','synonymous','UG'), las=3)
  
}

#==========chromosome =========
gnmtf<-var[,c(1,6)]
gnmtf<-unique(gnmtf)
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))

{
  gnmtf$Freq<-gnmtf$Freq/nrow(y)
  gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
  gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
  gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
  gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
  gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
}
{
  gnmtf$Freq<-gnmtf$Freq/sum(gnmtf$Freq)
  gnmtf1$Freq<-gnmtf1$Freq/sum(gnmtf1$Freq)
  gnmtf2$Freq<-gnmtf2$Freq/sum(gnmtf2$Freq)
  gnmtf3$Freq<-gnmtf3$Freq/sum(gnmtf3$Freq)
  gnmtf4$Freq<-gnmtf4$Freq/sum(gnmtf4$Freq)
  gnmtf5$Freq<-gnmtf5$Freq/sum(gnmtf5$Freq)
}

par(mfrow=c(2,3 ))
barplot(gnmtf$Freq, main = 'Whole patients',ylim=c(0, 1))
barplot(gnmtf1$Freq, main = 'SC1',ylim=c(0, 1))
barplot(gnmtf2$Freq, main = 'SC2',ylim=c(0, 1))
barplot(gnmtf3$Freq, main = 'SC3',ylim=c(0, 1))
barplot(gnmtf4$Freq, main = 'SC4',ylim=c(0, 1))
barplot(gnmtf5$Freq, main = 'SC5',ylim=c(0, 1))

barplot(gnmtf$Freq,col = c(1:24), main = 'Whole patients',ylim=c(0, 0.2))
barplot(gnmtf1$Freq,col = c(1:24), main = 'SC1',ylim=c(0, 0.2))
barplot(gnmtf2$Freq,col = c(1:24), main = 'SC2',ylim=c(0, 0.2))
barplot(gnmtf3$Freq,col = c(1:24), main = 'SC3',ylim=c(0, 0.2))
barplot(gnmtf4$Freq,col = c(1:24), main = 'SC4',ylim=c(0, 0.2))
barplot(gnmtf5$Freq,col = c(1:24), main = 'SC5',ylim=c(0, 0.2))

#========= project ========
sample<-read.csv('./Desktop/allllll/Colorectal/sample.tsv', sep = '\t', header = TRUE)
specimen<-read.csv('./Desktop/allllll/Colorectal/specimen.tsv', sep = '\t', header = TRUE)
donor<-read.csv('./Desktop/allllll/Colorectal/donor.tsv', sep = '\t', header = TRUE)

gnmtf<-subset(sample,sample$icgc_sample_id %in% y$X)
gnmtf1<-subset(sample,sample$icgc_sample_id %in% y1$X)
gnmtf2<-subset(sample,sample$icgc_sample_id %in% y2$X)
gnmtf3<-subset(sample,sample$icgc_sample_id %in% y3$X)
gnmtf4<-subset(sample,sample$icgc_sample_id %in% y4$X)
gnmtf5<-subset(sample,sample$icgc_sample_id %in% y5$X)
gnmtf6<-subset(sample,sample$icgc_sample_id %in% y6$X)
gnmtf7<-subset(sample,sample$icgc_sample_id %in% y7$X)

gnmtf<-as.data.frame(table(gnmtf$))
gnmtf1<-as.data.frame(table(gnmtf1$project_code))
gnmtf2<-as.data.frame(table(gnmtf2$project_code))
gnmtf3<-as.data.frame(table(gnmtf3$project_code))
gnmtf4<-as.data.frame(table(gnmtf4$project_code))
gnmtf5<-as.data.frame(table(gnmtf5$project_code))
gnmtf6<-as.data.frame(table(gnmtf6$project_code))
gnmtf7<-as.data.frame(table(gnmtf7$project_code))

{
  gnmtf$Freq<-gnmtf$Freq/nrow(y)
  gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
  gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
  gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
  gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
  gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
  gnmtf6$Freq<-gnmtf6$Freq/nrow(y6)
  gnmtf7$Freq<-gnmtf7$Freq/nrow(y7)
}
t<-gnmtf$Freq
{
  gnmtf$Freq<-gnmtf$Freq/t
  gnmtf1$Freq<-gnmtf1$Freq/t
  gnmtf2$Freq<-gnmtf2$Freq/t
  gnmtf3$Freq<-gnmtf3$Freq/t
  gnmtf4$Freq<-gnmtf4$Freq/t
  gnmtf5$Freq<-gnmtf5$Freq/t
  gnmtf6$Freq<-gnmtf6$Freq/t
  gnmtf7$Freq<-gnmtf7$Freq/t
}
{
  gnmtf$Freq<-gnmtf$Freq/sum(gnmtf$Freq)
  gnmtf1$Freq<-gnmtf1$Freq/sum(gnmtf1$Freq)
  gnmtf2$Freq<-gnmtf2$Freq/sum(gnmtf2$Freq)
  gnmtf3$Freq<-gnmtf3$Freq/sum(gnmtf3$Freq)
  gnmtf4$Freq<-gnmtf4$Freq/sum(gnmtf4$Freq)
  gnmtf5$Freq<-gnmtf5$Freq/sum(gnmtf5$Freq)
  gnmtf6$Freq<-gnmtf6$Freq/sum(gnmtf6$Freq)
  gnmtf7$Freq<-gnmtf7$Freq/sum(gnmtf7$Freq)
}
gnmtf1$subtype<-"SC1"
gnmtf2$subtype<-"SC2"
gnmtf3$subtype<-"SC3"
gnmtf4$subtype<-"SC4"
gnmtf5$subtype<-"SC5"
gnmtf6$subtype<-"SC6"
gnmtf7$subtype<-"SC7"
gnmtf<-rbind(gnmtf1,gnmtf2,gnmtf3,gnmtf4,gnmtf5,gnmtf6,gnmtf7)

colnames(gnmtf)<- c("Project","Percentage","Subtype")
pr<-ggplot(gnmtf, aes(fill=Project, y=Percentage, x=Subtype)) + 
  geom_bar( stat="identity", position="fill")

colnames(gnmtf)<- c("Gender","Percentage","Subtype")
sex<-ggplot(gnmtf, aes(fill=Gender, y=Percentage, x=Subtype)) + 
  geom_bar( stat="identity", position="fill")

ggarrange(sex,pr,nrow = 2)

par(mfrow=c(2,4 ))

mids <- barplot(gnmtf$Freq, xlab="",col = c(1:3),ylim=c(0, 1), main = 'Whole patients')
axis(1, at=mids, labels=c('COAD-US','COCA-CN','READ-US'), las=3)
mids <- barplot(gnmtf1$Freq, xlab="",col = c(1:3),ylim=c(0, 1), main = 'SC1')
axis(1, at=mids, labels=c('COAD-US','COCA-CN','READ-US'), las=3)
mids <- barplot(gnmtf2$Freq, xlab="",col = c(1:3),ylim=c(0, 1), main = 'SC2')
axis(1, at=mids, labels=c('COAD-US','COCA-CN','READ-US'), las=3)
mids <- barplot(gnmtf3$Freq, xlab="",col = c(1:3),ylim=c(0, 1), main = 'SC3')
axis(1, at=mids, labels=c('COAD-US','COCA-CN','READ-US'), las=3)
mids <- barplot(gnmtf4$Freq, xlab="",col = c(1:3),ylim=c(0, 1), main = 'SC4')
axis(1, at=mids, labels=c('COAD-US','COCA-CN','READ-US'), las=3)
mids <- barplot(gnmtf5$Freq, xlab="",col = c(1:3),ylim=c(0, 1), main = 'SC5')
axis(1, at=mids, labels=c('COAD-US','COCA-CN','READ-US'), las=3)
mids <- barplot(gnmtf6$Freq, xlab="",col = c(1:3),ylim=c(0, 1), main = 'SC6')
axis(1, at=mids, labels=c('COAD-US','COCA-CN','READ-US'), las=3)
mids <- barplot(gnmtf7$Freq, xlab="",col = c(1:3),ylim=c(0, 1), main = 'SC7')
axis(1, at=mids, labels=c('COAD-US','COCA-CN','READ-US'), las=3)

#=======furtur =======
d<-as.data.frame(gnmtf$icgc_donor_id)
d1<-as.data.frame(gnmtf1$icgc_donor_id)
d2<-as.data.frame(gnmtf2$icgc_donor_id)
d3<-as.data.frame(gnmtf3$icgc_donor_id)
d4<-as.data.frame(gnmtf4$icgc_donor_id)
d5<-as.data.frame(gnmtf5$icgc_donor_id)
d6<-as.data.frame(gnmtf6$icgc_donor_id)
d7<-as.data.frame(gnmtf7$icgc_donor_id)
s<-as.data.frame(gnmtf$icgc_specimen_id)
s1<-as.data.frame(gnmtf1$icgc_specimen_id)
s2<-as.data.frame(gnmtf2$icgc_specimen_id)
s3<-as.data.frame(gnmtf3$icgc_specimen_id)
s4<-as.data.frame(gnmtf4$icgc_specimen_id)
s5<-as.data.frame(gnmtf5$icgc_specimen_id)
s6<-as.data.frame(gnmtf6$icgc_specimen_id)
s7<-as.data.frame(gnmtf7$icgc_specimen_id)

gnmtf<-subset(donor, donor$icgc_donor_id %in% d$`gnmtf$icgc_donor_id`)
gnmtf1<-subset(donor, donor$icgc_donor_id  %in% d1$`gnmtf1$icgc_donor_id`)
gnmtf2<-subset(donor, donor$icgc_donor_id  %in% d2$`gnmtf2$icgc_donor_id`)
gnmtf3<-subset(donor, donor$icgc_donor_id  %in% d3$`gnmtf3$icgc_donor_id`)
gnmtf4<-subset(donor, donor$icgc_donor_id  %in% d4$`gnmtf4$icgc_donor_id`)
gnmtf5<-subset(donor, donor$icgc_donor_id  %in% d5$`gnmtf5$icgc_donor_id`)
gnmtf6<-subset(donor, donor$icgc_donor_id  %in% d6$`gnmtf6$icgc_donor_id`)
gnmtf7<-subset(donor, donor$icgc_donor_id  %in% d7$`gnmtf7$icgc_donor_id`)

gnmtf<-subset(specimen, specimen$icgc_specimen_id %in% s$`gnmtf$icgc_specimen_id`)
gnmtf1<-subset(specimen, specimen$icgc_specimen_id  %in% s1$`gnmtf1$icgc_specimen_id`)
gnmtf2<-subset(specimen, specimen$icgc_specimen_id  %in% s2$`gnmtf2$icgc_specimen_id`)
gnmtf3<-subset(specimen, specimen$icgc_specimen_id  %in% s3$`gnmtf3$icgc_specimen_id`)
gnmtf4<-subset(specimen, specimen$icgc_specimen_id  %in% s4$`gnmtf4$icgc_specimen_id`)
gnmtf5<-subset(specimen, specimen$icgc_specimen_id  %in% s5$`gnmtf5$icgc_specimen_id`)
gnmtf6<-subset(specimen, specimen$icgc_specimen_id  %in% s6$`gnmtf6$icgc_specimen_id`)
gnmtf7<-subset(specimen, specimen$icgc_specimen_id  %in% s7$`gnmtf7$icgc_specimen_id`)

gnmtf<-as.data.frame(table(gnmtf$donor_age_at_diagnosis))
gnmtf1<-as.data.frame(table(gnmtf1$donor_age_at_diagnosis))
gnmtf2<-as.data.frame(table(gnmtf2$donor_age_at_diagnosis))
gnmtf3<-as.data.frame(table(gnmtf3$donor_age_at_diagnosis))
gnmtf4<-as.data.frame(table(gnmtf4$donor_age_at_diagnosis))
gnmtf5<-as.data.frame(table(gnmtf5$donor_age_at_diagnosis))
gnmtf6<-as.data.frame(table(gnmtf6$donor_age_at_diagnosis))
gnmtf7<-as.data.frame(table(gnmtf7$donor_age_at_diagnosis))
{
  gnmtf<-gnmtf[-1,]
  gnmtf1<-gnmtf1[-1,]
  gnmtf2<-gnmtf2[-1,]
  gnmtf3<-gnmtf3[-1,]
  gnmtf4<-gnmtf4[-1,]
  gnmtf5<-gnmtf5[-1,]
}
{
  gnmtf$Freq<-gnmtf$Freq/sum(gnmtf$Freq)
  gnmtf1$Freq<-gnmtf1$Freq/sum(gnmtf1$Freq)
  gnmtf2$Freq<-gnmtf2$Freq/sum(gnmtf2$Freq)
  gnmtf3$Freq<-gnmtf3$Freq/sum(gnmtf3$Freq)
  gnmtf4$Freq<-gnmtf4$Freq/sum(gnmtf4$Freq)
  gnmtf5$Freq<-gnmtf5$Freq/sum(gnmtf5$Freq)
  gnmtf6$Freq<-gnmtf6$Freq/sum(gnmtf6$Freq)
  gnmtf7$Freq<-gnmtf7$Freq/sum(gnmtf7$Freq)
}
barplot(gnmtf$Freq,col = c(1:43), main = 'Whole patients',ylim=c(0, 0.5))
barplot(gnmtf1$Freq,col = c(1:43), main = 'SC1',ylim=c(0, 0.5))
barplot(gnmtf2$Freq,col = c(1:43), main = 'SC2',ylim=c(0, 0.5))
barplot(gnmtf3$Freq,col = c(1:43), main = 'SC3',ylim=c(0, 0.5))
barplot(gnmtf4$Freq,col = c(1:43), main = 'SC4',ylim=c(0, 0.5))
barplot(gnmtf5$Freq,col = c(1:43), main = 'SC5',ylim=c(0, 0.5))

barplot(gnmtf$Freq,col = c(1:24), main = 'Whole patients',ylim=c(0, 1))
barplot(gnmtf1$Freq,col = c(1:24), main = 'SC1',ylim=c(0, 1))
barplot(gnmtf2$Freq,col = c(1:24), main = 'SC2',ylim=c(0, 1))
barplot(gnmtf3$Freq,col = c(1:24), main = 'SC3',ylim=c(0, 1))
barplot(gnmtf4$Freq,col = c(1:24), main = 'SC4',ylim=c(0, 1))
barplot(gnmtf5$Freq,col = c(1:24), main = 'SC5',ylim=c(0, 1))


mids <- barplot(gnmtf$Freq,col = c(1:43), main = 'Whole patients',ylim=c(0, 0.5))
axis(1, at=mids, labels=c('Metastatic','blood derived','other','solid tissue',' tissue adjacent to primary','Primary tumour','Recurrent tumour'), las=3)
mids <- barplot(gnmtf1$Freq,col = c(1:43), main = 'SC1',ylim=c(0, 0.5))
axis(1, at=mids, labels=c('Metastatic','blood derived','other','solid tissue',' tissue adjacent to primary','Primary tumour','Recurrent tumour'), las=3)
mids <- barplot(gnmtf2$Freq,col = c(1:43), main = 'SC2',ylim=c(0, 0.5))
axis(1, at=mids, labels=c('Metastatic','blood derived','other','solid tissue',' tissue adjacent to primary','Primary tumour','Recurrent tumour'), las=3)
mids <- barplot(gnmtf3$Freq,col = c(1:43), main = 'SC3',ylim=c(0, 0.5))
axis(1, at=mids, labels=c('Metastatic','blood derived','other','solid tissue',' tissue adjacent to primary','Primary tumour','Recurrent tumour'), las=3)
mids <- barplot(gnmtf4$Freq,col = c(1:43), main = 'SC4',ylim=c(0, 0.5))
axis(1, at=mids, labels=c('Metastatic','blood derived','other','solid tissue',' tissue adjacent to primary','Primary tumour','Recurrent tumour'), las=3)
mids <- barplot(gnmtf5$Freq,col = c(1:43), main = 'SC5',ylim=c(0, 0.5))
axis(1, at=mids, labels=c('Metastatic','blood derived','other','solid tissue',' tissue adjacent to primary','Primary tumour','Recurrent tumour'), las=3)

mids <- barplot(gnmtf$Freq,col = c(1:43), main = 'Whole patients',ylim=c(0, 1))
axis(1, at=mids, labels=c('female','male'), las=3)
mids <- barplot(gnmtf1$Freq,col = c(1:43), main = 'SC1',ylim=c(0, 1))
axis(1, at=mids, labels=c('female','male'), las=3)
mids <- barplot(gnmtf2$Freq,col = c(1:43), main = 'SC2',ylim=c(0, 1))
axis(1, at=mids, labels=c('female','male'), las=3)
mids <- barplot(gnmtf3$Freq,col = c(1:43), main = 'SC3',ylim=c(0, 1))
axis(1, at=mids, labels=c('female','male'), las=3)
mids <- barplot(gnmtf4$Freq,col = c(1:43), main = 'SC4',ylim=c(0, 1))
axis(1, at=mids, labels=c('female','male'), las=3)
mids <- barplot(gnmtf5$Freq,col = c(1:43), main = 'SC5',ylim=c(0, 1))
axis(1, at=mids, labels=c('female','male'), las=3)
mids <- barplot(gnmtf6$Freq,col = c(1:43), main = 'SC6',ylim=c(0, 1))
axis(1, at=mids, labels=c('female','male'), las=3)
mids <- barplot(gnmtf7$Freq,col = c(1:43), main = 'SC7',ylim=c(0, 1))
axis(1, at=mids, labels=c('female','male'), las=3)

mids <- barplot(gnmtf$Freq,col = c(1:43), main = 'Whole patients',ylim=c(0, 1))
axis(1, at=mids, labels=c('alive','deceased'), las=3)
mids <- barplot(gnmtf1$Freq,col = c(1:43), main = 'SC1',ylim=c(0, 1))
axis(1, at=mids, labels=c('alive','deceased'), las=3)
mids <- barplot(gnmtf2$Freq,col = c(1:43), main = 'SC2',ylim=c(0, 1))
axis(1, at=mids, labels=c('alive','deceased'), las=3)
mids <- barplot(gnmtf3$Freq,col = c(1:43), main = 'SC3',ylim=c(0, 1))
axis(1, at=mids, labels=c('alive','deceased'), las=3)
mids <- barplot(gnmtf4$Freq,col = c(1:43), main = 'SC4',ylim=c(0, 1))
axis(1, at=mids, labels=c('alive','deceased'), las=3)
mids <- barplot(gnmtf5$Freq,col = c(1:43), main = 'SC5',ylim=c(0, 1))
axis(1, at=mids, labels=c('alive','deceased'), las=3)

par(mfrow=c(2,4))
hist(gnmtf$donor_age_at_diagnosis,breaks=seq(26,90,1),prob=TRUE, main="Whole patients",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="")
hist(gnmtf1$donor_age_at_diagnosis,breaks=seq(26,90,1),prob=TRUE, main="SC1",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="")
hist(gnmtf2$donor_age_at_diagnosis,breaks=seq(26,90,1),prob=TRUE, main="SC2",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="")
hist(gnmtf3$donor_age_at_diagnosis,breaks=seq(26,90,1),prob=TRUE, main="SC3",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="")
hist(gnmtf4$donor_age_at_diagnosis,breaks=seq(26,90,1),prob=TRUE, main="SC4",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="age")
hist(gnmtf5$donor_age_at_diagnosis,breaks=seq(26,90,1),prob=TRUE, main="SC5",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="age")
hist(gnmtf6$donor_age_at_diagnosis,breaks=seq(26,90,1),prob=TRUE, main="SC6",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="age")
hist(gnmtf7$donor_age_at_diagnosis,breaks=seq(26,90,1),prob=TRUE, main="SC7",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="age")


g = gnmtf$donor_age_at_diagnosis
m<-mean(g)
std<-sqrt(var(g))
hist(g,breaks=seq(26,90,1),prob=TRUE, main="Whole patients",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = median(g), col="green", lwd=3, lty=2)
abline(v = mean(g), col="red", lwd=3, lty=2)
g = gnmtf1$donor_age_at_diagnosis
m<-mean(g)
std<-sqrt(var(g))
hist(g,breaks=seq(26,90,1),prob=TRUE, main="SC1",xlim = c(25,91),ylim = c(0,0.15),ylab="",xlab="")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = median(g), col="green", lwd=3, lty=2)
abline(v = mean(g), col="red", lwd=3, lty=2)
g = gnmtf2$donor_age_at_diagnosis
m<-mean(g)
std<-sqrt(var(g))
hist(g,breaks=seq(26,90,1),prob=TRUE, main="SC2",xlim = c(25,91),ylim = c(0,0.15),ylab="",xlab="")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = median(g), col="green", lwd=3, lty=2)
abline(v = mean(g), col="red", lwd=3, lty=2)
g = gnmtf3$donor_age_at_diagnosis
m<-mean(g)
std<-sqrt(var(g))
hist(g,breaks=seq(26,90,1),prob=TRUE, main="SC3",xlim = c(25,91),ylim = c(0,0.15),ylab="",xlab="")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = median(g), col="green", lwd=3, lty=2)
abline(v = mean(g), col="red", lwd=3, lty=2)
g = gnmtf4$donor_age_at_diagnosis
m<-mean(g)
std<-sqrt(var(g))
hist(g,breaks=seq(26,90,1),prob=TRUE, main="SC4",xlim = c(25,91),ylim = c(0,0.15),ylab="Probability",xlab="age")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = median(g), col="green", lwd=3, lty=2)
abline(v = mean(g), col="red", lwd=3, lty=2)
g = gnmtf5$donor_age_at_diagnosis
m<-mean(g)
std<-sqrt(var(g))
hist(g,breaks=seq(26,90,1),prob=TRUE, main="SC5",xlim = c(25,91),ylim = c(0,0.15),ylab="",xlab="age")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = median(g), col="green", lwd=3, lty=2)
abline(v = mean(g), col="red", lwd=3, lty=2)
g = gnmtf6$donor_age_at_diagnosis
m<-mean(g)
std<-sqrt(var(g))
hist(g,breaks=seq(26,90,1),prob=TRUE, main="SC6",xlim = c(25,91),ylim = c(0,0.15),ylab="",xlab="age")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = median(g), col="green", lwd=3, lty=2)
abline(v = mean(g), col="red", lwd=3, lty=2)
g = gnmtf7$donor_age_at_diagnosis
m<-mean(g)
std<-sqrt(var(g))
hist(g,breaks=seq(26,90,1),prob=TRUE, main="SC7",xlim = c(25,91),ylim = c(0,0.15),ylab="",xlab="age")
curve(dnorm(x, mean=m, sd=std), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
#abline(v = median(g), col="green", lwd=3, lty=2)
abline(v = mean(g), col="red", lwd=3, lty=2)
#===========Transcript factor=====================
cg<-read.csv('../soudbeh/Desktop/Colorectal/Colorectal_gene.csv')
cg<-cg[,-1]
cg$gm<-paste0(cg$geneID,'_',cg$motif)
#tf <- as.data.frame(read.table("./Desktop/wgEncodeRegTfbsClusteredV3.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
tf<-read.csv('./Desktop/Autism/TF.csv')
tf<-as.data.frame(read.table('./Desktop/Autism/TF.csv', ))
tf<-tf[,c(1:4)]
colnames(tf)<-c('chromosome','start','end','TF')
df1<-cg[,c(1,2,3,4)]
df2<-tf[,c(4,1,2,3)]
df1<-unique(df1)
df2<-unique(df2)
tf<-df2

for (i in 1:4380444) {
  df2<-tf[i,]
  factor<- df1 %>% inner_join(df2, "chromosome") %>% 
    mutate(n = if_else(position >= start & position <= end, 1, 0)) 
  factor<-subset(factor,factor$n!=0)
  write.table(factor, "./Desktop/TF.csv", sep = ",", col.names = F, append = T)
}



gtf<-read.csv('./Desktop/TF.csv', header=T,fill=T,  col.names=c('number',"sample_id","motif","chromosome",'position','TF','start','end','n'),row.names=NULL)


gnmtf<-gtf[,-1]
gnmtf<-unique(gnmtf)
gnmtf<-gtf[,c(2,3,4,6)]
gnmtf$n<-paste0(gnmtf$motif,'_',gnmtf$TF)
gnmtf$n<-paste0(gnmtf$chromosome,'_',gnmtf$TF)
gnmtf<-unique(gnmtf)
gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))

gnmtf<-as.data.frame(table(gnmtf$TF))
gnmtf1<-as.data.frame(table(gnmtf1$TF))
gnmtf2<-as.data.frame(table(gnmtf2$TF))
gnmtf3<-as.data.frame(table(gnmtf3$TF))
gnmtf4<-as.data.frame(table(gnmtf4$TF))
gnmtf5<-as.data.frame(table(gnmtf5$TF))

gnmtf<-as.data.frame(table(gnmtf$n))
gnmtf1<-as.data.frame(table(gnmtf1$n))
gnmtf2<-as.data.frame(table(gnmtf2$n))
gnmtf3<-as.data.frame(table(gnmtf3$n))
gnmtf4<-as.data.frame(table(gnmtf4$n))
gnmtf5<-as.data.frame(table(gnmtf5$n))

gnmtf<-as.data.frame(table(gnmtf$motif))
gnmtf1<-as.data.frame(table(gnmtf1$motif))
gnmtf2<-as.data.frame(table(gnmtf2$motif))
gnmtf3<-as.data.frame(table(gnmtf3$motif))
gnmtf4<-as.data.frame(table(gnmtf4$motif))
gnmtf5<-as.data.frame(table(gnmtf5$motif))

{
  gnmtf$Freq<-gnmtf$Freq/sum(gnmtf$Freq)
  gnmtf1$Freq<-gnmtf1$Freq/sum(gnmtf1$Freq)
  gnmtf2$Freq<-gnmtf2$Freq/sum(gnmtf2$Freq)
  gnmtf3$Freq<-gnmtf3$Freq/sum(gnmtf3$Freq)
  gnmtf4$Freq<-gnmtf4$Freq/sum(gnmtf4$Freq)
  gnmtf5$Freq<-gnmtf5$Freq/sum(gnmtf5$Freq)
}
par(mfrow=c(2,3 ))
barplot(gnmtf$Freq, col = c(1:12),main = 'Whole patients',ylim=c(0, 0.1))
barplot(gnmtf1$Freq, col = c(1:12),main = 'SC1',ylim=c(0, 0.1))
barplot(gnmtf2$Freq, col = c(1:12),main = 'SC2',ylim=c(0, 0.1))
barplot(gnmtf3$Freq, col = c(1:12),main = 'SC3',ylim=c(0, 0.1))
barplot(gnmtf4$Freq, col = c(1:12),main = 'SC4',ylim=c(0, 0.1))
barplot(gnmtf5$Freq, col = c(1:12),main = 'SC5',ylim=c(0, 0.1))
#================freq ============
freq<- data.frame(subtype=c('SC1','SC2','SC3','SC4','SC5','SC6','SC7'),all=c(1:7),gene=c(1:7),lncRNA = c(1:7))
t<-c(nrow(gnmtf1),nrow(gnmtf2),nrow(gnmtf3),nrow(gnmtf4),nrow(gnmtf5),nrow(gnmtf6),nrow(gnmtf7))
freq$all<-t

#================ Gene associate significant ============

g1<- read.csv('./Desktop/table/gene-sc1.csv')
g2<- read.csv('./Desktop/table/gene-sc2.csv')
g3<- read.csv('./Desktop/table/gene-sc3.csv')
g4<- read.csv('./Desktop/table/gene-sc4.csv')
g5<- read.csv('./Desktop/table/gene-sc5.csv')
g6<- read.csv('./Desktop/table/gene-sc6.csv')
g7<- read.csv('./Desktop/table/gene-sc7.csv')

tg1<-subset(g1,g1$Var1 %in% sgngm$Var1)
tg2<-subset(g2,g2$Var1 %in% sgngm$Var1)
tg3<-subset(g3,g3$Var1 %in% sgngm$Var1)
tg4<-subset(g4,g4$Var1 %in% sgngm$Var1)
tg5<-subset(g5,g5$Var1 %in% sgngm$Var1)
tg6<-subset(g6,g6$Var1 %in% sgngm$Var1)
tg7<-subset(g7,g7$Var1 %in% sgngm$Var1)

{
  write.csv(tg1,'./Desktop/table/gene-sc1-sig.csv')
  write.csv(tg2,'./Desktop/table/gene-sc2-sig.csv')
  write.csv(tg3,'./Desktop/table/gene-sc3-sig.csv')
  write.csv(tg4,'./Desktop/table/gene-sc4-sig.csv')
  write.csv(tg5,'./Desktop/table/gene-sc5-sig.csv')
  write.csv(tg6,'./Desktop/table/gene-sc6-sig.csv')
  write.csv(tg7,'./Desktop/table/gene-sc7-sig.csv')
}
#================gene ontology==========
go1<-read.csv('./Desktop/allllll/table/Ontolgy/50 top/1/enrichment_results_wg_result1541682632.txt',sep = '\t')
go2<-read.csv('./Desktop/allllll/table/Ontolgy/50 top/2/enrichment_results_wg_result1541682662.txt',sep = '\t')
go3<-read.csv('./Desktop/allllll/table/Ontolgy/50 top/3/enrichment_results_wg_result1541682689.txt',sep = '\t')
go4<-read.csv('./Desktop/allllll/table/Ontolgy/50 top/4/enrichment_results_wg_result1541682709.txt',sep = '\t')
go5<-read.csv('./Desktop/allllll/table/Ontolgy/50 top/5/enrichment_results_wg_result1541682736.txt',sep = '\t')
go7<-read.csv('./Desktop/allllll/table/Ontolgy/50 top/7/enrichment_results_wg_result1541682774.txt',sep = '\t')

go1$Subtype <-'SC1'
go2$Subtype <-'SC2'
go3$Subtype <-'SC3'
go4$Subtype <-'SC4'
go5$Subtype <-'SC5'
go7$Subtype <-'SC7'

go<-rbind(go1,go2,go3,go4,go5,go7)

g<-go[,c(12,2)]
colnames(g)<-c('Sample','motif')


gene_sample<-g  %>% mutate(value=1) %>% complete(motif,Sample,fill=list(value=0)) %>%
  mutate(key=paste0(motif)) %>%
  group_by(Sample,key) %>%
  summarize(value = sum(value)) %>%
  spread(key,value) %>% 
  as.data.frame
write.csv(gene_sample,'./Desktop/allllll/table/Ontolgy/GO/heatmap-all-50.csv')
t<-gene_sample[,-1]
row.names(t)<-gene_sample$Sample
t<-as.matrix(t)


hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')
pp1<-ggplot(data = go, aes(x = description, y = Subtype)) +
  geom_tile(aes(fill = PValue)) + coord_equal() +
  scale_fill_gradientn(colours = hm.palette(100))

pp1<-ggplot(data = Go, aes(x = description, y = Subtype),) +
  geom_tile(aes(fill = PValue)) + coord_equal() +
  scale_fill_gradientn(colours = hm.palette(100))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Gene Ontology")+
  theme(plot.title = element_text(hjust = 0.5))


ggarrange(pp1,pp2,heights = c(2, 2),
          ncol = 2, nrow = 1, align = "v")

library("cowplot")
ggdraw() +
  draw_plot(pp1, x = 0, y =0, width = .5) +
  draw_plot(pp2, x = .5, y = 0, width = .5) 



library(highcharter)
library(tidyr)
library(dplyr)

hchart(t, "heatmap", hcaes(x = variable, y = name, value = value)) %>% 
  hc_colorAxis(stops = color_stops(2, c("red","green")))
#================pathway================

p3<-read.csv('./Desktop/allllll/table/Pathway/top 50/3/enrichment_results_wg_result1541696200.txt',sep = '\t')
p4<-read.csv('./Desktop/allllll/table/Pathway/top 50/4/enrichment_results_wg_result1541696222.txt',sep = '\t')
p6<-read.csv('./Desktop/allllll/table/Pathway/top 50/6/enrichment_results_wg_result1541696278.txt',sep = '\t')


p3$Subtype <-'SC3'
p4$Subtype <-'SC4'
p6$Subtype <-'SC6'


p<-rbind(p3,p4,p6)

hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')
pp2<-ggplot(data = p, aes(x = description, y = Subtype),) +
  geom_tile(aes(fill = PValue)) + coord_equal() +
  scale_fill_gradientn(colours = hm.palette(100))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("Pathway")+
  theme(plot.title = element_text(hjust = 0.5))


g<-p[,c(12,2)]
colnames(g)<-c('Sample','motif')


gene_sample<-g  %>% mutate(value=1) %>% complete(motif,Sample,fill=list(value=0)) %>%
  mutate(key=paste0(motif)) %>%
  group_by(Sample,key) %>%
  summarize(value = sum(value)) %>%
  spread(key,value) %>% 
  as.data.frame
write.csv(gene_sample,'./Desktop/allllll/table/Pathway/heatmap-all-50.csv')
t<-gene_sample[,-1]
row.names(t)<-gene_sample$Sample
t<-as.matrix(t)

library(highcharter)
library(tidyr)
library(dplyr)

hchart(t, "heatmap", hcaes(x = variable, y = name, value = value)) %>% 
  hc_colorAxis(stops = color_stops(2, c("red","green")))
#================disease================
p3<-read.csv('./Desktop/allllll/table/Disease/50-GLAD4U/1/enrichment_results_wg_result1543150130.txt',sep = '\t')
p4<-read.csv('./Desktop/allllll/table/Disease/50-GLAD4U/3/enrichment_results_wg_result1543150275.txt',sep = '\t')
p6<-read.csv('./Desktop/allllll/table/Disease/50-GLAD4U/4/enrichment_results_wg_result1543150335.txt',sep = '\t')


p3$sc <-'SC1'
p4$sc <-'SC3'
p6$sc <-'SC4'


p<-rbind(p3,p4,p6)

g<-p[,c(12,2)]
colnames(g)<-c('Sample','motif')


gene_sample<-g  %>% mutate(value=1) %>% complete(motif,Sample,fill=list(value=0)) %>%
  mutate(key=paste0(motif)) %>%
  group_by(Sample,key) %>%
  summarize(value = sum(value)) %>%
  spread(key,value) %>% 
  as.data.frame
write.csv(gene_sample,'./Desktop/allllll/table/Disease/heatmap-all-50.csv')
t<-gene_sample[,-1]
row.names(t)<-gene_sample$Sample
t<-as.matrix(t)

library(highcharter)
library(tidyr)
library(dplyr)

hchart(t, "heatmap", hcaes(x = variable, y = name, value = value)) %>% 
  hc_colorAxis(stops = color_stops(2, c("red","green")))

#====== Venn Diagram ========
library(VennDiagram)
SC1<- read.csv('./Desktop/allllll/table/gene-sc1.csv')
SC2<- read.csv('./Desktop/allllll/table/gene-sc2.csv')
SC3<- read.csv('./Desktop/allllll/table/gene-sc3.csv')
SC4<- read.csv('./Desktop/allllll/table/gene-sc4.csv')
SC5<- read.csv('./Desktop/allllll/table/gene-sc5.csv')
SC6<- read.csv('./Desktop/allllll/table/gene-sc6.csv')
SC7<- read.csv('./Desktop/allllll/table/gene-sc7.csv')


a1<-as.character(SC1[1:100,2])
a2<-as.character(SC2[1:100,2])
a3<-as.character(SC3[1:100,2])
a4<-as.character(SC4[1:100,2])
a5<-as.character(SC5[1:100,2])
a6<-as.character(SC6[1:100,2])
a7<-as.character(SC7[1:100,2])
tt<- cbind(a1,a2,a3,a4,a5,a6,a7)
#-----------Heat map
t<- matrix(0, nrow = 7, ncol = 7)

for (i in 1:7) {
  for (j in 1:7) {
    t[i,j] = nrow(as.data.frame(intersect(tt[,i],tt[,j])))
  }
}
row.names(t)<-c('SC1','SC2','SC3','SC4','SC5','SC6','SC7')
colnames(t)<-c('SC1','SC2','SC3','SC4','SC5','SC6','SC7')

pheatmap(t, cluster_rows = F, cluster_cols = F,display_numbers = T)

####Venn Diagramm

a1<-as.vector(SC1[1:100,]$Var1)
a2<-as.vector(SC2[1:100,]$Var1)
a3<-as.vector(SC3[1:100,]$Var1)
a4<-as.vector(SC4[1:100,]$Var1)
a5<-as.vector(SC5[1:100,]$Var1)
a6<-as.vector(SC6[1:100,]$Var1)
a7<-as.vector(SC7[1:100,]$Var1)

overlap <- calculate.overlap(
  x = list(
    "SC1" = a1,
    "SC2" = a2,
    "SC3" = a3,
    "SC4" = a4,
    "SC5" = a5
  )
)
x = list( "SC6" = a6, "SC7" = a7)
venn.plot <- venn.diagram(x,"Venn_3set_simple.tiff", fill = color[c(1:5)],cat.cex = 1.75,cex = 3)

color=c("blue", "green",'red','purple','pink','orange','brown')  

t <- t[,match(y$X, colnames(t))] # reorder colnames

tt<-as.matrix(t)




t<-read.csv('./Desktop/allllll/table/signature/exposures.tsv', sep = '\t')
library(pheatmap)
m <- as.matrix(t)
m <- t(t(m)/colSums(m))
rownames(m) <- paste0('Signature',c(1:dim(m)[1]))

library(grDevices)
d1 <- seq(0,1,by = 0.01)
c1 <- colorRampPalette(colors = c('blue','green','yellow','orange','red'))(length(d1))
pheatmap(m, cluster_rows = F, cluster_cols = F, gaps_row = c(1:7), show_colnames = F,
         gaps_col = c(44, 184, 248, 453, 466, 515), breaks = d1, color = c1)

t1<- m[,1:44]
t2<- m[,45:184]
t3<- m[,185:248]
t4<- m[,249:453]
t5<- m[,454:466]
t6<- m[,467:515]
t7<- m[,516:536]
fviz_cluster(list(data = t(m), cluster = y1),
             geom = ('point'),ellipse.type = 'euclid',stand=T,
             show.clust.cent=FALSE,
             main = "Clustering results in TWO first PCAs vector")  ## from ‘factoextra’ package

m <- t(gene_sample)
pres <- as.data.frame(prcomp(m)[[2]])[,c(1,2)]
plot(pres$PC1, pres$PC2)
library(ggplot2)
identical(as.character(pcatt$X), rownames(pres))

pcat<- pres
pcat$name<-rownames(pcat)
pcatt<- merge(y,pcat)
pcatt<-subset(pcatt, pcatt$name == pcatt$X)
pcatt$color <- pcatt$X1

pres$color <- factor(pcatt$X1)
ggplot(pres, aes(PC1,PC2), fill = color)+
  geom_point(aes(colour = color))

t<-t(m)

t2<-subset(t, rownames(t) %in% as.character(y2$X))
t3<-subset(t, rownames(t) %in% as.character(y3$X))
t4<-subset(t, rownames(t) %in% as.character(y4$X))
t2<-t(t2)
t3<-t(t3)
t4<-t(t4)

m24<- cbind(t2,t4)
y24<- rbind(y2,y4)

pres <- as.data.frame(prcomp(m24)[[2]])[,c(1,2)]
pcat<- pres
pcat$name<-rownames(pcat)
pcatt<- merge(y24,pcat)
pcatt<-subset(pcatt, pcatt$X == pcatt$name)

pres$color <- factor(pcatt$X1)
ggplot(pres, aes(PC1,PC2), fill = color)+ ylim(-0.07,0.07)+ xlim(-0.001,0.09)+
  geom_point(aes(colour = color))


p1<-ggplot(pres, aes(PC1,PC2), colour = pres$color)+ ylim(-0.05,0.05)+ xlim(-0.025,0.025)+
  geom_point(aes(colour = color)) + ggtitle("All patients")
p2<-ggplot(pres, aes(PC1,PC2), colour = pres$color)+ ylim(-0.025,0.025)+ xlim(-0.025,0.05)+
  geom_point(aes(colour = color)) + ggtitle("SC2, SC3, and SC4 patients")
p3<-ggplot(pres, aes(PC1,PC2), colour = pres$color)+ ylim(-0.05,0.05)+ xlim(-0.0005,0.025)+
  geom_point(aes(colour = color)) + ggtitle("SC2 and SC4 patients")
grid.arrange(p1, p2,p3, nrow = 1)
grid.arrange( p1, p2, p3, nrow=1, ncol=3)

par(mfrow=c(1,3 ))
boxplot(t(m), main="all")
boxplot(t(t1), main="SC1")
boxplot(t(t2), main="SC2")
boxplot(t(t3), main="SC3")
boxplot(t(t4), main="SC4")
boxplot(t(t5), main="SC5")
boxplot(t(t6), main="SC6")
boxplot(t(t7), main="SC7")

pca<-prcomp(m)
pca<-prcomp(m24)
eigs <- pca$sdev^2
eigs[1] / sum(eigs)
r<-rbind(
  SD = sqrt(eigs),
  Proportion = eigs/sum(eigs),
  Cumulative = cumsum(eigs)/sum(eigs))
tm<-t(m)
m24<-subset(tm,rownames(tm) %in% y24$X)
library(rgl)
pca <- princomp(m24, cor=TRUE, scores=TRUE)
plot3d(pca$scores[,1:3], col=y24$X1)



library(factoextra)
res.pca <- prcomp(gene_sample, scale = T)

par(mfrow=c(2,1 ))
fviz_eig(res.pca)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
) 
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
library(ggbiplot)

ggbiplot(res.pca)

library(ggfortify)
autoplot(prcomp(gene_sample, scale. = F), data = gene_sample, 
         loadings = TRUE, loadings.colour = 'blue', main = " Selected features in two most significant PCAs before scaling",
         loadings.label = FALSE, loadings.label.size = 3 )

autoplot(prcomp(gene_sample, scale. = T), data = gene_sample, 
         loadings = TRUE, loadings.colour = 'blue', main = " Selected features in two most significant PCAs after sacling",
         loadings.label = FALSE, loadings.label.size = 3 )

#### Motif rate in specific gene####
var<-read.csv('./Desktop/allllll/Colorectal/var.csv')
var<-var[,-1]
var<-unique(var)
var1<-var
var<-subset(var1, var1$gene %in% c("ENSG00000134982","ENSG00000157764","ENSG00000165323","ENSG00000133703","ENSG00000121879","ENSG00000196090","ENSG00000185008","ENSG00000141510"))
gnmtf<-var[,c(1,9)]
colnames(gnmtf)<-c('sample_id','gm')

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))
gnmtf6<-as.data.frame(table(gnmtf6$gm))
gnmtf7<-as.data.frame(table(gnmtf7$gm))

{
  gnmtf$Freq<-gnmtf$Freq/sum(gnmtf$gm)
  gnmtf1$Freq<-gnmtf1$Freq/sum(gnmtf1$Freq)
  gnmtf2$Freq<-gnmtf2$Freq/sum(gnmtf2$Freq)
  gnmtf3$Freq<-gnmtf3$Freq/sum(gnmtf3$Freq)
  gnmtf4$Freq<-gnmtf4$Freq/sum(gnmtf4$Freq)
  gnmtf5$Freq<-gnmtf5$Freq/sum(gnmtf5$Freq)
  gnmtf6$Freq<-gnmtf6$Freq/sum(gnmtf6$Freq)
  gnmtf7$Freq<-gnmtf7$Freq/sum(gnmtf7$Freq)
}

clr<-c(1:96)
clr[1:16]<-1
clr[17:32]<-2
clr[33:48]<-3
clr[49:64]<-4
clr[65:80]<-5
clr[81:96]<-6
par(mfrow=c(2,4 ))
barplot(gnmtf$Freq,col = clr, main = 'Whole patients motif rate in PCDHA1',ylim=c(0, 0.5))
barplot(gnmtf1$Freq,col = clr, main = 'SC1 motif rate in PCDHA1',ylim=c(0, 0.5))
barplot(gnmtf2$Freq,col = clr, main = 'SC2 motif rate in PCDHA1',ylim=c(0, 0.5))
barplot(gnmtf3$Freq,col = clr, main = 'SC3 motif rate in PCDHA1',ylim=c(0, 0.5))
barplot(gnmtf4$Freq,col = clr, main = 'SC4 motif rate in PCDHA1',ylim=c(0, 0.5))
barplot(gnmtf5$Freq,col = clr, main = 'SC5 motif rate in PCDHA1',ylim=c(0, 0.5))
barplot(gnmtf6$Freq,col = clr, main = 'SC6 motif rate in PCDHA1',ylim=c(0, 0.5))
barplot(gnmtf7$Freq,col = clr, main = 'SC7 motif rate in PCDHA1',ylim=c(0, 0.5))
#------ sankey diagram----
library(networkD3)
library(tidyverse)
library(ggalluvial)

{
  var<-subset(var1, var1$gene == "ENSG00000134982")
  gnmtf<-var[,c(1,9)]
  colnames(gnmtf)<-c('sample_id','gm')
  
  gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
  gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
  gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
  gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
  gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
  gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
  gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))
 
  gnmtf1<-as.data.frame(table(gnmtf1$gm))
  gnmtf2<-as.data.frame(table(gnmtf2$gm))
  gnmtf3<-as.data.frame(table(gnmtf3$gm))
  gnmtf4<-as.data.frame(table(gnmtf4$gm))
  gnmtf5<-as.data.frame(table(gnmtf5$gm))
  gnmtf6<-as.data.frame(table(gnmtf6$gm))
  gnmtf7<-as.data.frame(table(gnmtf7$gm))
  
  gnmtf1$subtype<-"SC1"
  gnmtf2$subtype<-"SC2"
  gnmtf3$subtype<-"SC3"
  gnmtf4$subtype<-"SC4"
  gnmtf5$subtype<-"SC5"
  gnmtf6$subtype<-"SC6"
  gnmtf7$subtype<-"SC7"
  
  gnmtf<-rbind(gnmtf1,gnmtf2,gnmtf3,gnmtf4,gnmtf5,gnmtf6,gnmtf7)
  gnmtf$Gene<-"APC"
}
df<-rbind(df,gnmtf)
colnames(df)<-c("Motif","Freq","Subtype","Gene")
df<-df[,c(3,4,1,2)]
alluvial(df[,c(1:3)], freq=df$Freq, col  = c(1:7),
         hide = df$Freq == 0,
         cex = 0.7
)
ggplot(as.data.frame(df1),
       aes(y = Freq, axis1 = Subtype, axis2 = Gene, axis3 = Motif)) +
  geom_alluvium(aes(fill = Subtype), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("Gene", "Motif"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Frequency of motifs in significant CRC related gene")

var$gm<-paste0(var$)
df<-as.data.frame(table(var,var$gene,var$motif))

ggplot(as.data.frame(var),
       aes(y = Freq, axis1 = icgc_samlpe_id, axis2 = gene, axis3 = motif)) +
  geom_alluvium(aes(fill = Subtype), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("Gene", "Motif"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Frequency of motifs in significant CRC related gene")

df1<-df

df1<-subset(df1, df1$Motif %in% t[c(96:90),1] & df1$Freq>=10 & df1$Gene %in% c("TP53","APC"))

for (i in 1:5280) {
  if(df1[i,1]=="SC1"){
    df1[i,4] <-  df1[i,4]/884 
  }
  if(df1[i,1]=="SC2"){
    df1[i,4] <-  df1[i,4]/ 2339
  }
  if(df1[i,1]=="SC3"){
    df1[i,4] <-  df1[i,4]/ 1230
  }
  if(df1[i,1]=="SC4"){
    df1[i,4] <-  df1[i,4]/ 4334
  }
  if(df1[i,1]=="SC5"){
    df1[i,4] <-  df1[i,4]/ 1529
  }
  if(df1[i,1]=="SC6"){
    df1[i,4] <-  df1[i,4]/ 1375
  }
  if(df1[i,1]=="SC7"){
    df1[i,4] <-  df1[i,4]/ 404
  }
    
}

#### Transcript in gene ####

gnmtf<-var[,c(1,3,4)]
gnmtf<-unique(gnmtf)
colnames(gnmtf)<-c('sample_id',"gene",'gm')
idnamefindgene <- subset(signame, signame$external_gene_name %in% findgene$Gene.Symbol)
gnmtf<-subset(gnmtf,gnmtf$gene %in% idnamefindgene$ensembl_gene_id)

gnmtf1<-subset(gnmtf,gnmtf$sample_id %in% row.names(y1))
gnmtf2<-subset(gnmtf,gnmtf$sample_id %in% row.names(y2))
gnmtf3<-subset(gnmtf,gnmtf$sample_id %in% row.names(y3))
gnmtf4<-subset(gnmtf,gnmtf$sample_id %in% row.names(y4))
gnmtf5<-subset(gnmtf,gnmtf$sample_id %in% row.names(y5))
gnmtf6<-subset(gnmtf,gnmtf$sample_id %in% row.names(y6))
gnmtf7<-subset(gnmtf,gnmtf$sample_id %in% row.names(y7))

gnmtf<-as.data.frame(table(gnmtf$gm))
gnmtf1<-as.data.frame(table(gnmtf1$gm))
gnmtf2<-as.data.frame(table(gnmtf2$gm))
gnmtf3<-as.data.frame(table(gnmtf3$gm))
gnmtf4<-as.data.frame(table(gnmtf4$gm))
gnmtf5<-as.data.frame(table(gnmtf5$gm))
gnmtf6<-as.data.frame(table(gnmtf6$gm))
gnmtf7<-as.data.frame(table(gnmtf7$gm))

{
  gnmtf$Freq<-gnmtf$Freq/nrow(y)
  gnmtf1$Freq<-gnmtf1$Freq/nrow(y1)
  gnmtf2$Freq<-gnmtf2$Freq/nrow(y2)
  gnmtf3$Freq<-gnmtf3$Freq/nrow(y3)
  gnmtf4$Freq<-gnmtf4$Freq/nrow(y4)
  gnmtf5$Freq<-gnmtf5$Freq/nrow(y5)
  gnmtf6$Freq<-gnmtf6$Freq/nrow(y6)
  gnmtf7$Freq<-gnmtf7$Freq/nrow(y7)
  gnmtf1$subtype<-"SC1"
  gnmtf2$subtype<-"SC2"
  gnmtf3$subtype<-"SC3"
  gnmtf4$subtype<-"SC4"
  gnmtf5$subtype<-"SC5"
  gnmtf6$subtype<-"SC6"
  gnmtf7$subtype<-"SC7"
}
gnmtf<-rbind(gnmtf1,gnmtf2,gnmtf3,gnmtf4,gnmtf5,gnmtf6,gnmtf7)
colnames(genename)<- c("gene","gene_name")
trans<-merge(genename,tt)

df<-subset(trans,trans$gene_name %in% c("APC"))

ggplot(as.data.frame(df),
       aes(y = Freq, axis1 = subtype, axis2 = gene_name, axis3 = transcript)) +
  geom_alluvium(aes(fill = subtype), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("Gene", "Motif"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Frequency of motifs in significant CRC related gene")


#--- circular bar---
selectedtrans<- subset(trans, trans$gene_name=="APC")
data<-selectedtrans[,c(3,5,4)]
colnames(data)<-c("individual","group","value")




#### Survival ####

library(survival)
library(survminer)
library(dplyr)

donor<-read.csv("~/Desktop/allllll/Colorectal/donor.tsv",sep = '\t')
sample<-read.csv("~/Desktop/allllll/Colorectal/sample.tsv",sep = '\t')
surv<- donor[,c("icgc_donor_id", "donor_age_at_diagnosis", "donor_survival_time")]
t<- subset(sample,sample$icgc_sample_id %in% y$X)
yt<-y
colnames(yt)<-c("icgc_sample_id","Subtype")
t<-merge(t,yt)
t<-k
k<-merge(t,surv)
k$fustat<-0
for (i in 1:530) {
  if(k[i,5] == 0){
    k[i,6]=0
  }
  if(k[i,5] > 0){
    k[i,6]=1
  }
}

surv_object <- Surv(time = k$donor_survival_time, event = k$fustat)
fit1 <- survfit(surv_object ~ Subtype, data = k)
ggsurvplot(fit1, data = k, pval = TRUE)

fit.coxph <- coxph(surv_object ~ Subtype + donor_age_at_diagnosis , 
                   data = k)
ggforest(fit.coxph, data = k)

data(ovarian)
glimpse(ovarian)

ovarian$rx <- factor(ovarian$rx, 
                     levels = c("1", "2"), 
                     labels = c("A", "B"))
ovarian$resid.ds <- factor(ovarian$resid.ds, 
                           levels = c("1", "2"), 
                           labels = c("no", "yes"))
ovarian$ecog.ps <- factor(ovarian$ecog.ps, 
                          levels = c("1", "2"), 
                          labels = c("good", "bad"))

hist(ovarian$age) 
ovarian <- ovarian %>% mutate(age_group = ifelse(age >=50, "old", "young"))
ovarian$age_group <- factor(ovarian$age_group)

surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat)
surv_object 

fit1 <- survfit(surv_object ~ rx, data = ovarian)
summary(fit1)
ggsurvplot(fit1, data = ovarian, pval = TRUE)

fit2 <- survfit(surv_object ~ resid.ds, data = ovarian)
ggsurvplot(fit2, data = ovarian, pval = TRUE)

fit.coxph <- coxph(surv_object ~ rx + resid.ds + age_group + ecog.ps, 
                   data = ovarian)
ggforest(fit.coxph, data = ovarian)



#####Gene name######

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",GRCh = 37)
chr1_genes <- getBM(attributes=c('ensembl_gene_id',
                                 'external_gene_name'), mart = ensembl)
significant<-read.csv("~/Desktop/allllll/table/significant_gene.csv")
cancer<-read.csv('~/Desktop/Census_allWed Oct 24 08_14_27 2018.csv')

signame<- subset(chr1_genes,chr1_genes$ensembl_gene_id %in% significant$Var1)
findgene<- subset(cancer, cancer$Gene.Symbol %in% signame$external_gene_name)



#####Figure 2 Exposure of signatures#####
{
  t<-read.csv('./Desktop/allllll/table/signature/exposures.tsv', sep = '\t')
  library(pheatmap)
  m <- as.matrix(t)
  m <- t(t(m)/colSums(m))
  rownames(m) <- paste0('Signature',c(1:dim(m)[1]))
  expo<- as.data.frame(t(m))
  
  y<-read.csv('./Desktop/allllll/table/cluster7.csv')
  row.names(y)<-y$X
  y1<-subset(y,y$X1==1)
  y2<-subset(y,y$X1==2)
  y3<-subset(y,y$X1==3)
  y4<-subset(y,y$X1==4)
  y5<-subset(y,y$X1==5)
  y6<-subset(y,y$X1==6)
  y7<-subset(y,y$X1==7)
  
  expo$X<-row.names(expo)
  t<-merge(expo,y)
  colnames(t)<-c("X","Signature1","Signature2","Signature3","Signature4","Signature5","Signature6","Signature7",
                 "Subtype")
  t$Subtype<-paste0("SC",t$Subtype)
  {
    s1<- ggplot(t, aes(Subtype,Signature1 )) +
      geom_violin(aes(fill = Subtype))
    s2<- ggplot(t, aes(Subtype,Signature2 )) +
      geom_violin(aes(fill = Subtype))
    s3<- ggplot(t, aes(Subtype,Signature3 )) +
      geom_violin(aes(fill = Subtype))
    s4<- ggplot(t, aes(Subtype,Signature4 )) +
      geom_violin(aes(fill = Subtype))
    s5<- ggplot(t, aes(Subtype,Signature5 )) +
      geom_violin(aes(fill = Subtype))
    s6<- ggplot(t, aes(Subtype,Signature6 )) +
      geom_violin(aes(fill = Subtype))
    s7<- ggplot(t, aes(Subtype,Signature7 )) +
      geom_violin(aes(fill = Subtype))
    
  }
  grid.arrange(s1,s2,s3,s4,s5,s6,s7, nrow = 7)
  
}
pdf(file="./Desktop/Figures/exposure.pdf") 
plot(s2)
dev.off() 


