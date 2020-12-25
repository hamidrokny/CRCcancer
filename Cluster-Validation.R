library('clValid')
library('kohonen')
library('clv')
data(mouse)
## internal validation
express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID[1:25]
intern <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                  validation="internal")
intern <- clValid(gene_sample, 3:15, clMethods=c( "kmeans",  "pam", "diana", "fanny", "sota", "pam", "clara","agnes", "model"),
                  validation="internal")
slotNames(intern)
## view results
intern

summary(intern)
optimalScores(intern)
plot(intern)
## Extract objects from slots
measures(intern)
hierClust <- clusters(intern,"model")
plot(hierClust)
measNames(intern)
nClusters(intern)

## stability measures
stab <- clValid(express, 4:7, clMethods=c("hierarchical", "kmeans",  "pam", "diana", "fanny", "sota", "pam", "clara", "model"),
                validation="stability")
optimalScores(stab)
plot(stab)
## biological measures
## first way - functional classes predetermined
fc <- tapply(rownames(express),mouse$FC[1:25], c)
fc <- fc[-match( c("EST","Unknown"), names(fc))]
bio <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
               validation="biological", annotation=fc)
optimalScores(bio)
plot(bio)


data(mouse)
express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID[1:25]
## hierarchical clustering
Dist <- dist(gene_sample,method="euclidean")
clusterObj <- hclust(Dist, method="average")
nc <- 2 ## number of clusters
cluster <- cutree(clusterObj,nc)

#===============================================================================
connectivity(Dist, s)
dunn(Dist, s)
ss<-silhouette(s,Dist)
plot(ss)
summary(ss)
t<-clv.Davies.Bouldin(cls.scatt.data(y,s),intracls = c("complete","average","centroid"),intercls = c("single", "complete", "average","centroid", "aveToCent", "hausdorff"))
#===============================================================================

data(mouse)
express <- mouse[1:25,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID[1:25]
## hierarchical clustering
Dist <- dist(gene_sample,method="canberra")
clusterObj <- hclust(Dist, method="average")
nc <- 4 ## number of clusters
cluster <- cutree(clusterObj,nc)
stab <- matrix(0,nrow=ncol(express),ncol=4)
colnames(stab) <- c("APN","AD","ADM","FOM")
## Need loop over all removed samples
for (del in 1:ncol(express)) {
  matDel <- express[,-del]
  DistDel <- dist(matDel,method="euclidean")
  clusterObjDel <- hclust(DistDel, method="model")
  clusterDel <- cutree(clusterObjDel,nc)
  stab[del,] <- stability(gene_sample, Dist, del, s, clusterDel)
}
colMeans(stab)
ts<-as.data.frame(t(a))

#### Ensmeble ID to gene name ####
library('biomaRt')
grch37 = useEnsembl(biomart="ensembl",GRCh=37)
listDatasets(grch37)[31:35,]
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",GRCh = 37)
chr1_genes <- getBM(attributes=c('ensembl_gene_id',
                                   'external_gene_name'), mart = ensembl)
significant<-read.csv("~/Desktop/allllll/table/significant_gene.csv")
cancer<-read.csv('~/Desktop/Census_allWed Oct 24 08_14_27 2018.csv')

signame<- subset(chr1_genes,chr1_genes$ensembl_gene_id %in% significant$Var1)
findgene<- subset(cancer, cancer$Gene.Symbol %in% signame$external_gene_name)

