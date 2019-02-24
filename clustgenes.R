#Program for spectral clustering of genes
tfs_expression_colombos_sh <- read.delim("tfs_expression_colombos-sh.txt")
View(tfs_expression_colombos_sh)
summary(tfs_expression_colombos_sh)
z=as.matrix(tfs_expression_colombos_sh[,c(4:20)])
image(z)
z=as.matrix(tfs_expression_colombos_sh[,4:ncol(tfs_expression_colombos_sh)])
#PCA 
pca1=prcomp(na.omit(z))
plot(pca1)
biplot(pca1)
length(pca1$rotation[pca1$rotation[,1]>0.05,1])
length(pca1$rotation[pca1$rotation[,2]>0.05,2])
pca1$rotation[pca1$rotation[,1]>0.05,1]
pca1$rotation[pca1$rotation[,2]>0.05,2]
pca1$rotation[pca1$rotation[,3]>0.05,3]
names(tfs_expression_colombos.sh)
ind=c(6,7,10,13,37,48,57,62,76,90,116,163,227,230,215,294)
z.clus=as.matrix(tfs_expression_colombos.sh[,ind])

#Kmeans
image(z)
nc=16
kmed1=kmeans(na.omit(z.clus),centers=nc)
kmed1
resulkmed=sapply(4:64, function(j) {km=kmeans(na.omit(z.clus),centers=j);return=c(km[5],km[6]) })
plot(4:64,resulkmed[1,])
#el minimo es nclus=16

par(mfrow=c(2,2))
 image(na.omit(z.clus)[kmed1$cluster==1,])
 image(na.omit(z.clus)[kmed1$cluster==2,])
 image(na.omit(z.clus)[kmed1$cluster==3,])
 image(na.omit(z.clus)[kmed1$cluster==4,])
image(na.omit(z.clus)[kmed1$cluster==5,])
 image(na.omit(z.clus)[kmed1$cluster==6,])
 image(na.omit(z.clus)[kmed1$cluster==7,])
 image(na.omit(z.clus)[kmed1$cluster==8,])
 image(na.omit(z.clus)[kmed1$cluster==9,])
 image(na.omit(z.clus)[kmed1$cluster==10,])
 image(na.omit(z.clus)[kmed1$cluster==11,])
 image(na.omit(z.clus)[kmed1$cluster==12,])
 image(na.omit(z.clus)[kmed1$cluster==13,])
 image(na.omit(z.clus)[kmed1$cluster==14,])
 image(na.omit(z.clus)[kmed1$cluster==15,])
 image(na.omit(z.clus)[kmed1$cluster==16,])
 
 randomkm=lapply(1:20,function(j)kmeans(na.omit(z.clus),centers=16))
 
 
 kmedrm1=randomkm[[10]]#el mejor detalles guardados en archivo
 
 clusgenes=sapply(1:16,function(j) na.omit(as.matrix(tfs_expression_colombos.sh[,c(2,ind)]))[kmedrm1$cluster==j,1])
 
 names(tfs_expression_colombos.sh)[ind]#variables
 
 #Imagenes de la mejor agrupaci?n
 
 par(mfrow=c(2,2))
 image(na.omit(z.clus)[kmedrm1$cluster==1,])
 image(na.omit(z.clus)[kmedrm1$cluster==2,])
 image(na.omit(z.clus)[kmedrm1$cluster==3,])
 image(na.omit(z.clus)[kmedrm1$cluster==4,])
 image(na.omit(z.clus)[kmedrm1$cluster==5,])
 image(na.omit(z.clus)[kmedrm1$cluster==6,])
 image(na.omit(z.clus)[kmedrm1$cluster==7,])
 image(na.omit(z.clus)[kmedrm1$cluster==8,])
 image(na.omit(z.clus)[kmedrm1$cluster==9,])
 image(na.omit(z.clus)[kmedrm1$cluster==10,])
 image(na.omit(z.clus)[kmedrm1$cluster==11,])
 image(na.omit(z.clus)[kmedrm1$cluster==12,])
 image(na.omit(z.clus)[kmedrm1$cluster==13,])
 image(na.omit(z.clus)[kmedrm1$cluster==14,])
 image(na.omit(z.clus)[kmedrm1$cluster==15,])
 image(na.omit(z.clus)[kmedrm1$cluster==16,])
 

 #clusterizacion espectral basado en knn
 #specclust(data, centers=NULL, nn = 7, method = "symmetric", gmax=NULL, ...)
 
 ## S3 method for class 'specclust':
 #plot((x, ...))
 
 install.packages("kknn")
 library(kknn)
 
 nclsp=12
 nn=3
 clrd <- specClust(na.omit(z.clus), nclsp, nn)
 
 pcol <- as.character(as.numeric(clrd$cluster))
 colme=colors(1)[sample(30:500,12)]
pairs(na.omit(z.clus)[,1:4], pch = pcol, col = colme[clrd$cluster])
pairs(na.omit(z.clus)[,5:8], pch = pcol, col = colme[clrd$cluster])
pairs(na.omit(z.clus)[,9:12], pch = pcol, col = colme[clrd$cluster])
 pairs(na.omit(z.clus)[,13:16], pch = pcol, col = colme[clrd$cluster])
 
 resultcl=sapply(4:64, function(j) {km=specClust(na.omit(z.clus),j,nn=10);return=c(km[5],km[6]) })
 
 par(mfrow=c(1,2))
 plot(4:64,resultcl[1,])
 plot(4:64,resulkmed[1,])

 
  
 
 randomspc=lapply(1:20, function(j) {specClust(na.omit(z.clus),nclsp,nn=10) })
 
 clrd=randomspc[[1]]
 #Genes in each cluster
 clusgenesspectr=sapply(1:nclsp,function(j) na.omit(as.matrix(tfs_expression_colombos.sh[,c(2,ind)]))[clrd$cluster==j,1])
 #Image plots for genes over 16 selected attributes for each cluster
 #Rows stand for 16 attributes and columns for genes in each group
 par(mfrow=c(3,2))
 image(na.omit(z.clus)[clrd$cluster==1,])
 image(na.omit(z.clus)[clrd$cluster==2,])
 image(na.omit(z.clus)[clrd$cluster==3,])
 image(na.omit(z.clus)[clrd$cluster==4,])
 image(na.omit(z.clus)[clrd$cluster==5,])
 image(na.omit(z.clus)[clrd$cluster==6,])
 image(na.omit(z.clus)[clrd$cluster==7,])
 image(na.omit(z.clus)[clrd$cluster==8,])
 image(na.omit(z.clus)[clrd$cluster==9,])
 image(na.omit(z.clus)[clrd$cluster==10,])
 image(na.omit(z.clus)[clrd$cluster==11,])
 image(na.omit(z.clus)[clrd$cluster==12,])
 
 
 
  
 
 