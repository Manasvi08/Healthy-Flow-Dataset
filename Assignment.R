install.packages("corrr")
install.packages("corrplot")
install.packages("Rtsne")
install.packages("ggplot2")
library(corrr)
library(Rtsne)
library(corrplot)
install.packages("purrr")
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(purrr)
library(ggplot2)
library(ggpubr)

hfd_csv <- read.csv("hfd_sub.csv")
plot(hfd_csv)
dim(hfd_csv)
summary(hfd_csv)
correlation<-cor(hfd_csv)
corrplot(correlation,diag = TRUE,outline = TRUE,tl.col = "black", tl.offset = 1,tl.srt = 45,tl.cex=0.5,addCoefasPercent=TRUE, main = "Correlation Plot")
corrplot()
cov(hfd_csv[,1:4])
diag(cov(hfd_csv[, 1:4]))

plot(density(hfd_csv$CD4,main="CD4"))
hist(hfd_csv$CD4,main="CD4")
boxplot(hfd_csv$CD4, main="CD4")

plot(density(hfd_csv$CD3))
hist(hfd_csv$CD3, main="CD3")
boxplot(hfd_csv$CD3,main="CD3")

plot(density(hfd_csv$CD8))
hist(hfd_csv$CD8,main="CD8")
boxplot(hfd_csv$CD8,main="CD8")

plot(density(hfd_csv$CD19))
hist(hfd_csv$CD19)
boxplot(hfd_csv$CD19,main="CD19")

plot(density(hfd_csv$gate))
hist(hfd_csv$gate,main="gate")
boxplot(hfd_csv$gate,main="gate")


plot(hfd_csv[, 1], hfd_csv[, 2], col=hfd_csv[, 5])
pairs(hfd_csv[,1:4],col=hfd_csv[,5])


#PCA
#We can use prcomp to directly apply PCA to a dataset, instead of using 
#the cov and eigen functions separately, like we did in last week’s lab.
#We’ve discussed how we obtain principal components (PCs) via an eigendecomposition 
#of the covariance matrix. The function prcomp uses a different but related technique that is more stable numerically.
#The output Standard deviations refers to the square root of the eigenvalues of the covariance matrix. Rotation are the eigenvectors of the covariance matrix, i.e., the PCs.
pca_analysis <- prcomp(hfd_csv[, 1:4],scale=TRUE)
pca_analysis
summary(pca_analysis)
round(pca_analysis$rotation, 2)
round(pca_analysis$sdev,2)
plot(pca_analysis, main="Healthy Flow Dataset")
healthy_var_explain <- (pca_analysis$sdev^2) / (sum(pca_analysis$sdev^2))
plot(healthy_var_explain, type = "b", main = "Healthy Flow Data", 
     xlab = "No. of components", ylab = "Proportion of variance explained", xaxt = "n")
axis(1, at = 1:4)




pairs(new_pca_analysis,col=hfd_csv[,5])
#From the second plot you should be able to see that the values of the leading principal component PC1 provide a good separation of the observations by species.
#plot(new_pca_analysis[,1], new_pca_analysis[,2], type="n", xlab="PC1", ylab="PC2")
#text(new_pca_analysis[,1], new_pca_analysis[,2], labels=substr(hfd_csv[,5],1,2), col=as.integer(hfd_csv[,5]))
#write why standardisation, why plotting what plotting is done, PCA using correlation,
#Just use the lab doc to write what is done for the analysis and interpretation see from internet. 
biplot(pca_analysis)
#A biplot has been plotted which gives how the two principal Components are aligned with their observations.
#It is seen that CD3 contributes most to PC1 while CD8 to PC2. CD4 is tend to contribute more or less equally to both  PC2 and PC1 but CD19 is more inclined parallelly to PC1. Also, more variability is seen in CD8 and CD19 by the length of the vectors. 


#https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/Principal-Component-Analysis/principal-components-basics/Interpretation-and-visualization/index.html

#Further, prediction of how each of the observation maps on to Principal Components is done for visualisation of how Principal Components prediction go in line with the different groups using the gate variable from the dataset.It is seen that there are few groups formed with the Principal Components 1 and 2. And the chart gives that there are some groups in aligning with the gate variables. 
new_pca_analysis <- predict(pca_analysis)
#head(new_pca_analysis, n = 5)
pairs(new_pca_analysis,col=hfd_csv[,5])
#it is seen that the leading pC1 provides a good separation between different gating numbers. 
#plot(new_pca_analysis[,1], new_pca_analysis[,2], type="n", xlab="PC1", ylab="PC2")
#text(new_pca_analysis[,1], new_pca_analysis[,2], labels=substr(hfd_csv[,5],1,2), col=as.integer(hfd_csv[,5]))


#Hierarchial Clustering
hfd_matrix <- as.matrix(hfd_csv[,1:4])
hfd_dis <- dist(hfd_matrix, method="euclidean")
hfd_dis_mat <- as.matrix(hfd_dis)
clust1 <- hclust(hfd_dis, method = "complete")
plot(clust1)
head(clust1$height)
h_mean<-mean(clust1$height)
sd_height<-sd(clust1$height)
recommended_height<-h_mean+(3*sd_height)
recommended_height
plot(clust1)
abline(h = recommended_height, lty=2, col=2)

cuttree_k <- cutree(clust1, k=4)
cuttree_height <- cutree(clust1, h=recommended_height)
palette(rainbow(10))
#plot(hfd_matrix[,1], hfd_matrix[,2], col = acids_label1)
pairs(hfd_matrix, col = cuttree_k)
pairs(hfd_matrix, col = cuttree_height)


table(cuttree_height, hfd_csv$gate)
table(cuttree_k, hfd_csv$gate)

#k-means clustering
totalDistance1 <- function(k1)
{
  #getting the medoids
  kmedoids12 <- pam(gow1, k1)
  #print(kmedoids12$clusinfo)
  ##summing to find the dissimilarity
  totalDis1 <- sum(kmedoids12$clusinfo[,1]*kmedoids1$clusinfo[,3])
  return(totalDis1)
}

#setting the potential k values to analyse from 1 to 15
k_values <- 1:10

#using map function to apply a function to each element and returning a vector of the same length 
distance_values1 <- map_dbl(k_values, totalDistance1)

library(purrr)

##using the elbow method we can see the different values of k and how the dissmilarity between medoids decreases as k increases. Where there is a 'kink' we can see a suitable choice for k
plot(k_values, distance_values1,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Average dissimilarity",
     main = "Elbow Method")




gow1 <- daisy(hfd_matrix[, 1:4], metric="gower")
kmedoids12 <-pam(gow1, 4)
kmedoids12$clusinfo
table(kmedoids12$clustering, hfd_csv$gate)

require(dplyr)
kmed_results <- hfd_csv %>%
  mutate(cluster = kmedoids12$clustering) %>%
  group_by(cluster)
kmed_summary <-kmed_results %>%
  do(the_summary = summary(kmed_results[1:4]))
kmed_summary$the_summary

require(reshape2)
melt_res <- melt(kmed_results[-5], id.var = "cluster")
melt_res$cluster <- as.factor(melt_res$cluster)
ggplot(data = melt_res, aes(x=variable, y=value)) + geom_boxplot(aes(fill=cluster))
p <- ggplot(data = melt_res, aes(x=cluster, y=value)) +
  geom_boxplot(aes(fill=cluster))
p + facet_wrap( ~ variable, scales="free")


sil <- silhouette(kmedoids12$cluster, dist(hfd_csv))
fviz_silhouette(sil)


fviz_cluster(res.km, data = hfd_csv[, -5],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#00AFBB"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

res.km <- kmeans(scale(hfd_csv[, -5]), 4, nstart = 25)

res.pca <- prcomp(hfd_csv[, -5],  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res.km$cluster)
# Add Species groups from the original data sett
ind.coord$gate <- hfd_csv$gate
# Data inspection
head(ind.coord)


eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)

ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
   size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)

clust1_label <- cutree(clust1, k=4)
clust2_label <- cutree(clust1, k=5)

install.packages("flexclust")
library(flexclust)
randIndex(clust1_label, clust2_label, correct = FALSE)

#https://www.datanovia.com/en/blog/k-means-clustering-visualization-in-r-step-by-step-guide/


#Supervised
#knn



#k-means
WSS <- rep(0,10)
WSS[2] <- sum(kmeans(hfd_csv, centers = 2)$withinss)


plot(scale.hfd.csv, col = cl2$cluster,main="Centroid of the clusters")
points(cl2$centers, col=1:k, pch=8, cex=5)
cl2 <- kmeans(scale.hfd.csv,centers = 4)


table(cl2$cluster,hfd_csv$gate)

g1 <- hfd_csv[which(cl2$cluster==1),]
ng1 <- cl2$size[1]
total1 <- sum(as.matrix(dist(rbind(g1, cl2$centers[1,])))[ng1+1,])
ave1 <- total1/ng1
library("cluster")
dist_mat <- dist(hfd_csv)*dist(hfd_csv)
silhouette_cl2 <- silhouette(cl2$cluster, dist_mat)
summary(silhouette_cl2)
plot(silhouette_cl2,border=NA)





