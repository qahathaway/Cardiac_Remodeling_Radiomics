
library(ez.combat)

cb <- ez.combat(df = df0, batch.var = "Group")

pca_res7 <- prcomp(cb$df, scale. = TRUE)

autoplot(pca_res7, data = df0, colour = "Group")

df = iris
library(ggfortify)

df0 <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/Radiomics_Revised_Values.csv'))

df <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/Radiomics_Values.csv'))

df2 <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/Radiomics_Normalized.csv'))

df3 <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/Radiomics_Normalized_All.csv'))

df4 <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/Radiomics_Normalized_Reference.csv'))

df6 <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/ASE-REWARD_Revised2.csv'))

df6 <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/Radiomics_Revised_Normalized.csv'))


pca_res <- prcomp(df, scale. = TRUE)
pca_res2 <- prcomp(df2, scale. = TRUE)
pca_res3 <- prcomp(df3, scale. = TRUE)
pca_res4 <- prcomp(df4, scale. = TRUE)
pca_res6 <- prcomp(df6, scale. = TRUE)


autoplot(pca_res, data = df, colour = "Group")
autoplot(pca_res2, data = df2, colour = "Group")
autoplot(pca_res3, data = df3, colour = "Group")
autoplot(pca_res4, data = df4, colour = "Group")
autoplot(pca_res6, data = df6, colour = "Group")


library(cluster)
autoplot(clara(df2, 3))
df <- autoplot(clara(df2, 3))

df$data
write.csv(df$data, file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/Radiomics_Normalized_Cluster.csv')


library(randomForest)
library(rfUtilities)

df5 <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/Radiomics_Normalized_ASE.csv'))

n = 2
clust.df <- rf.unsupervised(df5, n=n, proximity = TRUE,
                              silhouettes = TRUE)

mds <- stats:::cmdscale(clust.df$distances, eig=TRUE, k=n)
colnames(mds$points) <- paste("Dim", 1:n)
mds.col <- ifelse(clust.df$k == 1, rainbow(4)[1],
                  ifelse(clust.df$k == 2, rainbow(4)[2],
                         ifelse(clust.df$k == 3, rainbow(4)[3],
                                ifelse(clust.df$k == 4, rainbow(4)[4], NA))))
plot(mds$points[,1:2],col=mds.col, pch=20)
pairs(mds$points, col=mds.col, pch=20)
