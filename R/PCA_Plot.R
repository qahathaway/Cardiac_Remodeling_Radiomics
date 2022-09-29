
library(ez.combat)

cb <- ez.combat(df = df0, batch.var = "Group")

pca_res7 <- prcomp(cb$df, scale. = TRUE)

autoplot(pca_res7, data = df0, colour = "Group")


library(ggfortify)

df <- data.frame(read.csv(
  file = 'path/to/file.csv'))

pca_res <- prcomp(df, scale. = TRUE)
autoplot(pca_res, data = df, colour = "Group")

library(cluster)

autoplot(clara(df2, 3))
df <- autoplot(clara(df2, 3))

df$data
write.csv(df$data, file = 'path/to/file.csv')


library(randomForest)
library(rfUtilities)

n = 2
clust.df <- rf.unsupervised(df, n=n, proximity = TRUE,
                              silhouettes = TRUE)

mds <- stats:::cmdscale(clust.df$distances, eig=TRUE, k=n)
colnames(mds$points) <- paste("Dim", 1:n)
mds.col <- ifelse(clust.df$k == 1, rainbow(4)[1],
                  ifelse(clust.df$k == 2, rainbow(4)[2],
                         ifelse(clust.df$k == 3, rainbow(4)[3],
                                ifelse(clust.df$k == 4, rainbow(4)[4], NA))))
plot(mds$points[,1:2],col=mds.col, pch=20)
pairs(mds$points, col=mds.col, pch=20)
