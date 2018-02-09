#function to cut heirarchical clusterization results at given clusters numer and plot all pfrofiles from each

cuplcl <- function(data, method, k) {
library(fastcluster)
library(dendextend)
hclust <- hclust.vector(data, method = method)
cut_hclust<-cutree(hclust, k=k)
m <- rbind(rep(1, k), c(2:(k+1)))
layout(m)
par(mar = c(1, 1, 0, 0))

dend <- color_branches(as.dendrogram(hclust), col = cut_hclust[order.dendrogram(as.dendrogram(hclust))])
plot(dend%>%set('labels_cex', NA))

for (i in unique(cut_hclust[order.dendrogram(as.dendrogram(hclust))])) {
matplot(t(data[which(cut_hclust==i),]), type='l', col=i)

print(c(i, names(which(cut_hclust==i))))
}
}
