source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("cummeRbund")
library(cummeRbund)
browseVignettes("cummeRbund")
setwd("~/Desktop/cuffdiff_out")

#rebuild cuffdiff data
cuff<-readCufflinks(rebuild=T)
cuff

#squared coefficient of variation, SCV, for evaluating the quality your RNA-seq data
fpkmSCVPlot(genes(cuff))

#plot dispersion
dispersionPlot(genes(cuff))

#density plot
csDensity(genes(cuff))
csDensity(genes(cuff),replicates=T)

#scatterplot
csScatter(genes(cuff), 'j17', 'J37',smooth=T)
csScatterMatrix(genes(cuff))

#boxplot
csBoxplot(genes(cuff))
csBoxplot(genes(cuff),replicates=T)

#volcano plot
csVolcanoMatrix(genes(cuff))

#significant genes
sig <- getSig(cuff, alpha=0.01, level='genes')
length(sig)

#get data from cuff database
sigGenes <- getGenes(cuff,sig)
sigGenes

#get last 100 significant genes, heatmap
sigGenes100 <- getGenes(cuff,tail(sig,100))
csHeatmap(sigGenes50, cluster='row')
csHeatmap(sigGenes, cluster='both')

#dendrogram
csDendro(genes(cuff),replicates=T)
csDendro(genes(cuff))

# PCA
PCAplot(genes(cuff),"PC1","PC2")
PCAplot(genes(cuff),"PC1","PC2",replicates=T)

#MDS
MDSplot(genes(cuff))
MDSplot(genes(cuff),replicates=T)


#More information could be found in 
#https://github.com/solgenomics/bioinfocourse/blob/master/scripts/cummeRbund.R
#http://www.rpubs.com/aemoore62/Cuff_prelim_analysis
#https://www.jianshu.com/p/51400369594c
