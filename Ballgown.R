# set the working directory
setwd(...)

# install ballgown R package along with its dependencies
source("http://bioconductor.org/biocLite.R")
biocLite("ballgown")
library(ballgown)
library(ggplot2)
library(gplots)
library(genefilter)
library(GenomicRanges)
install.packages("plyr")
library(plyr)

# Read the design_matrix file
pheno_data = read.table(file ="Design_MatrixminusJ3.txt", header = TRUE, sep = "\t")

# full path to the sample directories
sample_full_path <- paste("ballgown_input_files",pheno_data[,1], sep = '/')

# Load ballgown data structure and save it to a variable “bg”
bg = ballgown(samples=as.vector(sample_full_path),pData=pheno_data)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)

# Load gene names for lookup later in the tutorial
bg_table = texpr(bg_filt, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Pull the gene_expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg_filt))

# View the column names; change the column names to your sample names
colnames(gene_expression) <- c("J1_2_1_S1_L007_R1_001","J1_2_2_S2_L007_R1_001",
                               "J1_2_3_S3_L007_R1_001","J1_7_1_S10_L007_R1_001",
                               "J1_7_2_S11_L007_R1_001","J1_7_3_S12_L007_R1_001",
                               "J3_2_1_S4_L007_R1_001","J3_2_2_S5_L007_R1_001",
                               "J3_2_3_S6_L007_R1_001","J3_7_2_S14_L007_R1_001",
                               "J3_7_3_S15_L007_R1_001","J4_2_1_S7_L007_R1_001",
                               "J4_2_2_S8_L007_R1_001","J4_2_3_S9_L007_R1_001",
                               "J4_7_1_S16_L007_R1_001","J4_7_2_S17_L007_R1_001",
                               "J4_7_3_S18_L007_R1_001")

# or
colnames(gene_expression) <- c("J1_2_1","J1_2_2","J1_2_3","J1_7_1","J1_7_2","J1_7_3",
                               "J3_2_1","J3_2_2","J3_2_3","J3_7_1","J3_7_2","J3_7_3","J4_2_1",
                               "J4_2_2","J4_2_3","J4_7_1","J4_7_2","J4_7_3")
                               
                              
# View the row names
row.names(gene_expression)

# Determine the dimensions of the dataframe. ‘dim()’ will return the number of rows and columns
dim(gene_expression)

# Assign colors to each. You can specify color by RGB, Hex code, or name To get a list of color names:
data_colors=c("121","122","123","171","172","173","321","322","323","371","372","373",
              "421","422","423","471","472","473")

# View expression values for the transcripts of a particular gene e.g “MSTRG.12187”, 
  # then display only those rows of the data.frame
i = row.names(gene_expression) == "MSTRG.12187"
gene_expression[i,]

# if we want to view values for a list of genes of interest all at once? e,g: 
        #“MSTRG.12017” “MSTRG.1204” “MSTRG.11942” “MSTRG.11963”
genes_of_interest = c("MSTRG.12017", "MSTRG.1204", "MSTRG.11942", "MSTRG.11963")
i = which(row.names(gene_expression) %in% genes_of_interest)
gene_expression[i,]

# Load the transcript to gene index from the ballgown object
transcript_gene_table = indexes(bg)$t2g
head(transcript_gene_table)

#Each row of data represents a transcript. Many of these transcripts represent 
  #the same gene. Determine the numbers of transcripts and unique genes
length(row.names(transcript_gene_table)) #Transcript count
length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count

#Plot #1 - the number of transcripts per gene.
#Many genes will have only 1 transcript, some genes will have several transcripts
#Use the ‘table()’ command to count the number of times each gene symbol occurs 
#(i.e. the # of transcripts that have each gene symbol) Then use the ‘hist’ command 
        #to create a histogram of these counts How many genes have 1 transcript? 
                #More than one transcript? What is the maximum number of transcripts 
        #for a single gene?
counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)

#To plot in old version R
dev.new()
dev.off()

#Plot #2 - the distribution of transcript sizes as a histogram In this analysis 
#we supplied StringTie with transcript models so the lengths will be those of known 
#transcripts However, if we had used a de novo transcript discovery mode, this step 
#would give us some idea of how well transcripts were being assembled If we had a low 
#coverage library, or other problems, we might get short ‘transcripts’ that are actually 
#only pieces of real transcripts
full_table <- texpr(bg , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

# Summarize FPKM values for all samples What are the minimum and maximum FPKM values for a particular library
min(gene_expression[,"J1_2_1_S1_L007_R1_001"])
max(gene_expression[,"J1_2_1_S1_L007_R1_001"])

# Set the minimum non-zero FPKM values for use later. Do this by grabbing a copy of 
# all data values, coverting 0’s to NA, and calculating the minimum or all non NA values
min_nonzero=1

# Set the columns for finding FPKM and create shorter names for figures
data_columns=c(1:18)
short_names=c("J12_1","J12_2","J12_3","J17_1","J17_2","J17_3","J32_1","J32_2",
              "J32_3","J37_1","J37_2","J37_3","J42_1","J42_2","J42_3","J47_1","J47_2","J47_3")

#Plot #3 - View the range of values and general distribution of FPKM values for all 
#libraries Create boxplots for this purpose Display on a log2 scale and add the minimum non-zero value to avoid log2(0)
boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, 
        names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 17 libraries")

#Compare the correlation ‘distance’ between all replicates Do we see the expected 
#pattern for all libraries (i.e. replicates most similar, then DS vs. WW)? Calculate 
#the FPKM sum for all 17 libraries
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)
i = which(gene_expression[,"sum"] > 5)
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
r

#Plot #8 - Convert correlation to ‘distance’, and use ‘multi-dimensional scaling’ 
#to display the relative differences between libraries This step calculates 2-dimensional 
#coordinates to plot points for each library Libraries with similar expression patterns 
#(highly correlated to each other) should group together What pattern do we expect to see, 
#given the types of libraries we have (technical replicates, biologal replicates, DS/WW)?
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes) for all libraries", xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

# Calculate the differential expression results including significance
results_genes = stattest(bg_filt, feature="gene", covariate="condition", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

# Plot #9 - View the distribution of differential expression values as a histogram 
# Display only those that are significant according to Ballgown
sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) Ctrl vs Lso", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

# Plot #10 - Display the grand expression values from Ctrl and Lso and mark those that are significantly differentially expressed
gene_expression[,"CT"]=apply(gene_expression[,c(1:3)], 1, mean)
gene_expression[,"Lso"]=apply(gene_expression[,c(3:6)], 1, mean)
x=log2(gene_expression[,"CT"]+min_nonzero)
y=log2(gene_expression[,"Lso"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="CT (log2)", ylab="Lso FPKM (log2)", main="CT vs Lso FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

# Run the PCA
bg = ballgown(samples=as.vector(sample_full_path),pData=pheno_data)
y = log2(texpr(bg)+1)
pca = prcomp(t(y))
pcmat = pca$x
# make a new data frame with the data we want
install.packages("stringr")
library(stringr)
pcmat = data.frame(Name=str_replace(rownames(pcmat), 'FPKM.', ''), PC1=pcmat[,1], 
                   PC2=pcmat[,2], group=pheno_data$condition)
install.packages("ggrepel")
install.packages("ggplot2")
library(ggrepel)
library(ggplot2)
# plot the PCA with nice labels
ggplot(data=pcmat, aes(x=PC1, y=PC2, colour=group, label=Name)) + geom_point() + geom_label()

# download table from R
write.table(pcmat, file= "/home/tangxt1_1990/ballgown_analysis/pcmat1.csv") 
# and then change the column name


# Write a simple table of differentially expressed transcripts to an output file 
# Each should be significant with a log2 fold-change >= 2
sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
sigde = which(abs(sigp[,"de"]) >= 2)
sig_tn_de = sigp[sigde,]

# Order the output by or p-value and then break ties using fold-change
o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)
output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
write.table(output, file="SigDE.txt", sep="\t", row.names=FALSE, quote=FALSE)
#View selected columns of the first 25 lines of output
output[1:25,c(1,4,5)]

# Create a heatmap to vizualize expression differences between the eight samples

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

main_title="sig DE Transcripts"
par(cex.main=0.8)
sig_genes_de=sig_tn_de[,"id"]
sig_gene_names_de=sig_tn_de[,"gene_name"]

data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
install.packages("gplots")
library(gplots)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,col=rev(heat.colors(75)))

