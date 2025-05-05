BiocManager::install('WGCNA')
BiocManager::install('flashClust')

library(WGCNA)
library(flashClust)
library(curl)
library(DESeq2)

setwd("~/Documents/Graduate School/Functional Genomics/Final Project")

################### Load in gene cout matrix#########################################
dog.data <- as.matrix(read.csv("dog_gene_count_matrix.csv", row.names="gene_id"))
coldata <-(read.table("~/Documents/Graduate School/Functional Genomics/Final Project/Dog_PHENO_DATA .txt", header=TRUE, row.names=1))
#Deseq2 to normalize counts
dds <- DESeqDataSetFromMatrix(countData = dog.data, colData=coldata,  design = ~treatment)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="dog.data.normalized.txt", sep="\t")

#Data Tranformation
norm.data<- read.delim("dog.data.normalized.csv", sep="\t")
expression.data <- as.data.frame(t(norm.data)) #transforming the data.frame so columns now represent genes and rows represent samples

#Outlier gene removal 
gsg <-goodSamplesGenes(expression.data)
summary(gsg)
gsg$allOK
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
abline(h = 15, col = "red");

#cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10) #returns numeric vector
#Remove outlier
#expression.data <- expression.data[cut.sampleTree==1, ]

spt <- pickSoftThreshold(expression.data) 
spt

par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")

par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")

softPower <- 10
adjacency <- adjacency(expression.data, power = softPower)


TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1-TOM


#creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 
#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module. 

ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)
#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75


merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs


plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")


############# Load in phenotypic trait data ###########################################
traitData <- read.csv("Dog_Traits.csv",header = TRUE, stringsAsFactors = FALSE)
head(traitData)
#allTraits <- traitData[, -c(31, 16)] #removing notes and comments sections 
allTraits <- traitData[, c(1, 5, 6) ] #pulling out only continuous traits 


Samples <- rownames(expression.data)
traitRows <- match(Samples, allTraits$SRR)
datTraits <- allTraits[traitRows, -1]
#rownames(datTraits) <- allTraits[traitRows]



# Define numbers of genes and samples
nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
par(mar = c(6, 8.5, 3, 1))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.correlation,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


################ define module membership and gene trait significance for specific trait######
# Define variable weight containing the size column of datTrait
Size = as.data.frame(datTraits$Size)
names(Size) = "Size"

modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expression.data, Size, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Size), sep="")
names(GSPvalue) = paste("p.GS.", names(Size), sep="")
head(GSPvalue)


##################### Gene intramoular connectivity Analysis correlation of the eigengene and the gene expression profile###########
par(mar=c(4,4,4,4))
module = "indianred4"
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Dog Size",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


par(mar=c(4,4,4,4))
module = "firebrick4"
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


############### Plot dendrogram and heatmap for specific trait ################
# Isolate size from the clinical traits
Size = as.data.frame(datTraits$Size);
names(Size) = "size"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, Size))
# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


# Plot the dendrogram
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)



################## Pull out list of genes in modules of interest ####################

library(tidyverse)

submod = geneModuleMembership %>%
  select(c("MMindianred4", "MMfirebrick4"))
submod<- rownames_to_column(submod, var = "geneIDs")



indianred4_exp<- names(expression.data)[ModuleColors=="indianred4"]
firebrick4_exp<- names(expression.data)[ModuleColors=="firebrick4"]

indian<- submod%>% 
  filter(geneIDs %in% indianred4_exp) %>% select(c(1,2))

fire<-submod%>% 
  filter(geneIDs %in% firebrick4_exp)%>% select(c(1,3))

dog.data.df<- as.data.frame(normalized_counts)
dog.data.named<- rownames_to_column(dog.data.df, var = "geneIDs")

indian_dge<- dog.data.named%>% 
  filter(geneIDs %in% indianred4_exp)

fire_dge<- dog.data.named%>% 
  filter(geneIDs %in% firebrick4_exp)


############## Annotate DEGs for possible use in GSEA #######################
Anno <- read.csv("dog_annotation.csv", stringsAsFactors = FALSE, na.strings=c(""))
summary(Anno)
dim(Anno)

DGEresults <- read.csv("DGE_results_dog.csv", stringsAsFactors = FALSE)
summary(DGEresults)
dim(DGEresults)

## Rename first column so it matches "gene_id" in annotation file
names(DGEresults)[1]<- "gene_id" 

#Merge anno with DGE results
DGE_Anno <- merge(Anno,DGEresults,by="gene_id",all.y = FALSE)
dim(DGE_Anno)
summary(DGE_Anno)



############################# Make ranked list for GSEA ####################
## example aa <- within(resOrdered, z <- x + y - 2)
DGE_Anno_Rank <-  within(DGE_Anno, rank <- sign(log2FoldChange) * -log10(pvalue))
DGE_Anno_Rank 

indian_rank<- DGE_Anno_Rank%>% 
  filter(gene_id %in% indianred4_exp)

fire_rank<- DGE_Anno_Rank%>% 
  filter(gene_id %in% firebrick4_exp)


write.table(as.data.frame(indian_dge), file="indian_DGErankName.rnk", quote=FALSE, row.names=FALSE, sep = "\t")  
write.table(as.data.frame(fire_rank), file="fire_DGErankName.rnk", quote=FALSE, row.names=FALSE, sep = "\t") 

###########################################################################
#GO enrichment analysis of module genes
BiocManager::install("org.Cf.eg.db")
library(clusterProfiler)
library(AnnotationDbi)
library(org.Cf.eg.db)

#select only gene names of log2foldchange of more than 0.5, positive only.
genes_go<- fire_rank[fire_rank$log2FoldChange>0.5, ]

genes_go<-genes_go$gene_name #select only gene_name symbols

#run go enrichment analysis on all go terms
Go_result<- enrichGO(gene=genes_go, OrgDb= org.Cf.eg.db, keyType = "SYMBOL", ont="ALL")
data.frame(Go_result)

#barplot top 20 enriched terms
fit<- plot(barplot(Go_result, showCategory = 20))        

################ WGCNA Code Tutorial Aknowledgment #########################
#Thank you to Victoria French, CeCe Gerstenbacher, Warrenkevin Henderson, and Elizabeth Varghese 11/29/2021 for WGCNA tutorial https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html#Weighted_Gene_Correlation_Network_Analysis

