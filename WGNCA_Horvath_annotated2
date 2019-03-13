#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
    #enableWGCNAThreads() #replace with below, updated function
allowWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
lnames = load(file = "Urchin-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames

#So start with a df for expr data (samples x treatment) and trait data (based on Wong ea,
# appears to relate which samples were in which treatment group?)


#=====================================================================================
#
#  Code chunk 2 
#
#=====================================================================================

# Thresholding: hard threshold is: absolute value of correlation matrix, choose a cutoff 
# (e.g. 0.85), anything above is considered connected in the network. But then there is 
# soft thresholding which is when you exponentiate the correlation matrix and that 
# accentuates larger connections
# The soft thresholding, is a value used to power the correlation of the genes to that 
# threshold. The assumption on that by raising the correlation to a power will reduce the 
# noise of the correlations in the adjacency matrix. To pick up one threshold use the 
# pickSoftThreshold function, which calculates for each power if the network resembles to 
# a scale-free graph. The power which produce a higher similarity with a scale-free network 
# is the one you should use.

# Summary - some thresholding power is set to determine whether nodes (genes) are notably
# connected in a network; ie if the expression values across all samples is similar for a given gene?
# the dataExpr df relates gene expression values to samples; so a network could be cosntructed by 
# observing genes whose expression patterns seem to correlate/anticorrelate. Is this averaged across 
# experimental treatments or samples, or done for each sample indpendently?

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#powers = (1:10, 12,14,16,18,20)
# Call the network topology analysis function
# pickSoftThreshold - Analysis of scale free topology for multiple soft thresholding 
# powers. The aim is to help the user pick an appropriate soft-thresholding power 
# for network construction.
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


softPower = 6; # what part of the above graphs made 6 an obvious choice? What do those
# graphs convey?

  softPower=18;#used by Wong et al
adjacency = adjacency(datExpr, power = softPower);
# adjacency - Calculates (correlation or distance) network adjacency from given expression 
# data or from a similarity.
# adjacency refers to the number of edges connected to that node (how many nodes it is 
# connected to?) 

# seeking to know how connected different genes are, which ones create "hubs" or networks of 
# high interconnectivity


#=====================================================================================
# 
#  Code chunk 4
#
#=====================================================================================

# from the name, seems to compare the similarities in topologies of multiple networks.
# likely means the adjacencies/ networks are being constructed for each sample, and then
# this function is actually comparing the shapes of these networks to one another?
# TOMsimilarity - Calculation of the topological overlap matrix, and the corresponding
# dissimilarity, from a given adjacency matrix.

# adjmat input must be a square symmetric adjacency matrix with values between 0 and 1. So 
# adjacency can't be the number of other nodes a node is connected to

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM #dissimalrity matrix - where prev each value is percent agreement of topologies,
# dissTOM is percent difference in topology


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

# clustering the dissimilatiry matrix as a distance matrix? Degree of dissimilarity in the 
# entries (topological networks, 22) is quantified and compared, so the networks with the 
# closest topologies are clustered together


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
# No, dendrogram shows individual genes (many more than 22 tips)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30; #
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Calculate eigengenes
  #WONG sets softPower=18
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
  #WONG USES flashclust instead of hclust, says works the same?
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

### WONG USES dynamicMergeCut first, which Calculates a suitable threshold for module 
# merging based on the number of samples and a desired Z quantile

  #WONG USES many more sig figs in MEDISSThres...
MEDissThres = 0.25
0.2559876
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")

save(MEs, moduleLabels, moduleColors, geneTree, file = "Urchin-networkConstruction-stepByStep.RData")





