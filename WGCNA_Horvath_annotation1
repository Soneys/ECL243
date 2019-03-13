#=====================================================================================
#
#  Code chunk 1 - Data grab
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
setwd("~/Documents/UC Davis/Winter 2019/ECL 243/Project/WGCNA_tutorial/");
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv");

# Take a quick look at what is in the data set:
dim(femData);
names(femData);


#=====================================================================================
#
#  Code chunk 2 - Place in a new working df (datExpr0)
#
#=====================================================================================


datExpr0 = as.data.frame(t(femData[, -c(1:8)])); #all rows for all columns except 1-8
names(datExpr0) = femData$substanceBXH; #create a vector out of column 1
rownames(datExpr0) = names(femData)[-c(1:8)]; #set the names in the rows of the new data frame to column1 vector
# output is a df with samples in each row and "BXH compound" (treatment) in each column

#=====================================================================================
#
#  Code chunk 3 - check for missing/poor data
#
#=====================================================================================

  # goodSamplesGenes - checks data for missing entries, entries with weights below a threshold, and 
  # zero-variance genes, and returns a list of samples and genes that pass criteria on maximum number 
  # of missing or low weight values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK #this is a logical of length 1; binary - all ok or no?
# goodGenes - logical for each gene, good or no
# goodSamples - logical for each, good or no
# allOK - binary logical response - TRUE if all genes/samples good, FALSE if not


#=====================================================================================
#
#  Code chunk 4 - remove missing/poor data
#
#=====================================================================================


if (!gsg$allOK) #if gsg$allOK = FALSE, ie there are genes/samples that don't meet criteria:
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) # if more than 0 genes did not make it into goodGenes, 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    # remove them verbosely, but printFlush flushes the statement after erasure complete (regarding verbose)
  if (sum(!gsg$goodSamples)>0)
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # same for bad samples, if more that 0 didn't make it into goodSamples remove from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  # set datExpr0 equal to itself with only good genes/samples
}


#=====================================================================================
#
#  Code chunk 5 - cluster rows (genes) using dissimilarity matrix
#
#=====================================================================================

# hclust from fastcluster package
# agglomerative (bottom up) method - 
    # dist(data) computes distance between rows of a data matrix, returns a "distance matrix"
    # hclust(matrix, method) - average method is UPGMA (ultrametric tree - the distance from the 
    # root to the branch tips is the same?)
  # UPGMA - at each step, two nearest clusters are combined into a cluster (lowest dist in distmatrix).
  # The dist between two clusters a and b is taken to be the average of all the distances d(x,y)
  # between the pairs of objects x in a and y in b, that is the mean distance between any two elements in a cluster
    # In this case, the distance matrix is constructed as the 
sampleTree = hclust(dist(datExpr0), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut
#  abline just adds a horizontal line with and intercept h; how is a cutoff value of 15 here decided?
#  must relate to the height of the dendrogram, has it been visually determined?
abline(h = 15, col = "red");
# Determine cluster under the line
# cutTreeStatic - Module detection in hierarchical dendrograms using a constant-height tree cut.
# Only branches whose size is at least minSize are retained
  # returns a vector with binary value for each sample where height >15 = 1, not >15 = 0
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
# table returns how many smaples are 0,1
table(clust)
# clust 1 contains the samples we want to keep.
# clust==1 returns a TF vector of whether a sample is 1 or not (above cutoff or not)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ] # new df, with only samples that made the cutoff
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

ncol(datExpr0) == ncol(datExpr) # should be true, no columns (treatments) were culled
nrow(datExpr0) == nrow(datExpr) # should be false, removed 1 sample here


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

# Wong et al. create datTraits in the beginning of their script

traitData = read.csv("ClinicalTraits.csv");
# names(traitData) = colnames(traitData)
# 361 rows - 1-361; there are nrow(datExpr)
# cells contain values of the traits defined in the column names
# what do the rows mean then??
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)]; # all rows in columns that are not 31 0r 16; 
# "Notes" and "comments"
allTraits = allTraits[, c(2, 11:36) ]; # all rows in columns 2, 11-36 (mice, some of the measures)
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr); #vector of the samples
traitRows = match(femaleSamples, allTraits$Mice);
  # match returns a vector of the positions of (first) matches of its first argument 
  # in its second. x %in% table
datTraits = allTraits[traitRows, -1]; #remove first column
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage(); # Performs garbage collection until free memory indicators show no change



#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples on the df without the sample that didn't make the cutoff
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
# red = high value of trait is observed in a sample, white = low value. Does this apply only to cont numerical? 
# did we crop the ones that were discrete/non-numerical values?
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")
save(datExpr0, datTraits, file = "Urchin-dataInput.RData")


