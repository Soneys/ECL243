# Juliet M Wong
# limma (Ritchie et al. 2015)
# Wong et al. 2017
# Transcriptomics reveal transgenerational effects in purple sea urchins exposed to upwelling conditions, and the response of their progeny to differential pCO2 levels

# Download limma and load required libraries
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR") # analysis performed in limma, but uses a few of the functions in edgeR
library("edgeR")
biocLite("limma")
biocLite("statmod")
library("limma")
library("statmod")
library("ggplot2")

# Read in counts data (.matrix files)
GeneCounts <- read.delim("Gene_Results_Single.counts.matrix", header=TRUE, row.names=1)

All<-GeneCounts[ ,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)] 

# Filtering for more that 0.5 cpm across all samples
keep <- rowSums(cpm(All)>0.5) >=24
y<- All[keep,] 
# Check how many sequences remain after filtering
dim(y)

# Read in the model (contains Sample name, Parent treatment, Offspring treatment, Biological replicate #s)
model<-read.delim("model.txt", header=TRUE)
# Create a treatment factor (Parental+Offspring treatment)
f <- paste(model$ParTreat, model$OffTreat,sep=".")
f <- factor(f)

# Create a model matrix of samples x treatment
design <-model.matrix(~0+f)
colnames(design) <- levels(f)
design

# Apply TMM normalization to RNA-seq read counts
# Create a DGEList object
dge <- DGEList(counts=y, group=f)
dge <- calcNormFactors(dge)
dge$samples
# Create MDS plot (not yet voom-transformed)
colors <- rep(c("purple","blue","red","orange"))
colorsECS <- rep(c("pink","green","black","brown")) #alteration of original code
plotMDS(dge, labels=1:24,col=colors[f])

# voom transformation using voomWithQualityWeights function
# combine observational-level weighting from voom with sample-specific quality weights
v1 <- voomWithQualityWeights(dge, design=design, lib.size=dge$samples$lib.size, plot = TRUE)
labelMDS <- rep(c("NH","NL","UH","UL"))
plotMDS(v1, labels=labelMDS[f],col=colors[f])

# Run voom again, this time adding a blocking variable and estimating correlation
# Samples are blocked by biological replicate
corfit <- duplicateCorrelation(v1,design,block=model$BioRep) 
corfit$consensus
v2 <- voomWithQualityWeights(dge,design,plot=TRUE,lib.size=dge$samples$lib.size,block=model$BioRep,correlation=corfit$consensus)

# MDS plot of voom-corrected data
par(xpd=T, mar=par()$mar+c(0,0,0,4))
mdsplot <- plotMDS(v2,pch=16,col=colorsECS[f])
plot(mdsplot, pch=16, col=colorsECS[f],xlab="Leading logFC dim 1", ylab="Leading logFC dim 2")
legend(2, 0 ,legend = c("NL","NH","UL","UH"), col=colorsECS , pch=16, bty="n")
dev.off()

# PCA across voom-corrected data
pcavoom<- prcomp(t(na.omit(v2$E)))
summary(pcavoom)
pca <- as.data.frame(prcomp(t(na.omit(v2$E)))$x)
pca$f<-f
shapes <- c(rep("N",3), rep("U",3),rep("N",3), rep("U",3),rep("N",3), rep("U",3),rep("N",3), rep("U",3))
pcaplot <- qplot(x=PC1, y=PC2, data=pca,colour=pca$f, size=I(4), shape=shapes)
pcaplot <- pcaplot + scale_colour_manual(values=colorECS), breaks=c("N.L","N.H","U.L","U.H"), name=NULL, labels=c("NL","NH","UL","UH"))
pcaplot <- pcaplot + theme_bw()
pcaplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=14),legend.text=element_text(size=14),
               panel.background = element_blank(), axis.line = element_line(color = "black"), axis.text.y = element_text(angle = 90), legend.key = element_blank())

# Apply limma pipelines for differential expression
fit <- lmFit(v2,design,block=model$BioRep,correlation=corfit$consensus) 
fit <- eBayes(fit) 
v2
All_data_output <- topTable(fit,coef=ncol (design), adjust.method="BH", number=2000, p.value=0.05) 
write.table(All_data_output, "All_data_output.txt")

# Define a contrast matrix to make all pairwise comparisons of interest
cont.matrix <- makeContrasts(U.L-N.L, U.H-N.H, N.H-N.L, U.H-U.L, levels=design)
cont.matrix

# Extract the linear model fit for the contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2)

# Assess differential expression
colnames(fit2)
Allcomp <- topTable(fit2,adjust="BH", number=Inf)

# UL vs NL
UL.vs.NL <- topTable(fit2,coef=1,adjust="BH", number=Inf, p.value=0.05)
write.table(UL.vs.NL, "UL.vs.NL.txt")

# UH vs NH
UH.vs.NH <- topTable(fit2,coef=2,adjust="BH", number=Inf, p.value=0.05)
write.table(UH.vs.NH, "UH.vs.NH.txt")

# NH vs NL
NH.vs.NL <- topTable(fit2,coef=3,adjust="BH", number=Inf, p.value=0.05)
write.table(NH.vs.NL, "NH.vs.NL.txt")

# UH vs UL
UH.vs.UL <- topTable(fit2,coef=4,adjust="BH", number=Inf, p.value=0.05)
write.table(UH.vs.UL, "UH.vs.UL.txt")

# Store results and summarize
results <- decideTests(fit2)
summary(results)
vennDiagram(results, include=c("down","up"),counts.col=c("red","dark blue"), circle.col=c(colorsECS), cex=1.5, names=c("UL vs. NL","UH vs. NH","NH vs. NL","UH vs. UL"))

# To make a Grouped Bar Plot of count data
DEcounts <- structure(list(UL.vs.NL= c(3000,2005), UH.vs.NH = c(159,238), NH.vs.NL = c(504,231), UH.vs.UL = c(1355,1322)),
                  .Names = c("UL vs. NL", "UH vs. NH", "NH vs.NL", "UH vs. UL"), class = "data.frame", row.names = c(NA,-2L))
attach(data)
print(data)
barcolors <- c("pink", "green")
barplot(as.matrix(DEcounts), ylab = "Gene counts", ylim = c(0,3500),cex.lab = 1.2, beside=TRUE, col=barcolors)
legend(7,3300, c("Down-regulated","Up-regulated"), cex=1, bty="n", fill=barcolors)

# LogFC comparison between NH vs. NL and UH vs. UL
NH.vs.NL <- topTable(fit2,coef=3,adjust="BH", number=Inf, sort.by = "none")
UH.vs.UL <- topTable(fit2,coef=4,adjust="BH", number=Inf, sort.by = "none")
# Create a factor for no significant expression for either (0), sig de for NH.vs.NL (1),
# sig de for UH.vs.UL (2), and sig de in both NH.vs.NL and UH.vs.UL (3)
a <- as.factor((NH.vs.NL$adj.P.Val < 0.05) + (UH.vs.UL$adj.P.Val < 0.05)*2)
# plot the data
plot(NH.vs.NL$logFC, UH.vs.UL$logFC, xlab = "NH vs. NL Log-fold Change in Expression", ylab = "UH vs. UL Log-fold Change in Expression", 
     pch = 20, col='gray', type='n')
#points(NH.vs.NL$logFC [a==0], UH.vs.UL$logFC[a==0], col = 'gray', pch = 20)
points(NH.vs.NL$logFC [a==1], UH.vs.UL$logFC[a==1], col = 'green', pch = 20)
points(NH.vs.NL$logFC [a==2], UH.vs.UL$logFC[a==2], col = 'pink', pch = 20)
points(NH.vs.NL$logFC [a==3], UH.vs.UL$logFC[a==3], col = 'yellow', pch = 20)
abline(h=0, v=0)
abline(0,1)
legend(1,-1.5, c("UH vs. UL only","NH vs. NL only","Both comparisons"), col=c("pink","green", "yellow"), cex = 0.5, pch = c(20,20,20), bty="n")

