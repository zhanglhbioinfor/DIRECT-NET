rm(list=ls())
suppressPackageStartupMessages(library(ArchR))
addArchRGenome("hg38")
addArchRThreads(1)
#Get Input Fragment Files
inputFiles <- getInputFiles("/Users/lihuazhang/Documents/DIRECTNET_revise/compare/PBMC")
inputFiles
names(inputFiles) <- c("PBMC")
#reformatFragmentFiles(inputFiles)
# remove original and change reformat name
ArrowFiles <- createArrowFiles(inputFiles,force = TRUE)
proj1 <- ArchRProject(ArrowFiles)
proj1
getAvailableMatrices(proj1)
head(proj1$cellNames)
quantile(proj1$TSSEnrichment)
df <- getCellColData(proj1, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
p
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj1, addDOC = FALSE)

p1 <- plotGroups(
  ArchRProj = proj1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1

p2 <- plotGroups(
  ArchRProj = proj1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2


p1 <- plotFragmentSizes(ArchRProj = proj1)
p1
p2 <- plotTSSEnrichment(ArchRProj = proj1)
p2
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


## filter remove doublets
proj2 <- addDoubletScores(proj1)
proj2 <- filterDoublets(proj2)

proj2 <- addIterativeLSI(
  ArchRProj = proj2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2,
  # clusterParams = list( #See Seurat::FindClusters
  #  resolution = c(0.2), 
  #  sampleCells = 10000, 
  #  n.start = 10
  #), 
  varFeatures = 25000, # 25000
  dimsToUse = 1:30,
  force = TRUE
)

### clustering
proj2 <- addClusters(
  input = proj2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1.2, force = TRUE
)
table(proj2$Clusters)

### Umap
proj2 <- addUMAP(
  ArchRProj = proj2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)



# pseudo-bulk data
proj3 <- addGroupCoverages(ArchRProj = proj2, groupBy = "Clusters") # based on scATAC-seq cluster
pathToMacs2 <- "/Users/lihuazhang/opt/anaconda3/bin/MACS2"
proj3 <- addReproduciblePeakSet(
  ArchRProj = proj3, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2
)
getPeakSet(proj3)

proj4 <- addPeakMatrix(proj3)

proj4 <- addCoAccessibility(
  ArchRProj = proj4,
  reducedDims = "IterativeLSI",
  maxDist = 500000
)
cA <- getCoAccessibility(
  ArchRProj = proj4,
  corCutOff = -1,
  resolution = 1,
  returnLoops = FALSE
)

# integrate with scRNA-seq data
pbmc <- readRDS(file = "pbmc_rna.rds")

proj5 <- addGeneIntegrationMatrix(
  ArchRProj = proj4, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = pbmc,
  addToArrow = TRUE,
  force= TRUE,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)
########### here need recall peaks????

proj5 <- addPeak2GeneLinks(
  ArchRProj = proj5,
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = proj5,
  corCutOff = -1,
  resolution = 1,
  returnLoops = FALSE
)


peaks <- as.character(metadata(cA)[[1]])
peaks <- as.vector(peaks)
peaks <- gsub("-", "_", peaks)
peaks <- gsub(":", "_", peaks)
# construct dataframe of the links
links <- data.frame(Peak1 = peaks[cA@listData$queryHits],Peak2 = peaks[cA@listData$subjectHits],coaccess = cA@listData$correlation)
save(links,file = "PBMC_ArchR_links.RData")

peaks <- as.character(metadata(p2g)[[1]])
peaks <- as.vector(peaks)
peaks <- gsub("-", "_", peaks)
peaks <- gsub(":", "_", peaks)
genes_all <- proj5@geneAnnotation$genes@elementMetadata@listData$symbol
genes_all <- as.vector(genes_all)

# construct dataframe of the links
links <- data.frame(Peak1 = genes_all[p2g@listData$idxRNA],Peak2 = peaks[p2g@listData$idxATAC],coaccess = p2g@listData$Correlation)
save(links,file = "PBMC_ArchR_links_genes.RData")
