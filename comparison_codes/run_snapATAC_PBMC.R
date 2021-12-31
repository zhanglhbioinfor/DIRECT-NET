# run SnapATAC
rm(list=ls())
library(SnapATAC);
library(ggplot2)
snap.file = "PBMC_data.snap"
x.sp = createSnap(file=snap.file, sample="PBMC")
# barcodes 
library(data.table)
fragments <- data.table::fread(file = "pbmc_granulocyte_sorted_10k_atac_fragments.tsv")
fragments[1:5,]
barcodes_name <- unique(fragments$V4)
showBinSizes("PBMC_data.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)
x.sp = makeBinary(x.sp, mat="bmat")#binary
library(GenomicRanges)
# remove balacklist
black_list = read.table("hg38.blacklist.bed");
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
)
idy = queryHits(
  findOverlaps(x.sp@feature, black_list.gr)
)
if(length(idy) > 0){
  x.sp = x.sp[,-idy, mat="bmat"];
}
x.sp

# remove unwanted 
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

# dimension reduction
x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat", 
  num.eigs=50
)
plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
)
# graph based clustering
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:20,
  k=15
)
x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=10
)
x.sp@metaData$cluster = x.sp@cluster

# plot reuslts
x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20, 
  method="Rtsne",
  seed.use=10
)



plotViz(
  obj=x.sp,
  method="tsne", 
  main="GM12878",
  point.color=x.sp@cluster, 
  point.size=1, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
)

# scRNA-seq based annotation
library(Seurat)
library(dplyr)
library(ggplot2)
inputdata.10x <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
pbmc.rna <- CreateSeuratObject(counts = rna_counts)
pbmc.rna[["percent.mt"]] <- PercentageFeatureSet(pbmc.rna, pattern = "^MT-")
plot1 <- FeatureScatter(pbmc.rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc.rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
pbmc.rna <- subset(pbmc.rna, subset = nFeature_RNA > 1000 & nFeature_RNA < 25000 & percent.mt < 20)
pbmc.rna <- NormalizeData(pbmc.rna, normalization.method = "LogNormalize", scale.factor = 10000)
#Variable Genes for Patient A (identify, plot, and use to run PCA)
pbmc.rna <- FindVariableFeatures(pbmc.rna, selection.method = "vst", nfeatures = 2000)

variable.genes = VariableFeatures(object = pbmc.rna)
genes.df = data.table::fread("gencode.v38.annotation.gene.bed")
gene_names <- unlist(genes.df[,10])
c <- strsplit(gene_names,split = ";")
# extract gene names
library(stringr)
names <- rep(NA)
for (i in 1:length(c)) {
  print(i)
  c1 <- str_subset(c[[i]], "gene_name")
  names[i] <- substr(c1,13,(str_length(c1)-1))
}

genes.gr = GRanges(as.vector(unlist(genes.df[,1])), IRanges(as.vector(unlist(genes.df[,2])), as.vector(unlist(genes.df[,3]))), name=names)
genes.sel.gr = genes.gr[which(genes.gr$name %in% variable.genes)]
## reload the bmat, this is optional but highly recommanded
x.sp = addBmatToSnap(x.sp)
x.sp = createGmatFromMat(
  obj=x.sp, 
  input.mat="bmat",
  genes=genes.sel.gr,
  do.par=TRUE,
  num.cores=6
)
#### prepare for integration
pbmc.atac <- snapToSeurat(
  obj=x.sp, 
  eigs.dims=1:20, 
  norm=TRUE,
  scale=FALSE
)
transfer.anchors <- FindTransferAnchors(
  reference = pbmc.rna, 
  query = pbmc.atac, 
  features = variable.genes, 
  reference.assay = "RNA", 
  query.assay = "ACTIVITY", 
  reduction = "cca"
)

pbmc.rna <- ScaleData(object = pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna, pc.genes = pbmc.rna@var.genes, pcs.compute = 40, do.print = FALSE)
PCAPlot(object = pbmc.rna)
numPC <- 30
pbmc.rna <- FindNeighbors(pbmc.rna, reduction = "pca", dims = 1:numPC)
pbmc.rna <- FindClusters(object = pbmc.rna, resolution = 0.9,algorithm = 1)

pbmc.rna$celltype <- Idents(pbmc.rna)
celltype.predictions <- TransferData(
  anchorset = transfer.anchors, 
  refdata = pbmc.rna$celltype,
  weight.reduction = pbmc.atac[["SnapATAC"]],
  dims = 1:20
)
x.sp@metaData$predicted.id = celltype.predictions$predicted.id
x.sp@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max)
x.sp@cluster = as.factor(x.sp@metaData$predicted.id)

refdata <- GetAssayData(
  object = pbmc.rna, 
  assay = "RNA", 
  slot = "data"
)
imputation <- TransferData(
  anchorset = transfer.anchors, 
  refdata = refdata, 
  weight.reduction = pbmc.atac[["SnapATAC"]], 
  dims = 1:20
)
x.sp@gmat = t(imputation@data)
rm(imputation) # free memory
rm(refdata)    # free memory
rm(pbmc.rna)   # free memory
rm(pbmc.atac) # free memory

hist(
  x.sp@metaData$predict.max.score, 
  xlab="prediction score", 
  col="lightblue", 
  xlim=c(0, 1),
  main="PBMC 10X"
)
abline(v=0.5, col="red", lwd=2, lty=2)
table(x.sp@metaData$predict.max.score > 0.5)
x.sp = x.sp[x.sp@metaData$predict.max.score > 0.5,]
x.sp

# call peaks for all cluster with more than 100 cells
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)]
peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i]);
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("PBMC.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/Users/lihuazhang/opt/anaconda3/bin/snaptools",
    path.to.macs="/Users/lihuazhang/opt/anaconda3/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=1,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder=tempdir()
  );
  peaks
}, mc.cores=5)
peaks.names = system("ls | grep narrowPeak", intern=TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls))

peaks.df = as.data.frame(peak.gr)[,1:3]
write.table(peaks.df,file = "peaks.combined.bed",append=FALSE,
              quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")

x.sp = addPmatToSnap(x.sp)
### predict gene-enhancer interaction for each gene
# here we focus on the genes used in DIRECTNET

x.sp1 = makeBinary(x.sp, mat="pmat")#binary
#### note: the chromsome name include b

peak.use = x.sp1@peak
c <- as.character(peak.use)
c_new <- gsub("b","",c)
c_new <- gsub("'","",c_new)
library(Signac)
peak.use_new <- StringToGRanges(c_new, sep = c(":", "-"))
x.sp1@peak <- peak.use_new

links_record <- list()
for (i in 1:length(TSS_used$genes)) {
  print(i)
  TSS.loci = GRanges(TSS_used$Chrom[i], IRanges(TSS_used$Starts[i], TSS_used$Ends[i]))
  tryCatch({
    pairs = predictGenePeakPair(
      x.sp1, 
      input.mat="pmat",
      gene.name=TSS_used$genes[i], 
      gene.loci=resize(TSS.loci, width=500000, fix="center"),
      do.par=FALSE
    )
    
    # convert the pair to genome browser arc plot format
    pairs.df = as.data.frame(pairs)
    pairs.df = data.frame(
      chr1=pairs.df[,"seqnames"],
      start1=pairs.df[,"start"],
      end1=pairs.df[,"end"],
      chr2=TSS_used$Chrom[i],
      start2=TSS_used$Starts[i],
      end2=TSS_used$Ends[i],
      Pval=pairs.df[,"logPval"]
    )
    links_record[[i]] <- pairs.df
  }, error = function(e){
  })
}
save(links_record, file = "PBMC_snapATAC_pmat_500k.RData")
