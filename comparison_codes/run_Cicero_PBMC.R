rm(list=ls())
library(cicero)
library(monocle)
library(Matrix)
library(data.table)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
options(stringsAsFactors = FALSE)
inputdata.10x <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
# extract ATAC data
atac_counts <- inputdata.10x$Peak
### randomly select cells 
cell_sample <- sample.int(ncol(atac_counts))
counts <- atac_counts[,cell_sample[1:2000]]
s <- rowSums(counts)
index <- which(s > 100)
counts <- counts[index,]
# convert to hg19 first
peaks <- rownames(counts)
peaks <- gsub("-",":",peaks)
write.table(peaks, file = "peaks_for_cicero.txt",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
#cat peaks_for_cicero.txt | cut -d':' -f1 > cicero_peaks_chr.txt
#cat peaks_for_cicero.txt | cut -d':' -f2 > cicero_peaks_starts.txt
#cat peaks_for_cicero.txt | cut -d':' -f3 > cicero_peaks_ends.txt
chrs <- read.table(file = "cicero_peaks_chr.txt",sep = "\t")
starts <- read.table(file = "cicero_peaks_starts.txt",sep = "\t")
ends <- read.table(file = "cicero_peaks_ends.txt",sep = "\t")
chep <- data.frame(R1.chrom = chrs$V1,R1.start = starts$V1,R1.end = ends$V1)

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(magrittr)
options(stringAsFactors=FALSE)

colnames(chep) <-c("R1.chrom","R1.start","R1.end")
head(chep)
# convert to hg19
## will need a liftover
#f <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
#system(paste("wget",f))
#system(paste("gunzip",basename(f)))
#f <- sub(".gz","",f)
f <- "hg38ToHg19.over.chain.gz"
f <- sub(".gz","",f)
liftover.chain<-basename(f) %>% import.chain() 
basename(f) %>% unlink()

chep$line<-1:nrow(chep)
chepr1<-with(chep,GRanges(seqnames=Rle(R1.chrom),ranges=IRanges(start=R1.start,end=R1.end),id=line))
chepr1.19<-unlist(liftOver(chepr1,liftover.chain))  
df1 <-  data.frame(iranges = chepr1.19)
# export bed file
head(df1)
# just take out the corresponding peak
index <- which(!duplicated(df1$iranges.id))
hg38_id <- df1$iranges.id[index]
counts_hg19 <- counts[hg38_id,]
new_peaks <- paste0(df1$iranges.seqnames,"-",df1$iranges.start,"-",df1$iranges.end)
rownames(counts_hg19) <- new_peaks[index]
counts_frag <- as.data.frame(as.table(as.matrix(counts_hg19)))
counts_frag$Var1 <- as.character(counts_frag$Var1)
counts_frag$Var2 <- as.character(counts_frag$Var2)
input_cds <- make_atac_cds(counts_frag, binarize = TRUE)

set.seed(2017)
peaks_hg19 <- rownames(counts_hg19)
rm(atac_counts)
library(monocle3)
input_cds <- monocle::detectGenes(input_cds)
#input_cds <- estimate_size_factors(input_cds)
#input_cds <- preprocess_cds(input_cds, method = "LSI")
#input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
#                              preprocess_method = "LSI")
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none")

tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- colnames(input_cds)
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

#D <- assay(cicero_cds)
#D <- as.matrix(D)
#D[1:5,1:5]
#write.table(D,file = "Cicero_processed_data.txt",sep = "\t")
data("human.hg19.genome")

conns_cicero <- run_cicero(cicero_cds,human.hg19.genome)
head(conns_cicero)
conns_cicero$coaccess[is.na(conns_cicero$coaccess)] <- 0
write.table(conns_cicero, file = "PBMC_cicero_connections_hg19.txt",sep = "\t",quote = FALSE)

### convert to hg38
# based on the original corresponing relationships
conns_cicero$Peak2 <- as.character(conns_cicero$Peak2)

peaks_hg38 <- peaks[hg38_id]
conns_cicero_hg38 <- conns_cicero
id1 <- match(conns_cicero$Peak1,peaks_hg19)
f1 <- which(is.na(id1))
length(f1)

id2 <- match(conns_cicero$Peak2,peaks_hg19)
f2 <- which(is.na(id2))
length(f2)
conns_cicero_hg38$Peak1 <- peaks_hg38[id1]
conns_cicero_hg38$Peak2 <- peaks_hg38[id2]
write.table(conns_cicero_hg38, file = "PBMC_cicero_connections.txt",sep = "\t",quote = FALSE)
save(conns_cicero_hg38, file = "PBMC_cicero_connections.RData")

