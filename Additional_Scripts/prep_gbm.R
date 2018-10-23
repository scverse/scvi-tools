rm(list = ls())
library(compositions)
library(Rtsne)
library(data.table)
# library(cellrangerRkit)
library(Matrix)


## ------ RNA UMI -----
RNA_umi <-
  read.table(
    file = "../cite/GSE100866_PBMC_vs_flow_10X-RNA_umi.csv.gz",
    sep = ",",
    header = TRUE,
    row.names = 1
  )

# Calling sample type:
RNA_SPECIES_TYPE <- factor(gsub("_.*", "", rownames(RNA_umi)))
RNA_umi_DT <- data.table(as.matrix(RNA_umi))
RNA_species_sums_DT <- setDT(RNA_umi_DT)[, lapply(.SD, sum), by = RNA_SPECIES_TYPE]
RNA_species_sums = matrix(
  as.matrix(RNA_species_sums_DT[, -1]),
  nrow = nrow(RNA_species_sums_DT),
  dimnames = list(levels(RNA_species_sums_DT$RNA_SPECIES_TYPE),
                  colnames(RNA_species_sums_DT)[-1])
)

RNA_species_comp <- ccomp(t(RNA_species_sums))
sample_species <- factor(c("mouse","mixed","human")[1 + (clo(RNA_species_comp)[,"HUMAN"] > 0.9) + (clo(RNA_species_comp)[,"HUMAN"] > 0.1)])
save(sample_species, file = "../cite/sample_species.rda")

## ------ ADT UMI -----
ADT_umi <- read.table(file = "../cite/GSE100866_PBMC_vs_flow_10X-ADT_umi.csv.gz",
                      sep = ",", header = TRUE,row.names = 1)

stopifnot(all(colnames(ADT_umi) == colnames(RNA_umi)))

# CLR ADT (0 values ignored for g calculation and mapped to zero)
ADT_comp <- ccomp(t(as.matrix(ADT_umi)))
plot(missingSummary(ADT_comp))
ADT_clr <- t(clo(clr(ADT_comp)))

# Mouse-derived Cutoffs
ADT_cut <- apply(ADT_clr[,sample_species == "mouse"],1,mean) + 
  apply(ADT_clr[,sample_species == "mouse"],1,sd)
ADT_cut_clr = ADT_clr - ADT_cut

boxplot(t(ADT_clr[,sample_species == "human"]), main = "CLR")
print(paste("CLR IQR before thresh:",IQR(ADT_clr[,sample_species == "human"])))

boxplot(t(ADT_cut_clr[,sample_species == "human"]),  main = "CLR w/ thresh")
print(paste("CLR IQR after subtracting mouse cutoff:",IQR(ADT_cut_clr[,sample_species == "human"])))

# tSNE
set.seed(101)
tsne_obj <- Rtsne(t(ADT_cut_clr),
                  pca = FALSE, verbose = TRUE)
plot(tsne_obj$Y, pch = 16, cex = .5, col = sample_species)

tsne_obj_prot_all <- tsne_obj
save(tsne_obj_prot_all, file = "../cite/tsne_obj_prot_all.rda")

# ----- Save Data -----

# RNA UMI Matrix
filt_dat_mat0 <- RNA_umi[RNA_SPECIES_TYPE == "HUMAN",][,sample_species == "human"]
rownames(filt_dat_mat0)<- gsub("HUMAN_","",rownames(filt_dat_mat0))
filt_dat_mat <- as(as.matrix(filt_dat_mat0), "dgTMatrix")
dimnames(filt_dat_mat) <- dimnames(filt_dat_mat0)

# QC for Cells
qc = RNA_species_comp[sample_species == "human",]
colnames(qc) <- paste0("qc_num_",colnames(qc),"_umi")

# New GBM Object
pd = data.frame(cellid = colnames(filt_dat_mat),t(ADT_cut_clr[,sample_species == "human"]))
                
fd = data.frame(rownames(filt_dat_mat), row.names = rownames(filt_dat_mat))
colnames(fd) <- c("gene_symbol")

write.csv(pd,file='../cite/ADT_cut_clr.csv',row.names=F,quote=F)
write.csv(fd,file='../cite/genenames.csv',row.names=F,quote=F)
writeMM(filt_dat_mat,file='../cite/count.mtx')
gbm = newGeneBCMatrix(filt_dat_mat, fd = fd, pd = pd)
save(gbm, file = "data/gbm.rda")
print("Done!")
