library(ape)
library(adegenet)
library(parallel)
# === Specify working paths before running
run.name <- "test"  # basename of BED file
bed <- "/path/to/data/"
plot.path <- "/path/to/save/plots/"
# == Read sample info from FAM
sample.group <- read.table(paste(ped,run.name, ".fam",sep=""))[,1]
samples <- read.table(paste(bed, run.name, ".fam", sep=""))[,2]
# === Renamed copy of BIM file for compatibility
file.copy(paste(ped,run.name, ".bim",sep=""),
          paste(ped,run.name, ".map",sep=""))
# === Reading data
snp.data <-read.PLINK(paste(bed,run.name, "_A.raw",sep=""),
                      map.file=paste(bed,run.name, ".map",sep=""),
                      parallel=TRUE)
snp.data.df <- as.data.frame(snp.data)
row.names(snp.data.df) <- samples
# === Distance matrix and tree calculation
snp.dist <- dist.gene(snp.data.df)
snp.nj <- nj(snp.dist)
write.nexus(snp.nj, file=paste(plot.path, "ovine_tree.nexus")
