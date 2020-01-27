# === Run PCA and create plots.
# === Plink1.9 must be present in $PATH
library(adegenet)
library(parallel)
library(pca3d)
library(heatmap3)

# === Specify working paths before running
run.name <- "test"  # basename of BED file
bed <- "/path/to/data/"
plot.path <- "/path/to/save/plots/"
stat <- "/path/to/save/results/"
# === Convert data for adegenet input
system(paste("plink1.9 --chr-set 26 --bfile ", 
             bed, run.name, 
             " --recodeA",
             " --out ", bed, run.name, "_A", 
             sep=""))
# === Extract breed and population info from FAM
sample.group <- read.table(paste(ped,run.name, ".fam",sep=""))[,1]
breed <- unlist(lapply(strsplit(as.character(sample.group), split="/"), FUN=function(x){return(x[1]}))
# === Renamed copy of BIM file for compatibility
file.copy(paste(ped,run.name, ".bim",sep=""),
          paste(ped,run.name, ".map",sep=""))
# === Read data and perform PCA
snp.data <-read.PLINK(paste(ped,run.name, "_A.raw",sep=""),
                      map.file=paste(ped,run.name, "_SNP_filt.map",sep=""),
                      parallel=TRUE)
snp.pca <- glPca(snp.data, parallel=TRUE, nf=100)

# === Calculate explained and cumulative variance
var.explained <- apply(snp.pca$scores,2,var)/sum(apply(snp.pca$scores,2,var))
var.cumul <- numeric()
for ( i in c(1:length(var.explained))) cumul <- c(cumul, sum(var.explained[1:i]))
write.table(data.frame(var.explained, var.cumul),
            paste(stat,run.name,"_variance.txt"), 
            row.names=FALSE, col.names=TRUE, sep="\t")

# === Plots
# === PCA heatmap
png(paste(plot.path, "PCA_heatmap.png",sep=""), width=4500, height=4500, res=300)
heatmap3(snp.pca$scores[order(sample.group),1:50], Colv=NA, Rowv=NA, 
        labRow=sort(sample.group), col=heat.colors(100),
        margins = c(20,10)
        )
dev.off()
# === Four panel plot of variance and PCA for first 3 PC
# === with one free panel 
png(paste(plot.path, "pca2d.png",sep=""),width=2000,height=2000,res=300)
layout(matrix(c(1:4),nrow=2, ncol=2,byrow=TRUE))
# + Barplot of variance
palette(c("white", "darkgrey"))
barplot(var.cumul, col=rep(c(rep(1,10),rep(2,10)),5))
barplot(var.explained, las=2, col="black",add=TRUE,names.arg=NA, axes=NULL)
abline(h=0.95,col="red")
palette(c("green3", "darkgreen", "mediumseagreen",
          "blue", "deepskyblue2",
          "olivedrab4",
          "darkorange", "darkred","firebrick" , "deeppink",
          "purple"))
pca2d(snp.pca$scores[,1:2],group=breed,
      col=as.numeric(sample.group),
      show.ellipses=TRUE, show.group.labels = FALSE,
      axes.color="black",
      axe.titles = c(paste("PC 1 (", round(var.explained[1]*100,digits=1), "% of variance explained)"),
                     paste("PC 2 (", round(var.explained[2]*100,digits=1), "% of variance exlained)")))
pca2d(snp.pca$scores[,2:3],group=breed,
      col=as.numeric(sample.group),
      show.ellipses=TRUE, show.group.labels = FALSE,
      axes.color="black",
      axe.titles = c(paste("PC 2 (", round(var.explained[2]*100,digits=1), "% of variance explained)"),
                     paste("PC 3 (", round(var.explained[3]*100,digits=1), "% of variance exlained)")))
dev.off()

# === 3D plot of first 3 PC
palette(c("green3", "lawngreen", "darkgreen",
          "blue", "deepskyblue2",
          "saddlebrown",
          "darkorange", "darkred","firebrick" , "deeppink",
          "purple"))
pca3d(snp.pca$scores[,1:3], group=breed,
      col=as.numeric(sample.group),
      show.ellipses=TRUE, show.group.labels = FALSE,
      axes.color="black")





