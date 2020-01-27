# === Process LD tables produced by PLINK1.9

# === Specify path with LD data
# === LD tables are supposed to be in folders
# === according groups (breed and population)
# === For example: LD/group1/pop1
ld <- "/path/to/ld"
plot.path <- "/path/to/save/plots/"
# === Read data to list
ld.tabs <- list()
for (folder in list.dirs(ld))
{
  if (file.exists(paste(folder, "/ld_table.ld", sep="")))
  {
   ld.tabs[[basename(folder)]] <- read.table(paste(folder, "/ld_table.ld", sep=""),
                                             header=TRUE, sep="")
  }
}
# === Calculate LD decay for each table in list
ld.decay <- list()
for (tab in names(ld.tabs))
{
  mean_ld <- numeric(0) 
  for (i in c(1:30)) {
  attach(ld.tabs[[tab]]) 
  R2n <- (R2 - mean(R2)) / sd(R2)
  mean_ld <- c(mean_ld, mean(R2n[abs(BP_A-BP_B)<=i*20000]))
  ld.decay[[tab]] <- mean_ld
  detach(ld.tabs[[tab]])
}}
# === Plot LD decay
png(paste(plot.path, "LD_decay.png"), width=2000, height=2000, res=300)
layout(matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE))
# LD by breeds
plot(ld.decay$LD, type="l", col="black", lwd=3, lty=9,
     main=" LD decay by breeds",
     xlab="Distance between SNP, x20 Kb",
     ylab="Mean R2 by distance inervals")
points(ld.decay$Akzhaiyk, col="darkgreen", pch=16)
lines(ld.decay$Akzhaiyk, col="darkgreen",lwd=3)
points(ld.decay$Edilbay, col="red", pch=16)
lines(ld.decay$Edilbay, col="red",lwd=3)
points(ld.decay$KazTon, col="blue", pch=16)
lines(ld.decay$KazTon, col="blue",lwd=3)
points(ld.decay$Saryarka, col="purple", pch=16)
lines(ld.decay$Saryarka, col="purple",lwd=3)
points(ld.decay$PG, col="orange", pch=16)
lines(ld.decay$PG, col="orange",lwd=3)
legend(x=10, y=0.4, pch=19, bty="n",
       fill="white", border="white",
       legend=c("Akzhayik", "Edilbay", "Kazakh Semi-coarse", "Kazakh Fine-wool", "Saryarka"),
       col=c("darkgreen", "red", "orange", "blue", "purple"),
       y.intersp = 1)
# LD by Akzhayuk populations
plot(ld.decay$LD, type="l", col="black", lwd=3, lty=9,
     main=" LD decay in Akzhayik breed",
     xlab="Distance between SNP, x20 Kb",
     ylab="Mean R2 by distance inervals")
points(ld.decay$Atameken, col="darkgreen", pch=1)
lines(ld.decay$Atameken, col="darkgreen",lwd=3)
points(ld.decay$Kuanysh, col="darkgreen", pch=2)
lines(ld.decay$Kuanysh, col="darkgreen",lwd=3)
points(ld.decay$Saltanat, col="darkgreen", pch=3)
lines(ld.decay$Saltanat, col="darkgreen",lwd=3)
legend(x=10, y=0.4, bty="n",
       fill="white", border="white",
       legend=c("Atameken", "Kuanysh", "Saltanat"),
       pch=c(1,2,3),
       y.intersp = 1)
# LD by Edilbay populations
plot(ld.decay$LD, type="l", col="black", lwd=3, lty=9,
     main=" LD decay in Edilbay breed",
     xlab="Distance between SNP, x20 Kb",
     ylab="Mean R2 by distance inervals")
points(ld.decay$Azhar, col="red", pch=1)
lines(ld.decay$Azhar, col="red",lwd=3)
points(ld.decay$Birlik, col="red", pch=2)
lines(ld.decay$Birlik, col="red",lwd=3)
legend(x=10, y=0.4, bty="n",
       fill="white", border="white",
       legend=c("Azhar","Birlik"),
       pch=c(1,2),
       y.intersp = 1)
# LD by Kazakh Semi-coarse popukations
plot(ld.decay$LD, type="l", col="black", lwd=3, lty=9,
     main=" LD decay in Kazakh Semi-coarse breed",
     xlab="Distance between SNP, x20 Kb",
     ylab="Mean R2 by distance inervals")
points(ld.decay$AltynAsel, col="orange", pch=1)
lines(ld.decay$AltynAsel, col="orange",lwd=3)
points(ld.decay$KaraAdyr, col="orange", pch=2)
lines(ld.decay$KaraAdyr, col="orange",lwd=3)
points(ld.decay$Khasiev, col="orange", pch=3)
lines(ld.decay$Khasiev, col="orange",lwd=3)
points(ld.decay$Otkanzhar, col="orange", pch=4)
lines(ld.decay$Otkanzhar, col="orange",lwd=3)
legend(x=10, y=0.4, bty="n",
       fill="white", border="white",
       legend=c("Altyn Asel", "Kara Adyr", "Khasiev", "Otkanzhar"),
       pch=c(1,2,3,4),
       y.intersp = 1)
dev.off()







