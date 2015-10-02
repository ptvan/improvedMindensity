# example workflow of using improvedMindensity() and associated functions
# to improve 1D gating of cyTOF samples
# detailed usage are in comments of improvedMindensity.R

library(openCyto)
library(data.table)
library(latticeExtra)
library(plyr)
library(ggplot2)

dataPath <- "/shared/silo_researcher/Gottardo_R/pvan_working/ragon/CAVD2015Jul/"
path <- "/shared/silo_researcher/Gottardo_R/pvan_working/improvedMindensity/"
setwd(path)

# load the data
gs <- load_gs(paste(dataPath, "output/gs_auto/", sep="/"))
pd <- pData(gs)

# plot the markers
#  png("output/improvedMindensity_live_marker_Sep17.png", width=1000, height=1000)
#  plotGate(gs, "live", type="densityplot", main="live(Rh103Di) marker, improvedMindensity Sep17" )
#  dev.off()
#  
# png("output/improvedMindensity_CD19_marker_Sep17.png", width=1000, height=1000)
# plotGate(gs, "live/CD19neg", type="densityplot", main="CD19(Nd142Di) marker, improvedMindensity Sep17")
# dev.off()

# Rh103Di = live
# Nd142Di = CD19
# Sm154Di = CD3
# Nd145Di = CD4
# Nd146Di = CD8a
# Er167Di = CD27
# Nd144Di = CD38

fs <- getData(gs, "live")
markers <- pData(parameters(fs[[1]]))

chnls <- markers[!is.na(markers$desc),]
chnls <- chnls[chnls$desc != "EQ_Beads",]
chnls <- chnls[!grepl("DNA", chnls$desc),]
chnls <- chnls[!grepl("IL", chnls$desc),]
chnls <- chnls[!grepl("Granzyme", chnls$desc),]
chnls <- chnls[!grepl("PD_1", chnls$desc),]
chnls <- chnls[!grepl("IFN_g", chnls$desc),]
chnls <- chnls[!grepl("TNF_a", chnls$desc),]
chnls <- chnls[!grepl("CD107a", chnls$desc),]
chnls <- chnls[!grepl("Granzyme", chnls$desc),]
chnls <- chnls[!grepl("Perforin", chnls$desc),]
chnls <- as.vector(chnls$name)
names(chnls) <- NULL

#chnls <- c("Rh103Di", "Nd142Di", "Sm154Di")

source("improvedMindensity.R")

gates <- getChannelsPops(gs)
gates <- gates[2]

out <- getSampleStats(gs, gates)
sampleStats <- out[[1]]
featureList <- out[[2]]

chnls <- unlist(gates)
chnl <- chnls

# png("output/improvedMindensity_live_marker_before_Oct01.png", width=1000,height=1000)
plotGate(gs,"live", default.y = "Rh103Di", type="densityplot", main = "live, before")
# dev.off()

png("output/improvedMindensity_CD19_marker_before_Oct01.png", width=1000,height=1000)
plotGate(gs,"live", default.y = "Nd142Di", type="densityplot", main = "live, before")
dev.off()



sampleStats <- flagBadSamples(sampleStats, chnl)
sampleStats[channel == "Rh103Di" & flags != ""]

# Rh103Di = live
# Nd142Di = CD19
# Sm154Di = CD3
# Nd145Di = CD4
# Nd146Di = CD8a
# Er167Di = CD27
# Nd144Di = CD38

# save stats
save(sampleStats, file="output/sampleStats.RData")
save(featureList, file="output/featureList.RData")

# regate problematic samples
regateBadSamples(gs, sampleStats, chnl, plot=F, verbose=T)

png("output/improvedMindensity_live_marker_after_Oct01.png", width=1000,height=1000)
plotGate(gs,"live", default.y = "Rh103Di", type="densityplot", main = "ragonCAVD2015Jul live marker, after regating")
dev.off()

png("output/improvedMindensity_CD19_marker_after_Oct01.png", width=1000,height=1000)
plotGate(gs,"CD19neg", default.y = "Nd142Di", type="densityplot", main = "CD19, after")
dev.off()


