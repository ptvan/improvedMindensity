# improvedMindensity.R
# Phu T. Van & Greg Finak, FredHutch
# September 2015
# 
# this file contains functions for the improved Mindensity 1D gating workflow
# it requires openCyto and data.table 
# 
# To use: 
# 1. call getSampleStats() on a gatingSet and the channel to work on. 
# This will perform 1D gating using improvedMindensity() and generate
# sampleStats and featureList tables
#
# 2. call flagBadSamples() on generated sampleStats. This will flag samples
# with cutpoints outside the threshold (default = 3 Mean Adjusted Deviates) 
# as Far Left or Far Right within their respective channel 
# 
# 3. call regateBadSamples() on gatingSet and sampleStats, which will generate
# new cutpoints for bad samples, set `execute=T` to write new cutpoints to the gatingSet


improvedMindensity <- function(D,adjust=2,gate_range=NA, plot = FALSE, ...){
  # construct the density from data and adjust params we were given
  dens <- density(D,adjust=adjust)
  
  # restrict data to range
  if (!is.null(gate_range)) {
    if (length(gate_range) == 2 & gate_range[1] < gate_range[2]) {
      filter <- dens$x > gate_range[1] & dens$x < gate_range[2]
      dens$x <- dens$x[filter]
      dens$y <- dens$y[filter]
    }else{
      #no range provided, do nothing
    }
  }
  
  sp <- smooth.spline(dens$x,dens$y)
  pred <- predict(sp)
  
  d1 <- predict(sp,deriv = 1)
  d2 <- predict(sp,deriv = 2)
  d3 <- predict(sp,deriv = 3)
  
  # find features
  inf1 <- sign(d2$y[-length(d2$y)])>0&sign(d2$y[-1L])<0
  inf2 <- sign(d2$y[-length(d2$y)])<0&sign(d2$y[-1L])>0
  minima <- sign(d1$y[-1])>0&sign(d1$y[-length(d1$y)])<0
  maxima <- sign(d1$y[-1])<0&sign(d1$y[-length(d1$y)])>0
  shoulders <- sign(d3$y[-1])<0&sign(d3$y[-length(d3$y)])>0
  
  minima_xcoords <- sp$x[which(minima)] # x-coords of minima
  maxima_xcoords <- sp$x[which(maxima)] # x-coords of maxima
  
  minima_ycoords <- sp$y[which(sp$x %in% minima_xcoords)] # y-coords of minima
  maxima_ycoords <- sp$y[which(sp$x %in% maxima_xcoords)] # y-coords of maxima
  
    
  if (length(which(minima == TRUE)) == 0){ # no minima found, look through shoulders
    # if there is a peak, pick first shoulder to the right of peak  
    if (length(which(maxima == TRUE)) == 1  ){ 
      pkidx <- which(maxima == TRUE)
      pt <- sp$x[median(which(shoulders)[which(shoulders) > pkidx])]
      if (is.na(pt)){
        # there is only one peak, and there is no shoulder to the right of it
        # pick an inflection point as a last resort
        pt <- sp$x[median(which(inf2)[which(inf2) > pkidx])]
      }
      
    }  else {
      # no peak (or multiple peaks), select overall median shoulder as cutpoint (for now...)
      pt <- sp$x[median(which(shoulders))]    
    }                                         
  } else if (length(which(minima == TRUE)) > 1) { # multiple minima
    
    m <- min(sp$y[which(sp$x %in% minima_xcoords)])  # pick the minima with lowest y
    pt <- sp$x[which(sp$y == m)]
    
    
  } else if (length(which(minima == TRUE)) == 1){ # only 1 minima, use it as cut point
    pt <- sp$x[which(minima)]
    
  }
  
  .plots = function(){
#     abline(v = d2$x[which(inf1)],col="red")       # red    == inflection point
#     abline(v = d2$x[which(inf2)],col="green")     # green  == inflection point
    abline(v = d2$x[which(maxima)],col="blue")    # blue   == maxima
    abline(v = d2$x[which(minima)],col="orange")  # orange == minima
    abline(v = d2$x[which(shoulders)],col="pink") # pink   == shoulders
  }
  if(plot){
    par(mfrow=c(2,1))
    plot(sp,type="l",main="features")
    .plots()
    plot(sp,type="l",main="final_cut")
    abline(v = pt, col="black", lwd=2)
    par(mfrow=c(1,1)) # be nice, restore plot settings
  }    
 
  return(list(density = sp,
              inf_rising = sp$x[which(inf1)],
              inf_falling = sp$x[which(inf2)],
              maxima = maxima_xcoords,
              minima = minima_xcoords,
              maxima_heights = maxima_ycoords,
              minima_heights = minima_ycoords,
              shoulders = sp$x[which(shoulders)], 
              final_cut = pt))
  
}

getSampleStats <- function(gs, chnl, g_range=c(1,4), adj=2) {
  # NOTE: WILL NEED SOME INPUT CHECKING HERE
  require(openCyto)
  sampleCount <- length(gs) 
  samples <- sampleNames(gs)
#   gates <- getChannelsPops(gs)
  
  # dat <- openCyto:::.truncate_flowframe(dat, channels = "Rh103Di", min=1)
  featureList <- data.table(sample = "",
                            channel = "",
                            feature_type = "",
                            feature_coord =""
  )
  sampleStats <- data.table(sample = rep(samples, length(gates)),
                            channel = unlist(lapply(gates, function(x) rep(x,sampleCount))),
                            pop="",
                            cut_x=0,
                            cut_y=0,
                            num_min=0, 
                            num_max=0, 
                            num_shoulders=0,
                            peaks_xdist=0,
                            flags=" "
  )
  
  keys <- c("sample","channel")
  setkeyv(featureList, keys)
  setkeyv(sampleStats, keys)
  gateName <- names(chnl)
  cat("   ",  chnl, gateName, " \n")

  for (i in 1:sampleCount){
    
    s <- samples[i]
    cat(s, "\n")
    
    # get current gate cutpoints: cut_x is the gate boundary
    g <- getGate(gs[[s]], gateName)
    if (is.infinite(g@min)) cut_x <- g@max
    if (is.infinite(g@max)) cut_x <- g@min
      
    tmp <- getData(gs[[s]], gateName)
    tmp <- openCyto:::.truncate_flowframe(tmp, channels = chnl, min=1)
    tmp <- exprs(tmp[,chnl])
    
    # cut_y has to be calculated from the same density originally used to generate the gate
    d <- density(tmp, adjust=adj)
    cut_y <- d$density$y[which(d$density$x == cut_x)]
    
    x <- improvedMindensity(tmp, plot=F, gate_range=g_range, adjust=adj)
    
    num_min <- length(x$minima)
    num_max <- length(x$maxima)
    num_shoulders <- length(x$shoulders)
    
    # save sample stats in sampleStats
    sampleStats[sample == s & channel == chnl]$pop <- names(which(gates == chnl))
    sampleStats[sample == s & channel == chnl]$cut_y <- cut_y
    sampleStats[sample == s & channel == chnl]$cut_x <- cut_x
    sampleStats[sample == s & channel == chnl]$num_min <- num_min
    sampleStats[sample == s & channel == chnl]$num_max <- num_max
    sampleStats[sample == s & channel == chnl]$num_shoulders <- num_shoulders
    
    if (num_max == 2){
      # there are two peaks, record their X-distance
      sampleStats[sample == s & channel == chnl]$peaks_xdist <- abs(x$maxima[1] - x$maxima[2])
    } 
    if (num_max == 1) { 
      # just 1 peak, set X-distance to distance between peak & cut point
      sampleStats[sample == s & channel == chnl]$peaks_xdist <- abs(x$maxima[1] - x$final_cut)
    }
    
    # save all features into featureList
    if (num_shoulders > 0){
      l <- lapply(as.list(x$shoulders), function(x) c(s,chnl,"shoulder",x))
      l <- as.data.table(t(sapply(l, rbind)))
      colnames(l) <- colnames(featureList)
      featureList <-  rbindlist(list(l, featureList), use.names = F, fill = F)  
    }
    
    if (num_max > 0){
      l <- lapply(as.list(x$maxima), function(x) c(s,chnl,"maxima",x))
      l <- as.data.table(t(sapply(l, rbind)))
      colnames(l) <- colnames(featureList)
      featureList <-  rbindlist(list(l, featureList), use.names = F, fill = F) 
    }
    
    if (num_min > 0){
      l <- lapply(as.list(x$minima), function(x) c(s,chnl,"minima",x))
      l <- as.data.table(t(sapply(l, rbind)))
      colnames(l) <- colnames(featureList)
      featureList <-  rbindlist(list(l, featureList), use.names = F, fill = F)    
    }
    
  
  cat("\n")
  }
  featureList$feature_coord <- as.numeric(featureList$feature_coord)
  return(list(sampleStats, featureList))  
} 
  
flagBadSamples <- function(sampleStats, chnl, mad_thresh = 3){
  # NEEDS INPUT CHECKING HERE !!!
  cat(chnl)
  sampleStats$flags <- " "
  g <- sampleStats[channel == chnl]
  ut <- (mad(g$cut_x) * mad_thresh) + median(g$cut_x)
  lt <- median(g$cut_x) - (mad(g$cut_x) * mad_thresh)  
  far_right <- g$sample[which(g$cut_x > ut)]
  far_left <- g$sample[which(g$cut_x < lt)]
  sampleStats[sample %in% far_left & channel == chnl]$flags <- "FL"
  sampleStats[sample %in% far_right & channel == chnl]$flags <- "FR"
  
  cat(" median =", median(g$cut_x), "[", lt, ",", ut, "] @ MAD =", mad_thresh , "\n")
  cat(length(far_left), "far_left, ", length(far_right), "far_right")
  cat("\n")
  
  return(sampleStats)
}
  
regateBadSamples <- function(gs, sampleStats, chnl, plot=F, negative=F, execute=F){
   cat("gs has", length(gs), "samples...\n")
  # for each channel, get good samples and bad samples
    gate_name <- names(chnl)
#     isParentRoot <- FALSE
#     if (getParent(gs, gate_name)) { isParentRoot <- TRUE}
    
    bad_samples <- subset(sampleStats, flags=="FL"|flags=="FR")
    bad_samples <- subset(bad_samples, channel == chnl)$sample

    cat(gate_name, ":", length(bad_samples), "outlier samples \n")
    
    if (length(bad_samples) > 0){
      good_samples <- subset(sampleStats, flags!="FR"|flags!="FL")
      good_samples <- subset(good_samples, channel == chnl)$sample
      
      cuts_x <- sampleStats[sample %in% good_samples & channel == chnl & peaks_xdist != 0]$cut_x
      cuts_xmedian <- median(cuts_x)
      peaks_xdist_median <- median(sampleStats[sample %in% good_samples & channel == chnl & peaks_xdist != 0]$peaks_xdist)
      for (j in 1:length(bad_samples)){ 
        s <- bad_samples[j]
        p <- getParent(gs, gate_name)
        tmp <- getData(gs[[s]], p)
        
        tmp <- openCyto:::.truncate_flowframe(tmp, channels = chnl, min=1)
        tmp <- exprs(tmp[,chnl])
        x <- improvedMindensity(tmp, plot=F, gate_range=c(1,4))
        # get cut candidates (features detected by improvedMindensity() excluding current cut point)
        # and set feature closest in X-coordinate to the sample's median
        cut_candidates <- c(x$minima, x$shoulders)
        new_cut <- cut_candidates[which.min(abs(cut_candidates - cuts_xmedian))]
        cat("   ", s, ":",x$final_cut," -> ", new_cut)
        
        if (execute){
          cat(" Updating gate....")
          if (negative){
            g_coords <- list(c(-Inf, new_cut))
            cat(" negative gate...")
          } else {
            g_coords <- list(c(new_cut, Inf))
          }
          names(g_coords) <- chnl
          g_filterId <- getGate(gs[[s]], gate_name)@filterId
          g <- rectangleGate(g_coords, filterId = g_filterId)
          
          # remove old gate and add new one
          flowWorkspace::Rm(gate_name, gs[[s]])
          add(gs[[s]], g, parent=getParent(gs, gate_name))
          cat("DONE !")
        }
          
        if (plot){
          par(mfrow=c(3,1))
          plot(x$density, type="l", main=paste(s, chnl, "features"), xlab="", ylab="density")
          #   abline(v = x$inf_rising, col="red")       # red    == inflection point
          #   abline(v = x$inf_falling, col="green")     # green  == inflection point
          abline(v = x$maxima ,col="blue")    # blue   == maxima
          abline(v = x$minima,col="orange")  # orange == minima
          abline(v = x$shoulders,col="pink") # pink   == shoulders
          
          plot(x$density, type="l", main="old_cut", xlab="", ylab="density")
          abline(v = x$final_cut,col="gray")
          
          plot(x$density, type="l", main="new_cut", xlab="", ylab="density")
          rug(cuts_x, col="red")
          abline(v = x$final_cut,col="gray")
          abline(v = new_cut,col="black")
          
          par(mfrow=c(1,1))
          li <- readline()
        }
        cat("\n")
      }
    }  
  
  # update the gatingSet with new cell counts
  if (execute) {
    cat("recomputing event stats...")
    recompute(gs, gate_name)
  }
}  

getChannelsPops <- function(gs){
  # translate channel names to gate names, only works with rectangleGates for now
  nodes <- getNodes(gs[[1]])[-1]
  gates <- lapply(nodes, function(x) getGate(gs[[1]],x) )
  names(gates) <- nodes
  gates <- gates[unlist(lapply(gates, function(x)class(x)) ) == "rectangleGate"]
  gates <- lapply(gates, function(x)names(x@min) )
  return(gates)
}