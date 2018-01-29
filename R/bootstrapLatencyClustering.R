library(cluster)
library(mclust)


# Function to compute probability of LL from an Mclust object. Requires two Gaussian clusters
# Probabilities returned by calling emProbLL(clusters, clusters$data) should be very similar to
# LL class prob from clusters$z
emProbLL <- function(clusters, latency) {

    if (clusters$parameters$variance$G != 2) {
        stop("This function requires exactly two clusters")
    }

    sigmaSq <- clusters$parameters$variance$sigmasq
    
    if (length(sigmaSq) == 1) {
        # Equal variance, copy into array
        sigmaSq <- c(sigmaSq, sigmaSq)
    }
    
    LL_index <- 2
    
    if ( clusters$parameters$mean[2] < clusters$parameters$mean[1]) {
        LL_index <- 1
    }
    
    SL_index <- 3 - LL_index
    
    l <- clusters$parameters$pro[LL_index] * dnorm(latency, mean = clusters$parameters$mean[LL_index], sd = sqrt(sigmaSq[LL_index]))
    
    s <- clusters$parameters$pro[SL_index] * dnorm(latency, mean = clusters$parameters$mean[SL_index], sd = sqrt(sigmaSq[SL_index]))
    
    r <- l / (l + s)
    
    return(r)
    
}



bootstrapEM_Class <- function(latency, its=1000, stratification="initial", modelNames = "E") {
    
    # The prior groups are used to stratify the bootstrap, which maintains a consistent fraction of LL and SL
    # in each bootstrap.
    #
    # Optionally, we can resample the groups by using the prior probabilities, which are derived from the EM
    # classification of the original data.
    #
    # This allows us to vary the stratification according to uncertainty in the prior classification.
    #
    # For example, consider the subjects (A,B,C,D,E) with prior probabilities (0, 0.05, 0.45, 0.95, 1).
    # This contains 3 SL and 2 LL. But if we resample the groups using these probabilities, we will classify A
    # as LL with probability 0, B with probability 0.05, C with probability 0.45, and so on.
    # 

    numSubj <- length(latency)

    # Assume two clusters, equal variance
    
    clusters <- Mclust(latency, G = 2, modelNames)

    priorProbLL <- emProbLL(clusters, latency)
    
    # Record mix, mean, and std dev of Gaussians of all bootstraps
    em_mix <- matrix(nrow = 2, ncol = its)
    em_mean <- matrix(nrow = 2, ncol = its)
    em_std <- matrix(nrow = 1, ncol = its)
    
    if (length(modelNames) > 1 || modelNames == "V") {
        em_std <- matrix(nrow = 2, ncol = its)
    }

    # Probability of LL classification for each subject at each iteration
    r_boot <- matrix(nrow = numSubj, ncol = its)

    # Data for plotting bootstrap confidence bounds of prob(LL) curve
    r_curve_boot <- matrix(nrow = 900, ncol = its)


    priorLL <- which(priorProbLL > 0.5)
    priorSL <- which(priorProbLL <= 0.5)


    stratifyInitial <- TRUE
    stratifyProbabilistic <- FALSE

    if (tolower(stratification) == "none") {
        stratifyInitial <- FALSE
        stratifyProbabilistic <- FALSE
    }
    else if (tolower(stratification) == "initial") {
        stratifyInitial <- TRUE
        stratifyProbabilistic <- FALSE
    }
    else if (tolower(stratification) == "prob") {
        stratifyInitial <- FALSE
        stratifyProbabilistic <- TRUE
    }
    else {
        stop(paste("Unrecognized stratification option ", stratification, sep = ""))
    }
    
    
    for (i in 1:its) {


        # boot is a collection of indices to sample, latencies are latencies[boot]      
        boot <- c()

        if (stratifyInitial || stratifyProbabilistic) {
            bootPriorLL <- priorLL
            bootPriorSL <- priorSL
            
            if (stratifyProbabilistic) {
                isLL <- vector("numeric", numSubj)
                
                for (n in 1:numSubj) {
                    isLL[n] <- rbinom(1,1, priorProbLL[n])
                }
                
                isLL <- as.logical(isLL)
                
                bootPriorLL <- which(isLL)
                bootPriorSL <- which(!isLL)
                
            }
            
            bootLL <- sample(bootPriorLL, replace = TRUE)
            bootSL <- sample(bootPriorSL, replace = TRUE)
            
            boot <- c(bootLL, bootSL)
        }
        else {
            # Sample the full array of input latencies without regard to initial classification
            boot <- sample(1:numSubj, replace = TRUE)
        }
            
        clusters_boot <- Mclust(latency[boot], G = 2, modelNames = modelNames)

        LL_index <- 2

        # Make sure clusters are correctly ordered
        if ( clusters_boot$parameters$mean[2] < clusters_boot$parameters$mean[1]) {
            LL_index <- 1
        }

        SL_index <- 3 - LL_index

        r <- emProbLL(clusters_boot, latency)

        r_boot[,i] <- r

        em_mix[,i] <- clusters_boot$parameters$pro[c(SL_index, LL_index)]
        em_mean[,i] <- clusters_boot$parameters$mean[c(SL_index, LL_index)]

        if (length(clusters_boot$parameters$variance$sigmasq) > 1) {
            em_std[,i] <- sqrt(clusters_boot$parameters$variance$sigmasq[c(SL_index, LL_index)])
        }
        else {
            em_std[,i] <- sqrt(clusters_boot$parameters$variance$sigmasq)
        }
        

        r_curve_boot[,i] <-  emProbLL(clusters_boot, c(1:900))

    }


    # Now return some interesting things
    bootProbLL <- rowSums(r_boot > 0.5) / its


    # Compute approximate LL classification threshold
    #
    # Lazily do an O(its / 2) search here.
    #
    # Could do a binary search or even solve analytically

    LL_thresh = vector("numeric", its)

    for (i in 1:its) {
        n <- 1
        while(r_curve_boot[n, i] < 0.5) {
            if (r_curve_boot[n+1,i] > 0.5) {
                LL_thresh[i] = n + 0.5 
            }
            n <- n + 1
        }
    }


    stuff <- list(bootProbLL = bootProbLL, bootThreshLL = LL_thresh, clusters = clusters, em_mean = em_mean, em_mix = em_mix, em_std = em_std, its = its, latency = latency, priorProbLL = priorProbLL, r_boot = r_boot, r_curve_boot = r_curve_boot)

    return(stuff)
    
}


plotBootstrapEM_Class <- function(bootClass, bootCurveAlpha = 0.05, subjectID = NULL) {

    # Plotting stuff
    
    # Plot one after the other with dev.new(), could also do lm style with
    # devAskNewPage(ask = TRUE) but this is annoying because the plot window
    # steals focus
    
    its <- bootClass$its

    numSubj <- length(bootClass$latency)

    if (is.null(subjectID)) {
        subjectID <- as.character(c(1:numSubj))
    }
    else {
        if (length(subjectID) != numSubj) {
            stop("Length of subjectID vector does not match number of entries in latency data")
        }
    }
    
    
    percentileLower <- round(its * bootCurveAlpha / 2)
    percentileUpper <- round((1.0 - bootCurveAlpha / 2) * its)

    r_curve_boot_lower <- vector("numeric", length = 900)
    r_curve_boot_upper <- vector("numeric", length = 900)

    r_curve_boot = bootClass$r_curve_boot
    
    for (n in 1:900) {
        r_curve_boot_lower[n] <- r_curve_boot[n,order(r_curve_boot[n,])[percentileLower]]
        r_curve_boot_upper[n] <- r_curve_boot[n,order(r_curve_boot[n,])[percentileUpper]]
    }

    dev.new(width=10, height=7)
    
    par(mar = c(5,5,2,3))

    plot(x = bootClass$latency, y = bootClass$priorProbLL, cex = 2, cex.axis = 1.5, cex.lab = 1.5, xlim = c(1, 900), ylim = c(0.0, 1.0), xlab = "Average latency (s)", ylab = "Probability of stress resilience", lwd = 2, pch = 1, main = "Defeat latency and probability of classification as stress resilient")

    ## # Annotate subjects where classification is not clear, which we'll say is between 0.4 and 0.6
    ## labelSubj = (priorProbLL > 0.4) & (priorProbLL < 0.6)

    ## text(bootClass$latency[labelSubj], bootClass$priorProbLL[labelSubj], subjectID[labelSubj], pos = 2)
    
    # Probability curve computed with all data
    r_curve <- emProbLL(bootClass$clusters, c(1:900))
    
    points(1:900, y = r_curve, lwd = 2, lty = 1, type = "l")

    points(c(1:900), r_curve_boot_lower, lty = 3, lwd = 2, type = "l")
    points(c(1:900), r_curve_boot_upper, lty = 3, lwd = 2, type = "l")

    abline(h = 0, lty = 2, lwd = 1)
    abline(h = 0.25, lty = 2, lwd = 1)
    abline(h = 0.5, lty = 2, lwd = 1)
    abline(h = 0.75, lty = 2, lwd = 1)
    abline(h = 1, lty = 2, lwd = 1)
    
    legend(0,0.99, c("Subject","EM Model", "Bootstrap CI"),  pch = c(1, NA, NA), lty = c(NA, 1, 3), lwd = 2, cex = 1 )

    
    dev.new()
    par(mar = c(5,5,2,3))
    
    # Plot histogram of LL threshold
    hist(bootClass$bootThreshLL, 25, xlim = c(min(bootClass$bootThreshLL) - 10, max(bootClass$bootThreshLL) + 10), xlab = "Threshold average latency (s)", ylab = paste("Frequency (", its, " iterations)", sep = ""), lwd = 2, cex = 1.25, cex.axis = 1.5, cex.lab = 1.5, main = "Stress resilience classification threshold over all bootstraps")

    
    # Count number of times a subject is classified as LL, compare to prior EM probs
    
    dev.new(width=14, height=7)
    par(mar = c(5,7,2,2), las = 2)
    

    diff <- bootClass$priorProbLL - bootClass$bootProbLL

    barCols <- matrix(nrow=2, ncol=numSubj)
    
    barCols[1,] = bootClass$priorProbLL
    barCols[2,] = bootClass$bootProbLL

    ord = order(bootClass$bootProbLL, bootClass$priorProbLL)
    
    pos = barplot(barCols[,ord], names.arg = subjectID[ord], args.legend = list(x = 20, y = 0.95), legend.text = c("EM", "Boot"), main = "Stress-resilience classification probabilities", ylab = "Probability of stress resilience", lwd = 2, cex.main = 1.5, cex.axis = 1.25, cex.lab = 1.25, beside = TRUE, ylim = c(-0.01, 1.01))

    abline(h = 0, lty = 2, lwd = 1)
    abline(h = 0.25, lty = 2, lwd = 1)
    abline(h = 0.5, lty = 2, lwd = 1)
    abline(h = 0.75, lty = 2, lwd = 1)
    abline(h = 1, lty = 2, lwd = 1)
    
}




bootstrapClusterClass <- function(latency, its=1000, stratification="initial", algorithm=c("kmeans", "hc", "pam"), emModelNames = "E") {

    # Does clustering with either the original kmeans (stats), hierarchical clustering (mclust) or pam (cluster)
    
    # The prior groups are derived from the original data, or optionally derived from EM and used to probabilistically stratify the bootstrap
    #
    
    algorithm <- algorithm[1]

    numSubj <- length(latency)

    priorProbLL <- NULL
    
    # Centers (or medoids) and the clusters from clustering the original data
    centers <- NULL
    clusters <- NULL
    boundary <- NULL

    if (algorithm == "kmeans") {
        
        kmeans <- kmeans(latency, centers = 2, iter.max = 20, nstart = 10)
        
        centers <- kmeans$centers
        
        clusters <- kmeans$cluster
        
    }
    else if (algorithm == "pam") {
        pam <- pam(latency, 2)
        
        centers <- pam$medoids
        
        clusters <- pam$clustering
        
    }
    else if (algorithm == "hc") {
        hc <- hc("E", latency)
        clusters <- hclass(hc, 2)
        
        centers = c( mean(latency[clusters == 1]), mean(latency[clusters == 2]) )
        
    }
    else {
        stop(paste("Unrecognized algorithm ", algorithm[1], sep = ""))
    }
    
    
    LL_index <- 2
    
    # Make sure clusters are correctly ordered
    if ( centers[2] < centers[1]) {
        LL_index <- 1
    }
    
    SL_index <- 3 - LL_index
    
    # Re-order clusters if needed
    
    centers <- centers[c(SL_index, LL_index)]
    
    LL <- which(clusters == LL_index)
    SL <- which(clusters == SL_index)
    
    clusters[LL] <- 2
    clusters[SL] <- 1
    
    minLL <- min(latency[LL])
    maxSL <- max(latency[SL])

    boundary <- (centers[1] + centers[2]) / 2

    stratifyInitial <- TRUE
    stratifyProbabilistic <- FALSE

    if (tolower(stratification) == "none") {
        stratifyInitial <- FALSE
        stratifyProbabilistic <- FALSE
    }
    else if (tolower(stratification) == "initial") {
        stratifyInitial <- TRUE
        stratifyProbabilistic <- FALSE
    }
    else if (tolower(stratification) == "prob") {
        stratifyInitial <- FALSE
        stratifyProbabilistic <- TRUE
    }
    else {
        stop(paste("Unrecognized stratification option ", stratification, sep = ""))
    }
    
    if (stratifyProbabilistic) {
        clusters <- Mclust(latency, G = 2, modelNames = emModelNames)
        priorProbLL <- emProbLL(clusters, latency)
    }
    else {
        priorProbLL[LL] <- 1
        priorProbLL[SL] <- 0
    }
    
    # Record classification for samples used in each bootstrap
    class_boot <- matrix(nrow = numSubj, ncol = its)

    # Decision boundary for each bootstrap
    boundary_boot <- matrix(nrow = 1, ncol = its)
    
    priorLL <- which(priorProbLL > 0.5)
    priorSL <- which(priorProbLL <= 0.5)

    for (i in 1:its) {

        # boot is a collection of indices to sample, latencies are latencies[boot]      
        boot <- c()

        if (stratifyInitial || stratifyProbabilistic) {
            
            bootPriorLL <- priorLL
            bootPriorSL <- priorSL
            
            if (stratifyProbabilistic) {
                isLL <- vector("numeric", numSubj)
                
                for (n in 1:numSubj) {
                isLL[n] <- rbinom(1,1, priorProbLL[n])
                }
                
                isLL <- as.logical(isLL)
                
                bootPriorLL <- which(isLL)
                bootPriorSL <- which(!isLL)
            }
        
            bootLL <- sample(bootPriorLL, replace = TRUE)
            bootSL <- sample(bootPriorSL, replace = TRUE)
            
            boot <- c(bootLL, bootSL)
            
        }
        else {
            # Sample the full array of input latencies without regard to initial classification
            boot <- sample(1:numSubj, replace = TRUE)
        }
            
        centersBoot <- NULL

        clustersBoot <- NULL
        
        if (algorithm == "kmeans") {
            kmeansBoot <- kmeans(latency[boot], centers = 2, iter.max = 20, nstart = 10)

            centersBoot <- kmeansBoot$centers

            clustersBoot <- kmeansBoot$cluster
            
        }
        else if (algorithm == "pam") {
            pamBoot <- pam(latency[boot], 2)

            centersBoot <- pamBoot$medoids

            clustersBoot <- pamBoot$clustering
            
        }
        else if (algorithm == "hc") {

            latBoot = latency[boot]
            
            hcBoot <- hc("E", latBoot)
            clustersBoot <- hclass(hcBoot, 2)
            
            centersBoot = c( mean(latBoot[clustersBoot == 1]), mean(latBoot[clustersBoot == 2]) )
            
        }
        else {
            stop(paste("Unrecognized algorithm ", algorithm[1], sep = ""))
        } 
        
        
        LL_index <- 2
        
        # Make sure clusters are correctly ordered, as classifier might put either one first
        if ( centersBoot[2] < centersBoot[1]) {
            LL_index <- 1
        }
        
        SL_index <- 3 - LL_index
        
        class_boot[,i] = vector("numeric", numSubj)
        
        bootLL_Subj <- boot[which(clustersBoot == LL_index)]
        bootSL_Subj <- boot[which(clustersBoot == SL_index)]

        # Boundary is halfway between cluster means   
        boundary_boot[i] <- (centersBoot[LL_index] + centersBoot[SL_index]) / 2

        # Here we enforce SL = 1, LL = 2
        class_boot[,i] = as.numeric(latency > boundary_boot[i]) + 1
        
        
    }

    SL_Count <- rowSums(class_boot == 1)
    LL_Count <- rowSums(class_boot == 2)

    bootProbLL <- LL_Count / (LL_Count + SL_Count)

    # Return interesting stuff
    
    stuff <- list(bootProbLL = bootProbLL, its = its, latency = latency, boundary = boundary, boundary_boot = boundary_boot, centers = centers, class_boot = class_boot, clusters = clusters)

    return(stuff)
    
}



plotBootstrapClusterClass <- function(bootClass, bootBoundaryAlpha = 0.05, subjectID = NULL) {
    
    # Plotting stuff
    
    # Plot one after the other with dev.new(), could also do lm style with
    # devAskNewPage(ask = TRUE) but this is annoying because the plot window
    # steals focus
    
    its <- bootClass$its

    numSubj <- length(bootClass$latency)

    if (is.null(subjectID)) {
        subjectID <- as.character(c(1:numSubj))
    }
    else {
        if (length(subjectID) != numSubj) {
            stop("Length of subjectID vector does not match number of entries in latency data")
        }
    }
    
    
    percentileLower <- round(its * bootBoundaryAlpha / 2)
    percentileUpper <- round((1.0 - bootBoundaryAlpha / 2) * its)

    boundaryConf <- bootClass$boundary_boot[ order(bootClass$boundary_boot)[c(percentileLower, percentileUpper)] ]

    dev.new(width=10, height=7)
    
    par(mar = c(5,5,2,3))

    plot(x = bootClass$latency, y = bootClass$bootProbLL, cex = 2, cex.axis = 1.5, cex.lab = 1.5, xlim = c(1, 900), ylim = c(0.0, 1.0), xlab = "Average latency (s)", ylab = "Bootstrap probability of stress resilience", lwd = 2, pch = bootClass$clusters, main = "Initial clustering and bootstrap probability of stress resilience")

    abline(h = 0, lty = 2, lwd = 1)
    abline(h = 0.25, lty = 2, lwd = 1)
    abline(h = 0.5, lty = 2, lwd = 1)
    abline(h = 0.75, lty = 2, lwd = 1)
    abline(h = 1, lty = 2, lwd = 1)
    
    abline(v = bootClass$boundary, lty = 1, lwd = 2)
    abline(v = boundaryConf[1], lty = 3, lwd = 2)
    abline(v = boundaryConf[2], lty = 3, lwd = 2)

    legend(0,0.99, c("SL subject", "LL subject", "Cluster boundary", "Bootstrap CI"),  pch = c(1, 2, NA, NA), lty = c(NA, NA, 1, 3), lwd = 2, cex = 1 )
    
    # Plot histogram of LL classification threshold
    
    dev.new()
    par(mar = c(5,5,2,3))
    
    # Plot results
    hist(bootClass$boundary_boot, 25, xlim = c(min(bootClass$boundary_boot) - 10, max(bootClass$boundary_boot) + 10), xlab = "Threshold average latency (s)", ylab = paste("Frequency (", its, " iterations)", sep = ""), lwd = 2, cex = 1.25, cex.axis = 1.5, cex.lab = 1.5, main = "Stress resilient cluster boundary over all bootstraps")

    
    
    dev.new(width=14, height=7)
    par(mar = c(5,7,2,2), las = 2)
    
    ord <- order(bootClass$bootProbLL)
    
    pos = barplot(bootClass$bootProbLL[ord], names.arg = subjectID[ord], args.legend = list(title = "Original cluster", x = 10, y = 0.95, fill = c("gray20", "gray65")), legend.text = c("SL subject", "LL subject"), main = "Bootstrap stress-resilience probabilities", ylab = "Probability of stress resilience", lwd = 2, cex.main = 1.5, cex.axis = 1.25, cex.lab = 1.25, col = c("gray20", "gray65")[bootClass$clusters[ord]], ylim = c(-0.01, 1.01))

    abline(h = 0, lty = 2, lwd = 1)
    abline(h = 0.25, lty = 2, lwd = 1)
    abline(h = 0.5, lty = 2, lwd = 1)
    abline(h = 0.75, lty = 2, lwd = 1)
    abline(h = 1, lty = 2, lwd = 1)
    
}
