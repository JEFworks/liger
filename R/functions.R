#' Gene set enrichment analysis
#'
#' @param values vector of values with associated gene names; values must be named, according to names appearing in set elements
#' @param geneset vector of genes in the gene set
#' @param power an exponent to control the weight of the step (default: 1)
#' @param rank whether to use ranks as opposed to values (default: FALSE)
#' @param weight additional weights associated with each value (default: rep(1,length(values)))
#' @param n.rand number of random permutations used to assess significance (default: 1e4)
#' @param plot whether to plot (default: TRUE)
#' @param main plot title (default: "")
#' @param return.details whether to return extended details (default: FALSE)
#' @param quantile.threshold threshold used (default: min(100/n.rand,0.1))
#' @param random.seed random seed (default: 1)
#' @param mc.cores number of cores for parallel processing (default: 1)
#'
#' @examples
#' data("org.Hs.GO2Symbol.list")  
#' universe <- unique(unlist(org.Hs.GO2Symbol.list))  # get universe
#' gs <- org.Hs.GO2Symbol.list[[1]]  # get a gene set
#' # fake dummy example where everything in gene set is perfectly enriched
#' vals <- rnorm(length(universe), 0, 10)  
#' names(vals) <- universe
#' vals[gs] <- rnorm(length(gs), 100, 10)
#' # test obviously enriched set, reduce n.rand for speed
#' gsea(values=vals, geneset=gs, mc.cores=1, n.rand=100, main="GO:Random") 
#' 
#' @export
#' 
gsea <- function(values, geneset, power=1, rank=FALSE, weight=rep(1,length(values)), n.rand=1e4, plot=TRUE, main="", return.details=FALSE, quantile.threshold=min(100/n.rand,0.1), random.seed=1, mc.cores=1) {

    # Binary vector indicating presence in the gene set
    set <- names(values) %in% geneset

    # Order the vectors by the values
    vo <- order(values, decreasing = TRUE)
    values <- values[vo]
    set <- set[vo]
    weight <- weight[vo]
    
    # Calculate enrichment score
    if (rank) {
      # Use the gene weights
      es <- weight
    } else {
      # Use the values raised to a power and weight them
      es <- (abs(values) ^ power) * weight
    }
    
    eso <- es
    
    # Evaluate the fraction of hits and misses
    es[set] <- es[set] / sum(es[set])  # Phit
    es[!set] <- -1 * es[!set] / sum(es[!set])  # -Pmiss
    es <- es * length(values)
    sv <- cumsum(es) 
    # From original paper: "The ES is the maximum deviation from zero of Phit - Pmiss. For
    # a randomly distributed S, ES(S) will be relatively small, but if it is
    # concentrated at the top or bottom of the list, or otherwise nonrandomly
    # distributed, then ES(S) will be correspondingly high."
    
    sv.p <- max(0, max(sv))
    sv.n <- min(0, min(sv))
    p.x.p <- which.max(sv)
    p.x.n <- which.min(sv)

    # Randomization
    # From original paper: "We assess the significance of an observed
    # ES by comparing it with the set of scores ESNULL computed with
    # randomly assigned phenotypes."
    # Though here, we randomly permute set labels as opposed to recomputing full phenotype recalculation
    if (mc.cores > 1) {
        rvll <- parallel::mclapply(1:mc.cores, function(i) {
            set.seed(random.seed + i)
            gseaRandCore(set, eso, nsamples = ceiling(n.rand / mc.cores))
        }, mc.preschedule=TRUE, mc.cores = mc.cores)
        rvl <- list(
          p = unlist(lapply(rvll, function(x) x$p)),
          n = unlist(lapply(rvll, function(x) x$n))
        )
    } else {
        set.seed(random.seed)
        rvl <- gseaRandCore(set, eso, nsamples = n.rand)
    }
    # From original paper: "Estimate nominal P value for S from ESNULL by using the
    # positive or negative portion of the distribution corresponding to
    # the sign of the observed ES(S)."
    p.v.p <- (sum(rvl$p >= sv.p)+1)/(length(rvl$p+1))
    p.v.n <- (sum(rvl$n <= sv.n)+1)/(length(rvl$n+1))
    p.x <- ifelse(p.v.p < p.v.n, p.x.p, p.x.n)  # index of edge value
    p.v <- ifelse(p.v.p < p.v.n, p.v.p, p.v.n)  # P-value

    if(p.v.p < p.v.n) {
        q.thr <- stats::quantile(rvl$p, p = (1-quantile.threshold))
    } else {
        q.thr <- stats::quantile(rvl$n, p = quantile.threshold)
    }

    if(plot) {
        l <- graphics::layout(matrix(c(1,2,3),3, 1, byrow=T), c(1,1,1), c(1.4,0.3,0.9), FALSE)
        graphics::par(mar = c(0.1,3.5, 1.5, 0.5), mgp = c(2,0.65,0), cex = 0.9)
        # Plot scores along gene list rank
        plot( c(1:length(sv)), sv, type = 'l', 
              ylab = "score", xaxt="n", xlab="", main=main, xaxs="i", 
              xlim=c(1,length(sv)), col="darkblue")
        # Plot the maximum deviation from zero ie. the enrichment score
        graphics::segments(p.x, 0, p.x, sv[p.x], col=2, lty=2) 
        graphics::segments(0, sv[p.x], p.x, sv[p.x], col=2, lty=2)

        # Plot quantile thresholds
        qv.p <- stats::quantile(rvl$p,p=1-10^(-1*seq(1,round(log10(length(rvl$p))))))
        qv.n <- -1*stats::quantile(-1*rvl$n,p=1-10^(-1*seq(1,round(log10(length(rvl$n))))))
        graphics::abline(h=c(qv.p,qv.n),lty=3,col=8)
        graphics::abline(h=0,lty=2,col=8)
        
        # Plot P-value
        xpos <- "right"; if(p.x>round(length(sv)/2)) {  xpos <- "left"; }
        graphics::legend( x = ifelse( sv[round(length(sv)/2)]>sv[p.x], paste("bottom",xpos,sep=""), paste("top",xpos,sep="")), 
                bty="n",
                legend = paste("P-value <", format(p.v, digits=3))
                )
        
        # Plot gene set
        graphics::par(mar = c(0.1,3.5,0.1,0.5))
        mset <- set
        mset[!set] <- NA
        plot( mset, xaxt="n", ylim=c(0,1), xlab="", ylab="", xaxs="i", yaxt="n", type='h', col="blue")
        graphics::abline( v = p.x, col=2, lty=2)
        graphics::box()
        
        # Plot values 
        graphics::par(mar = c(0.5,3.5,0.1,0.5))
        plot( values, ylab="values", xaxs="i", type='h', xaxt="n", col="darkblue")
        graphics::legend( x = "topright",
                bty="n", 
                legend=paste("edge value = ", format(values[p.x], digits=2))
                )
        graphics::abline( v = p.x, col=2, lty=2)
        graphics::box()
    }
    
    if(return.details) {
        rv <- c(p.val = p.v, edge.score = as.numeric(sv[p.x]), edge.value = as.numeric(values[p.x]), scaled.score = as.numeric(sv[p.x]/q.thr))
        #if(pareto.estimate) { rv <- c(rv, pareto.pval = pareto.tail.estimate(abs(rvl), as.numeric(abs(sv[p.x])))) }
        return(rv)
    } else {
        return(p.v)
    }
}


#' Bulk gene set enrichment analysis
#'
#' @param values vector of values with associated gene names; values must be named, according to names appearing in set.list elements
#' @param set.list list of gene sets
#' @param power an exponent to control the weight of the step (default: 1)
#' @param rank whether to use ranks as opposed to values (default: FALSE)
#' @param weight additional weights associated with each value (default: rep(1,length(values)))
#' @param n.rand number of random permutations used to assess significance (default: 1e4)
#' @param return.details whether to return extended details (default: FALSE)
#' @param quantile.threshold threshold used (default: min(100/n.rand,0.1))
#' @param mc.cores number of cores for parallel processing (default: 1)
#' @param skip.qval.estimation whether to skip q-value estimation for multiple testing (default: FALSE)
#'
#' @examples
#' data("org.Hs.GO2Symbol.list")  
#' universe <- unique(unlist(org.Hs.GO2Symbol.list))  # get universe
#' gs <- org.Hs.GO2Symbol.list[[1]]  # get a gene set
#' vals <- rnorm(length(universe), 0, 10)  # simulate values
#' names(vals) <- universe
#' vals[gs] <- rnorm(length(gs), 100, 10)  
#' gs.list <- org.Hs.GO2Symbol.list # get gene sets
#' # reduce n.rand for speed
#' bulk.gsea(values = vals, set.list = gs.list[1:3], mc.cores = 1, n.rand=100)
#' 
#' @export
#' 
bulk.gsea <- function(values, set.list, power=1, rank=FALSE, weight=rep(1,length(values)), n.rand=1e4, mc.cores=1, quantile.threshold=min(100/n.rand,0.1), return.details=FALSE, skip.qval.estimation=FALSE) {

    # Determine set matrix
    setm <- do.call(rbind, parallel::mclapply(set.list, function(set) names(values) %in% set, mc.cores=mc.cores))
    # Only bother testing if more than 2 genes present
    setm <- setm[rowSums(setm) > 2, , drop = FALSE]
    
    # Order the vectors by the values
    vo <- order(values, decreasing=TRUE)
    values <- values[vo]
    setm <- setm[, vo, drop=FALSE]
    weight <- weight[vo]
    if(rank) { 
      es <- weight 
    } else { 
      es <- (abs(values)^power) * weight
    }
    
    eso <- es

    ies <- t(t(setm)*es)
    ies <- ies/rowSums(ies)
    ies[is.nan(ies)] <- 0  # when values null out membership
    pes <- t(t(!setm)*es)
    pes <- pes/rowSums(pes)
    pes[is.nan(pes)] <- 0
    esm <- ies-pes
    svm <- t(apply(esm,1,cumsum))*length(values)
    
    svm.maxp <- apply(svm,1,function(x) x[which.max(x)])
    svm.maxp[svm.maxp<0] <- 0
    svm.maxn <- apply(svm,1,function(x) x[which.min(x)])
    svm.maxn[svm.maxn>0] <- 0
    p.xmp <- apply(svm, 1, which.max)
    p.xmn <- apply(svm, 1, which.min)
    
    # Randomization
    if(mc.cores>1) {
        rvlp <- parallel::mclapply(1:mc.cores,function(i) {
            set.seed(i)
            gseaBulkCore(setm,eso,ceiling(n.rand/mc.cores))
        }, mc.preschedule=TRUE, mc.cores=mc.cores)
        rvl <- list(p=do.call(cbind,lapply(rvlp,function(x) x$p)),
                    n=do.call(cbind,lapply(rvlp,function(x) x$n)));
        rm(rvlp); gc();
    } else {
        set.seed(1)
        rvl <- gseaBulkCore(setm,eso,n.rand)
    }

    # Raw p-values
    p.v.p <- (rowSums(rvl$p - svm.maxp >= 0) +1)/(ncol(rvl$p)+1)
    p.v.n <- (rowSums(rvl$n - svm.maxn <= 0) +1)/(ncol(rvl$n)+1)
    # Upper quantiles for scaling thresholds
    q.thr.p <- matrixStats::rowQuantiles(rvl$p,probs=(1-quantile.threshold),drop=T)
    q.thr.n <- matrixStats::rowQuantiles(rvl$n,probs=quantile.threshold,drop=T)

    p.val <- ifelse(p.v.p<p.v.n,p.v.p,p.v.n)
    sscore <- ifelse(p.v.p<p.v.n,svm.maxp/q.thr.p,-1*svm.maxn/q.thr.n)
    x.val <- values[ifelse(p.v.p<p.v.n,p.xmp,p.xmn)]

    if(!skip.qval.estimation) {
        # Sign-aware mean-scaling
        s.pm <- apply(rvl$p,1,mean)
        s.nm <- apply(rvl$n,1,mean)
        s.pm[is.na(s.pm)] <- 0; s.nm[is.na(s.nm)] <- 0;

        s.pss <- stats::ecdf(as.numeric(rvl$p/s.pm))
        s.nss <- stats::ecdf(as.numeric(rvl$n/s.nm))

        # Scale the actual Smax
        s.svm.maxp <- svm.maxp/s.pm;
        s.svm.maxn <- svm.maxn/s.nm;

        vo <- order(s.svm.maxp,decreasing=F)
        vr <- rank(s.svm.maxp)
        qv <- (1-s.pss(s.svm.maxp-1e-10*max(s.svm.maxp)))/(1-stats::ecdf(s.svm.maxp)(s.svm.maxp-1e-10))
        q.v.p <- cummin(qv[vo])[vr];

        vo <- order(s.svm.maxn,decreasing=F)
        vr <- rank(s.svm.maxn)
        qv <- (1-s.pss(s.svm.maxn-1e-10*max(s.svm.maxn)))/(1-stats::ecdf(s.svm.maxn)(s.svm.maxn-1e-10))
        q.v.n <- cummin(qv[vo])[vr];

        q.val <- ifelse(q.v.p<q.v.n,q.v.p,q.v.n)
    } else {
        q.val <- stats::p.adjust(p.val)
    }

    # P-value, Q-value table
    df <- data.frame(p.val=p.val,q.val=q.val,sscore=sscore,edge=x.val)
    rownames(df) <- rownames(setm)

    if(return.details) {
        ddf <- data.frame(svm.maxp=svm.maxp,svm.maxn=svm.maxn,p.xmp=p.xmp,p.xmn=p.xmn,p.v.p=p.v.p,p.v.n=p.v.n,q.thr.p=q.thr.p,q.thr.n=q.thr.n,q.v.p=q.v.p,q.v.n=q.v.n)
        return(list(df=df,ddf=ddf,svm=svm,rvl=rvl))
    } else {
        rm(rvl); gc();
        return(df);
    }
}


#' Iterative bulk gene set enrichment analysis
#'
#' @param set.list list of gene sets
#' @param threshold.eval threshold for applying additional permutations (default: 10)
#' @param n.rand list of number of random permutations used to assess significance (default: c(1e2,1e3,1e4))
#' @param verbose whether to use high verbosity level (default: TRUE)
#' @param ... arguments to be passed to \code{\link{bulk.gsea}}
#'
#' @examples
#' data("org.Hs.GO2Symbol.list")  
#' universe <- unique(unlist(org.Hs.GO2Symbol.list))  # get universe
#' gs <- org.Hs.GO2Symbol.list[[1]]  # get a gene set
#' vals <- rnorm(length(universe), 0, 10)  # simulate values
#' names(vals) <- universe
#' vals[gs] <- rnorm(length(gs), 100, 10)  
#' gs.list <- org.Hs.GO2Symbol.list # get gene sets
#' # reduce n.rand for speed
#' iterative.bulk.gsea(values = vals, set.list = gs.list[1:3], mc.cores = 1, n.rand=100) 
#' 
#' @export
#' 
iterative.bulk.gsea <- function(..., set.list, threshold.eval=10, n.rand=c(1e2,1e3,1e4), verbose=TRUE) {
  
    # initial screen
    if(verbose) { cat(paste("initial: [",format(n.rand[1],scientific=T)," - ",sep="")) }
  
    df <- bulk.gsea(..., set.list = set.list, n.rand=n.rand[1])
    vs <- rownames(df)[df$p.val<=(threshold.eval+1)/(n.rand[1]+1)]
    if(verbose) { cat(paste(length(vs),"] ",sep=""));  }
    
    for(nr in n.rand[-1]) {
        if(length(vs)>0) {
          
            if(verbose) { cat(paste("[",format(nr,scientific=T)," - ",sep=""))  }
          
            dfr <- bulk.gsea(..., set.list = set.list[vs], n.rand=nr, skip.qval.estimation=TRUE)
            df[match(rownames(dfr),rownames(df)),] <- dfr
            vs <- rownames(df)[df$p.val<=(threshold.eval+1)/(nr+1)]
            
            if(verbose) { cat(paste(length(vs),"] ",sep=""))  }
        }
    }
    
    df$q.val <- stats::p.adjust(df$p.val,method="BH")
    
    if(verbose) { cat("done\n")  }
    
    return(df)
}