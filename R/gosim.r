# common.r
sn <- function(x) {
    names(x) <- x
    return(x)
}

# Rcpp implementation of the GSEA resampling procedure
# Used by \code{\link{gsea}}
gseaRandCore <- cxxfunction(signature(Set="integer",Eso="numeric",Nsamples="integer",Seed="integer"),body='
    std::vector<int> sset=Rcpp::as< std::vector<int> >(Set);
    std::vector<double> eso=Rcpp::as< std::vector<double> >(Eso);
    int nsamples=Rcpp::as<int>(Nsamples);
    std::vector<double> pscores(nsamples);
    std::vector<double> nscores(nsamples);
    int nelem=sset.size();

    int seed=Rcpp::as<int>(Seed);
    srand(seed);

    for(int i=0;i<nsamples;i++) {
      // shuffle the set
      std::random_shuffle(sset.begin(),sset.end());
      // determine normalizing factors
      double insum=0; double outsum=0;
      for(int j=0;j<nelem;j++) {
         if(sset[j]) {
            insum+=eso[j];
         } else {
            outsum+=eso[j];
         }
      }
      insum/=((double)nelem);
      outsum=(-1.0)*outsum/((double)nelem);
      // calculate score
      double smaxn=0; double smaxp=0; double cs=0;
      for(int j=0;j<nelem;j++) {
         if(sset[j]) {
            cs+=eso[j]/insum;
         } else {
            cs+=eso[j]/outsum;
         }
         if(cs>smaxp) { smaxp=cs; } else if(cs<smaxn) { smaxn=cs; }
      }
      pscores[i]=smaxp; nscores[i]=smaxn;
    }
    return Rcpp::List::create(Rcpp::Named( "p" ) = wrap(pscores),
                              Rcpp::Named( "n" ) = wrap(nscores));
   ',plugin="Rcpp")

#' Gene set enrichment analysis
#'
#' @param values - vector of values with associated gene names
#' @param set - vector of genes in the gene set
#' @param power - power used in auto-weight calculations
#' @param rank - do not weigh by values (default: FALSE)
#' @param weight - additional weights associated with each value (default: rep(1,length(values)))
#' @param n.rand - number of randomization iterations (default: 1e4)
#' @param plot - whether to plot (default: TRUE)
#'
#' @examples
gsea <- function(values, geneset, power=1, rank=FALSE, weight=rep(1,length(values)), n.rand=1e4, plot=TRUE, return.details=FALSE, quantile.threshold=min(100/n.rand,0.1), random.seed=1, mc.cores=10) {

    # set former options
    decreasing=T
    values.lab="values"
    score.lab="score"
    main.lab=""
    cex=0.9
    body.col="darkblue"
    randomize.order=T

    set <- names(vals) %in% geneset

    # randomize first
    ro <- sample(1:length(values))
    values <- values[ro]; set <- set[ro]; weight <- weight[ro]

    # order values
    vo <- order(values,decreasing=decreasing)
    values <- values[vo]; set <- set[vo]; weight <- weight[vo];
    if(rank) { es <- weight } else { es <- (abs(values)^power)*weight; }
    eso <- es;
    es[set] <- es[set]/sum(es[set]);
    es[!set] <- -1*es[!set]/sum(es[!set]);
    es <- es*length(values);
    sv <- cumsum(es)
    sv.p <- max(0,max(sv));
    sv.n <- min(0,min(sv));
    p.x.p <- which.max(sv);
    p.x.n <- which.min(sv);

    # randomizations
    if(mc.cores>1) {
        rvll <- mclapply(1:mc.cores,function(i) {
            gseaRandCore(set,eso,ceiling(n.rand/mc.cores),random.seed+i)
        },mc.preschedule=T,mc.cores=mc.cores)
        rvl <- list(p=unlist(lapply(rvll,function(x) x$p)),n=unlist(lapply(rvll,function(x) x$n)));
        rm(rvll); gc();
    } else {
        rvl <- gseaRandCore(set,eso,n.rand,random.seed)
    }
    p.v.p <- (sum(rvl$p >= sv.p)+1)/(length(rvl$p+1))
    p.v.n <- (sum(rvl$n <= sv.n)+1)/(length(rvl$n+1))

    p.x <- ifelse(p.v.p<p.v.n,p.x.p,p.x.n)
    p.v <- ifelse(p.v.p<p.v.n,p.v.p,p.v.n)

    if(p.v.p<p.v.n) {
        q.thr <- quantile(rvl$p,p=(1-quantile.threshold))
    } else {
        q.thr <- quantile(rvl$n,p=quantile.threshold)
    }

    if(plot) {
        l <- layout(matrix(c(1,2,3),3,1,byrow=T),c(1,1,1),c(1.4,0.3,0.9),FALSE);
        par(mar = c(0.1,3.5,ifelse(main.lab=="",0.5,3.5),0.5), mgp = c(2,0.65,0), cex = cex)
        plot(c(1:length(sv)),sv,type='l',ylab=score.lab,xaxt="n",xlab="",main=main.lab,xaxs="i",xlim=c(1,length(sv)),col=body.col)
        segments(p.x,0,p.x,sv[p.x],col=2,lty=2)
        segments(0,sv[p.x],p.x,sv[p.x],col=2,lty=2)

        qv.p <- quantile(rvl$p,p=1-10^(-1*seq(1,round(log10(length(rvl$p))))))
        qv.n <- -1*quantile(-1*rvl$n,p=1-10^(-1*seq(1,round(log10(length(rvl$n))))))
        abline(h=c(qv.p,qv.n),lty=3,col=8)
        abline(h=0,lty=2,col=8)
        xpos <- "right"; if(p.x>round(length(sv)/2)) {  xpos <- "left"; }
        legend(x=ifelse(sv[round(length(sv)/2)]>sv[p.x],paste("bottom",xpos,sep=""),paste("top",xpos,sep="")),legend=paste("P-value <",format(p.v,digits=3)),bty="n")
        par(mar = c(0.1,3.5,0.1,0.5))
        #plot(set,xaxt="n",xlab="",ylab="",xaxs="i",yaxt="n",type='h',col=body.col)
        mset <- set; mset[!set] <- NA;
        #den <- densum(which(set),from=0,to=length(set))$y[which(set)]; den <- round((den)*255/max(den)+1);
        #plot(mset,xaxt="n",ylim=c(0,1),xlab="",ylab="",xaxs="i",yaxt="n",type='h',col=colorRampPalette(brewer.pal(9, "Blues")[-(1:2)])(256)[den])
        #plot(mset,xaxt="n",ylim=c(0,1),xlab="",ylab="",xaxs="i",yaxt="n",type='h',col=densCols(which(set),nbin=256,colramp=colorRampPalette(brewer.pal(9, "Blues")[-(1:2)])))
        plot(mset,xaxt="n",ylim=c(0,1),xlab="",ylab="",xaxs="i",yaxt="n",type='h',col="blue")
        abline(v=p.x,col=2,lty=2)
        box()
        par(mar = c(0.5,3.5,0.1,0.5))
        plot(values,ylab=values.lab,xaxs="i",type='h',xaxt="n",col=body.col)
        legend(x=ifelse(decreasing,"topright","topleft"),bty="n",legend=paste("edge value = ",format(values[p.x],digits=2)))
        abline(v=p.x,col=2,lty=2)
        box()
    }
    if(return.details) {
        rv <- c(p.val=p.v,edge.score=as.numeric(sv[p.x]),edge.value=as.numeric(values[p.x]),scaled.score=as.numeric(sv[p.x]/q.thr));
        #if(pareto.estimate) { rv <- c(rv,pareto.pval=pareto.tail.estimate(abs(rvl),as.numeric(abs(sv[p.x])))) }
        return(rv);
    } else {
        return(p.v);
    }
}


# gsea randomization core method for multiple sets
gseaBulkCore <- cxxfunction(signature(SetM="integer",Eso="numeric",Nsamples="integer",Seed="integer"),includes="#include <iostream>",body='
    arma::mat setm=Rcpp::as<arma::mat>(SetM);
    arma::vec eso=Rcpp::as< arma::vec >(Eso);
    int nsamples=Rcpp::as<int>(Nsamples);
    int seed=Rcpp::as<int>(Seed);
    srand(seed);
    int nelem=setm.n_cols;
    int nsets=setm.n_rows;
    arma::mat pscores(nsets,nsamples);
    arma::mat nscores(nsets,nsamples);

    std::vector<int> colord(nelem);
    for(int i=0;i<nelem;i++) { colord[i]=i; }

    arma::vec innv(nsets);
    arma::vec outnv(nsets);
    arma::vec smax(nsets);
    arma::vec smin(nsets);
    arma::vec cs(nsets);

    for(int i=0;i<nsamples;i++) {
      // shuffle the order
      std::random_shuffle(colord.begin(),colord.end());

      innv.zeros(); outnv.zeros();

      // determine normalizing factors
      for(int j=0;j<nelem;j++) {
        int rj=colord[j];
        for(int k=0;k<nsets;k++) {
          if(setm(k,rj)) {
            innv[k]+=eso[rj];
          } else {
            outnv[k]+=eso[rj];
          }

        }
      }
      innv*=(1.0)/((double)nelem);
      outnv*=(-1.0)/((double)nelem);
      // calculate score

      smin.zeros(); smax.zeros(); cs.zeros();
      for(int j=0;j<nelem;j++) {
        int rj=colord[j];
        // update cumulative
        for(int k=0;k<nsets;k++) {
          if(setm(k,rj)) {
            cs[k]+=eso[rj]/innv[k];
         } else {
            cs[k]+=eso[rj]/outnv[k];
         }
         if(cs[k]>smax[k]) { smax[k]=cs[k]; } else if(cs[k]<smin[k]) { smin[k]=cs[k]; }
        }
      }
      pscores.col(i)=smax;
      nscores.col(i)=smin;
    }
    return Rcpp::List::create(Rcpp::Named( "p" ) = wrap(pscores),
                              Rcpp::Named( "n" ) = wrap(nscores));
   ',plugin="RcppArmadillo")

# set.list - a list of character vectors corresponding to sets to be tested
# values must be named, according to names appearing in set.list elements
bulk.gsea <- function(values, set.list, power=1, rank=FALSE, weight=rep(1,length(values)), n.rand=1e4, mc.cores=10, quantile.threshold=min(100/n.rand,0.1), return.details=FALSE, skip.qval.estimation=FALSE) {

    # old options
    cex=0.9
    use.Rcpp=TRUE
    decreasing=TRUE

    # determine set matrix
    setm <- do.call(rbind, mclapply(set.list, function(set) names(values) %in% set, mc.cores=mc.cores))
    # only bother testing if more than 2 genes present
    setm <- setm[rowSums(setm)>2,,drop=F]

    # randomize first
    ro <- sample(seq_along(values))
    values <- values[ro]
    setm <- setm[, ro, drop=F]
    weight <- weight[ro]

    # order
    vo <- order(values,decreasing=decreasing)
    values <- values[vo]; setm <- setm[,vo,drop=F]; weight <- weight[vo];
    if(rank) { es <- weight } else { es <- (abs(values)^power)*weight; }
    eso <- es;

    ies <- t(t(setm)*es); ies <- ies/rowSums(ies);
    ies[is.nan(ies)] <- 0; # when values null out membership
    pes <- t(t(!setm)*es); pes <- pes/rowSums(pes);
    pes[is.nan(pes)] <- 0;
    esm <- ies-pes;
    svm <- t(apply(esm,1,cumsum))*length(values);
    svm.maxp <- apply(svm,1,function(x) x[which.max(x)]); svm.maxp[svm.maxp<0] <- 0;
    svm.maxn <- apply(svm,1,function(x) x[which.min(x)]); svm.maxn[svm.maxn>0] <- 0;
    p.xmp <- apply(svm,1,which.max)
    p.xmn <- apply(svm,1,which.min)

    # randomizations
    if(mc.cores>1) {
        rvlp <- mclapply(1:mc.cores,function(i) {
            gseaBulkCore(setm,eso,ceiling(n.rand/mc.cores),i)
        },mc.preschedule=T,mc.cores=mc.cores)
        rvl <- list(p=do.call(cbind,lapply(rvlp,function(x) x$p)),
                    n=do.call(cbind,lapply(rvlp,function(x) x$n)));
        rm(rvlp); gc();
    } else {
        rvl <- gseaBulkCore(setm,eso,n.rand,1)
    }

    # raw p-values
    p.v.p <- (rowSums(rvl$p - svm.maxp >= 0) +1)/(ncol(rvl$p)+1)
    p.v.n <- (rowSums(rvl$n - svm.maxn <= 0) +1)/(ncol(rvl$n)+1)
    # upper quantiles for scaling thresholds
    q.thr.p <- rowQuantiles(rvl$p,probs=(1-quantile.threshold),drop=T)
    q.thr.n <- rowQuantiles(rvl$n,probs=quantile.threshold,drop=T)

    p.val <- ifelse(p.v.p<p.v.n,p.v.p,p.v.n)
    sscore <- ifelse(p.v.p<p.v.n,svm.maxp/q.thr.p,-1*svm.maxn/q.thr.n)
    x.val <- values[ifelse(p.v.p<p.v.n,p.xmp,p.xmn)]

    if(!skip.qval.estimation) {
        # sign-aware mean-scaling
        s.pm <- apply(rvl$p,1,mean)
        s.nm <- apply(rvl$n,1,mean)
        s.pm[is.na(s.pm)] <- 0; s.nm[is.na(s.nm)] <- 0;

        s.pss <- ecdf(as.numeric(rvl$p/s.pm))
        s.nss <- ecdf(as.numeric(rvl$n/s.nm))

        # scale the actual Smax
        s.svm.maxp <- svm.maxp/s.pm;
        s.svm.maxn <- svm.maxn/s.nm;

        vo <- order(s.svm.maxp,decreasing=F)
        vr <- rank(s.svm.maxp)
        qv <- (1-s.pss(s.svm.maxp-1e-10*max(s.svm.maxp)))/(1-ecdf(s.svm.maxp)(s.svm.maxp-1e-10))
        q.v.p <- cummin(qv[vo])[vr];

        vo <- order(s.svm.maxn,decreasing=F)
        vr <- rank(s.svm.maxn)
        qv <- (1-s.pss(s.svm.maxn-1e-10*max(s.svm.maxn)))/(1-ecdf(s.svm.maxn)(s.svm.maxn-1e-10))
        q.v.n <- cummin(qv[vo])[vr];

        q.val <- ifelse(q.v.p<q.v.n,q.v.p,q.v.n)
    } else {
        q.val <- p.adjust(p.val);
    }

    # p-value, q-value table
    df <- data.frame(p.val=p.val,q.val=q.val,sscore=sscore,edge=x.val)
    rownames(df) <- rownames(setm);

    if(return.details) {
        ddf <- data.frame(svm.maxp=svm.maxp,svm.maxn=svm.maxn,p.xmp=p.xmp,p.xmn=p.xmn,p.v.p=p.v.p,p.v.n=p.v.n,q.thr.p=q.thr.p,q.thr.n=q.thr.n,q.v.p=q.v.p,q.v.n=q.v.n)
        return(list(df=df,ddf=ddf,svm=svm,rvl=rvl))
    } else {
        rm(rvl); gc();
        return(df);
    }
}

# calls bulk.gsea
iterative.bulk.gsea <- function(..., set.list, threshold.eval=10, n.rand=c(1e3,1e4,1e5,1e6),verbose=TRUE) {
    # initial screen
    if(verbose) { cat(paste("initial: [",format(n.rand[1],scientific=T)," - ",sep=""));  }
    df <- bulk.gsea(...,set.list=set.list, n.rand=n.rand[1]);
    #df <- bulk.gsea(values=values,set.list=set.list, power=power,mc.cores=28,n.rand=n.rand[1]);
    vs <- rownames(df)[df$p.val<=(threshold.eval+1)/(n.rand[1]+1)];
    if(verbose) { cat(paste(length(vs),"] ",sep=""));  }
    for(nr in n.rand[-1]) {
        if(length(vs)>0) {
            if(verbose) { cat(paste("[",format(nr,scientific=T)," - ",sep=""));  }
            dfr <- bulk.gsea(..., set.list=set.list[vs],n.rand=nr,skip.qval.estimation=T);
            #dfr <- bulk.gsea(values=values,set.list=set.list[vs], power=power,mc.cores=28,n.rand=nr);
            df[match(rownames(dfr),rownames(df)),] <- dfr;
            vs <- rownames(df)[df$p.val<=(threshold.eval+1)/(nr+1)];
            if(verbose) { cat(paste(length(vs),"] ",sep=""));  }
        }
    }
    df$q.val <- p.adjust(df$p.val,method="BH")
    if(verbose) { cat("done\n");  }
    # update qvalues
    return(df);
}






# move into examples or saved datasets

make.mouse.go.tables <- function() {
    # make org.Mm.symbol2GO
    library(org.Mm.eg.db)
    library(GO.db);

    etr <- as.list(org.Mm.egENSEMBL[mappedkeys(org.Mm.egENSEMBL)])
    go.valid.codes <- c("IDA","IPI","IMP","IGI","IEP","ISS","TAS")

    # ENSEMBL2GO
    x <- as.list(org.Mm.egGO2ALLEGS)
    gl <- mclapply(x,function(d) unique(unlist(etr[as.character(d)[names(d) %in% go.valid.codes]])),mc.cores=10)
    ggs <- do.call(rbind,lapply(sn(names(gl)),function(go) cbind(gene=as.character(gl[[go]]),go=rep(go,length(gl[[go]])))))
    gl <- tapply(ggs[,2],as.factor(ggs[,1]),I)
    org.Mm.ENSEMBL2GO <- new.env(parent=globalenv());
    x <- lapply(names(gl),function(g) assign(g,gl[[g]],envir=org.Mm.ENSEMBL2GO))

    # symbol2GO
    etr <- as.list(org.Mm.egSYMBOL[mappedkeys(org.Mm.egSYMBOL)])
    x <- as.list(org.Mm.egGO2ALLEGS)
    gl <- mclapply(x,function(d) unique(unlist(etr[as.character(d)[names(d) %in% go.valid.codes]])),mc.cores=10)
    ggs <- do.call(rbind,lapply(sn(names(gl)),function(go) cbind(gene=as.character(gl[[go]]),go=rep(go,length(gl[[go]])))))
    gl <- tapply(ggs[,2],as.factor(ggs[,1]),I)
    org.Mm.Symbol2GO <- new.env(parent=globalenv());
    x <- lapply(names(gl),function(g) assign(g,gl[[g]],envir=org.Mm.Symbol2GO))
    rm(x); gc();

    # GO2ENSEMBL
    org.Mm.GO2ENSEMBL <- new.env(parent=globalenv());
    xn <- ls(env=org.Mm.ENSEMBL2GO)
    xl <- mget(xn,env=org.Mm.ENSEMBL2GO)
    gel <- tapply(rep(xn,unlist(lapply(xl,length))),unlist(xl),I)
    x <- lapply(names(gel),function(n) assign(n,gel[[n]],envir=org.Mm.GO2ENSEMBL));
    rm(xn,xl,x,gel); gc();

    # GO2Symbol
    org.Mm.GO2Symbol <- new.env(parent=globalenv());
    xn <- ls(env=org.Mm.Symbol2GO)
    xl <- mget(xn,env=org.Mm.Symbol2GO)
    gel <- tapply(rep(xn,unlist(lapply(xl,length))),unlist(xl),I)
    x <- lapply(names(gel),function(n) assign(n,gel[[n]],envir=org.Mm.GO2Symbol));
    rm(xn,xl,x,gel); gc();

    #save(org.Mm.ENSEMBL2GO,org.Mm.GO2ENSEMBL,org.Mm.Symbol2GO,org.Mm.GO2Symbol,file="~/keith/me3/org.Mm.GOenvs.RData")

}

make.human.go.tables <- function() {
    library(org.Hs.eg.db)
    library(GO.db);

    etr <- as.list(org.Hs.egENSEMBL[mappedkeys(org.Hs.egENSEMBL)])
    go.valid.codes <- c("IDA","IPI","IMP","IGI","IEP","ISS","TAS")

    # ENSEMBL2GO
    x <- as.list(org.Hs.egGO2ALLEGS)
    gl <- mclapply(x,function(d) unique(unlist(etr[as.character(d)[names(d) %in% go.valid.codes]])),mc.cores=10)
    #gl <- mclapply(x,function(d) unique(unlist(etr[as.character(d)])),mc.cores=10) # all codes
    gl <- gl[!(names(gl) %in% get.top.go.terms(level=2))]
    gl <- gl[unlist(lapply(gl,length))>5]

    ggs <- do.call(rbind,lapply(sn(names(gl)),function(go) cbind(gene=as.character(gl[[go]]),go=rep(go,length(gl[[go]])))))
    gl <- tapply(ggs[,2],as.factor(ggs[,1]),I)
    str(ggs)

    org.Hs.ENSEMBL2GO <- new.env(parent=globalenv());
    x <- lapply(names(gl),function(g) assign(g,gl[[g]],envir=org.Hs.ENSEMBL2GO))


    # symbol2GO
    etr <- as.list(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])
    x <- as.list(org.Hs.egGO2ALLEGS)
    gl <- mclapply(x,function(d) unique(unlist(etr[as.character(d)[names(d) %in% go.valid.codes]])),mc.cores=10)
    #gl <- mclapply(x,function(d) unique(unlist(etr[as.character(d)])),mc.cores=10) # all codes
    #gl <- gl[!(names(gl) %in% get.top.go.terms(level=2))]
    #gl <- gl[unlist(lapply(gl,length))>5];
    ggs <- do.call(rbind,lapply(sn(names(gl)),function(go) cbind(gene=as.character(gl[[go]]),go=rep(go,length(gl[[go]])))))
    gl <- tapply(ggs[,2],as.factor(ggs[,1]),I)
    org.Hs.Symbol2GO <- new.env(parent=globalenv());
    x <- lapply(names(gl),function(g) assign(g,gl[[g]],envir=org.Hs.Symbol2GO))
    rm(x); gc();


    # GO2ENSEMBL
    org.Hs.GO2ENSEMBL <- new.env(parent=globalenv());
    xn <- ls(env=org.Hs.ENSEMBL2GO)
    xl <- mget(xn,env=org.Hs.ENSEMBL2GO)
    gel <- tapply(rep(xn,unlist(lapply(xl,length))),unlist(xl),I)
    x <- lapply(names(gel),function(n) assign(n,gel[[n]],envir=org.Hs.GO2ENSEMBL));
    rm(xn,xl,x,gel); gc();


    # GO2Symbol
    org.Hs.GO2Symbol <- new.env(parent=globalenv());
    xn <- ls(env=org.Hs.Symbol2GO)
    xl <- mget(xn,env=org.Hs.Symbol2GO)
    gel <- tapply(rep(xn,unlist(lapply(xl,length))),unlist(xl),I)
    x <- lapply(names(gel),function(n) assign(n,gel[[n]],envir=org.Hs.GO2Symbol));
    rm(xn,xl,x,gel); gc();

    #save(org.Hs.ENSEMBL2GO,org.Hs.GO2ENSEMBL,org.Hs.Symbol2GO,org.Hs.GO2Symbol,file="~/keith/me3/org.Hs.GOenvs.RData")

}

# sample runner
# need to make org.Hs.GO2Symbol first
main <- function() {

    # get universe
    universe <- unique(unlist(as.list(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])))
    # get a gene set
    go <- ls(org.Hs.GO2Symbol)[[1]]
    gs <- get(go, org.Hs.GO2Symbol)
    # fake dummy example where everything in gene set is perfectly enriched
    vals <- rnorm(length(universe), 0, 10)
    names(vals) <- universe
    vals[gs] <- rnorm(length(gs), 100, 10)
    barplot(sort(vals, decreasing=TRUE))
    # test obviously enriched set
    gstest <- gs
    gsea(values=vals, geneset=gstest)
    # hard to see...too many genes in universe

    # test obviously not enriched set
    go <- ls(org.Hs.GO2Symbol)[[2]]
    gs <- get(go, org.Hs.GO2Symbol)
    gstest <- gs
    gsea(values=vals, geneset=gstest)

    # universe is huge...make smaller universe example
    # get universe
    universe <- unique(unlist(as.list(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])))[1:1000]
    # get a gene set
    gs <- universe[1:10]
    # fake dummy example where everything in gene set is perfectly enriched
    vals <- rnorm(length(universe), 0, 10)
    names(vals) <- universe
    vals[gs] <- rnorm(length(gs), 100, 10)
    # add some noise
    vals[sample(1:length(universe), 100)] <-  rnorm(length(gs), 100, 10)
    barplot(sort(vals, decreasing=TRUE))
    gstest <- gs
    gsea(values=vals, geneset=gstest)

    #pvl <- iterative.bulk.gsea(
    #    values=vals,
    #    set.list=gsl)

    bulk.gsea(values, gel[1:10])

}
