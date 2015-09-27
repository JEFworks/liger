# How rda data files were generated
# (later removed many to save space)

# common.r
sn <- function(x) {
  names(x) <- x
  return(x)
}

make.mouse.go.lists <- function() {
  library(org.Mm.eg.db)
  library(GO.db)
  
  etr <- as.list(org.Mm.egENSEMBL[mappedkeys(org.Mm.egENSEMBL)])
  go.valid.codes <- c("IDA","IPI","IMP","IGI","IEP","ISS","TAS")
  
  # ENSEMBL2GO
  x <- as.list(org.Mm.egGO2ALLEGS)
  gl <- mclapply(x,function(d) unique(unlist(etr[as.character(d)[names(d) %in% go.valid.codes]])),mc.cores=10)
  ggs <- do.call(rbind,lapply(sn(names(gl)),function(go) cbind(gene=as.character(gl[[go]]),go=rep(go,length(gl[[go]])))))
  gl <- tapply(ggs[,2],as.factor(ggs[,1]),I)
  org.Mm.ENSEMBL2GO <- new.env(parent=globalenv());
  x <- lapply(names(gl),function(g) assign(g,gl[[g]],envir=org.Mm.ENSEMBL2GO))
  rm(x); gc();
  
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
  xn <- ls(env=org.Mm.ENSEMBL2GO)
  xl <- mget(xn,env=org.Mm.ENSEMBL2GO)
  gel <- tapply(rep(xn,unlist(lapply(xl,length))),unlist(xl),I)
  org.Mm.GO2ENSEMBL.list <- gel
  
  devtools::use_data(org.Mm.GO2ENSEMBL.list)
  
  # GO2Symbol
  xn <- ls(env=org.Mm.Symbol2GO)
  xl <- mget(xn,env=org.Mm.Symbol2GO)
  gel <- tapply(rep(xn,unlist(lapply(xl,length))),unlist(xl),I)
  org.Mm.GO2Symbol <- gel
  
  devtools::use_data(org.Mm.GO2Symbol.list)
  
}

make.human.go.lists <- function() {
  library(org.Hs.eg.db)
  library(GO.db)
  
  etr <- as.list(org.Hs.egENSEMBL[mappedkeys(org.Hs.egENSEMBL)])
  go.valid.codes <- c("IDA","IPI","IMP","IGI","IEP","ISS","TAS")
  
  # ENSEMBL2GO
  x <- as.list(org.Hs.egGO2ALLEGS)
  gl <- mclapply(x,function(d) unique(unlist(etr[as.character(d)[names(d) %in% go.valid.codes]])),mc.cores=10)
  ggs <- do.call(rbind,lapply(sn(names(gl)),function(go) cbind(gene=as.character(gl[[go]]),go=rep(go,length(gl[[go]])))))
  gl <- tapply(ggs[,2],as.factor(ggs[,1]),I)
  org.Hs.ENSEMBL2GO <- new.env(parent=globalenv());
  x <- lapply(names(gl),function(g) assign(g,gl[[g]],envir=org.Hs.ENSEMBL2GO))
  rm(x); gc();
  
  # symbol2GO
  etr <- as.list(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)])
  x <- as.list(org.Hs.egGO2ALLEGS)
  gl <- mclapply(x,function(d) unique(unlist(etr[as.character(d)[names(d) %in% go.valid.codes]])),mc.cores=10)
  ggs <- do.call(rbind,lapply(sn(names(gl)),function(go) cbind(gene=as.character(gl[[go]]),go=rep(go,length(gl[[go]])))))
  gl <- tapply(ggs[,2],as.factor(ggs[,1]),I)
  org.Hs.Symbol2GO <- new.env(parent=globalenv());
  x <- lapply(names(gl),function(g) assign(g,gl[[g]],envir=org.Hs.Symbol2GO))
  rm(x); gc();
  
  # GO2ENSEMBL
  xn <- ls(env=org.Hs.ENSEMBL2GO)
  xl <- mget(xn,env=org.Hs.ENSEMBL2GO)
  gel <- tapply(rep(xn,unlist(lapply(xl,length))),unlist(xl),I)
  org.Hs.GO2ENSEMBL.list <- gel
  
  devtools::use_data(org.Hs.GO2ENSEMBL.list)
  
  # GO2Symbol
  org.Hs.GO2Symbol <- new.env(parent=globalenv());
  xn <- ls(env=org.Hs.Symbol2GO)
  xl <- mget(xn,env=org.Hs.Symbol2GO)
  gel <- tapply(rep(xn,unlist(lapply(xl,length))),unlist(xl),I)
  org.Hs.GO2Symbol.list <- gel
  
  devtools::use_data(org.Hs.GO2Symbol.list)
  
}