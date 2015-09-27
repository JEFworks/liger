# test
require(Rcpp); require(RcppArmadillo); require(parallel); require(matrixStats)
set.seed(0)

# load data
data("org.Hs.GO2Symbol.list")  

# get universe
universe <- unique(unlist(org.Hs.GO2Symbol.list))
# get a gene set
gs <- org.Hs.GO2Symbol.list[[1]]
# fake dummy example where everything in gene set is perfectly enriched
vals <- rnorm(length(universe), 0, 10)
names(vals) <- universe
vals[gs] <- rnorm(length(gs), 100, 10)
# test obviously enriched set
gsea(values=vals, geneset=gs, mc.cores=1)
  
# test obviously not enriched set
gs <- org.Hs.GO2Symbol.list[[2]]
gsea(values=vals, geneset=gs)
  
# add some noise
vals[sample(1:length(universe), 1000)] <-  rnorm(1000, 100, 10)
# test previously perfectly enriched gene set again
gs <- org.Hs.GO2Symbol.list[[1]]
gsea(values=vals, geneset=gs)
  
# bulk test
ptm <- proc.time()
bulk.gsea(values=vals, set.list=org.Hs.GO2Symbol.list[1:10])
proc.time() - ptm
# iterative test
ptm <- proc.time()
iterative.bulk.gsea(values=vals, set.list=org.Hs.GO2Symbol.list[1:10])
proc.time() - ptm
