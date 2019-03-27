test_that(context("Test GSEA"), {
  library(liger)
  data("org.Hs.GO2Symbol.list")  
  universe <- unique(unlist(org.Hs.GO2Symbol.list))  # get universe
  gs <- org.Hs.GO2Symbol.list[[1]]  # get a gene set
  vals <- rnorm(length(universe), 0, 10)  # simulate values
  names(vals) <- universe
  vals[gs] <- rnorm(length(gs), 100, 10)  
  gs.list <- org.Hs.GO2Symbol.list # get gene sets
  
  # test obviously enriched set, reduce n.rand for speed
  pv1 <- gsea(values=vals, geneset=gs, mc.cores=1, n.rand=100) 
  expect_equal(pv1 < 0.05, TRUE)
  
  pv2 <- gsea(values=vals, geneset=gs, mc.cores=1, n.rand=100, rank=TRUE) 
  expect_equal(pv2 < 0.05, TRUE)
})

test_that(context("Test speed"), {
  library(liger)
  data("org.Hs.GO2Symbol.list")  
  universe <- unique(unlist(org.Hs.GO2Symbol.list))  # get universe
  gs <- org.Hs.GO2Symbol.list[[1]]  # get a gene set
  vals <- rnorm(length(universe), 0, 10)  # simulate values
  names(vals) <- universe
  vals[gs] <- rnorm(length(gs), 100, 10)  
  gs.list <- org.Hs.GO2Symbol.list # get gene sets
 
  start_time <- Sys.time()
  bulk.gsea(values = vals, set.list = gs.list[1:3], mc.cores = 1, n.rand=100)
  end_time <- Sys.time()
  t1 <- end_time - start_time
 
  start_time <- Sys.time()
  iterative.bulk.gsea(values = vals, set.list = gs.list[1:3], mc.cores = 1, n.rand=100) 
  end_time <- Sys.time()
  t2 <- end_time - start_time
 
  expect_equal(t2 <= t1, TRUE)
})