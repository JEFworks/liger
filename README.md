# LIGER - Lightweight Iterative Gene set Enrichment in R

**[See vignette for tutorial ☞ Gene Set Enrichment Analysis with LIGER](vignettes/gsea.pdf)**

Gene Set Enrichment Analysis (GSEA) is a computational method that determines whether an a priori defined set of genes shows statistically significant, concordant differences between two biological states. The original algorithm is detailed in [Subramanian, Tamayo, et al.](http://www.pnas.org/content/102/43/15545.full.pdf) with Java implementations available through the [Broad Institute](http://www.broadinstitute.org/gsea/index.jsp).

The `liger` package provides a lightweight R implementation of this enrichment test on a list of values. Given a list of values, such as p-values or log-fold changes derived from differential expression analysis or other analyses comparing biological states, this package enables you to test a priori defined set of genes for enrichment to enable interpretability of highly significant or high fold-change genes.

# Sample plots

## Testing individual gene sets
```
gsea(values=vals, geneset=gs, mc.cores=1, plot=TRUE)
```

![](images/gsea_tp.png)
![](images/gsea_tn.png)

## Testing multiple gene sets

```
bulk.gsea(values=vals, set.list=org.Hs.GO2Symbol.list[1:10])
```

```
iterative.bulk.gsea(vals, set.list=org.Hs.GO2Symbol.list[1:10])
```

```
> initial: [1e+02 - 1] [1e+03 - 1] [1e+04 - 1] done
                p.val      q.val     sscore       edge
		GO:0000002 0.00009999 0.00059994  2.6054741  70.912194
		GO:0000003 0.22772277 0.34158416  0.9048369  13.170093
		GO:0000012 0.44554455 0.44554455  0.7918879   8.392397
		GO:0000014 0.18811881 0.34158416  0.8532878  -4.458762
		GO:0000018 0.18811881 0.34158416 -0.9206904  11.111976
		GO:0000022 0.42574257 0.44554455 -0.6432124 -11.015244
```
