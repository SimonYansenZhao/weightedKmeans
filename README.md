# weightedKmeans: An R Package for Weighted k-means Clustering

The [weightedKmeans](http://cran.r-project.org/package=weightedKmeans)
is an R package for weighted k-means clustering.  And this repo is
used for the next final version, because all the work afterwards will
be continued on the new package
[wksm](http://cran.r-project.org/package=wskm).

Entropy weighted kmeans (ewkm) is a weighted subspace clustering
algorithm that is well suited to very high dimensional data. Weights
are calculated as the importance of a variable with regard to cluster
membership. The feature group weighted kmenas (fgkm) extends this
concept by grouping features and weighting the group in addition to
weihgting individual features.

## Plan

As more functions and utilities are put into the package, we need a
better name for it.  The _weightedKmeans_ package will be left as a
transitional one for those R packages who depend on it.  And no
further updates will be made after we finish this transition.  

## Related links

* [weightedKmeans on CRAN](http://cran.r-project.org/web/packages/weightedKmeans/index.html)
* [wskm on CRAN](http://cran.r-project.org/web/packages/wskm/index.html)
* [wskm on GitHub](https://github.com/SimonYansenZhao/wskm)


## License

GPL (>= 3)
