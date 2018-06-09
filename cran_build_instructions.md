1. Build a source package:
```
R CMD build .
```

2. Look in the package to make sure it doesn't contain spurious files. Modify .Rbuildignore as needed.
```
tar tvf liger_*.tar.gz
```

3. Check the package:
```
R CMD check --as-cran liger_*.tar.gz
```

