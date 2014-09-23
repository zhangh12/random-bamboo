Random Bamboo
=============


#Introduction

`Random Bamboo` (`RB`) is an ultra-fast C++ implement of the [`Random Forest`](https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm) algorithm. The program is designed for analyzing the large-scale [genome-wide association study](http://en.wikipedia.org/wiki/Genome-wide_association_study) (GWAS) and building predictive model on hundreds of thousands of SNPs for personalized medicine. `RB` adopts fancy [data structure](http://bioinformatics.oxfordjournals.org/content/30/15/2171) and [Boolean operation method](http://bioinformatics.ust.hk/BOOST.html) to speedup the algorithm dramatically. The single-thread version is 40 times faster than the `R` package [`randomForest`](http://cran.r-project.org/web/packages/randomForest/index.html). [`OpenMP`](http://openmp.org/wp/) is used to enable parallelization feature. `RB` is thus appliable to large study with thousands or tens of thousands of individuals. 

#Usage

`RB` requires input data with format defined by [`PLINK`](http://pngu.mgh.harvard.edu/~purcell/plink/). To accelerate loading the training data, the binary files should be generated via `PLINK` before applying `RB`. Suppose we have large `train.ped` file containing genotypes information and a mapping file `train.map`, both are described [here](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml). The following command generates files `train.bed`, `train.bim` and `train.fam`
```
plink --file train --out train --make-bed --noweb
```

Similarly we generate the testing data `test.bed`, `test.bim` and `test.fam` from `test.ped` and `test.map`
```
plink --file test --out test --make-bed --noweb
```

Run the following command to train a model with SNPs only
```
bamboo --file train --mtry 5000 --ntree 10000 --out rb
```

files `rb.cof`, `rb.err`, `rb.imp` and `rb.bam` will be saved, which contain the confusion matrix, oob error, Gini importance and trained model, respectively.

If continuous and categorical covariates are available in files `train.con` and `train.cat`, run the following command to train a model with covariates
```
bamboo --file train --cont train --cate train --out rb --mtry 5000 --ntree 10000
```
Please note that the first line of `train.con` and `train.cat` is the header and the first column is the individual IDs used to align the genotypes information.

We can predict testing data after training a model
```
bamboo --file train --cont train --cate train --out rb --mtry 5000 --ntree 10000 --pred test
```
Five files `test.bed`, `test.bim`, `test.fam`, `test.con` and `test.cat` are needed in prediction. Please note that the files' names of training data can be different, but they are required to be the same in testing data. The predicton results are saved in `rb.prd`.

We can use single bam file (`rb.bam`) saved in local to predict
```
bamboo --pred test --bam rb --out rb
```
If multiple bam files are save in a directory `./path`, we can use all of them to make a prediction
```
bamboo --pred test --bam ./path --out rb
```
`RB` will parse the meaning of `--bam` automatically.

`RB` can print information of given bam file(s)
```
bamboo --bam rb
bamboo --bam ./model
```

`RB` parses the input files according to their extensions. Only the file name without extension is needed when invoking `RB`.


#Options


* `-f`, `--file`. File names of three PLINK-formated datasets, including `bed`, `bim` and `fam`. Character. No default.
* `-c`, `--cont`. File name of continuous or ordered covariates, with extension `con`. The first row is header. The first column contains individual IDs. Character. No default.
* `-a`, `--cate`. File name of categorical covariates, with extension `cat`. The first row is header. The first column contains individual IDs. Character. No default.
* `-o`, `--out`. File name of all output files. Character. If unspecified, it is set as `--file` or `--pred`.
* `-p`, `--pred`. File names of datasets used in prediction, including `bed`, `bim`, `fam`, `con` and `cat`. The last two are optional if no covariates are used in training model. Character. No default.
* `-b`, `--bam`. Used to specify trained model. It can be the name of a single `bam` file or a directory contains multiple `bam` files. Character. No default.
* `-y`, `--trainid`. File name of individual IDs used in training model. The file contains one column and no header. With extension `iid`. If unspecified, all individuals in the datasets are used.
* `-z`, `--testid`. File name of individual IDs used in prediction. If `--file` is specified, the individuals are predicted with model trained from the specified datasets. If `--bam` is specified, the individuals are predicted with specified model(s). The file contains one column and no header. With extension `iid`. Character. If unspecified, all individuals in the datasets are used.
* `-S`, `--snpid`. File name of SNPs' names used in training model. The file contains one column and no header. With extension `sid`. Character. No default.
* `-t`, `--ntree`. Number of trees in the forest. Positive integer. Default: 1.
* `-m`, `--mtry`. Number of candidate variables in determining best split in each node. Positive integer. Default: sqrt of total number of variables.
* `-s`, `--seed`. Random seed. Positive integer. Default: 1.
* `-l`, `--maxnleaf`. Maximum leaves in a single tree. Positive integer. Default: 1000000.
* `-e`, `--minleafsize`. Minimum sample size in a single tree. Positive integer. Default: 1.
* `-i`, `--imp`. Types of variable importance. 1 - Gini; 2 - Breiman & Cutler's; 3 - Liaw & Wiener's; 4 - Raw; 5 - Meng's; 6 - All. Default: 1.
* `-d`, `--nthread`. Number of threads in parallelization. Positive integer. Default: maximum available number of CPUs allowed by the hardware.
* `-w`, `--classwt`. Class weight assigned to the case group. Positive double. Default: 1.0.
* `-u`, `--cutoff`. Unimplemented.
* `-g`, `--noflip`. Keep the genotypes of SNPs with MAF > 0.5 unchanged. Default: switched off.
* `-x`, `--prox`. Save the proximity matrix to local file. Default: switched off.
* `-n`, `--noimp`. Don't save variable importances to local file. Default: switched off.
* `-B`, `--balance`. Balance the ratio of cases and controls in training data to 1. Default: switched off.
* `-r`, `--trace`. Print debug information. Valid in single-thread mode. Default: switched off.
* `-N`, `--nobam`. Don't save trained model to local file. Default: switched off.
* `-h`, `--help`. Print help information. Unimplemented. Default: switched off.
* `-v`, `--version`. Print version information. Unimplemented. Default: switched off.


#Todo

* [ ] allow specifying must-included variables (e.g., STUDY, AGE, etc.)
* [ ] optimize memory and speed further
* [ ] mute output except `bam` file
* [ ] print data information (e.g., MAF, etc.)
* [ ] filter genotypes by MAF
* [ ] test AUC function with tie
* [ ] allow specifying cutoff in prediction
* [ ] allow individual weights in sampling
* [ ] print help information


