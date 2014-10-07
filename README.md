Random Bamboo
=============


#Introduction

`Random Bamboo` (`RB`) is an ultra-fast C++ implement of the [`Random Forest`](https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm) algorithm. The program is designed for analyzing the large-scale [genome-wide association study](http://en.wikipedia.org/wiki/Genome-wide_association_study) (GWAS) and building predictive model on hundreds of thousands of SNPs for personalized medicine. `RB` adopts novel [data structure](http://bioinformatics.oxfordjournals.org/content/30/15/2171) and [Boolean operation method](http://bioinformatics.ust.hk/BOOST.html) to speedup the algorithm dramatically. The single-thread version is 40 times faster than the `R` package [`randomForest`](http://cran.r-project.org/web/packages/randomForest/index.html). [`OpenMP`](http://openmp.org/wp/) is used to enable parallelization feature. `RB` is thus appliable to large study with thousands or tens of thousands of individuals. I name this fast program as "Bamboo" because bamboos are some of the [fastest-growing plants](http://en.wikipedia.org/wiki/Bamboo) in the world.

#Usage

`RB` requires input data with format defined by [`PLINK`](http://pngu.mgh.harvard.edu/~purcell/plink/). To accelerate loading the genotype data, binary files should be generated via `PLINK` before applying `RB`. Suppose we have large `train.ped` file containing genotypes information and a mapping file `train.map`, both are described [here](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml). The following command generates three *genotype* files `train.bed`, `train.bim` and `train.fam`
```
plink --file train --out train --make-bed --noweb
```
Similarly we generate the testing data `test.bed`, `test.bim` and `test.fam` from `test.ped` and `test.map`
```
plink --file test --out test --make-bed --noweb
```

**NOTE** The three genotype files MUST share the same name.

**NOTE** The sixth column in `fam` file is the phenotype, with 1 as control and 2 as case.

**NOTE** Although six columns are in `fam` file, `RB` only uses the second (INDIVIDUAL ID) and sixth (PHENOTYPE) columns. The SEX in the fifth column is ignored. If SEX is required to be used in training, place it in covariate file `cat` (character) or covariate file `con` (dummy variable).

**NOTE** `RB` parses the input files according to their extensions. Only the file name without extension is needed when invoking `RB`.

Run the following command to train a model with SNPs only
```
bamboo --file train --out rb --mtry 5000 --ntree 10000
```
Four files `rb.cof`, `rb.err`, `rb.imp` and `rb.bam` will be saved to local, which contain the confusion matrix, oob error, Gini importance and trained model, respectively. The `bam` file is generally large, which can be muted by switching on the option `--nobam` if the trained model will not be used in future prediction.

If continuous and categorical covariates are available in *covariate* files `train.con` and `train.cat`, run the following command to train a model with covariates
```
bamboo --file train --cont train --cate train --out rb [...]
```
`[...]` can be any other available options to customize the training stage. The first line of the covariate files is the header and the first column is the individual IDs used for aligning the genotypes information given in genotype files. Below is an example of covariate files

```{r}
> fam <- read.table("train.fam", header = FALSE, as.is = TRUE)
> cat <- read.table("train.cat", header = TRUE, as.is = TRUE)
> con <- read.table("train.con", header = TRUE, as.is = TRUE)
> head(fam, 3)
            V1           V2 V3 V4 V5 V6
1 CG-L03-81390 CG-L03-82390  0  0  2  2
2 CG-L03-24262 CG-L03-22262  0  0  2  2
3 CG-L03-69034 CG-L03-69734  0  0  2  2
> head(cat, 3)
       GWAS_ID AGE_CAT GENDER STUDY CIGDAY_CAT CIG_CAT          QUIT_CAT
1 CG-L02-83435  61to65   MALE  PLCO     21to30 CURRENT CURRENT_NEVER_lt1
2 CG-L02-90230  61to65   MALE  PLCO     21to30 CURRENT CURRENT_NEVER_lt1
3 CG-L02-67314  66to70   MALE  PLCO     11to20 CURRENT CURRENT_NEVER_lt1
> head(con, 3)
       GWAS_ID EAGLE_EV2     PLCO_EV4     PLCO_EV5 ATBC_EV2
1 CG-L02-83435         0 -0.005957644  0.010665938        0
2 CG-L02-90230         0 -0.017674167 -0.006934786        0
3 CG-L02-67314         0  0.003294853 -0.010221391        0
```

**NOTE** We can create dummy variables for categorical covariates and treat them as ordinal covariates. In that case, only `con` file is needed.

**NOTE** If the individual contained in the three files are not the same, `RB` will find the intersection and train the the model on it.

We can predict testing data after training a model
```
bamboo --file train --cont train --cate train --out rb --pred test [...]
```
We can also use single bam file (`rb.bam`) saved in local to predict
```
bamboo --pred test --bam rb --out rb
```
If multiple `bam` files are save in a directory `./path`, we can use all of them to make a prediction
```
bamboo --pred test --bam ./path --out rb
```
`RB` will parse the meaning of `--bam` automatically. The feature is useful as we can run `RB` on multiple unparallelizable nodes, each with option `--nthread` enabled and different random seed set by `--seed`, and then all outputed `bam` files are used together in prediction.

**NOTE** The prediction results are saved in `prd` file.

**NOTE** Genotype files are always needed in prediction. If the model is trained with covariates, then `con` and `cat` files are needed as well.

**NOTE** The names of covariate files can be different from each other in training stage. They can even be different from the names of genotype files. For example, the following command works fine
```
bamboo --file geno_train --cont cont_covar --cate cate_covar --out rb [...]
```
However, when we are using the outputed `rb.bam` to predict new data, the covariate files and genotype files MUST share the same name. In the example below, all the data files share the name `test` as only the option `--pred` can be used to specify the data in predicting stage
```
bamboo --bam rb --pred test --out rb
```

Sometimes we prefer to use a subset of data in training or predicting. For example, all available genotype data are imputed and saved in the same genotype files; or all covariates are put in the same covariate files. We can use option `--trainid` to specify which individuals are going to be used in training, and use option `--testid` to specify those used in predicting. For example, we have data files `all_data.fam`, `all_data.bim`, `all_data.bed`, `all_data.con`, `all_data.cat` and two files `train_id.iid` and `test_id.iid`. Both of the `iid` files contain one column of individual IDs. They look like
```
CG-L02-61545
CG-L02-03491
CG-L02-67300
... ...
```
The following command builds model with individuals listed in `train_id.iid` only
```
bamboo --file all_data --cont all_data --cate all_data --trainid train_id --out rb [...]
```
With the model (`rb.bam`) fitted by the command above, we can predict the data specified in `test_id.iid`
```
bamboo --bam rb --pred all_data --testid test_id
```
Or we can do the above in one line
```
bamboo --file all_data --cont all_data --cate all_data --trainid train_id --pred all_data --testid test_id --out rb [...]
```

**NOTE** The file specified by the options `--trainid` and `--testid` must have extention `iid`.

As suggested by Leo Breiman in [here](https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm#varimp)

> ... If the number of variables is very large, forests can be run once with all the variables, then run again using only the most important variables from the first run ... 

We can train a forest using specified markers only by using the option `--snpid`. A file with extension `sid` is needed, in which all SNPs to be used in training a model are listed in one column.
```
bamboo --file train --snpid snp [...]
```
The `snp.sid` looks like
```
rs20112
rs320078
rs213980
...
```

To print summary information of given bam file(s)
```
bamboo --bam rb
bamboo --bam ./path
```


#Options


* `-f`, `--file`. File names of three PLINK-formated genotype files, with extensions `bed`, `bim` and `fam`. Character. No default.
* `-c`, `--cont`. File name of continuous or ordinal covariates, with extension `con`. The first row is header. The first column contains individual IDs. Character. No default.
* `-a`, `--cate`. File name of categorical covariates, with extension `cat`. The first row is header. The first column contains individual IDs. Character. No default.
* `-o`, `--out`. File name of all output files. Character. If unspecified, it is set by `--file` or `--pred`.
* `-p`, `--pred`. File names of datasets used in prediction, including `bed`, `bim`, `fam`, `con` and `cat`. The last two are optional and is activated if covariates are used in training model. Character. No default.
* `-b`, `--bam`. Used to specify model(s) used in prediction. It can be the name of a single `bam` file or a directory contains multiple `bam` files. Character. No default.
* `-y`, `--trainid`. File name of individual IDs used in training model. The file, with extension `iid`, contains one column and no header. Character. If unspecified, all individuals in the datasets are used.
* `-z`, `--testid`. File name of individual IDs used in prediction. If `--file` is specified, the individuals are predicted with model trained from the specified datasets. If `--bam` is specified, the individuals are predicted with specified model(s). The file, with extension `iid`, contains one column and no header. Character. If unspecified, all individuals in the datasets are used.
* `-S`, `--snpid`. File name of SNPs' names used in training model. The file, with extension `sid`, contains one column and no header. Character. If unspecified, all SNPs in the training dataset are used.
* `-t`, `--ntree`. Number of trees in the forest. Positive integer. Default: 1.
* `-m`, `--mtry`. Number of candidate variables in determining best split in each node. Positive integer. Default: sqrt of total number of variables.
* `-s`, `--seed`. Random seed. Positive integer. Default: 1.
* `-l`, `--maxnleaf`. Maximum leaves in a single tree. Positive integer. Default: 1000000.
* `-e`, `--minleafsize`. Minimum sample size in a single tree. Positive integer. Default: 1.
* `-i`, `--imp`. Types of variable importance. 1 - Gini; 2 - Permuted (Breiman & Cutler's fortran verion); 3 - Permuted (Liaw & Wiener's `R` package `randomForest`); 4 - Permuted (Raw, unrescaled); 5 - Permuted (Meng's rescaled method); 6 - print Gini and all permuted importance. Default: 1.
* `-d`, `--nthread`. Number of threads in parallelization. Positive integer. Default: maximum available number of CPUs allowed by the hardware.
* `-w`, `--classwt`. Class weight assigned to the case group. Positive double. Default: 1.0.
* `-u`, `--cutoff`. Unimplemented.
* `-g`, `--flip`. Flip the genotypes of SNPs with MAF > 0.5. Default: switched off. Switching on this option is tricky and thus should be used with caution. 
* `-x`, `--prox`. Save the proximity matrix to local file. Default: switched off. Computing proximity matrix is memory-consuming. 
* `-n`, `--noimp`. Don't save variable importances to local file. Default: switched off.
* `-B`, `--balance`. Balance the ratio of cases and controls in training data to 1. Default: switched off.
* `-r`, `--trace`. Print debug information. Valid in single-thread mode. Default: switched off.
* `-N`, `--nobam`. Don't save trained model to local file. Default: switched off.
* `-h`, `--help`. Print help information. Unimplemented. 
* `-v`, `--version`. Print version information. Unimplemented. 


#Todo

* [ ] allow stratified sampling (e.g., STUDY)
* [ ] allow specifying must-included variables (e.g., STUDY, AGE, etc.)
* [ ] optimize memory and speed further
* [ ] mute output except `bam` file
* [ ] print data information (e.g., MAF, etc.)
* [ ] filter genotypes by MAF
* [ ] test AUC function with tie
* [ ] allow specifying cutoff in prediction
* [ ] allow individual weights in sampling
* [ ] print help information


