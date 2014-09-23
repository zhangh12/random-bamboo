Random Bamboo
=============


#INTRODUCTION

Random Bamboo (RB) is a ultra-fast C++ implement of the Random Forest (RF) algorithm. The program is designed for analyzing the large-scale genome-wide association study (GWAS) and building predictive model on hundreds of thousands of SNPs for personalized medicine. RB adopts fancy data structure and Boolean operation method to speedup the algorithm dramatically. The single-thread version is 40 times faster than the R package `randomForest`. OpenMP is used to enable parallelization feature. RB is thus appliable to large study with thousands or tens of thousands of individuals. 

#OPTIONS


* -f, --file. File names of three plink-formated datasets, including bed, bim and fam. Character. No default.
* -c, --cont. File name of continuous or ordered covariates, with extension con. The first row is header. The first column contains individual IDs. Character. No default.
* -a, --cate. File name of categorical covariates, with extension cat. The first row is header. The first column contains individual IDs. Character. No default.
* -o, --out. File name of all output files. Character. If unspecified, it is set as --file or --pred.
* -p, --pred. File names of datasets used in prediction, including bed, bim, fam, con and cat. The last two are optional if no covariates are used in training model. Character. No default.
* -b, --bam. Used to specify trained model. It can be the name of a single bam file or a directory contains multiple bam files. Character. No default.
* -y, --trainid. File name of individual IDs used in training model. The file contains one column and no header. With extension iid. If unspecified, all individuals in the datasets are used.
* -z, --testid. File name of individual IDs used in prediction. If --file is specified, the individuals are predicted with model trained from the specified datasets. If --bam is specified, the individuals are predicted with specified model(s). The file contains one column and no header. With extension iid. Character. If unspecified, all individuals in the datasets are used.
* -S, --snpid. File name of SNPs' names used in training model. The file contains one column and no header. With extension sid. Character. No default.
* -t, --ntree. Number of trees in the forest. Positive integer. Default: 1.
* -m, --mtry. Number of candidate variables in determining best split in each node. Positive integer. Default: sqrt of total number of variables.
* -s, --seed. Random seed. Positive integer. Default: 1.
* -l, --maxnleaf. Maximum leaves in a single tree. Positive integer. Default: 1000000.
* -e, --minleafsize. Minimum sample size in a single tree. Positive integer. Default: 1.
* -i, --imp. Types of variable importance. 1 - Gini; 2 - Breiman & Cutler's; 3 - Liaw & Wiener's; 4 - Raw; 5 - Meng's; 6 - All. Default: 1.
* -d, --nthread. Number of threads in parallelization. Positive integer. Default: maximum allowed by the hardware.
* -w, --classwt. Class weight assigned to the case group. Positive double. Default: 1.0.
* -u, --cutoff. Unimplemented.
* -g, --noflip. If specified, the genotypes of SNPs with MAF > 0.5 are flipped.
* -x, --prox. If specified, save the proximity matrix to local file.
* -n, --noimp. If specified, don't save variable importances to local file.
* -B, --balance. If specified, balance the ratio of cases and controls in training data to 1.
* -r, --trace. If specified, print debug information. Valid in single-thread mode.
* -N, --nobam. If specified, don't save trained model to local file.
* -h, --help. If specified, print help information. Unimplemented.
* -v, --version. If specified, print version information. Unimplemented.




