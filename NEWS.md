# CHANGES IN RANDOM BAMBOO VERSION 0.5.x

## MAJOR CHANGES

- before 0.5.3, a covariate is randomly selected as a candidate split, which is not a good practice in GWAS since covariates are supposed to have significant effects on the risk of disease. In 0.5.3, all covariates become the must-included candidates by default. A new option `--selcovar` is added if we want to switch back to the strategy before 0.5.3, although it is not recommended. 
- an option, `--swt` is added to allow the users to specify sample weights for each sample. This is particular useful when the training dataset consists of samples from multiple studies with different sampling design. Only the individuals specified in both `--swt` (if any) and `--trainid` (if any) are involved in training. The individuals with zero weights are discarded, too
- two experimental options, `--neighbor` and `--searchback` introduced in v0.4.x are removed, since they have no significant effect on improving power due to the multiple-comparison issue as expected; all other major/minor changes and bug fixes in 0.4.x are inherited

## BUG FIXES

- fix a minor bug in function `PredictTestingSampleFromMultipleForest()`. It could lead to incorrect ntree

# CHANGES IN Random Bamboo VERSION 0.4.x

## MAJOR CHANGES

- discarded option `--noflip` by introducing a new one `--flip` instead; in the previous version, markers with MAF > 0.5 are flipped by default; this flag is saved in the local model and prompts `Random Bamboo` to flip markers with MAF > 0.5 in testing data. This is dangerous as the MAFs are calculated from the testing data and it can be misleading if sample size is small and MAF is close to 0.5
- two experimental options, `--neighbor` and `--searchback`, are added to explored LDs between markers

## BUG FIXES

- fixed a bug that the individuals may not be aligned correctly when loading testing data; all previous versions should not be used in predicting new dataset
- `RB` allows users to specify `--ntree` in prediction step; only specified number of bamboos are used in prediction; there is a minor bug when a single bam file is used as the trained model.
- if all markers specified by `--snpid` are not included in the data, the program may crash; it is fixed by printing a message before quiting if this happens 
