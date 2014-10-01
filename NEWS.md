# CHANGES IN Random Bamboo VERSION 0.4.0

## MAJOR CHANGES

- discarded option `--noflip` by introducing a new one `--flip` instead; in the previous version, markers with MAF > 0.5 are flipped by default; this flag is saved in the local model and prompts `Random Bamboo` to flip markers with MAF > 0.5 in testing data. This is dangerous as the MAFs are calculated from the testing data and it can be misleading if sample size is small and MAF is close to 0.5

## BUG FIXES

- fixed a bug that the individuals may not be aligned correctly when loading testing data; all previous versions should not be used in predicting new dataset
