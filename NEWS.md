# CHANGES IN Random Bamboo VERSION 0.5.x

## MAJOR CHANGES

- two experimental options, `--neighbor` and `--searchback` introduced in v0.4.x are removed, since they have no significant effect on improving power due to the multiple-comparison issue as expected; all other major/minor changes and bug fixes in v0.4.x are inherited


# CHANGES IN Random Bamboo VERSION 0.4.x

## MAJOR CHANGES

- discarded option `--noflip` by introducing a new one `--flip` instead; in the previous version, markers with MAF > 0.5 are flipped by default; this flag is saved in the local model and prompts `Random Bamboo` to flip markers with MAF > 0.5 in testing data. This is dangerous as the MAFs are calculated from the testing data and it can be misleading if sample size is small and MAF is close to 0.5
- two experimental options, `--neighbor` and `--searchback`, are added to explored LDs between markers

## BUG FIXES

- fixed a bug that the individuals may not be aligned correctly when loading testing data; all previous versions should not be used in predicting new dataset
- `RB` allows users to specify `--ntree` in prediction step; only specified number of bamboos are used in prediction; there is a minor bug when a single bam file is used as the trained model.
- if all markers specified by `--snpid` are not included in the data, the program may crash; it is fixed by printing a message before quiting if this happens 
