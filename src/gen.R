

# This script will create several files required by the random bamboo:
# (1) *.ped, *.map: These two files will be used with the program PLINK to generate files *.fam, *.bim and *.bed. The latter three files are used as inputs of random bamboo
# (2) *.con, *.cat: These two files will be used as (optional) inputs of random bamboo. They contain continuouse and categorical covariates
# (3) The files train.* will be used in training, while test.* will be used in testing

set.seed(1)

n<-2000 # total sample size
p<-15000 # number of SNPs in testing dataset
po<-10000 # number of SNPs in training dataset

# You can use smaller values for the three parameters above for convenience
# You don't need to modify others below

#sel<-sort(sample(1:p, po, FALSE))
sel<-1:po
sel.col<-sort(c(2*sel-1, 2*sel))

maf<-runif(p,.4,.4)
pr<-rep(rep(maf,each=2), times=n)

hap<-matrix(rbinom(2*n*p,1,pr),nrow=n,byrow=T)+1

#h<-apply(hap,2,mean)

fid<-paste0("FAM", 1:n)
iid<-paste0("IID", 1:n)
pid<-paste0("PID", 1:n)
mid<-paste0("MID", 1:n)

fid<-1:n
iid<-1:n
pid<-1:n
mid<-1:n

sex<-rep(1,n)
x1<-rnorm(n) #continuous var with effect
x2<-sample(1:5, n, TRUE) #ordered var with effect
x3<-sample(c("aa", "bb_b", "cc#cc", "dddd"), n, TRUE) #categorical var with effect
x4<-rnorm(n) #continuous var without effect
x5<-sample(1:6, n, TRUE) #ordered var without effect
x6<-sample(c("a1a", "b2b", "c3", "dd4"), n, TRUE) #categorical var without effect
x7<-rnorm(n) #continuous var with effect but not observed
x8<-runif(n)
x9<-runif(n)
x10<-sample(1:10, n, TRUE)
x11<-sample(c("a1","a2","a3","a4","a5","a6"), n, TRUE)

x12<-sample(11:20, n, TRUE)
x13<-runif(n,-1,1)
x14<-sample(c("haha", "heyhey", "xixi", "hoho"), n, TRUE)


hap<-data.frame(hap)
colnames(hap)<-paste0("hap", 1:ncol(hap))

p0<-min(p,200)
h1<-hap[,2*(1:p0)-1]-1
h2<-hap[,2*(1:p0)]-1
geno<-h1+h2


#bet<-seq(5,10,length.out=p0)
bet <- runif(p0,-.4,.4)
bet <- bet[order(abs(bet))]
eff<-as.matrix(geno) %*% bet + 1*x1+.5*x2+0.5*(x3=="aa")+1*(x3=="cc#cc")+3*x7
#eff<-as.matrix(geno) %*% bet
eff<-log(1/2/(1-1/2))+eff-mean(eff)
pr<-exp(eff)/(1+exp(eff))
phe<-rbinom(n,1,pr)+1
table(phe)


PED.test<-data.frame(fid,iid,pid,mid,sex,phe,hap)
PED.train<-data.frame(fid,iid,pid,mid,sex,phe,hap[, sel.col])
h1<-hap[,2*(1:p)-1]-1
h2<-hap[,2*(1:p)]-1
geno<-h1+h2
rm(hap)
for(i in 1:100){gc()}

PED.train<-as.matrix(PED.train)
PED.test<-as.matrix(PED.test)

write.table(PED.train, "train.ped",col.names=F,row.names=F,quote=F,append=F)
write.table(PED.test, "test.ped",col.names=F,row.names=F,quote=F,append=F)


chr<-rep(1,p)
snp<-paste0("rs",1:p)
dist<-rep(0,p)
pos<-1:p
MAP.test<-data.frame(chr,snp,dist,pos)
MAP.train<-MAP.test[sel,]
write.table(MAP.train, "train.map", col.names=F,row.names=F,quote=F,append=F)
write.table(MAP.test, "test.map", col.names=F,row.names=F,quote=F,append=F)

#write.table(data.frame(x1,x2,x4,x5,x8,x9), "test.con", col.names=F,row.names=F,quote=F,append=F,sep="\t")
#write.table(data.frame(x3,x6,x10,x11), "test.cat", col.names=F,row.names=F,quote=F,append=F,sep="\t")


con.train<-data.frame(iid,x1,x2,x4,x5,x8,x9)
cat.train<-data.frame(iid,x3,x6,x10,x11)

write.table(con.train, "train.con", col.names=T,row.names=F,quote=F,append=F,sep="\t")
write.table(cat.train, "train.cat", col.names=T,row.names=F,quote=F,append=F,sep="\t")

con.test<-data.frame(iid,x2, x12, x9, x4, x5, x13, x1, x8)
cat.test<-data.frame(iid, x6, x10, x3, x11, x14)
write.table(con.test, "test.con", col.names=T,row.names=F,quote=F,append=F,sep="\t")
write.table(cat.test, "test.cat", col.names=T,row.names=F,quote=F,append=F,sep="\t")

