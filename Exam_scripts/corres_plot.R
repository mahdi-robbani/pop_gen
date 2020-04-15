# read population labels and estimated admixture proportions
pop<-read.table("/home/xqr418/popgen/round2/data/pruned.fam")
q<-read.table("/home/xqr418/popgen/round2/admixture/filter.3.Q")


# order according to population and plot the ADMIXTURE reults
ord<-order(pop[,2])
barplot(t(q)[,ord],col=2:10,space=0,border=NA,xlab="Individuals",ylab="Demo2 Admixture proportions for K=3")
text(tapply(1:nrow(pop),pop[ord,2],mean),-0.05,unique(pop[ord,2]),xpd=T)
abline(v=cumsum(sapply(unique(pop[ord,2]),function(x){sum(pop[ord,2]==x)})),col=1,lwd=1.2)


r<-as.matrix(read.table("k3.corres"))

# Plot correlation of residuals
source("evalAdmix/NicePlotCorRes.R")
plotCorRes(cor_mat = r, pop = as.vector(pop[,2]), title="Evaluation of 1000G admixture proportions with K=3", max_z=0.1, min_z=-0.1)

