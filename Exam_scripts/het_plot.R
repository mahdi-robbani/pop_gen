# Read in each of the frequency files
cols <- c("CHR","SNP","A1","A2","MAF", "NCHROBS", "pos")
ceu<-read.table("~/popgen/round2/CEU/CEU.frq2",h=F)
colnames(ceu) <- cols
yri<-read.table("~/popgen/round2/YRI/YRI.frq2",h=F)
chb<-read.table("~/popgen/round2/CHB/CHB.frq2",h=F)
mxl<-read.table("~/popgen/round2/MXL/MXL.frq2",h=F)
# Function for estimating the expected heterozygosity
het<-function(x){2*x*(1-x)}
# Remove all fixed alleles in each population
yri <- yri[yri[,"MAF"]>0,]
ceu <- ceu[ceu[,"MAF"]>0,]
chb <- chb[chb[,"MAF"]>0,]
mxl <- mxl[mxl[,"MAF"]>0,]
# No obvious positions for yri and ceu plus a guestimate
# on the length of the included chromosome
yri <- cbind(yri, pi=het(yri$MAF)*(length(yri$MAF)/(2795722897)))
ceu <- cbind(ceu, pi=het(ceu$MAF)*(length(ceu$MAF)/(2795722897)))
chb <- cbind(chb, pi=het(chb$MAF)*(length(chb$MAF)/(2795722897)))
mxl <- cbind(mxl, pi=het(mxl$MAF)*(length(mxl$MAF)/(2795722897)))
# Making a barplot with the nucleotide diversity
val <- c(mean(ceu$pi), mean(yri$pi), mean(chb$pi), mean(mxl$pi))
val2 <- c(sd(ceu$pi), sd(yri$pi), sd(chb$pi), sd(mxl$pi))
par(mfrow=c(1,1))
pdf("Het_plot.pdf")
barplot(val, ylim = c(0, 0.0001), ylab="pi", xlab="Population", names.arg=c("CEU","YRI","CHB","MXL"))
# To shut down the plot window
dev.off()
#Make table of values for the plot
names <- c("CEU","YRI","CHB","MXL")
df <- data.frame("Group" = names, "Mean.Heterozygosity" = val, "SD.Heterozygosity" = val2)
write.table(df, "Het.txt", quote = F, row.names = F)



