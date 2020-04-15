#load data

chr_plot <- function(c){
  ceu<-read.table(paste0("~/popgen/round2/CEU/CEU.frq.", toString(c)),h=F)
  yri<-read.table(paste0("~/popgen/round2/YRI/YRI.frq.", toString(c)),h=F)
  chb<-read.table(paste0("~/popgen/round2/CHB/CHB.frq.", toString(c)),h=F)
  mxl<-read.table(paste0("~/popgen/round2/MXL/MXL.frq.", toString(c)),h=F)
  #change col names
  cols <- c("CHR","SNP","A1","A2","MAF", "NCHROBS", "position")
  colnames(ceu) <- cols
  colnames(yri) <- cols
  colnames(chb) <- cols
  colnames(mxl) <- cols
  #remove fixed alleles
  yri <- yri[yri[,"MAF"]>0,]
  ceu <- ceu[ceu[,"MAF"]>0,]
  chb <- chb[chb[,"MAF"]>0,]
  mxl <- mxl[mxl[,"MAF"]>0,]
  # Function for estimating the expected heterozygosity
  het<-function(x){2*x*(1-x)}
  # No obvious positions for yri and ceu plus a guestimate
  # on the length of the included chromosome
  yri <- cbind(yri, pi=het(yri$MAF)*(length(yri$MAF)/(2795722897)))
  ceu <- cbind(ceu, pi=het(ceu$MAF)*(length(ceu$MAF)/(2795722897)))
  chb <- cbind(chb, pi=het(chb$MAF)*(length(chb$MAF)/(2795722897)))
  mxl <- cbind(mxl, pi=het(mxl$MAF)*(length(mxl$MAF)/(2795722897)))
  
  ## Function for generating sliding windows
  slidingwindowplot <- function(mainv, xlabv, ylabv, ylimv, window.size,
                                step.size,input_x_data,input_y_data)
  {
    if (window.size > step.size)
      step.positions <- seq(window.size/2 + 1, length(input_x_data)- window.size/2,
                            by=step.size)
    else
      step.positions <- seq(step.size/2 + 1, length(input_x_data)- step.size,
                            by=step.size)
    n <- length(step.positions)
    means_x <- numeric(n)
    means_y <- numeric(n)
    for (i in 1:n) {
      chunk_x <- input_x_data[(step.positions[i]-
                                 window.size/2):(step.positions[i]+window.size-1)]
      means_x[i] <- mean(chunk_x,na.rem=TRUE)
      chunk_y <- input_y_data[(step.positions[i]-
                                 window.size/2):(step.positions[i]+window.size-1)]
      means_y[i] <- mean(chunk_y,na.rem=TRUE)
    }
    plot(means_x,means_y,type="b",main=mainv,xlab=xlabv,ylab=ylabv,ylim=ylimv,cex=0.25,
         pch=20,cex.main=0.75)
    vec <- c(0.025,0.5,0.975)
    zz <- means_y[!is.na(means_y)]
    abline(h=quantile(zz,0.025,na.rem=TRUE),col="blue")
    abline(h=quantile(zz,0.925,na.rem=TRUE),col="blue")
    abline(h=mean(input_y_data))
  }
  
  ## R is doing strange things on the graphics window; therefore, we plot it on
  ## a pdf file. You can view it with evince afterwards
  
  pdf(paste0("nucleotide_diversity_in_chr",toString(c), ".pdf"))
  par(mfrow=c(2,2))
  windowsize<- 3000
  steps<- 100
  start <- 0
  end <- 0.00001
  # CEU
  mainvv = paste("CEU pi = ",format(mean(ceu$pi,na.rem=TRUE), digits=3), "SNPs =",
                 length(ceu$pi), "Win: ", windowsize, "Step: ", steps)
  slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
                    ylab=expression(paste("pi")),ylimv=c(start,end), window.size=windowsize/4,
                    step.size=steps, input_x_data=ceu$position/1000000,input_y_data=ceu$pi)
  # yri
  mainvv = paste("YRI pi = ",format(mean(yri$pi,na.rem=TRUE), digits=3), "SNPs =",
                 length(yri$pi), "Win: ", windowsize, "Step: ", steps)
  slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
                    ylab=expression(paste("pi")),ylimv=c(start,end), window.size=windowsize/4,
                    step.size=steps, input_x_data=yri$position/1000000,input_y_data=yri$pi)
  # chb
  mainvv = paste("CHB pi = ",format(mean(chb$pi,na.rem=TRUE), digits=3), "SNPs =",
                 length(chb$pi), "Win: ", windowsize, "Step: ", steps)
  slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
                    ylab=expression(paste("pi")),ylimv=c(start,end), window.size=windowsize/4,
                    step.size=steps, input_x_data=chb$position/1000000,input_y_data=chb$pi)
  # mxl
  mainvv = paste("MXL pi = ",format(mean(mxl$pi,na.rem=TRUE), digits=3), "SNPs =",
                 length(mxl$pi), "Win: ", windowsize, "Step: ", steps)
  slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
                    ylab=expression(paste("pi")),ylimv=c(start,end), window.size=windowsize/4,
                    step.size=steps, input_x_data=mxl$position/1000000,input_y_data=mxl$pi)
  dev.off()
}

chr_plot(1)