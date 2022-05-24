# Plotting tvectors from MCMC

args=commandArgs(trailingOnly=TRUE)

theTable <- args[1]
thePrefix <- args[2]

myt <- read.table(theTable, header=T, stringsAsFactors=F, colClasses=c(rep('numeric', 10), 'character'))

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

barplot_tvectors_MCMC <- function(x, prefix){
    colvec <- c('#1f78b4', '#ff7f00')
    d <- matrix(c(unlist(strsplit(x, split=''))), nrow=length(x), ncol=nchar(x[1]), byrow=T)
    o <- data.frame(matrix(ncol=2, nrow=ncol(d)))
    for (i in 1:ncol(d)){
        o[i,] <- table(d[,i])/sum(table(d[,1]))
    }
    theName <- paste(c(prefix, '_barplot_tvector_MCMC.pdf'), sep='', collapse='')
    pdf(theName, width=ncol(d), height=5)
    barplot(t(o), col=colvec, beside=T, ylim=c(0,1.2))
    legend('topright', fill=colvec, legend=c('1 -> 2', '2 -> 1'))
    box()
    dev.off()
}

histogram_parameter_MCMC <- function(x, prefix){
    subtable <- x[x$state>0,]
    theName <- paste(c(prefix, '_histogram_parameter_distribution.pdf'), sep='', collapse='')
    # pdf('histogram_parameter_distribution.pdf', width=10, height=10)
    pdf(theName, width=12, height=12)
    par(mfrow=c(3,3), font=2)
    for (i in 3:(ncol(subtable)-1)){
        maintitle <- paste(c('histogram of ', colnames(subtable)[i]), sep='', collapse='')
        myhist.lim <- max(hist(subtable[,i], breaks=100, plot=F)$density)
        hist(subtable[,i], breaks=100, freq=F, main=maintitle, xlab=colnames(subtable)[i], col='#56B4E9', ylim=c(0,myhist.lim+myhist.lim/4))
        lines(density(subtable[,i]), col='#D55E00', lwd=3)
        box()
        print(i)
        print(estimate_mode(subtable[,i]))
        legend('topright', legend=paste(c("mode=" ,format(estimate_mode(subtable[,i]), scientific=F, digits=4)), sep='', collapse=''), bty='n')
    }
    myhist.lim <- max(hist(subtable[,'log.posterior'], breaks=100, plot=F)$density)
    hist(subtable[,'log.posterior'], breaks=100, freq=F, main='histogram of the log-posterior', xlab='log posterior', col='#56B4E9', ylim=c(0,myhist.lim+myhist.lim/4))
    lines(density(subtable[,'log.posterior']), col='#D55E00', lwd=3)
    box()
    legend('topright', legend=paste(c("mode=" ,format(estimate_mode(subtable[,'log.posterior']), scientific=F, digits=4)), sep='', collapse=''), bty='n')
    dev.off()
}

if( (!(is.na(myt$tvector[1]))) & (nchar(myt$tvector[1])>1)){barplot_tvectors_MCMC(myt$tvector, thePrefix)}

histogram_parameter_MCMC(myt, thePrefix)

print("#--- Mode of the distributions ---#")
out.df <- c(apply(myt[myt$state>0,seq(2,10)], MARGIN=2, FUN=estimate_mode))
print(out.df)
