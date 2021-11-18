# Plot utilities for epics/epocs
# I'm storing every possible plots for general use with epics/epocs.

args=commandArgs(trailingOnly=TRUE)

# 1. QQPLOT

qq_plot <- function(pval.vec, prefix){
    exp.vector <- (rank(pval.vec, ties.method='first')) / (length(pval.vec))
    # pdf(paste(c(prefix, '_qqplot.pdf'), sep='', collapse=''), width=5, height=5)
    png(paste(c(prefix, '_qqplot.png'), sep='', collapse=''), width=480, height=480)
    plot(-log10(exp.vector), -log10(pval.vec), pch=16, ylab='observed -log10(pval)', xlab='expected -log10(pval)')
    abline(0,1, lwd=2, col='#ef3b2c')
    dev.off()
    return
}

# 2. HEATMAP

heatmap_plot <- function(epocstable, prefix){
    library(gplots) # Best heatmap maker around...
    color_ia_vec <- colorRampPalette(c('grey100', 'navyblue'))
    color_ia <- color_ia_vec(100)

    color_ib_vec <- colorRampPalette(c('khaki1', '#d7191c'))
    color_ib <- color_ib_vec(100)

    color_vec <- c(color_ia, color_ib)

    # ia = 1->2, i.e., row -> column
    # ib = 2->1, i.e., column -> row

    # I'm showing only the top 50 according to the pvalue - to clarify the heatmap
    # Or I take everything if I have less than 50 pairs
    epocstable <- epocstable[order(-epocstable$logPVAL),][1:min(50, nrow(epocstable)),]

    condition=FALSE # If I have the event IDs in the table, will become True
    if ('name_e1' %in% colnames(epocstable)){
        uniqueNames <- c(unique(c(epocstable$name_e1, epocstable$name_e2)))
        condition=TRUE
    }
    else{uniqueNames <- c(unique(c(epocstable$Event1, epocstable$Event2)))}

    # Preparing the heatmap matrix
    ranks <- seq(1,length(uniqueNames))
    names(ranks) <- uniqueNames
    newtable <- matrix(rep(0, length(ranks)*length(ranks)), ncol=length(ranks), nrow=length(ranks))
    colnames(newtable) <- names(ranks)
    rownames(newtable) <- names(ranks)

    # And now filling the heatmap matrix
    for (i in 1:nrow(epocstable)){
        if (condition){
            pos1 <- as.numeric(ranks[names(ranks)==epocstable[i,]$name_e1])
            pos2 <- as.numeric(ranks[names(ranks)==epocstable[i,]$name_e2])
        }
        else{
            pos1 <- as.numeric(ranks[names(ranks)==epocstable[i,]$Event1])
            pos2 <- as.numeric(ranks[names(ranks)==epocstable[i,]$Event2])
        }

        if (epocstable[i,]$combinedModels == 'ib'){
            # E2 --> E1
            theValue <- (epocstable[i,]$lambda1)/10         # Normalize to 100 for display purposes
            if (theValue==Inf){theValue=100}                # Maximum lambda value
            theValue <- as.integer(theValue)
            if (theValue >= 1){theValue <- theValue+100}    # If I have a positive value, I'm drawing the color in the vector after 100 (such that I can show the direction in the heatmap)
        }

        else if (epocstable[i,]$combinedModels == 'ia'){
            # E1 --> E2
            theValue <- (epocstable[i,]$lambda2)/10
            if (theValue==Inf){theValue=100}
            theValue <- as.integer(theValue)
        }

        else{
            # E1 <--> E2
            # Logic: I'm retrieving here the highest lambda between the two!
            if ((epocstable[i,]$lambda2)/10 > (epocstable[i,]$lambda1)/10){
                theValue <- (epocstable[i,]$lambda2)/10
                if (theValue==Inf){theValue=100}
                theValue <- as.integer(theValue)
            }else{
                theValue <- (epocstable[i,]$lambda1)/10
                if (theValue==Inf){theValue=100}
                theValue <- as.integer(theValue)
                if (theValue >=1 ){theValue <- theValue+100}
            }
        }

        newtable[pos1,pos2] <- theValue

    }

    indRowZero <- as.numeric(which(apply(newtable, 1, FUN=sum) == 0))
    indColZero <- as.numeric(which(apply(newtable, 2, FUN=sum) == 0))

    # print(newtable)

    if (length(indRowZero)>0){
        if (is.null(nrow(newtable))){newtable <- newtable[-indRowZero]}
        else{newtable <- newtable[-indRowZero,]}
    }

    if (length(indColZero)>0){
        if (is.null(ncol(newtable))){newtable <- newtable[-indColZero]}
        else{newtable <- newtable[,-indColZero]}
    }

    if ((is.null(nrow(newtable))) || (is.null(ncol(newtable)))){
        stop("not possible to plot heatmaps! not enough data to show... Exiting.")
    }
    pdf(paste(c(prefix, 'HEATMAP.pdf'), sep='', collapse=''), width=15, height=15)
    myheatmap <- heatmap.2(newtable, col=color_vec, na.rm=T, breaks=seq(0,200), trace='none', key=T, colsep=1:ncol(newtable), rowsep=1:nrow(newtable), sepcolor='ivory2', sepwidth=c(0.01, 0.01), revC=F, margins=c(10,10), keysize=1, density.info='none', key.xtickfun=function(){
        return(list(at=c(0,.25, .75, 1), labels=c('', '\nstrength\ninduction\nrow->column', '\nstrength\ninduction\ncol->row', '')))
    }, key.xlab=NA, key.par=list(mar=c(12, 4, 4, 1) +0.1, mgp=c(3,3,0)))
    dev.off()
    # return
}

# 3. DOT FORMAT - WIP

convert_dot <- function(epocstable, prefix){
    return
}

# 4. CYTOSCAPE FORMAT - NOT NECESSARY RIGHT NOW

cytoscape_convert <- function(epocstable, prefix){
    return
}

# 5. MANHATTAN - WIP

manhattan_plot <- function(fulltable, prefix){
    return
}

# --- MAIN LOOP ---

if (length(args)!=2){stop('Rscript PlotUtilities.R [FolderPrefix] [kindOfPlot]
				 [FolderPrefix]: name of folder with results of epocs
				 [kindOfPlot]: QQ, HEATMAP')}

myFolder <- args[1]
thePlot <- args[2]

condition <- FALSE

if (grepl('/',myFolder))
{
    myFolder <- gsub('/', '', myFolder)
}

# If the user wants to run in current directory: need to adjust a bit the variable to account for it
if (myFolder == '.')
{
    myFolder=getwd()
    myOutputPrefix=tail(strsplit(myFolder, split='/')[[1]], n=1)
    # myPrefix=paste(c(myOutputPrefix, '_mat_epocs'), sep='', collapse='')
    myPrefix = myOutputPrefix
    condition=TRUE
}else{
    # myPrefix <- paste(c(myFolder, '_mat_epocs'), sep='', collapse='')
    myPrefix <- myFolder
}

# all output files from epocs - one per model!
allFiles <- list.files(myFolder, pattern=myPrefix)
if (length(allFiles)<2){stop('Error reading epocs files - not enough files present! Make sure you chose the correct folder')}

if(!condition){setwd(myPrefix)}

if (thePlot == 'QQ'){
    # myFileOfInterest <- allFiles[which(grepl('mat_epics_I', allFiles))]
    # myt <- read.table(myFileOfInterest, header=T, stringsAsFactors=F, sep='\t', check.names=F)
    # qq_plot(myt$Pval, 'EpicsMatI')
    # myFileOfInterest <- allFiles[which(grepl('mat_epics_S', allFiles))]
    # myt <- read.table(myFileOfInterest, header=T, stringsAsFactors=F, sep='\t', check.names=F)
    # qq_plot(myt$Pval, 'EpicsMatS')
    myFileOfInterest <- allFiles[which(grepl('LRT.tab', allFiles))]
    myt <- read.table(myFileOfInterest, header=T, stringsAsFactors=F, sep='\t', check.names=F)
    qq_plot(myt$Pvalue, 'EpocsBestModel')
}

if (thePlot == 'HEATMAP'){
    myFileOfInterest <- allFiles[which(grepl('BestModel.tab', allFiles))]
    myt <- read.table(myFileOfInterest, header=T, stringsAsFactors=F, sep='\t', check.names=F)
    heatmap_plot(myt, 'EpocsBestModel')
}

# if (thePlot == 'CYTOSCAPE'){
#     myFileOfInterest <- allFiles[which(grepl('BestModel.tab', allFiles))]
#     myt <- read.table(myFileOfInterest, header=T, stringsAsFactors=F, sep='\t', check.names=F)
#     cytoscape_convert(myt, 'EpocsBestModel')
# }
