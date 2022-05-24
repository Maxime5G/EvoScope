# Plot a 2D summary for epocs with either 4 or 8 models tested

args=commandArgs(trailingOnly=TRUE)

plotFromNormalizedTable <- function(myd, theSubtable.LRT, name, vectorOfYLabels, vectorOfXLabels){

    par(mgp=c(3,1.5,0))
    plot(x=seq(1,nrow(myd)+1), xlim=c(1,nrow(myd)+1), ylim=c(0,1), type="n", xaxt="n", yaxt="n", main=name, xlab="model", ylab="LRT values")
    # axis(side=1, at=seq(1,nrow(myd)+1), label=c("i",myd$m), family="serif")
    axis(side=1, at=seq(1,nrow(myd)+1), label=vectorOfXLabels, family="serif", cex.lab=.7)
    axis(side=2, at=c(0,vectorOfYLabels$n), label=c(0,vectorOfYLabels$v), las=2, family="serif")
    for (i in 1:nrow(myd)){
        # points(i+1, myd$normalized[i], lwd=2, pch=19, col=myCol[myd$normalized[i]*10000], cex=3)
        points(i+1, myd$normalized[i], lwd=2.5, pch=21, col="black", bg=myCol[myd$normalized[i]*10000], cex=3)
    }

  # Next plotting the segments

    for (i in 1:nrow(theSubtable.LRT)){
        theL <- theSubtable.LRT[i,]

        if(theL$Pvalue < 0.05){
            lwdForLine <- (length(unlist(gregexpr("\\*", theSubtable.LRT[1,]$LRT2)))+1)/2
            if (theL$Model1  ==  "i"){
                segments(x0=1,y0=0,x1=match(theL$Model2, myd$m)+1, y1=myd[which(myd$m == theL$Model2),]$normalized, lwd=lwdForLine)
            }else{
                segments(x0=match(theL$Model1, myd$m)+1, y0=myd[which(myd$m == theL$Model1),]$normalized, x1=match(theL$Model2, myd$m)+1, y1=myd[which(myd$m == theL$Model2),]$normalized, lwd=lwdForLine)
            }
        }

        else{
            if (theL$Model1  ==  "i"){
                segments(x0=1,y0=0,x1=match(theL$Model2, myd$m)+1, y1=myd[which(myd$m == theL$Model2),]$normalized, lty=2)
            }else{
                segments(x0=match(theL$Model1, myd$m)+1, y0=myd[which(myd$m == theL$Model1),]$normalized, x1=match(theL$Model2, myd$m)+1, y1=myd[which(myd$m == theL$Model2),]$normalized, lty=2)
            }
        }
    }
}

normalizeAllSubtables <- function(theTable){

    orderedTable <- theTable[order(theTable$Pvalue),]
    E1 <- c()
    E2 <- c()
    pairs <- c()
    theNormalizedTable <- NULL
    for (i in 1:nrow(orderedTable)){
        if ((length(pairs)<9) & (orderedTable[i,]$Pvalue < 0.05)){
            currentPair <- paste(c(orderedTable[i,]$Event1, orderedTable[i,]$Event2), sep="", collapse="_")

            if (!(currentPair %in% pairs)){
                pairs <- c(pairs, currentPair)
                presentTable <- orderedTable[which(orderedTable$Event1  ==  orderedTable[i,]$Event1 & orderedTable$Event2  ==  orderedTable[i,]$Event2),]
                posInPlot <- list()
                for (m in allowedModels){

                    if ((match(m, allowedModels)<6) && m!="i")
                    {
                        theCombinedModels <- paste(c("i", m), sep="", collapse="")
                        if(theCombinedModels %in% presentTable$combinedModels) {posInPlot[m] <- presentTable[which(presentTable$combinedModels  ==  theCombinedModels),]$LRT}
                    }

                    if (match(m, allowedModels)>5)
                    {
                        theCombinedModels <- paste(c("I", m), sep="", collapse="")
                        if(theCombinedModels %in% presentTable$combinedModels) {posInPlot[m] <- presentTable[which(presentTable$combinedModels  ==  theCombinedModels),]$LRT + posInPlot[["I"]]}
                    }
                }

                myTable <- data.frame(n=currentPair, v=unlist(posInPlot), m=names(unlist(posInPlot)), row.names=NULL)
                if(is.null(theNormalizedTable)){theNormalizedTable <- myTable}
                else{
                    theNormalizedTable <- rbind(theNormalizedTable, myTable)
                }
                # print(myTable)
                # stop()
            }
        }
    }
    theNormalizedTable$normalized <- theNormalizedTable$v/max(theNormalizedTable$v)
    return (theNormalizedTable)
}

if (length(args)>5){stop("Rscript PlotIn2D.R [Prefix] [LRT] [ConcatenatedFile] [WhichMode (BP or e1-e2 or all)] [eventIDs]")}

if (length(args)==4){eventIDs=NULL}
if (length(args)==5){
    eventIDs=args[5]
    evIDs <- read.table(eventIDs, header=F, stringsAsFactors=F, colClasses='character')
}

allowedModels <- c("i", "a", "b", "l", "I", "A", "B", "L")

options(scipen=999)       # Disabling scientific notations
options(warn=0)           # Original setting for warnings

myPrefix <- args[1]
myEpocsOutput.LRT <- args[2]
myEpocsOutput.Conc <- args[3]
myChoice <- args[4]

fl.m <- 0
fl.c <- 0

myTable.LRT <- read.table(myEpocsOutput.LRT, header=T, stringsAsFactors=F, sep="\t", check.names=F)
myTable.Conc <- read.table(myEpocsOutput.Conc, header=T, stringsAsFactors=F, sep="\t", check.names=F)
myTable.LRT$LRT2 <- rep(NA, nrow(myTable.LRT))


# some sanity checks
if (length(which(!(myTable.LRT$Model1 %in% allowedModels)))>0){stop()}
if (length(which(!(myTable.LRT$Model2 %in% allowedModels)))>0){stop()}

allLRT <- myTable.LRT$LRT

modelsInTable <- unique(c(which(allowedModels %in% myTable.LRT$Model2), which(allowedModels %in% myTable.LRT$Model1)))
myModels <- allowedModels[modelsInTable]
uppercaseModels <- unlist(gregexpr("[A-Z]", myModels))
lowercaseModels <- unlist(gregexpr("[a-z]", myModels))

if (myChoice  ==  "BP"){
    minPval <- myTable.LRT[myTable.LRT$Pvalue  ==  min(myTable.LRT$Pvalue),]
    E1 <- minPval[1,]$Event1
    E2 <- minPval[1,]$Event2
    theSubtable.LRT <- myTable.LRT[which(myTable.LRT$Event1  ==  E1 & myTable.LRT$Event2  ==  E2),]
    theSubtable.Conc <- myTable.Conc[which(myTable.Conc$Event1  ==  E1 & myTable.Conc$Event2  ==  E2),]
    fl.c <- 1
    for (i in 1:nrow(theSubtable.LRT)){
        theLRT <- format(as.numeric(theSubtable.LRT[i,]$LRT), digits=4)
        thePVAL <- theSubtable.LRT[i,]$Pvalue
        if (thePVAL < 0.05){
            theLRT <- paste(c(theLRT, "*"), sep="", collapse="")
        }
        if (thePVAL < 0.005){
            theLRT <- paste(c(theLRT, "*"), sep="", collapse="")
        }
        if (thePVAL < 0.0005){
            theLRT <- paste(c(theLRT, "*"), sep="", collapse="")
        }
        if (theLRT == 0){theLRT=0}
        theSubtable.LRT[i,]$LRT2 <- theLRT
    }
}

if (grepl("-", myChoice)){
    sublist <- as.numeric(strsplit(myChoice, split="-")[[1]])
    if (length(sublist) != 2){stop("please input only one pair please, exiting...")}
    if (length(which(is.na(sublist)))>0){stop("error in your event choice, exiting...")}

    sublist <- sort(sublist)
    theSubtable.LRT <- myTable.LRT[which(myTable.LRT$Event1  ==  sublist[1] & myTable.LRT$Event2  ==  sublist[2]),]
    if (nrow(theSubtable.LRT)<1){stop("error in your event choice, exiting...")}

    E1 <- theSubtable.LRT[1,]$Event1
    E2 <- theSubtable.LRT[1,]$Event2
    fl.c=1

    for (i in 1:nrow(theSubtable.LRT)){
        theLRT <- format(as.numeric(theSubtable.LRT[i,]$LRT), digits=4)
        thePVAL <- theSubtable.LRT[i,]$Pvalue
        if (thePVAL < 0.05){
            theLRT <- paste(c(theLRT, "*"), sep="", collapse="")
        }
        if (thePVAL < 0.005){
            theLRT <- paste(c(theLRT, "*"), sep="", collapse="")
        }
        if (thePVAL < 0.0005){
            theLRT <- paste(c(theLRT, "*"), sep="", collapse="")
        }
        if (theLRT == 0){theLRT=0}
        theSubtable.LRT[i,]$LRT2 <- theLRT
    }
}
if (myChoice  ==  "all"){fl.c=1}

if (fl.c  ==  0){stop("Error occurred - check your inputs!\n")}

colors <- colorRampPalette(c("blue", "red"))
myCol <- colors(10000)

listForXLabels <- list("i"="H0", "a"="E1 -> E2", "b"="E1 <- E2", "l"="E1 <-> E2","I"="H0\n+state_dep", "A"="E1 -> E2\n+state_dep", "B"="E1 <- E2\n+state_dep", "L"="E1 <-> E2\n+state_dep")

if ((length(which(uppercaseModels == 1)) == 0) && (length(which(lowercaseModels == 1))>0)){

    if (myChoice  ==  "BP" | grepl("-", myChoice)) {

        myNormTable <- normalizeAllSubtables(theSubtable.LRT)
        TheDivisionOfValues <- ceiling(max(myNormTable$v)/5)
        TheValues <- c(seq(from=TheDivisionOfValues, to=TheDivisionOfValues*4, by=TheDivisionOfValues), floor(max(myNormTable$v)))

        vectorOfYLabels <- data.frame(v=TheValues, n=c(TheValues[1]/max(myNormTable$v), TheValues[2]/max(myNormTable$v), TheValues[3]/max(myNormTable$v), TheValues[4]/max(myNormTable$v), TheValues[5]/max(myNormTable$v)))

        vectorOfXLabels <- c("H0")
        for (k in 1:nrow(myNormTable)){
            vectorOfXLabels <- c(vectorOfXLabels, listForXLabels[[myNormTable[k,]$m]])
        }

        main.name <- paste(c("event ", E1, " vs event ", E2), sep="", collapse="")
        filename <- paste(c(myPrefix, "_", paste(c(E1,"vs",E2), sep="", collapse="_"), ".pdf"),sep="", collapse="")

        pdf(filename, width=8, height=8)

        plotFromNormalizedTable(myNormTable, theSubtable.LRT, main.name, vectorOfYLabels, vectorOfXLabels)

        dev.off()
    }else if(myChoice == "all"){

        print("Extracting at most 9 pairs to plot! Decision is based on lowest p-values for each pair.")
        myNormTable <- normalizeAllSubtables(myTable.LRT)
        uniquePairs <- unique(myNormTable$n)
        TheDivisionOfValues <- ceiling(max(myNormTable$v)/5)
        TheValues <- c(seq(from=TheDivisionOfValues, to=TheDivisionOfValues*4, by=TheDivisionOfValues), floor(max(myNormTable$v)))

        vectorOfYLabels <- data.frame(v=TheValues, n=c(TheValues[1]/max(myNormTable$v), TheValues[2]/max(myNormTable$v), TheValues[3]/max(myNormTable$v), TheValues[4]/max(myNormTable$v), TheValues[5]/max(myNormTable$v)))

        filename <- paste(c(myPrefix, "_", "BestNinePairsInDataset.pdf"),sep="", collapse="")

        pdf(filename, width=24, height=24)
        par(mfrow=c(3,3))
        for (i in 1:length(uniquePairs)){
            spl.pairs <- strsplit(uniquePairs[i], split="_")[[1]]
            n1 <- spl.pairs[1]
            n2 <- spl.pairs[2]
            if (!(is.null(eventIDs))){
                n1 <- evIDs[n1,]
                n2 <- evIDs[n2,]
            }
            main.name <- paste(c("event ", n1, " vs event ", n2), sep="", collapse="")
            sub.normTable <- myNormTable[myNormTable$n  ==  uniquePairs[i],]
            sub.lrtTable <- myTable.LRT[which(myTable.LRT$Event1  ==  spl.pairs[1] & myTable.LRT$Event2  ==  spl.pairs[2]),]
            vectorOfXLabels <- c("H0")
            for (k in 1:nrow(sub.normTable)){
                vectorOfXLabels <- c(vectorOfXLabels, listForXLabels[[sub.normTable[k,]$m]])
            }
            plotFromNormalizedTable(sub.normTable, sub.lrtTable, main.name, vectorOfYLabels, vectorOfXLabels)
        }
        dev.off()
    }
    fl.m <- 1
}

if ((length(which(uppercaseModels == 1))>0) && (length(which(lowercaseModels == 1))>0)) {

    if (myChoice  ==  "BP" | grepl("-", myChoice)) {

        myNormTable <- normalizeAllSubtables(theSubtable.LRT)

        TheDivisionOfValues <- ceiling(max(myNormTable$v)/5)
        TheValues <- c(seq(from=TheDivisionOfValues, to=TheDivisionOfValues*4, by=TheDivisionOfValues), floor(max(myNormTable$v)))

        vectorOfYLabels <- data.frame(v=TheValues, n=c(TheValues[1]/max(myNormTable$v), TheValues[2]/max(myNormTable$v), TheValues[3]/max(myNormTable$v), TheValues[4]/max(myNormTable$v), TheValues[5]/max(myNormTable$v)))

        vectorOfXLabels <- c("H0")
        for (k in 1:nrow(myNormTable)){
            vectorOfXLabels <- c(vectorOfXLabels, listForXLabels[[myNormTable[k,]$m]])
        }

        main.name <- paste(c("event ", E1, " vs event ", E2), sep="", collapse="")
        filename <- paste(c(myPrefix, "_", paste(c(E1,"vs",E2), sep="", collapse="_"), ".pdf"),sep="", collapse="")

        pdf(filename, width=8, height=8)

        plotFromNormalizedTable(myNormTable, theSubtable.LRT, main.name, vectorOfYLabels, vectorOfXLabels)

        dev.off()

    }else if(myChoice == "all"){

        print("Extracting at most 9 pairs to plot! Decision is based on lowest p-values for each pair.")

        myNormTable <- normalizeAllSubtables(myTable.LRT)
        uniquePairs <- unique(myNormTable$n)

        TheDivisionOfValues <- ceiling(max(myNormTable$v)/5)
        TheValues <- c(seq(from=TheDivisionOfValues, to=TheDivisionOfValues*4, by=TheDivisionOfValues), floor(max(myNormTable$v)))

        vectorOfYLabels <- data.frame(v=TheValues, n=c(TheValues[1]/max(myNormTable$v), TheValues[2]/max(myNormTable$v), TheValues[3]/max(myNormTable$v), TheValues[4]/max(myNormTable$v), TheValues[5]/max(myNormTable$v)))

        filename <- paste(c(myPrefix, "_", "BestNinePairsInDataset.pdf"),sep="", collapse="")

        pdf(filename, width=24, height=24)
        par(mfrow=c(3,3))
        for (i in 1:length(uniquePairs)){
            spl.pairs <- strsplit(uniquePairs[i], split="_")[[1]]
            n1 <- spl.pairs[1]
            n2 <- spl.pairs[2]
            if (!(is.null(eventIDs))){
                n1 <- evIDs[n1,]
                n2 <- evIDs[n2,]
            }
            main.name <- paste(c("event ", n1, " vs event ", n2), sep="", collapse="")
            sub.normTable <- myNormTable[myNormTable$n  ==  uniquePairs[i],]
            sub.lrtTable <- myTable.LRT[which(myTable.LRT$Event1  ==  spl.pairs[1] & myTable.LRT$Event2  ==  spl.pairs[2]),]

            vectorOfXLabels <- c("H0")
            for (k in 1:nrow(sub.normTable)){
                vectorOfXLabels <- c(vectorOfXLabels, listForXLabels[[sub.normTable[k,]$m]])
            }
            plotFromNormalizedTable(sub.normTable, sub.lrtTable, main.name, vectorOfYLabels, vectorOfXLabels)
        }
        dev.off()
    }
    fl.m <- 1
}