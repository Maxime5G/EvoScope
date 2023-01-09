# Parsing the result from epocs
# Export the concatenated results table (i.e., all scores, likelihoods etc for each model)
# Export the LRT results for all comparisons
# Export the best selected model according to our criterion (i.e., significant likelihood gain should be 1.5x as high as the previous best model

args=commandArgs(trailingOnly=TRUE)

# --- FUNCTIONS ---

# --- Calculating the Likelihood ratio test ---
# Calculating the LRT using the formula: 2 * (log(LK_model2) - log(LK_model1))
# Then comparing the value to the appropriate chi-square distribution (df = diff in parameters between both models)
# Also formatting it nicely for the output table

pairwiseChiSquareAndFormat <- function(subtable, paramTable, possiblePairs, outTable)
{
	seqResult <- c()
	pvalueResults <- c()
	tmpTable <- data.frame(Event1=integer(), Event2=integer(), Model1=character(), Model2=character(), LRT=numeric(), Pvalue=numeric(),combinedModels=character(), df=numeric())

	for (i in seq(1,nrow(subtable)-1))
	{
		for (j in seq(i+1,nrow(subtable)))
		{
			criticalValue <- 2 * (subtable[j,]$lnML - subtable[i,]$lnML)
			m1 <- subtable[i,]$Model
			m2 <- subtable[j,]$Model
			combOfModels <- paste(c(m1, m2), sep='', collapse='')
			numbOfDF <- paramTable[names(paramTable)==m2]-paramTable[names(paramTable)==m1]
			# if (numbOfDF>=1){
			if (combOfModels %in% possiblePairs){

				pValue <- format(pchisq(criticalValue, df=numbOfDF, lower.tail=FALSE), scientific=TRUE, digits=4)

				seqResult <- c(seqResult, criticalValue)
				pvalueResults <- c(pvalueResults, pValue)
				if (nrow(tmpTable) == 0){
					tmpTable[1,] <- c(subtable[1,]$Event1, subtable[1,]$Event2, subtable[i,]$Model,subtable[j,]$Model,criticalValue, pValue, NA, NA)
				}else{
					tmpTable <- rbind.data.frame(tmpTable, c(subtable[1,]$Event1, subtable[1,]$Event2, subtable[i,]$Model,subtable[j,]$Model,criticalValue, pValue, NA, NA))
				}
			}
		}
	}

	if(nrow(tmpTable) == 0){return (outTable)}

	tmpTable$LRT <- as.numeric(tmpTable$LRT)
	tmpTable$Pvalue <- as.numeric(format(as.numeric(tmpTable$Pvalue), scientific=T, digits=4))

	tmpTable <- tmpTable %>% rowwise() %>% mutate(combinedModels=paste(Model1, Model2, sep='', collapse=''))

	vecOfDF <- c()
	for (i in 1:nrow(tmpTable)){
		vecOfDF <- c(vecOfDF, degreeOfFreedoms[as.character(tmpTable[i,4]), as.character(tmpTable[i,3])])
	}
	tmpTable$df <- vecOfDF
	# print(tmpTable)

	outTable <- rbind.data.frame(outTable, tmpTable)
	return(outTable)
}

# --- Extracting the best model ---
# Idea to select the appropriate model:
# The model with the lowest number of parameters and the highest gain of LRT
# Idea is to compare each df, from lower to higher, and everytime compare
# How much I gain vs the previous df
# If I gain more than 1.5x: I save it.
# Input table should have only significant LRT

extractTheBest2 <- function(mySubtable, paramTable){
	resLRT <- NULL
 	saveLRTdiff <- NULL
	if (nrow(mySubtable) == 0){return(resLRT)}
	if (nrow(mySubtable) == 1){
    	return(1)
	}
	for (i in 1:(nrow(mySubtable))){
    	if (is.null(resLRT)){				# If nothing compared yet, save the line
			resLRT <- i
		}else{								# Need to compare!

			l1 <- mySubtable[resLRT,]
			l2 <- mySubtable[i,]
			df1 <- l1$df
			df2 <- l2$df
			LRT1 <- l1$LRT
			LRT2 <- l2$LRT
			m11 <- paramTable[names(paramTable)==l1$Model1]
			m12 <- paramTable[names(paramTable)==l1$Model2]
			m21 <- paramTable[names(paramTable)==l2$Model1]
			m22 <- paramTable[names(paramTable)==l2$Model2]

			if (df2-df1 == 0){
				# There are two possibilities so far:
				# - either we have ia/ib vs al/bl
				# - either we have ia/ib vs ib/ia
				# in the first case, I would retain the second one only if the stat test is >>> (maybe 1.5?)
				# in the second case, I would retain the best one

				if (LRT2>LRT1){
					if ((m22-m12)>0){
						LRTdiff <- LRT2/LRT1
						if (LRTdiff>1.5){	# Putting 1.5 here - sounds great so far - see below
							resLRT <- i
						}
					}else{
						resLRT <- i
					}
				}
			}
			else if (df2-df1 == 1){
				LRTdiff <- LRT2/LRT1
				if (LRTdiff>1.5){	# Putting 1.5 here - sounds great so far
					resLRT <- i
					saveLRTdiff <- LRTdiff
				}
			}
		}
	}
	return(resLRT)
}

# --- INITIALIZATIONS ---

# Here I'm saving the order of all models in a single table
# With the associated difference between the number of parameters
# Practical when calculating the chisquare

orderOfAllModels <- c('1', 'i', 'x', 'u', 'U', 'X', 'r', 'a', 'b', 'l', 'I', 'R', 'A', 'B', 'L', '8')
numParamAllModels <- c(1,2,2,3,3,3,3,3,3,4,4,4,5,5,6,8)

degreeOfFreedoms <- data.frame(matrix(rep(0, length(orderOfAllModels)*length(orderOfAllModels)), ncol=length(orderOfAllModels), nrow=length(orderOfAllModels)))
for (i in 1:length(orderOfAllModels)){
  for (j in 1:length(orderOfAllModels)){
    degreeOfFreedoms[i,j] <- numParamAllModels[i]-numParamAllModels[j]
  }
}

rownames(degreeOfFreedoms) <- orderOfAllModels
colnames(degreeOfFreedoms) <- orderOfAllModels

# --- MAIN LOOP ---

if (length(args)>2){stop('Rscript ParseEpocsFinal.R [FolderPrefix] [eventIDsTable]
				 [FolderPrefix]: name of folder with results of epocs
				 [eventIDsTable]: if you have it, 1-column file with the name of your traits (full table)')}

if (length(args)<1){stop('Rscript ParseEpocsFinal.R [FolderPrefix] [eventIDsTable]
				 [FolderPrefix]: name of folder with results of epocs
				 [eventIDsTable]: if you have it, 1-column file with the name of your traits (full table)')}


# library(tidyverse)
suppressPackageStartupMessages(library(tidyverse))

myFolder <- args[1]
EvIds <- if (length(args)<2) NA else args[2]

if (!(is.na(EvIds))){myEventsIds <- tibble (read.table(EvIds, header=F, stringsAsFactors = F, colClasses = c("character")))}

# To evaluate if I have to run the script from inside the folder ('.' in myFolder variable)
condition <- FALSE

nullHypothesisModel <- 'i'

allPossibleCombinationsWith8Models <- c('iI', 'ia', 'ib', 'il', 'al', 'bl', 'IA', 'IB', 'IL', 'AL', 'BL', 'bB', 'aA', 'lL')

# If the user wants to run in current directory: need to adjust a bit the variable to account for it
if (myFolder == '.')
{
	myFolder=getwd()
	myOutputPrefix=tail(strsplit(myFolder, split='/')[[1]], n=1)
	basename=myOutputPrefix
	myPrefix=paste(c(myOutputPrefix, '_mat_epocs'), sep='', collapse='')
	condition=TRUE
}else{
	folderSplit <- strsplit(myFolder, '/')[[1]]
	basename <- folderSplit[length(folderSplit)]
	myFolder <- paste(folderSplit, sep='', collapse='/')

	myPrefix <- paste(c(basename, '_mat_epocs'), sep='', collapse='')
}

# all output files from epocs - one per model!
allFiles <- list.files(myFolder, pattern=myPrefix)

if (length(allFiles)<2){stop('Error reading epocs files - not enough files present! Make sure you chose the correct folder')}

df <- data.frame()

for (f in allFiles)
{
	# retrieving the content of each file
	myFile <- paste(c(myFolder, '/', f), sep='', collapse='')
	if (condition){myFile=f}
	myT <- read.table(myFile, header=T, stringsAsFactors=F, sep='\t', row.names=NULL, na.strings='-', colClasses=c('integer', 'integer','integer', 'character', rep('numeric', 10)))

	# and adding all of them in the same data frame
	df=rbind.data.frame(myT, df)
}

df$Model <- unlist(lapply(df[,4], trimws)) # Removing white spaces (because of C)

# I'm now reordering the data frame such that the events are sorted in ascending order.
# This is important because afterwards, I'll extract the subtables for each pair.
# Therefore, they need to be in the correct order.
df <- df[order(df$Event1, df$Event2),]

# Models in the analysis
inputModels <- names(table(df$Model))

# Sanity check - need to make sure there is at least two models chosen, differing in the number of df
subtableModels <- degreeOfFreedoms[inputModels,inputModels]
if(length(unique(unlist(subtableModels))) < 2){stop('Error reading models - maybe you dont have H0 - check your epocs run')}

# Number of models: needed to extract the subtables of all models for each pairs
numbEntriesPerPair <- length(inputModels)

# Number of pairwise comparisons possible - needed for the dimensions of the output table
# countComparisons <- (numbEntriesPerPair*(numbEntriesPerPair-1))/2

countComparisons <- NA

matResult <- data.frame(Event1=integer(), Event2=integer(), Model1=character(), Model2=character(), LRT=numeric(), Pvalue=numeric(),combinedModels=character(), df=numeric())

matBestModel <- NULL

count <- 1

for (k in seq(1,nrow(df), numbEntriesPerPair))
{
    # subtable for each pairs in the test analysis
    subtable <- df[seq(k,k+numbEntriesPerPair-1),]

    # Now I need to reorder the table such that the models are sorted in ascending order of complexity
    # This is essential since the LRT compares the more complex to the less complex model
    allModels <- subtable$Model
    orderOfMyModels <- match(orderOfAllModels, allModels)
    orderOfMyModels <- orderOfMyModels[which(!(is.na(orderOfMyModels)))]
    subtable <- subtable[c(orderOfMyModels),]
    df[seq(k,k+numbEntriesPerPair-1),] <- subtable

	names(orderOfMyModels) <- allModels
	names(numParamAllModels) <- orderOfAllModels

    # And now I'm doing the chi-square analysis (see function above)
	# Adds automatically the lines to the matResult object

    matResult <- pairwiseChiSquareAndFormat(subtable, numParamAllModels, allPossibleCombinationsWith8Models, matResult)

	# I'm extracting, for the first run, the number of comparisons performed.
	# This should work fine if the models in the epocs folder are correctly selected
	if (is.na(countComparisons)){countComparisons <- nrow(matResult)}

	# Extracting the best LRT for the pair considered - need to subset the full table
    test.subtable <- matResult[seq(count,count+countComparisons-1),] %>% filter(Pvalue < 0.05) %>% filter(df>0) %>% group_by(df)

	if ('iI' %in% test.subtable$combinedModels && test.subtable[which(test.subtable$combinedModels == 'iI'),]$Pvalue < 0.05){
		test.subtable <- test.subtable[c(which(test.subtable$Model1 %in% c('I','A','B','L')), which(test.subtable$combinedModels == 'iI')),]
	}

	# Retrieving the "best" pair
    myLRT <- extractTheBest2(test.subtable, numParamAllModels)
    if (!(is.null(myLRT))){
        if (is.null(matBestModel)){
          matBestModel <- test.subtable[myLRT,]
        }else{
          matBestModel <- rbind(matBestModel, test.subtable[myLRT,])
        }
    }

	count <- count+countComparisons
}

if (is.null(matBestModel)){
	# 1. Table of all LRTs
	outputTable <- paste(myFolder, '/', basename, '_LRT.tab', sep='', collapse='')
	if (condition){outputTable <- paste(myOutputPrefix, '_LRT.tab', sep='', collapse='')}
	write.table(matResult, outputTable, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

	# 2. Concatenated results (i.e., all pairs within the same table)
	outputTable2 <- paste(myFolder, '/', basename, '_ConcatenatedResults.tab', sep='', collapse='')
	if (condition){outputTable2 <- paste(myOutputPrefix, '_ConcatenatedResults.tab', sep='', collapse='')}
	write.table(df, outputTable2, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

	# 3. Warning
	print("No significant events found!")
	text_message = sprintf("Result files %s and %s are written in the output directory", outputTable, outputTable2)
	print(text_message)
	print("Exiting...")

	quit()
}

# Adding the trait ids (when provided)
if (!(is.na(EvIds))){
	e1 <- myEventsIds[as.numeric(matBestModel$Event1),]$V1
	e2 <- myEventsIds[as.numeric(matBestModel$Event2),]$V1
	matBestModel$name_e1 <- e1
	matBestModel$name_e2 <- e2
}

# Adding the -log(PVAL) for completeness
matBestModel$logPVAL <- -log10(matBestModel$Pvalue)

# Adding the lambda's
lambda1 <- c()
lambda2 <- c()

for (i in 1:nrow(matBestModel)){
	entryInDf <- df[which(df$Event1 == matBestModel[i,]$Event1 & df$Event2 == matBestModel[i,]$Event2 & df$Model == matBestModel[i,]$Model2),]
	lambda1 <- c(lambda1, entryInDf$mu1star/entryInDf$mu1)
	lambda2 <- c(lambda2, entryInDf$mu2star/entryInDf$mu2)
}

matBestModel$lambda1 <- lambda1
matBestModel$lambda2 <- lambda2

# --- END: Writing the outputs ---

# 1. Table of all LRTs
outputTable <- paste(myFolder, '/', basename, '_LRT.tab', sep='', collapse='')
if (condition){outputTable <- paste(myOutputPrefix, '_LRT.tab', sep='', collapse='')}
write.table(matResult, outputTable, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

# 2. Concatenated results (i.e., all pairs within the same table)
outputTable2 <- paste(myFolder, '/', basename, '_ConcatenatedResults.tab', sep='', collapse='')
if (condition){outputTable2 <- paste(myOutputPrefix, '_ConcatenatedResults.tab', sep='', collapse='')}
write.table(df, outputTable2, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

# 3. The best model table (i.e., for each pair, which model explains better the data observed)
outputTable3 <- paste(myFolder, '/', basename, '_BestModel.tab', sep='', collapse='')
if (condition){outputTable3 <- paste(myOutputPrefix, '_BestModel.tab', sep='', collapse='')}
write.table(matBestModel, outputTable3, sep='\t', quote=F, col.name=T, row.names=F)
