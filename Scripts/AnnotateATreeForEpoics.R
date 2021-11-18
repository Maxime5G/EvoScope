# Annotate a tree with node and tip labels of mutations
# reads the output of pastml and assigns mutation values to nodes and tips
# the pastml file should have one column per trait and one column for the tip and node ids

args=commandArgs(trailingOnly=TRUE)

if (length(args) == 4){

  options(scipen=999)       # Disabling scientific notations
  options(warn=0)           # Original setting for warnings
  library(ape)              # Library for manipulating phylogenetic trees

  myTree <- args[1]         # Newick tree
  myFile <- args[2]         # PastML file
  outputTree <- args[3]     # Output newick tree with mutations mapped
  outputMat <- args[4]      # Output table of mutations

  outlist = list()

  theTree <- read.tree(myTree)

  # CHECK 1: Tree has branch ids.

  if (is.null(theTree$node.label)){stop('Please input a phylogenetic tree with named nodes! Exiting...')}

  modifMat <- read.table(myFile, header=T, stringsAsFactors=F, check.names=F, sep='\t')

  # CHECK 2: Trait matrix is tab-delimited
  if (ncol(modifMat)==1){stop('Please make sure to have a tab-delimited file in input! Exiting...')}

  colnames(modifMat) <- c('tip', colnames(modifMat)[2:ncol(modifMat)])

  # Now I should parse the vectors to test the presence of mutations.

  allIDsInTheTree <- c(theTree$tip.label, theTree$node.label)

  # CHECK 3: Tab-delimited file has the correct format (first column = tip/node ids, columns 2-N are traits)
  if (length(allIDsInTheTree) != nrow(modifMat)){stop('Error reading ACR file - please make sure the format is correct! Exiting...')}

  outlist_mut = list()

  allPaths <- nodepath(theTree)

  for (k in 2:ncol(modifMat)){      # For each event
      for (m in 1:nrow(modifMat)){  # For each node

    		theVar <- modifMat[m,k]

            theID <- as.character(modifMat[m,1]) # current branch id
            if (!(theID %in% names(outlist_mut))){
                outlist_mut[[theID]] <- c()
            }

            # Here I'm converting the tip/node id to the position in the edge labelling from ape
            thePosInTheEdges <- match(theID,allIDsInTheTree)

            # Now I'm retreiving the ancestor in the edges tables from ape
            theEdgeID <- theTree$edge[which(theTree$edge[,2] == thePosInTheEdges),]

            # If length==0 -> root (i.e. no ancestor possible)
            # If length==2 -> intermediate edge
            if (length(theEdgeID) == 2){
                ancestorID <- allIDsInTheTree[theEdgeID[1]]
                ValueAncestor <- modifMat[which(modifMat$tip == ancestorID),k]

                if (ValueAncestor == theVar){
                    outlist_mut[[theID]] <- c(outlist_mut[[theID]], 0)
                }else{
                    outlist_mut[[theID]] <- c(outlist_mut[[theID]], 1)
                }
            }else{
                outlist_mut[[theID]] <- c(outlist_mut[[theID]], 0)
            }
      } # end of for m in 1:nrow
  }     # end of for k in 2:ncol

  for (i in 1:theTree$Nnode){
      theTree$node.label[i] = paste('[', paste(outlist_mut[[theTree$node.label[i]]], sep='', collapse='/'), ']', sep='', collapse='')
  }

  for (i in 1:length(theTree$tip.label)){
      tmp = theTree$tip.label[i]
      theTree$tip.label[i] = paste(tmp, paste('[', paste(outlist_mut[[tmp]], sep='', collapse='/'), ']', sep='', collapse=''), sep='', collapse='')
  }

  write.tree(theTree, file=outputTree, digits=8)

  for (n in names(outlist_mut)){
      modifMat[modifMat$tip == n, c(2:ncol(modifMat))] <- outlist_mut[[n]]
  }
  write.table(modifMat, outputMat, col.names=T, row.names=F, sep='\t', quote=F)

}else{writeLines("Rscript AnnotateATree_Correct_Uncertainty.R [Tree] [PastML_output] [OutputTree] [OutputMat]
                    [Tree]: Phylogenetic tree (NEWICK) with named nodes
                    [PastML_output]: Ancestral state reconstruction file by pastml. Rows are nodes/leaves. First column MUST be the node/leaf ID, columns 2->N are reconstructed traits
                    [OutputTree]: Output phylogenetic tree (NEWICK) with mutations mapped on branches - formatted for epics/epocs
                    [OutputMat]: Output tab-delimited file with mutations for each node/leaf. Rows are nodes/leaves, columns are traits")}
