# Annotate a tree with node and tip labels of mutations
# reads the output of pastml and assigns mutation values to nodes and tips
# the pastml file should have one column per trait and one column for the tip and node ids

args=commandArgs(trailingOnly=TRUE)

if ((length(args) < 6) & (length(args) > 3)){

  options(scipen=999)       # Disabling scientific notations
  options(warn=0)           # Original setting for warnings
  library(ape)              # Library for manipulating phylogenetic trees

  myTree <- args[1]         # Newick tree
  myFile <- args[2]         # PastML file
  outputTree <- args[3]     # Output newick tree with mutations mapped
  outputMat <- args[4]      # Output table of mutations
  mask <- 0
  eventOrder <- NULL

  if (length(args)==5) {
      mask <- args[5]
      # if (!(mask %in% c(1,2))){stop('Error: Please input either 1 or 2 for the mask! Exiting...')}
      if (!(mask %in% c(1,2))){
          mask <- 0
          eventOrderFile <- read.table(args[5])
          eventOrder <- eventOrderFile$V1
          if (nrow(eventOrderFile)==0){stop('Error: Please input a correct event ID file! Exiting...')}
      }
  }

  outlist = list()

  theTree <- read.tree(myTree)

  # CHECK 1: Tree has branch ids.

  if (is.null(theTree$node.label)){stop('Please input a phylogenetic tree with named nodes! Exiting...')}
  modifMat <- read.table(myFile, header=T, stringsAsFactors=F, check.names=F, sep='\t')
  if (!(is.null(eventOrder))){modifMat <- modifMat[,c('node', eventOrder)]}

  # CHECK 2: Tree has branch lengths>0
  if(length(theTree$edge.length[theTree$edge.length==0])>0){
      print("Setting branch lengths to minimum length 1e-10!")
      theTree$edge.length[theTree$edge.length==0] <- 1e-10
  }

  # CHECK 3: Trait matrix is tab-delimited
  if (ncol(modifMat)==1){stop('Please make sure to have a tab-delimited file in input! Exiting...')}

  colnames(modifMat) <- c('tip', colnames(modifMat)[2:ncol(modifMat)])

  # Now I should parse the vectors to test the presence of mutations.

  allIDsInTheTree <- c(theTree$tip.label, theTree$node.label)

  # CHECK 4: Tab-delimited file has the correct format (first column = tip/node ids, columns 2-N are traits)
  if (length(allIDsInTheTree) != nrow(modifMat)){stop('Error reading ACR file - please make sure the format is correct! Exiting...')}

  # Extracting the root ID
  theRoot <- theTree$node.label[1]

  # The plan going forward is to compute the difference between the daughter nodes and the parent nodes
  # I can therefore have different behaviours:
  # Null values - nothing happened (1-1, 0-0)
  # Positive values - gain/mutation (1-0)
  # Negative values - losses (0-1)
  # What I'm doing is to compute the differences between vectors of presence-absence of daughters and parent


  # Preparing the output table: same size as the ACR file - filled with NA's
  # Root node is always 0 (nothing happened yet)
  out_table <- modifMat[,]
  out_table[,2:ncol(out_table)] <- NA
  out_table[out_table$tip==theRoot,2:ncol(out_table)] <- rep(0,(ncol(out_table)-1))

  for (i in 1:nrow(theTree$edge)){

        name1 <- allIDsInTheTree[theTree$edge[i,1]]
        name2 <- allIDsInTheTree[theTree$edge[i,2]]
        # print(paste(c(name1, name2), sep='', collapse=' '))

        vec1 <- unlist(c(modifMat[modifMat$tip == name1,seq(2,ncol(out_table))]))
        vec2 <- unlist(c(modifMat[modifMat$tip == name2,seq(2,ncol(out_table))]))

        resVector <- as.numeric(vec2!=vec1)

        if (-1 %in% vec1){
            resVector[which(vec1 <0) ] <- -1
        }
        if (-1 %in% vec2){
            resVector[which(vec2 <0) ] <- -1
        }
        out_table[out_table$tip == allIDsInTheTree[theTree$edge[i,2]],seq(2,ncol(out_table))] <- resVector

  }

  for (i in 1:theTree$Nnode){
      theTree$node.label[i] = paste('[', paste(c(out_table[out_table$tip == theTree$node.label[i],seq(2,ncol(out_table))]),sep='', collapse='/'),']', sep='', collapse='')
  }
  for (i in 1:length(theTree$tip.label)){
      tmp = theTree$tip.label[i]
      theTree$tip.label[i] = paste(tmp, paste('[', paste(c(out_table[out_table$tip == tmp,seq(2,ncol(out_table))]),sep='', collapse='/'),']', sep='', collapse=''),sep='', collapse='')
  }
  write.tree(theTree, file=outputTree, digits=8)
  write.table(out_table, outputMat, col.names=T, row.names=F, sep='\t', quote=F)

}else{writeLines("Rscript AnnotateATree_Correct_Uncertainty.R [Tree] [PastML_output] [OutputTree] [OutputMat] [Mask?]
                    [Tree]: Phylogenetic tree (NEWICK) with named nodes
                    [PastML_output]: Ancestral state reconstruction file by pastml. Rows are nodes/leaves. First column MUST be the node/leaf ID, columns 2->N are reconstructed traits
                    [OutputTree]: Output phylogenetic tree (NEWICK) with mutations mapped on branches - formatted for epics/epocs
                    [OutputMat]: Output tab-delimited file with mutations for each node/leaf. Rows are nodes/leaves, columns are traits
                    [Mask?]:(Optional) Decide if you want to mask gains (1) or losses (2)")}
