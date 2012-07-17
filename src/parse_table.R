#!/usr/bin/Rscript

parse.table  <- function(file="", header=FALSE, sep="\t", ncol=8){ 
  # load annotation file and parse line by line
  annot <- readLines(file)

  # structure of one line as follows
  # index 1: probe ID CHAR
  # index 2: num of matched probes INT
  # index 3...n: various identifiers CHAR

  # define regexps for identifiers for eacjh columns.
  RefSeq3 <- '[X|N][M|R]_[[:digit:]]'
  UniGene4 <- 'Mm\\.[[:digit:]]'
  RIKEN5 <- '[[:alnum:]]Rik'
  EntrezGene6 <- '[[:digit:]]'
  Symbol7 <- '[[:alnum:]]'
  Descr8 <- '[[:graph:]]'
  Ensembl9 <- 'ENS[[:alnum:]]'
  CH_Entrez10 <- '[[:digit:]]'

  list.names <- c('ProbeID','NumMatched','RefSeq','UniGene','RIKEN','EntrezGene','Symbol','Ensembl','CH_Entrez')

  # helper function to parse one ID
  parseID <- function(charvec,regexp,col,c) {
    if (col == 8) {
      # parsing the description
      # basically, everything is allowed
      # yet make a sanity check
      if (grepl(regexp,charvec,fixed=FALSE)) {
	# assume evrything ok
	return (paste(charvec, collapse=","))
      }
      else {
	#print('Not sure about that description')
	#print(paste('Line number ',c,sep=''))
	#print(charvec)
	return (NA)
      }
    }
    else {
      parts <- unlist(strsplit(charvec," ",fixed=TRUE))
      if (all(grepl(regexp,parts,fixed=FALSE))) {
	return (paste(parts, collapse=","))
      }
      else {
	#print('Not sure about that entry')
	#print(paste('Line number ',c,' and col ',col,sep=''))
	#print(parts)
	return (NA)
      }
    }
    return (NA)
  }

  # create the entry for the final mapping
  createEntry <- function(c,line) {
    content <- unlist(strsplit(line,sep,fixed=TRUE))
    pid <- content[1]
    nummatch <- ifelse(content[2] == "",NA,content[2])
    refseq <- ifelse(content[3] == "",NA,parseID(content[3],RefSeq3,3,c))
    unigene <- ifelse(content[4] == "",NA,parseID(content[4],UniGene4,4,c))
    riken <- ifelse(content[5] == "",NA,parseID(content[5],RIKEN5,5,c))
    entrez <- ifelse(content[6] == "",NA,parseID(content[6],EntrezGene6,6,c))
    symbol <- ifelse(content[7] == "",NA,parseID(content[7],Symbol7,7,c))
    #descr <- ifelse(content[8] == "",NA,parseID(content[8],Descr8,8,c))
    ensembl <- ifelse(content[9] == "",NA,parseID(content[9],Ensembl9,9,c))
    chentrez <- list(NA)
    if (length(content) == 10) {
      chentrez <- ifelse(content[10] == "",NA,parseID(content[10],CH_Entrez10,10,c))
    }
    retlist <- unlist(c(pid,nummatch,refseq,unigene,riken,entrez,symbol,ensembl,chentrez))
    if (length(retlist) != 9) {
	    print("length of array not 9, error:")
	    print(length(retlist))
	    print(retlist)
	    stop()
    }
    #names(retlist) <- list.names
    return(retlist)
  }



  data  <- NULL
  for (c in c(2:length(annot))) {
    entry <- createEntry(c,annot[c])
    #print(entry)
    data <- rbind(data,entry)
    #assign(x=entry$ProbeID,value=entry,envir=gnf1m.map)
  }
  colnames(data) <- list.names
  rownames(data) <- data[["ProbeID"]]
  data = data.frame(data)

  return(data)

}
