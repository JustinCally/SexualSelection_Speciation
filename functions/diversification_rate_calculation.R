#Equal splits statistic calculation ("DR Statistic")
calculate.log.es <- function(phy) {
rootnode <- length(phy$tip.label) + 1
es <- numeric(length(phy$tip.label))
for (i in 1:length(es)){
  node <- i
  index <- 1
  qx <- 0
  while (node != rootnode){
    el <- phy$edge.length[phy$edge[,2] == node]
    node <- phy$edge[,1][phy$edge[,2] == node]			
    qx <- qx + el* (1 / 2^(index-1))			
    index <- index + 1
  }
  es[i] <- 1/qx
	
names(es) <- phy$tip.label
}

es <- log(es[phy$tip.label])
result <- as.vector(list(es))
names(result) <- "es"
return(result)
}

# node depth statistic calculation
calculate.nd <- function(phy) {
rootnode <- length(phy$tip.label) + 1
nd <- numeric(length(phy$tip.label)) # Number of tips
for (i in 1:length(nd)){ # For each tip
  node <- i
  nodecount <- 0
  while (node != rootnode){
    node <- phy$edge[,1][phy$edge[,2] == node] # Work back a node	
    nodecount <- nodecount + 1
  }
  nd[i] <- nodecount/max(branching.times(phy))
names(nd) <- phy$tip.label
}
result <- as.vector(list(nd))
names(result) <- "nd"
return(result)
}

#terminal branch statistic 
calculate.tb <- function(phy) {
  
  require(ape)
  require(mvtnorm)
  
  # Calculate terminal edge lengths
  n <- length(phy$tip.label)
  # based on post on Liam Revell's blog:
  invis <- setNames(phy$edge.length[sapply(1:n, function(x,y) which(y==x), y=phy$edge[,2])], phy$tip.label)
  tb <- 1/invis
  
  tb <- log(tb[phy$tip.label]) # log transform
  names(tb) <- phy$tip.label
  
  result <- as.vector(list(tb))
  names(result) <- "tb"
  return(result)
  
}
