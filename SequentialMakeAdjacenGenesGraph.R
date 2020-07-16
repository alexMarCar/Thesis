# Given a gene expression data frame, it calculates the squared matrix, by using a determined power
adjmatrix <- function(exprData) {
  library(WGCNA)
  paste("Adjusting beta value for the given matrix")
  b.study = WGCNA::pickSoftThreshold(exprData,
                                     powerVector=1:20,
                                     verbose=5,
                                     moreNetworkConcepts=TRUE,
                                     networkType="signed",
                                     corOptions=list(use='p'))
  r.sq.cutoff = 0.8
  max.k.cutoff = 150
  beta = min(b.study$fitIndices[as.numeric(b.study$fitIndices$SFT.R.sq) > r.sq.cutoff &
                                  as.numeric(b.study$fitIndices$slope) < 0 &
                                  as.numeric(b.study$fitIndices$max.k) > max.k.cutoff,"Power"])
  paste("Calculating adjacency matrix from expression data and a power equal to ",beta)
  #This is the squared matrix of adjacency between genes we will use as a meausure of "Friendlyness" between genes
  adjacency = WGCNA::adjacency(exprData, 
                               power = beta,
                               type = "signed", 
                               corOptions = list(use='p'))
}

# Given a gene expression data frame, the genesNumber are wanted to make the graph and the arcsSet to build the graph, it makes a graph.
graphmaker <- function(exprData, genesNumber,arcsSet) {
  library(igraph)
  
  adjacency = adjmatrix(exprData) ## Cacluclating the adjacency matrix by using the function we made previously
  
  alladjs = apply(adjacency,2,sum)
  alladjs = order(alladjs,decreasing=T)[1:genesNumber] ## Recalculating the adjacency matrix, for the genes selected
  adj = adjacency[alladjs,alladjs]
  adj[upper.tri(adj, diag = T)] = -1
  
  allArcsIndexes = arrayInd(order(adj, decreasing = T),dim(adj)) 
  
  arcsIndexes = NULL
  arcsSets = seq(from=1, to=arcsSet, length.out = 5)
  
  arcsIndexes = allArcsIndexes[1:arcsSet,]
  edges = NULL
  ws = NULL 
  
  for (i in 1:dim(arcsIndexes)[1]) {
    indexes = arcsIndexes[i, ]
    newWs = adj[indexes[1],indexes[2]] #adj
    
    if (newWs == -1) break
    ws = c(ws,newWs)
    edges = c(edges, c(colnames(adj)[indexes[1]], colnames(adj)[indexes[2]])) #adj
  }
  
  graphdf = NULL
  graphdf = make_graph(edges=edges, directed = F)
  E(graphdf)$weight = ws
  V(graphdf)$size = alladjs[match(names(V(graphdf)),names(alladjs))]
  V(graphdf)$size = V(graphdf)$size/max(V(graphdf)$size)
  
  saveRDS(graphdf, file = paste0("graph",arcsSet,".rds"))
}

