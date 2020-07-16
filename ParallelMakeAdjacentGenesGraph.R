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


unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


graphMakerParallel <- function(exprData, genesNumber,arcsSet) {
  library(igraph)
  library(doParallel)
  library(foreach)
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1]-3) #not to overload the server 
  registerDoParallel(cl)
  registerDoSEQ()
  
  adjacency = adjmatrix(expr.data) ## Cacluclating the adjacency matrix by using the function we made previously
  
  alladjs = apply(adjacency,2,sum)
  alladjs = order(alladjs,decreasing=T)[1:genesNumber]
  adj = adjacency[alladjs,alladjs]
  adj[upper.tri(adj, diag = T)] = -1
  
  allArcsIndexes = arrayInd(order(adj, decreasing = T),dim(adj)) 
  arcsIndexes = allArcsIndexes[1:arcsSet,]
  
  arcsSets = as.list(seq(from = 1, to =arcsSet, length.out = 6))
  
  previousSet =0
  nextSet = 0
  
  output = foreach(i = 1:(length(arcsSets)-1), .combine = list) %dopar% {
    
    previousSet = as.numeric(arcsSets[i])
    nextSet = as.numeric(arcsSets[i+1])
    edges = NULL
    ws = NULL
    for (arc in previousSet:nextSet) {
      indexes = arcsIndexes[arc, ]
      newWs = adj[indexes[1],indexes[2]] 
      
      if (newWs == -1) break
      ws = c(ws,newWs)
      edges = c(edges, c(colnames(adj)[indexes[1]], colnames(adj)[indexes[2]])) 
    }
    list(ws = ws,edges = edges)
  }
  
  edges = list(c(output[[1]][[1]][[1]][[1]]$edges,output[[1]][[1]][[1]][[2]]$edges,output[[1]][[1]][[2]]$edges,
                 output[[1]][[2]]$edges,output[[2]]$edges))
  ws = list(c(output[[1]][[1]][[1]][[1]]$ws, output[[1]][[1]][[1]][[2]]$ws, output[[1]][[1]][[2]]$ws, output[[1]][[2]]$ws,
              output[[2]]$ws))
  
  graphdf = NULL
  graphdf = make_graph(edges=edges[[1]], directed = F)
  E(graphdf)$weight = ws[[1]]
  V(graphdf)$size = alladjs[match(names(V(graphdf)),names(alladjs))]
  V(graphdf)$size = V(graphdf)$size/max(V(graphdf)$size)
  
  saveRDS(graphdf, file = paste0("graph",arcsSet,".rds"))
  
  stopCluster(cl)
  unregister()
}
