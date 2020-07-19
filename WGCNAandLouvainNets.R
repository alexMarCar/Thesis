
#We created WGCNA network as simple as we show here. 
#By using the function below from CoExpNets package, wecan create.
library(CoExpNets)
library(WGCNA)
netWGCNA = getDownstreamNetwork(tissue="WGCNAfiltCorrected",
                                 n.iterations=20,
                                 net.type = "signed",
                                 debug=F,
                                 expr.data=AMPPDdata,
                                 job.path=".")

# We created the Louvain network as we show here
# In this case, for creating the network, we first found communities over the graph made from expression datasets.
# Then, we could generate a network object as we can see below
library(igraph)
communitiesLouvain = cluster_louvain(graphAMP.PD)
netLouvain = NULL
netLouvain$moduleColors = WGCNA::labels2colors(communitiesLouvain$membership) # We add a color for each different community found
names(netLouvain$moduleColors) = communitiesLouvain$names # We name them
#And now the eigengenes
MEs = WGCNA::moduleEigengenes(AMPPDdata,
                                   colors = netLouvain$moduleColors,
                                   excludeGrey=F)
netLouvain$MEs = data.frame(MEs$eigengenes)

saveRDS(object = netLouvain, file = "LouvainNet20M.rds") # We save the net
