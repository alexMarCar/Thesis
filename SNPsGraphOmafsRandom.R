
graphSNPmaker4 <- function(SNPdata,SNPs,m,a) {
  library(igraph)
  
  # First, we seppate controls from cases
  CasesIndexes = which(SNPdata[,2] == 1)
  casesDF = SNPdata[CasesIndexes,]
  controlsDF = SNPdata[-CasesIndexes,]
  
  ## Then, we calculate the allele frequency for each SNP over cases and controls samples separately.
  totalCases = as.numeric(dim(casesDF)[1])
  totalControls = as.numeric(dim(controlsDF)[1])
  
  
  OMAFcases <- apply(casesDF[,3:length(colnames(casesDF))], 2, function(x){
    as.numeric(length(unique(which(x == 1 | x == 2)))) / totalCases
  })
  OMAFcontrols <- apply(controlsDF[,3:length(colnames(controlsDF))], 2, function(x){
    as.numeric(length(unique(which(x == 1 | x == 2)))) / totalControls  #sum( x !=0)
  })
  
  # This time we filter according to MOMAF (m)
  casesSNPSelect = apply(as.data.frame(OMAFcases), 2, function(x){
    which(x > m)
  })
  controlsSNPSelect = apply(as.data.frame(OMAFcontrols),2, function(x) {
    which(x > m)
  })
  
  OMAFcases = OMAFcases[casesSNPSelect]
  OMAFcontrols = OMAFcontrols[controlsSNPSelect]
  
  #Once we got rid of odd SNPs, we randomly select 1500 SNPs for both controls and cases dataset.
  SelectedControls = sample(OMAFcontrols, size = SNPs)
  SelectedCases = sample(OMAFcases, size = SNPs)
  
  controlsDF2 = controlsDF[,3:length(colnames(controlsDF))]
  casesDF2 = casesDF[,3:length(colnames(casesDF))]
  
  # Then, we can filter cases and controls dataframes, according to those SNPs whose OMAF was over m threshold
  controlsDF2 = controlsDF2[,names(SelectedControls)]
  rownames(controlsDF2) = controlsDF$ID
  casesDF2 = casesDF2[,names(SelectedCases)]
  rownames(casesDF2) = casesDF$ID
  
  
  
  #Now, we have to calculate the ratio for each pair of SNPs. We will have to compare all SNPs among them, and we will be measuring     the times each pair of SNPs appear simultuneously mutated.
  
  # I managed with indexes so the ones from the extremes do not suppose us an issue. I run "i" until the last SNP index -1.
  # As R indexes start from 1, we stablish the i range from 2.
  
  # For each i index, which is equal to a SNP, it is faced to any other SNP, and get the JOMAF values.
  # It it important to see that the j index, only takes SNPs which has not already been compared to the SNP stored in i.
  
  
  ## First we generate the graph for controls
  
  # Initializing variables
  edgesControls = NULL
  JOMAFcontrols = NULL
  GenesControls = length(colnames(controlsDF2))
  samplesControls = as.numeric(dim(controlsDF2)[1])
  
  
  # We run as I said already
  for (i in 1:(GenesControls-1)) {
    
    for (j in (i+1):GenesControls) {
      counts = 0
      #Then we compare i and j SNPs states for each sample. If both SNPs are not 0, it means they are not mutated.
      for (k in 1:samplesControls) {
        if ((controlsDF2[k,j] != 0) & (controlsDF2[k,i] !=0)) {
          counts = counts +1
        }
      }
      # Finally, if jomaf value is higher than omaf value, we add the SNPs and the jomaf value to the list we initialized previously
      jomafC = counts / samplesControls
      if (jomafC > a) {edgesControls = c(edgesControls,c(names(controlsDF2)[[i]],names(controlsDF2)[[j]]))
      JOMAFcontrols = c(JOMAFcontrols, jomafC)
      }
    }
    
  }
  #Then we can generate the graph and we save it into an RDS object
  graphdf = NULL
  graphdf = make_graph(edges=edgesControls, directed = F)
  E(graphdf)$weight = JOMAFcontrols
  saveRDS(object = graphdf, file = paste0("~/TFM/SNPs/",SNPs,"GraphControlsRandom.rds"))
  
  
  # We still have to generate the graph for cases
  
  # Initializing variables
  edgesCases = NULL
  JOMAFcases = NULL
  GenesCases = length(colnames(casesDF2))
  samplesCases = as.numeric(dim(casesDF2)[1])
  
  # We run as I said already
  for (i in 1:(GenesCases-1)) {
    
    for (j in (i+1):GenesCases) {
      counts = 0
      #Then we compare i and j SNPs states for each sample. If both SNPs are not 0, it means they are not mutated.
      for (k in 1:samplesCases) {
        if ((casesDF2[k,j] != 0) & (casesDF2[k,i] !=0)) {
          counts = counts +1
        }
      }
      # Finally, if jomaf value is higher than omaf value, we add the SNPs and the jomaf value to the list we initialized            previously
      jomafC = counts / samplesCases
      if (jomafC > a) {edgesCases = c(edgesCases,c(names(casesDF2)[[i]],names(casesDF2)[[j]]))
      JOMAFcases = c(JOMAFcases, jomafC)
      }
    }
    
  }
  #Then we can generate the graph and we save it into an RDS object
  graphdf = NULL
  graphdf = make_graph(edges=edgesCases, directed = F)
  E(graphdf)$weight = JOMAFcases
  saveRDS(object = graphdf, file = paste0("~/TFM/SNPs/",SNPs,"GraphCasesRandom.rds"))
  
}
