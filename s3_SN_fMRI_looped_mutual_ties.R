#Carolyn McNabb - January 2019 - functional connectivity and connectomics analysis as a function of social network structure

# years <- c("Lower_4_2017","Year_4_2017","Lower_4_2018")
years <- c("Cohort 1","Cohort 3","Cohort 2")

setwd("F:\\0_parcellation_analysis/scripts-data-sharing/")

networkfiles <- c("corr_all.csv","corrDMN_all.csv","corrsalience_all.csv","corrlFPN_all.csv","corrrFPN_all.csv","corr_s_star.csv","corr_h_star.csv")
networks <- lapply(X = networkfiles, FUN = read.csv, header = F)
nets <- c("wholebrain", "DMN", "salience", "lFPN", "rFPN","strength","diversity")
names(networks) <- nets

## list of girls in each year group who were scanned
scanned_list <- list()
scanned_list[[1]] <- c("N201", "N204", "N210", "N211", "N213", "N216", "N217", "N218", "N219", "N220", "N222", "N223", "N224", "N226", "N227", "N232", "N233", "N240", "N244", "N245", "N246", "N247", "N256")
scanned_list[[2]] <- c("N300", "N302", "N305", "N307", "N308", "N309", "N312", "N315", "N317", "N318", "N321", "N322", "N325", "N332", "N333", "N334", "N335", "N340", "N344", "N346", "N347", "N350", "N352", "N355", "N361", "N366", "N367", "N369")
scanned_list[[3]] <- c("N702", "N704", "N706", "N708", "N709", "N711", "N713", "N715", "N720", "N722", "N737", "N739", "N742", "N743", "N745", "N748", "N751")
names(scanned_list) <- years

same_community <- c("communities_L4.csv","communities_4.csv", "communities_700s.csv")
same_comms <- lapply(X = same_community, FUN = read.csv, header = F)
names(same_comms) <- years


#create variable called modularity difference that contains information about the similarities
#in brain modularity (determined using Brain Connectivity Toolbox) between members of each dyad.
#this may have to be added as a matrix to the list that contains whole brain, RSNs etc.
modularity_data <- c("modularity_sim_L4.csv","modularity_sim_4.csv","modularity_sim_700.csv")
brain_mod <- lapply(X = modularity_data, FUN = read.csv, header = F)
names(brain_mod) <- years
brain_modularity <- list()
for (y in years){
  colnames(brain_mod[[y]]) <- scanned_list[[y]]
  rownames(brain_mod[[y]]) <- colnames(brain_mod[[y]])
  brain_mod[[y]][upper.tri(brain_mod[[y]],diag=TRUE)] <- NA
  brain_modularity[[y]] <- reshape2::melt(brain_mod[[y]], na.rm=TRUE, value.name="modularity_difference")#melt the data and label the column "modularity_difference"
}



library(pacman);p_load(BBmisc);p_load(igraph);p_load(psych);
#p_load(SDMTools);
p_load(ggplot2);p_load(BioGeoBEARS);p_load(phylobase);p_load(dplyr);p_load(wesanderson);p_load(lme4);p_load("ggpubr")
setwd("F:\\0_parcellation_analysis/scripts-data-sharing/")
#source("motivationscript.R")
source("http://peterhaschke.com/Code/multiplot.R") # multiplot function allows you to plot multiple plots(from a list) in one pane. This sources the function from the internet.


## weighted matrices
filenames <- c("QAS_L4_Oct17_CM.csv","QAS_4_Oct17_CM.csv","QA_700s_CM.csv")
data <- lapply(X = filenames, FUN = read.csv, header = T)
names(data) <- years

# Eigenfiles <- c("centr_L4.csv","centr_4.csv", "centr_700s.csv")
# centrality <- lapply(X = Eigenfiles, FUN = read.csv, header = FALSE)
# names(centrality) <- years
# vec_centrality <- list() #create a list that will be used to input melted centrality similarity data
# ## Enter subject numbers - note there are only scanned girls in this dataset
# for (y in years){
#   colnames(centrality[[y]]) <- scanned_list[[y]]
#   rownames(centrality[[y]]) <- colnames(centrality[[y]])
#   centrality[[y]][upper.tri(centrality[[y]],diag=TRUE)] <- NA
#   vec_centrality[[y]] <- reshape2::melt(centrality[[y]], na.rm=TRUE, value.name="centrality_difference")#melt the centrality data and label the column "centrality_difference"
# }

data <- lapply(data, function(year_data){
  year_data <- as.matrix(year_data[,as.character(year_data$ID)])#will include only those columns that match the rows of girls who completed the survey
  diag(year_data) <- 0
  year_data
})
## Threshold the data 
threshold <- 4
data <- lapply(data, function(thresh_data){
  thresh_data <- thresh_data>=threshold
  rownames(thresh_data) <- colnames(thresh_data)
  thresh_data <- thresh_data*1
  thresh_data
})
vec_comms <- list() #create a list that will be used to input melted community similarity data
## copy the subject numbers from the data list into the same_comms list
for (y in years){
  colnames(same_comms[[y]]) <- colnames(data[[y]])
  rownames(same_comms[[y]]) <- colnames(same_comms[[y]])
  same_comms[[y]] <- same_comms[[y]][scanned_list[[y]],scanned_list[[y]]]#includes only those participants who were scanned
  same_comms[[y]][upper.tri(same_comms[[y]],diag=TRUE)] <- NA #remove upper triangle and replace with NA
  vec_comms[[y]] <- reshape2::melt(same_comms[[y]], na.rm=TRUE, value.name="communities_match")#melt the community data and label the column "communities_match"
}


## Create graph from adjacency matrix
adjacency_data <- lapply(data, function(adjacency){
  adjacency <- graph_from_adjacency_matrix(adjacency, mode = "directed")
  adjacency
})

## create undirected graph from adjacency graph created in previous step
undirected <- list()
for (y in years){
undirected[[y]] <- as.undirected(adjacency_data[[y]], mode = c("collapse", "each", "mutual"),
                      edge.attr.comb = igraph_opt("edge.attr.comb"))}

#if you wish to do the following steps using an undirected graph instead of a directed graph,
#uncomment the code below:
adjacency_data <- undirected

## create social network metrics data file
SNstats  <- list()
SNrows <- sapply(data, nrow)
for (y in years) {SN  <- as.data.frame(matrix(ncol=0, nrow=SNrows[[y]]))
rownames(SN) <- rownames(data[[y]])
#SN$indegree <- degree(adjacency_data[[y]], v=rownames(SN), mode=c("in"),loops = TRUE, normalized = FALSE)
SN$eigen_centrality <- eigen_centrality(adjacency_data[[y]], directed=F)$vector
SNstats[[y]] <-  SN
}
## community structure analysis - now uses cluster_louvain method as it is much quicker
community_cluster <- list()
for (y in years) {
  community_cluster[[y]] <- cluster_louvain(adjacency_data[[y]],weights = NULL)
  #community_cluster[[y]] <- cluster_optimal(adjacency_data[[y]],weights = NULL)
}

communities <- list()
for (y in years){
  communities[[y]] <- data.frame(student=character(),community=c())
}

for (y in years){
  modules <- community_cluster[[y]]
  for (i in 1:length(modules))
    communities[[y]] <- rbind(communities[[y]], data.frame(student=c(modules[[i]]), community=rep(i,length(modules[[i]]))))
  comm <- communities[[y]]
  communities[[y]] <- comm[order(as.character(comm$student)),]
  SNstats[[y]]$community <- as.double(communities[[y]]$community)
}



for (y in years) {
  rbPal <- c("chartreuse","olivedrab1","darkolivegreen1","Khaki2","bisque2","ivory3","thistle3","lavenderblush4","darkslategrey", "gray19")
  set.seed(12)
  plot.igraph(adjacency_data[[y]], vertex.color=rbPal[(cut(SNstats[[y]]$eigen_centrality,breaks = 10))], edge.arrow.size = 0.35, vertex.label="", main=y)
 # SDMTools::legend.gradient(cbind(x =c(0.9,1,1,0.9), y =c(1,1,1,0.4)), rbPal,limits=c("Low","High"),title="Eigenvector centrality score")
  rbPal <- c("chartreuse","darkolivegreen1","thistle3","darkslategrey")
  set.seed(12)
  if (y<3){
  plot.igraph(adjacency_data[[y]], vertex.color=rbPal[(cut(SNstats[[y]]$community,breaks = 4))], edge.arrow.size = 0.15, vertex.label="", main=y)
  legend("topright", title='Community structure - Module', legend=levels(as.factor(SNstats[[y]]$community)),col = rbPal, bty = "n", pch=20 , pt.cex = 3, cex = 1, text.col=c('black') , horiz = FALSE, inset = c(0.001, 0.001))
  }
  else {plot.igraph(adjacency_data[[y]], vertex.color=rbPal[(cut(SNstats[[y]]$community,breaks = 3))], edge.arrow.size = 0.15, vertex.label="", main=y)
    legend("topright", title='Community structure - Module', legend=levels(as.factor(SNstats[[y]]$community)),col = rbPal , bty = "n", pch=20 , pt.cex = 3, cex = 1, text.col=c('black') , horiz = FALSE, inset = c(0.001, 0.001))
  }
}

## determine social distance between dyads
#geodesic distance - removes upper triangle and replaces with NA
dist <- list()
dist_scanned <- list()#will eventually contain information about path length between members of each dyad
for (y in years){dist[[y]] <- distances(adjacency_data[[y]])
dist_scanned[[y]] <- dist[[y]][scanned_list[[y]],scanned_list[[y]]]
dist_scanned[[y]][upper.tri(dist_scanned[[y]],diag=TRUE)] <- NA
}

melted_dist <- list()#will eventually contain vectors of distance data
for (y in years){
  melted_dist[[y]] <- reshape2::melt(dist_scanned[[y]], na.rm=TRUE, value.name="distance")#Melt distance matrices into vectors
}

## correlations between members of each dyad 
#note that I am keeping the year groups separate from now on. 
#To look at all groups together, I will use a meta-analysis instead

#name columns and rows of dataframes in networks list
networks <- lapply(networks, function(net_type){
  rownames(net_type)<- c("N201", "N204", "N210", "N211", "N213", "N216", "N217", "N218", "N219", "N220", "N222", "N223", "N224", "N226", "N227", "N232", "N233", "N240", "N244", "N245", "N246", "N247", "N256", "N300", "N302", "N305", "N307", "N308", "N309", "N312", "N315", "N317", "N318", "N321", "N322", "N325", "N332", "N333", "N334", "N335", "N340", "N344", "N346", "N347", "N350", "N352", "N355", "N361", "N366", "N367", "N369", "N702", "N704", "N706", "N708", "N709", "N711", "N713", "N715", "N720", "N722", "N737", "N739", "N742", "N743", "N745", "N748", "N751")
  colnames(net_type) <- rownames(net_type)
  net_type
})

## make a new list called yrs_nets and copy the networks list into it for each year group 
#(so now you have three lists with 5 networks lists inside those lists)
yrs_nets <- list()
for (y in years){
  yrs_nets[[y]] <- networks
}

## now reduce each of the network lists within each of the year group lists into only those girls who were scanned
for (y in years){
  yrs_nets[[y]] <- lapply(yrs_nets[[y]], function(net_type){
    net_type <- net_type[scanned_list[[y]], scanned_list[[y]]]
    net_type
  })
}

## create a list called net_corr_vectors and put the melted data from all year groups and all networks into this list
net_corr_vectors <- list()
for (y in years){
  yrs_nets[[y]] <- lapply(yrs_nets[[y]], function(removediag){
    removediag[upper.tri(removediag,diag=TRUE)] <- NA #make all values in the upper triangle NA
    removediag <- as.matrix(removediag) #turn into a matrix
    removediag #return removediag as new value for yrs_nets[[y]]
  })
  net_corr_vectors[[y]] <- lapply(yrs_nets[[y]], function(melt_it){
    melted <- reshape2::melt(melt_it, na.rm=TRUE, value.name="correlation")#melt data
    melted
  })
  for (n in nets){
    net_corr_vectors[[y]][[n]] <- cbind(net_corr_vectors[[y]][[n]], melted_dist[[y]]$distance,vec_comms[[y]]$communities_match)#, vec_centrality[[y]]$centrality_difference)#bind the correlation, distance, same_communities and centrality difference lists
    colnames(net_corr_vectors[[y]][[n]]) <- c("ppt_1","ppt_2","correlation","distance","communities_match")#,"centrality_difference")#label the columns
  }
  net_corr_vectors[[y]][["wholebrain"]] <- cbind(net_corr_vectors[[y]][["wholebrain"]],brain_modularity[[y]]$modularity_difference)#add the modularity and total strength data to only the whole brain dataframes
  colnames(net_corr_vectors[[y]][["wholebrain"]]) <- c("ppt_1","ppt_2","correlation","distance","communities_match","modularity_difference")
}
#now you have a list (net_corr_vectors) that contains information about the correlation strength, distance, community similarities and centrality difference (and brain modularity difference, pos and neg strength) between each dyad

## create an empty list of lists of lists to put LME output into
model_type <-  c("community_structure", "social_distance")#, "Eigenvector_centrality_difference")
LMEmodels <- list()
models_yrs <- list()
for (y in years){
  for (n in nets){
    models_yrs[[y]][[n]] <- list()
  }
}
for (t in model_type){
  for (y in years){
    LMEmodels[[t]] <- models_yrs
  }
}
plot <- LMEmodels

## Produce linear mixed effects output for each condition
for (y in years){
  for (n in nets){
    model <- lmer(correlation ~ communities_match +(1|ppt_1)+(1|ppt_2), net_corr_vectors[[y]][[n]])
    LMEmodels[[1]][[y]][[n]] <- model
  }
}
for (y in years){
  for (n in nets){
    model <- lmer(correlation ~ distance +(1|ppt_1)+(1|ppt_2), net_corr_vectors[[y]][[n]])
    LMEmodels[[2]][[y]][[n]] <- model
  }
}
# for (y in years){
#   for (n in nets){
#     model <- lmer(correlation ~ centrality_difference +(1|ppt_1)+(1|ppt_2), net_corr_vectors[[y]][[n]])
#     LMEmodels[[3]][[y]][[n]] <- model
#   }
# }
### connectomics models (LME)

for (y in years){
  model <- lmer(modularity_difference ~ communities_match +(1|ppt_1)+(1|ppt_2), net_corr_vectors[[y]][["wholebrain"]])
  LMEmodels[[1]][[y]][["brain_modularity"]] <- model
}
for (y in years){
  model <- lmer(modularity_difference ~ distance +(1|ppt_1)+(1|ppt_2), net_corr_vectors[[y]][["wholebrain"]])
  LMEmodels[[2]][[y]][["brain_modularity"]] <- model
}
# for (y in years){
#   model <- lmer(modularity_difference ~ centrality_difference +(1|ppt_1)+(1|ppt_2), net_corr_vectors[[y]][["wholebrain"]])
#   LMEmodels[[3]][[y]][["brain_modularity"]] <- model
# }


## make distance and communities match columns as.factor
vecs_as_factors <- net_corr_vectors
for (y in years){
  for (n in nets){
    vecs_as_factors[[y]][[n]]$distance <- as.factor(vecs_as_factors[[y]][[n]]$distance)
    vecs_as_factors[[y]][[n]]$communities_match <- as.factor(vecs_as_factors[[y]][[n]]$communities_match)
  }
}

## Get plots of social network vs brain network correlations
for (t in model_type){
  for (y in years){
    for (n in nets){
      plot[[1]][[y]][[n]] <- ggplot(vecs_as_factors[[y]][[n]], aes(x=communities_match, y=correlation))+geom_boxplot()+labs(title=paste(y,n,"community structure"),x="Same module?", y = "Pearson's correlation between dyads")
      plot[[2]][[y]][[n]] <- ggplot(vecs_as_factors[[y]][[n]], aes(x=distance, y=correlation))+geom_boxplot()+labs(title=paste(y,n,"Social distance"),x="Social distance", y = "Pearson's correlation between dyads")
      #plot[[3]][[y]][[n]]<- ggscatter(vecs_as_factors[[y]][[n]], x = "centrality_difference", y = "correlation", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab=("centrality difference"),ylab=paste(y,n,"between-dyad correlation strength"))
    }
  }
}
for (y in years){
  #modularity difference plots
  plot[[1]][[y]][["brain_modularity"]] <- ggplot(vecs_as_factors[[y]][["wholebrain"]], aes(x=communities_match, y=modularity_difference))+geom_boxplot()+labs(title=paste(y,"brain modularity vs community structure"),x="Same module?", y = "Pearson's correlation between dyads")
  plot[[2]][[y]][["brain_modularity"]] <- ggplot(vecs_as_factors[[y]][["wholebrain"]], aes(x=distance, y=modularity_difference))+geom_boxplot()+labs(title=paste(y,"brain modularity vs Social distance"),x="Social distance", y = "Pearson's correlation between dyads")
  #plot[[3]][[y]][["brain_modularity"]]<- ggscatter(vecs_as_factors[[y]][["wholebrain"]], x = "centrality_difference", y = "modularity_difference", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab=("centrality difference"),ylab=paste(y,"brain modularity similarity"))
}

# here a subset of the plots are drawn in one pane with three cols. The matrix arg allows for more specification of plot order.
# for (t in model_type){
#   for (y in years){
#     multiplot(plotlist = plot[[t]][[y]],cols = 2)
#   }
# }


# save the relevant lists into an .RData file for easy loading
save(list=c("plot", "LMEmodels","net_corr_vectors"),file = 'shiny_lists.RData')


#meta-analysis outputs
library(devtools)
install_github("https://github.com/LilyFG/metaforlmer.git")

library(metaforlmer)
library(metafor)

MAplotlist <- list()
MetaAnalysis <- list()
analyses <- names(LMEmodels[[1]][[1]])
for (t in model_type){
  MetaAnalysis[[t]] <- list()
  for (a in analyses){
    model_list <- list(LMEmodels[[t]][["Cohort 1"]][[a]],LMEmodels[[t]][["Cohort 3"]][[a]],LMEmodels[[t]][["Cohort 2"]][[a]])
    names(model_list) <- years
    MetaAnalysis[[t]][[a]] <- meta_models(model_list = model_list)
    # MAplotlist[[t]][[a]] <- ggforest(MetaAnalysis)
    MAplotlist[[t]][[a]] <- metaforlmer::ggforest(MetaAnalysis[[t]][[a]], labels = paste(t,a))
    #MAplotlist[[t]][[a]] <- ggforest(data=MetaAnalysis[[t]][[a]], main = paste(t,a))#creates errors on my laptop
  }
}

whichPlot <- list(plot,MAplotlist)
names(whichPlot) <- c("plot","MAplotlist")
# save the relevant lists into an .RData file for easy loading
save(list=c("LMEmodels","whichPlot"),file = 'shiny_lists.RData')

##################################################

#making plots for CNS conference poster
rbPal <- c("chartreuse","olivedrab1","yellow","thistle3","darkslategrey")
yrs <- c("Cohort 1","Cohort 2","Cohort 3") #have changed the order of this as well as the order of LMEmodels below
MAplotlist <- list()
#plotforest <- list()
MetaAnalysis <- list()
colour=1
analyses <- names(LMEmodels[[1]][[1]])
for (t in model_type){
  MetaAnalysis[[t]] <- list()
  for (a in analyses){
    model_list <- list(LMEmodels[[t]][["Cohort 1"]][[a]],LMEmodels[[t]][["Cohort 2"]][[a]],LMEmodels[[t]][["Cohort 3"]][[a]])
    names(model_list) <- yrs
    MetaAnalysis[[t]][[a]] <- meta_models(model_list = model_list)
    # MAplotlist[[t]][[a]] <- ggforest(MetaAnalysis)
    #MAplotlist[[t]][[a]] <- ggforest(data=MetaAnalysis[[t]][[a]], main = paste(a))
    MAplotlist[[t]][[a]] <- metaforlmer::ggforest(MetaAnalysis[[t]][[a]], labels = paste(a), palette = rbPal[colour])
    colour <- colour+1
  if (colour>12){colour <- 1}
    }
}

plotlist_community <- list(MAplotlist[[1]][[1]],MAplotlist[[1]][[2]],MAplotlist[[1]][[3]],MAplotlist[[1]][[4]],MAplotlist[[1]][[5]])
multiplot(plotlist = plotlist_community, cols = 5)

plotlist_socialdist <- list(MAplotlist[[2]][[1]],MAplotlist[[2]][[2]],MAplotlist[[2]][[3]],MAplotlist[[2]][[4]],MAplotlist[[2]][[5]])
multiplot(plotlist = plotlist_socialdist, cols = 5)
###########################
setwd("F:\\0_parcellation_analysis/scripts-data-sharing/")
write.csv(net_corr_vectors[["Cohort 1"]][["wholebrain"]]$distance, file = "cohort1_distance.csv", row.names = F)
write.csv(net_corr_vectors[["Cohort 2"]][["wholebrain"]]$distance, file = "cohort2_distance.csv", row.names = F)
write.csv(net_corr_vectors[["Cohort 3"]][["wholebrain"]]$distance, file = "cohort3_distance.csv", row.names = F)
write.csv(net_corr_vectors[["Cohort 1"]][["wholebrain"]]$communities_match, file = "cohort1_comm.csv", row.names = F)
write.csv(net_corr_vectors[["Cohort 2"]][["wholebrain"]]$communities_match, file = "cohort2_comm.csv", row.names = F)
write.csv(net_corr_vectors[["Cohort 3"]][["wholebrain"]]$communities_match, file = "cohort3_comm.csv", row.names = F)
#these need to be edited to remove first column, then saved for later use in elastic lists file
#save as yL4_dist, y700_dist, y4_dist and yL4_comm, y700_comm and y4_comm, respectively
