#CAROLYN MCNABB 2019 community data analysis 

setwd("F:\\0_parcellation_analysis/scripts-data-sharing/")

# library(BBmisc)
# library(igraph)
# library(psych)
library(pacman)
p_load(BBmisc)
p_load(igraph)
p_load(psych)
p_load(SDMTools)

#source("motivationscript.R")


# weighted matrices

QA_L4 <- read.csv("QAS_L4_Oct17_CM.csv", header = T)#this has been updated to include the latest data
# Somehow L4 has the rows of zeros, which I exclude
#QA_L4 <- QA_L4[QA_L4[, 2] != 0, ]
QA_L4 <- QA_L4[, 1:60]
QA_L4 <- as.matrix(QA_L4[, as.character(QA_L4$ID)])


QA_4 <- read.csv("QAS_4_Oct17_CM.csv", header = T)#this has been updated to include the latest data
QA_4 <- QA_4[, 1:66]
QA_4 <- as.matrix(QA_4[, as.character(QA_4$ID)])

QA_700 <- read.csv("QA_700s_CM.csv", header = T)#this has been updated to include the latest data
QA_700 <- QA_700[, 1:50]
QA_700 <- as.matrix(QA_700[, as.character(QA_700$ID)])

# U5 is for now excluded
# check diagonals --- need doublecheck
years <- c("L4", "4", "700")#, "U4", "L5", "L6", "U6")
Qyears <- c("QA_L4", "QA_4", "QA_700")#, "QA_L5", "QA_L6", "QA_U6")

for (g in Qyears){
print(g)
eval(parse(text=paste("print(diag(", g, "))", sep = "")))
}

# all diagonals = 0
for (g in Qyears){
print(paste(g, " trasformation of diagonals done", sep = ""))
eval(parse(text=paste("diag(", g, ")<-0", sep = "")))
}


# thresholding
# 5 point scale. Equal or more than the threshold is 1.
# add "b" at the end of the object --- like QA_L4b
threshold <- 4
for (g in Qyears){
print(paste(g, " thresholding done", sep = ""))
eval(parse(text=paste(g, "b <-", g, ">=", threshold, sep = "")))
}


dim(QA_L4b)
dim(QA_4b)
dim(QA_700b)


rownames(QA_L4b) <- colnames(QA_L4b)
rownames(QA_4b) <- colnames(QA_4b)
rownames(QA_700b) <- colnames(QA_700b)

IDL4  <- as.factor(colnames(QA_L4b))
ID4  <- as.factor(colnames(QA_4b))
ID700  <- as.factor(colnames(QA_700b))

## Final data used

#QA_dfL4s <-  QA_dfL4[is.element(QA_dfL4$ID, IDL4), ]
QA_SocL4s <- QA_L4b[IDL4, IDL4]
binary_L4_thresh <- QA_SocL4s*1
write.csv(binary_L4_thresh, row.names = TRUE, file = "binary_L4_social_net.csv")


#QA_df4s <-  QA_df4[is.element(QA_df4$ID, ID4), ]
QA_Soc4s <- QA_4b[ID4, ID4]
binary_4_thresh <- QA_Soc4s*1
write.csv(binary_4_thresh, row.names = TRUE, file = "binary_4_social_net.csv")


QA_Soc700s <- QA_700b[ID700, ID700]
binary_700_thresh <- QA_Soc700s*1
write.csv(binary_700_thresh, row.names = TRUE, file = "binary_700_social_net.csv")

## Graph object - SELECT BASED ON WHICH VERSION OF R YOU ARE RUNNING
SocL4 <- graph_from_adjacency_matrix(QA_SocL4s, mode = "directed")
#SocL4 <- graph.adjacency(QA_SocL4s, mode = "directed")
Soc4 <- graph_from_adjacency_matrix(QA_Soc4s, mode = "directed")
#Soc4 <- graph.adjacency(QA_Soc4s, mode = "directed")
Soc700 <- graph_from_adjacency_matrix(QA_Soc700s, mode = "directed")


centr_4  <- as.data.frame(matrix(ncol=0, nrow=65))
rownames(centr_4) <- ID4
centr_L4  <- as.data.frame(matrix(ncol=0, nrow=59))
rownames(centr_L4) <- IDL4
centr_700  <- as.data.frame(matrix(ncol=0, nrow=49))
rownames(centr_700) <- ID700

#plot(Soc4)
centr_4$indegree <- degree(Soc4, v=ID4, mode=c("in"),loops = TRUE, normalized = FALSE)
centr_4$betweenness <- betweenness(Soc4)
centr_4$eccentricity <- eccentricity(Soc4)
centr_4$closeness <- closeness(Soc4)
centr_4$eigen_centrality <- eigen_centrality(Soc4, directed=TRUE)$vector#this was not directed before - have now changed
centr_4$transitivity <- transitivity(Soc4, "local")
centr_4$power <- power_centrality(Soc4, loops=FALSE, exponent=0.9, rescale=FALSE)

#plot(SocL4)
centr_L4$indegree <- degree(SocL4, v=IDL4, mode=c("in"),loops = TRUE, normalized = FALSE)
centr_L4$betweenness <- betweenness(SocL4)
centr_L4$eccentricity <- eccentricity(SocL4)
centr_L4$closeness <- closeness(SocL4)
centr_L4$eigen_centrality <- eigen_centrality(SocL4, directed=TRUE)$vector#this was not directed before - have now changed
centr_L4$transitivity <- transitivity(SocL4, "local")
centr_L4$power <- power_centrality(SocL4, loops=FALSE, exponent=0.9, rescale=FALSE)

#plot(Soc700)
centr_700$indegree <- degree(Soc700, v=ID700, mode=c("in"),loops = TRUE, normalized = FALSE)
centr_700$betweenness <- betweenness(Soc700)
centr_700$eccentricity <- eccentricity(Soc700)
centr_700$closeness <- closeness(Soc700)
centr_700$eigen_centrality <- eigen_centrality(Soc700, directed=TRUE)$vector#this was not directed before - have now changed
centr_700$transitivity <- transitivity(Soc700, "local")
centr_700$power <- power_centrality(Soc700, loops=FALSE, exponent=0.9, rescale=FALSE)

################################
#for subanalysis in only those who were scanned... get centrality scores for only scanned girls
rownames_L4 <- c("N201", "N204", "N210", "N211", "N213", "N216", "N217", "N218", "N219", "N220", "N222", "N223", "N224", "N226", "N227", "N232", "N233", "N240", "N244", "N245", "N246", "N247", "N256")
rownames_4 <- c("N300", "N302", "N305", "N307", "N308", "N309", "N312", "N315", "N317", "N318", "N321", "N322", "N325", "N332", "N333", "N334", "N335", "N340", "N344", "N346", "N347", "N350", "N352", "N355", "N361", "N366", "N367", "N369")
rownames_700s <- c("N702", "N704", "N706", "N708", "N709", "N711", "N713", "N715", "N720", "N722", "N737", "N739", "N742", "N743", "N745", "N748", "N751")
centr_L4_scanned <- as.data.frame(centr_L4[rownames_L4,"eigen_centrality"])
rownames(centr_L4_scanned) <- rownames_L4
colnames(centr_L4_scanned) <- c("centrality")
centr_4_scanned <- as.data.frame(centr_4[rownames_4,"eigen_centrality"])
rownames(centr_4_scanned) <- rownames_4
colnames(centr_4_scanned) <- c("centrality")
centr_700s_scanned <- as.data.frame(centr_700[rownames_700s,"eigen_centrality"])
rownames(centr_700s_scanned) <- rownames_700s
colnames(centr_700s_scanned) <- c("centrality")
scanned_both <- rbind(centr_L4_scanned,centr_4_scanned,centr_700s_scanned)
##uncomment the four lines below when you want to make these files for the first time
# write.csv(centr_4_scanned,file="centrality_4.csv",row.names=TRUE)
# write.csv(centr_L4_scanned,file="centrality_L4.csv",row.names=TRUE)
# write.csv(centr_700s_scanned,file="centrality_700s.csv",row.names=TRUE)
# write.csv(scanned_both,file="centrality_all.csv",row.names=TRUE)
#after this, remove "N"s from the first column in excel. then use the centrality_similarities.m script in Matlab to get the difference between dyads
#output from Matlab is centr_L4.csv, centr_4.csv and centr_all.csv
#use these as input in the R_correlation_analysis_CM.R script
################################
#indegree only

indegree_L4 <- as.data.frame(centr_L4[rownames_L4,"indegree"])
rownames(indegree_L4) <- rownames_L4
colnames(indegree_L4) <- c("indegree")
indegree_4 <- as.data.frame(centr_4[rownames_4,"indegree"])
rownames(indegree_4) <- rownames_4
colnames(indegree_4) <- c("indegree")
indegree_700 <- as.data.frame(centr_700[rownames_700s,"indegree"])
rownames(indegree_700) <- rownames_700s
colnames(indegree_700) <- c("indegree")
indegree_scanned <- rbind(indegree_L4,indegree_4,indegree_700)
################################

p_load(wesanderson)
rbPal <- colorRampPalette(c(wes_palette("Zissou1", 100, type = c("continuous"))))
# rbPal <- colorRampPalette(c(wes_palette("Darjeeling1", 100, type = c("continuous"))))
set.seed(1)
plot.igraph(Soc4, vertex.color=rbPal(10)[(cut(centr_4$eigen_centrality,breaks = 10))], edge.arrow.size = 0.35, vertex.label="", main="Year 4/U4 (2017) Eigenvector centrality")
#SDMTools::legend.gradient(cbind(x =c(0.8,0.9,1.9,2), y =c(1,1,0.3,0.3)), rbPal(100),c("Low","High"),title="Centrality score")
set.seed(1)
plot.igraph(SocL4, vertex.color=rbPal(10)[(cut(centr_L4$eigen_centrality,breaks = 10))], edge.arrow.size = 0.35, vertex.label="", main="Year L4/4 (2017) Eigenvector centrality")
#SDMTools::legend.gradient(cbind(x =c(0.8,0.9,1.9,2), y =c(1,1,0.3,0.3)), rbPal(100),c("Low","High"),title="Centrality score")
set.seed(1)
plot.igraph(Soc700, vertex.color=rbPal(10)[(cut(centr_700$eigen_centrality,breaks = 10))], edge.arrow.size = 0.35, vertex.label="", main="Year L4/4 (2018) Eigenvector centrality")
#SDMTools::legend.gradient(cbind(x =c(0.8,0.9,1.9,2), y =c(1,1,0.3,0.3)), rbPal(100),c("Low","High"),title="Centrality score")
############
#Plot indegree graphs
set.seed(1)
plot.igraph(Soc4, vertex.color=rbPal(10)[(cut(centr_4$indegree,breaks = 10))], edge.arrow.size = 0.35, vertex.label="", main="Year 4/U4 (2017) Indegree centrality")
#SDMTools::legend.gradient(cbind(x =c(0.8,0.9,1.9,2), y =c(1,1,0.3,0.3)), rbPal(100),c("Low","High"),title="Centrality score")
set.seed(1)
plot.igraph(SocL4, vertex.color=rbPal(10)[(cut(centr_L4$indegree,breaks = 10))], edge.arrow.size = 0.35, vertex.label="", main="Year L4/4 (2017) Indegree centrality")
#SDMTools::legend.gradient(cbind(x =c(0.8,0.9,1.9,2), y =c(1,1,0.3,0.3)), rbPal(100),c("Low","High"),title="Centrality score")
set.seed(1)
plot.igraph(Soc700, vertex.color=rbPal(10)[(cut(centr_700$indegree,breaks = 10))], edge.arrow.size = 0.35, vertex.label="", main="Year L4/4 (2018) Indegree centrality")
#SDMTools::legend.gradient(cbind(x =c(0.8,0.9,1.9,2), y =c(1,1,0.3,0.3)), rbPal(100),c("Low","High"),title="Centrality score")
##############


#GLM_L4_thresh5 <- centr_L4[c("eigen_centrality","degree", "betweenness", "transitivity", "power")]
#GLM_4_thresh5 <- centr_4[c("eigen_centrality", "degree", "betweenness", "transitivity", "power")]
# #GLM_L4_and_4_thresh5 <- rbind(GLM_L4_thresh5, GLM_4_thresh5)
# GLM_L4_thresh4 <- centr_L4[c("eigen_centrality", "betweenness", "transitivity")]
# GLM_4_thresh4 <- centr_4[c("eigen_centrality", "betweenness", "transitivity")]
# GLM_L4_and_4_thresh4 <- rbind(GLM_L4_thresh4, GLM_4_thresh4)
# 
# #write.csv(GLM_L4_and_4_thresh4, file = "centrality_measures_threshold_4.csv", row.names = TRUE)
# #write.csv(GLM_L4_and_4_thresh5, file = "centrality_measures_threshold_5.csv", row.names = TRUE)
# 
# 
# scanned <- GLM_L4_and_4_thresh4[c("N201", "N204", "N210", "N211", "N213", "N216", "N217", "N218", "N219", "N220", "N222", "N223", "N224", "N226", "N227", "N232", "N233", "N240", "N244", "N245", "N246", "N247", "N256", "N300", "N302", "N305", "N307", "N308", "N309", "N312", "N315", "N317", "N318", "N321", "N322", "N325", "N332", "N333", "N334", "N335", "N340", "N344", "N346", "N347", "N350", "N352", "N355", "N361", "N366", "N367", "N369"),]
# handedness <-  read.csv("QAS_handedness.csv", header = T)
# scanned <- cbind(scanned,handedness[2])
# 
#have changed to cluster louvain method to match mutual ties analysis
#NB THESE USED TO READ cluster_louvain(SocL4,weights = NULL) but had to be updated to work (2021)
communityL4 <- cluster_louvain(as.undirected(SocL4),weights = NULL)
#communityL4 <- cluster_optimal(SocL4,weights = NULL)
modularity_L4 <- modularity(communityL4)

community4 <- cluster_louvain(as.undirected(Soc4),weights = NULL)
#community4 <- cluster_optimal(Soc4,weights = NULL)
modularity_4 <- modularity(community4)

community700 <- cluster_louvain(as.undirected(Soc700),weights = NULL)
# community700 <- cluster_optimal(Soc700,weights = NULL)
modularity_700 <- modularity(community700)

char_path_length_L4 <- mean_distance(SocL4, directed = TRUE, unconnected = TRUE)
char_path_length_4 <- mean_distance(Soc4, directed = TRUE, unconnected = TRUE)
char_path_length_700 <- mean_distance(Soc700, directed = TRUE, unconnected = TRUE)
##############

comm4 <- data.frame(student=character(),community4=c())
for (i in 1:length(community4))
  comm4 <- rbind(comm4, data.frame(student=c(community4[[i]]), community4=rep(i,length(community4[[i]]))))
test1 <- comm4[order(as.character(comm4$student)),]
commL4 <- data.frame(student=character(),communityL4=c())
for (i in 1:length(communityL4))
  commL4 <- rbind(commL4, data.frame(student=c(communityL4[[i]]), communityL4=rep(i,length(communityL4[[i]]))))
test2 <- commL4[order(as.character(commL4$student)),]
comm700 <- data.frame(student=character(),community700=c())
for (i in 1:length(community700))
  comm700 <- rbind(comm700, data.frame(student=c(community700[[i]]), community700=rep(i,length(community700[[i]]))))
test3 <- comm700[order(as.character(comm700$student)),]

##Plot communities for each year group
# rbPal <- colorRampPalette(c(wes_palette("GrandBudapest2", 4, type = c("discrete"))))
# rbPal <- colorRampPalette(c(wes_palette("FantasticFox1", 4, type = c("discrete"))))
rbPal <- colorRampPalette(c(wes_palette("Darjeeling1", 4, type = c("discrete"))))
# rbPal <- colorRampPalette(c(wes_palette("Zissou1", 5, type = c("discrete"))))
set.seed(1)
plot.igraph(Soc4, vertex.color=rbPal(4)[(cut(test1$community4,breaks = 4))], edge.arrow.size = 0.15, vertex.label="", main="Year 4/U4 (2017) community structure")
legend("topright", title='Module', legend=levels(as.factor(test1$community4)),col = rbPal(4) , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col=c('grey') , horiz = FALSE, inset = c(0.001, 0.001))
set.seed(1)
plot.igraph(SocL4, vertex.color=rbPal(5)[(cut(test2$communityL4,breaks = 5))], edge.arrow.size = 0.15, vertex.label="", main="Year L4/4 (2017) community structure")
legend("topright", title='Module', legend=levels(as.factor(test2$communityL4)),col = rbPal(4) , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col=c('grey') , horiz = FALSE, inset = c(0.001, 0.001))
set.seed(1)
plot.igraph(Soc700, vertex.color=rbPal(5)[(cut(test3$community700,breaks = 5))], edge.arrow.size = 0.15, vertex.label="", main="Year L4/4 (2018) community structure")
legend("topright", title='Module',legend=levels(as.factor(test3$community700)),col = rbPal(4) , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col=c('grey') , horiz = FALSE, inset = c(0.001, 0.001))

################
write.csv(commL4,file='sorted_communities_L4.csv')
write.csv(comm4,file='sorted_communities_4.csv')
write.csv(comm700,file='sorted_communities_700.csv')
#now go into Matlab and create community similarity files!!!
#First you will need to remove the first column and first row and the N from the ppt numbers
