

# # keep the complete full tree ("RAxML_parsimonyTree.complete100") a separate variable since I'm using it for all the analyses below
# # I call it parseTreeComp
# parseTreeComp <- read.tree("RAxML_parsimonyTree.complete100")

# For Robinson-Foulds values, we can use foreach loop to iterate through each tree (There were 10 trees)

#I compare these 10 bias trees with the fully random complete full tree ("RAxML_bestTree.mgTree" in large dataset directory).  
# I call it fullyRandomparseTreeComp, a separate variable since I'm using it for all the analyses below.
# fullyRandomparseTreeComp <- read.tree("RAxML_bestTree.mgTree")

# fullyRandom20_1 <- read.tree("F20_1_NEW")
# fullyRandom40_1 <- read.tree("F40_1_NEW")
# fullyRandom60_1 <- read.tree("F60_1_NEW")
# fullyRandom80_1 <- read.tree("F80_1_NEW")

files <- list.files(path="/home/thanu/Desktop/FishData/large_dataset", pattern="*_[0-9]{1,2}_NEW", full.names=TRUE, recursive=FALSE)
treeList <- foreach(i=1:length(files)) %do% read.tree(files[i])

# Now I have a list with each element being its own separate tree, so treeList[[1]] is my first tree 20_1_new, treeList[[2]] is 20_2_new etc.

# Also we can name each list element by substituting out our path from the names of each file
names(treeList) <- gsub("/home/thanu/Desktop/FishData/large_dataset/", "", files)

# This will let us keep track of which list element corresponds to which tree file
names(treeList) 

# we can see the names of the files that correspond to each list element and we can reference any list element by name, for example typing:
#treeList$`20_1_new`

# Then I keep the complete full tree ("RAxML_parsimonyTree.complete100") a separate variable since I'm using it for all the analyses below
# I call it parseTreeComp
parseTreeComp <- read.tree("RAxML_bestTree.mgTree")

# For Robinson-Foulds values, we can use foreach loop to iterate through each tree (There were 40 trees)
RF_List <- foreach(i=1:length(treeList)) %do% RF.dist(parseTreeComp,treeList[[i]],normalize = TRUE)

# Same thing for path distances 
P_List <- foreach(i=1:length(treeList)) %do% path.dist(parseTreeComp,treeList[[i]])


metric <- function(L_List){
  #reordering
  L_List <- L_List[order]
  #group the values and assign values to variable names
  L_20 <- unlist(L_List[1:10])
  L_40 <- unlist(L_List[11:20])
  L_60 <- unlist(L_List[21:30])
  L_80 <- unlist(L_List[31:40])
  #generate a dataframe using backbone % and metric values
  results <- data.frame(L_20,L_40,L_60,L_80)
}

# Now we can use the above function to generate dataframes in different metrics
# First, Robinson-Foulds metric,
RF_metric<- metric(RF_List)
RF_random_stack <- stack(RF_metric)

P_metric<- metric(P_List)
P_random_stack <- stack(P_metric)

#for stratified

files <- list.files(path="/home/thanu/Desktop/FishData/Stratified", pattern="*_[0-9]{1,2}_new", full.names=TRUE, recursive=FALSE)

treeList <- foreach(i=1:length(files)) %do% read.tree(files[i])

# Now I have a list with each element being its own separate tree, so treeList[[1]] is my first tree 20_1_new, treeList[[2]] is 20_2_new etc.

# Also we can name each list element by substituting out our path from the names of each file

names(treeList) <- gsub("/home/thanu/Desktop/FishData/Stratified/", "", files)

# This will let us keep track of which list element corresponds to which tree file

names(treeList)

# we can see the names of the files that correspond to each list element and we can reference any list element by name, for example typing:

treeList$`S20_1_new`

# Then I keep the complete full tree ("RAxML_bestTree.complete100") a separate variable since I'm using it for all the analyses below

# I call it parseTreeComp

parseTreeComp_strat <- read.tree("RAxML_bestTree.StratTree")

# For Robinson-Foulds values, we can use foreach loop to iterate through each tree (There were 40 trees)

RF_List_S <- foreach(i=1:length(treeList)) %do% RF.dist(parseTreeComp_strat,treeList[[i]],normalize = TRUE)

# Same thing for path distances

P_List_S <- foreach(i=1:length(treeList)) %do% path.dist(parseTreeComp_strat,treeList[[i]])

RF_metric_S<- metric(RF_List_S)

P_metric_S<- metric(P_List_S)


# For Robinson-Foulds values, we can use foreach loop to iterate through each tree 
# RF_List_fullyRandomTree <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandomparseTreeComp,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_fullyRandomTree


# RF_List_RandomTree20 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom20_1,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree40 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom40_1,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree60 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom60_1,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree80 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyRandom80_1,treeList[[i]],normalize = TRUE,check.labels = FALSE)

# #playing around-----
# RF_List_RandomTree20 <- foreach(i=1:length(treeList)) %do% RF.dist(cF,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree40 <- foreach(i=1:length(treeList)) %do% RF.dist(mF,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree60 <- foreach(i=1:length(treeList)) %do% RF.dist(wF,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree80 <- foreach(i=1:length(treeList)) %do% RF.dist(ggF,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree20 <- RF.dist(aF,fullyRandomparseTreeComp,normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree40 <- RF.dist(lF,fullyRandomparseTreeComp,normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree60 <- RF.dist(vF,fullyRandomparseTreeComp,normalize = TRUE,check.labels = FALSE)
# RF_List_RandomTree80 <- RF.dist(ffF,fullyRandomparseTreeComp,normalize = TRUE,check.labels = FALSE)

#unlist the list for downstream analysis
# Unlist_RF_List_fullyRandomTree <- unlist(RF_List_fullyRandomTree)

# Unlist_RF_List_fullyRandomTree20 <- unlist(RF_List_RandomTree20)
# Unlist_RF_List_fullyRandomTree40 <- unlist(RF_List_RandomTree40)
# Unlist_RF_List_fullyRandomTree60 <- unlist(RF_List_RandomTree60)
# Unlist_RF_List_fullyRandomTree80 <- unlist(RF_List_RandomTree80)

# Same thing for path distances 
# P_List_fullyRandomTree <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandomparseTreeComp,treeList[[i]],check.labels = FALSE)


# p_List_RandomTree20 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom20_1,treeList[[i]],check.labels = FALSE)
# p_List_RandomTree40 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom40_1,treeList[[i]],check.labels = FALSE)
# p_List_RandomTree60 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom60_1,treeList[[i]],check.labels = FALSE)
# p_List_RandomTree80 <- foreach(i=1:length(treeList)) %do% path.dist(fullyRandom80_1,treeList[[i]],check.labels = FALSE)

# #playing around2---
# 
# p_List_RandomTree20 <- path.dist(fullyRandom20_1,fullyRandomparseTreeComp,check.labels = FALSE)
# p_List_RandomTree40 <- path.dist(fullyRandom40_1,fullyRandomparseTreeComp,check.labels = FALSE)
# p_List_RandomTree60 <- path.dist(fullyRandom60_1,fullyRandomparseTreeComp,check.labels = FALSE)
# p_List_RandomTree80 <- path.dist(fullyRandom80_1,fullyRandomparseTreeComp,check.labels = FALSE)


# #unlist the list for downstream analysis
# Unlist_P_List_fullyRandomTree<- unlist(P_List_fullyRandomTree)
# 
# Unlist_P_List_fullyRandomTree20<- unlist(p_List_RandomTree20)
# Unlist_P_List_fullyRandomTree40<- unlist(p_List_RandomTree40)
# Unlist_P_List_fullyRandomTree60<- unlist(p_List_RandomTree60)
# Unlist_P_List_fullyRandomTree80<- unlist(p_List_RandomTree80)
# # Then I compare these 10 bias trees with the stratified complete full tree ("RAxML_bestTree.StratTree").  
# # I call it fullyRandomparseTreeComp, a separate variable since I'm using it for all the analyses below
# #StratifiedparseTreeComp <- read.tree("RAxML_bestTree.StratTree")
# 
# fullyStrati20_1 <- read.tree("S20_1_new")
# fullyStrati40_1 <- read.tree("S40_1_new")
# fullyStrati60_1 <- read.tree("S60_1_new")
# fullyStrati80_1 <- read.tree("S80_1_new")
# 
# # For Robinson-Foulds values, we can use foreach loop to iterate through each tree (There were 40 trees)
# # RF_List_stratifiedTree <- foreach(i=1:length(treeList)) %do% RF.dist(StratifiedparseTreeComp,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# # RF_List_stratifiedTree 
# 
# RF_List_stratifiedTree20 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati20_1,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_stratifiedTree40 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati40_1,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_stratifiedTree60 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati60_1,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# RF_List_stratifiedTree80 <- foreach(i=1:length(treeList)) %do% RF.dist(fullyStrati80_1,treeList[[i]],normalize = TRUE,check.labels = FALSE)
# 
# 
# #unlist the list for downstream analysis
# # Unlist_RF_List_stratifiedTree<- unlist(RF_List_stratifiedTree)
# 
# Unlist_RF_List_stratifiedTree20<- unlist(RF_List_stratifiedTree20)
# Unlist_RF_List_stratifiedTree40<- unlist(RF_List_stratifiedTree40)
# Unlist_RF_List_stratifiedTree60<- unlist(RF_List_stratifiedTree60)
# Unlist_RF_List_stratifiedTree80<- unlist(RF_List_stratifiedTree80)
# 
# 
# # Same thing for path distances 
# # P_List_stratifiedTree <- foreach(i=1:length(treeList)) %do% path.dist(StratifiedparseTreeComp,treeList[[i]],check.labels = FALSE)
# 
# 
# P_List_stratifiedTree20 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati20_1,treeList[[i]],check.labels = FALSE)
# P_List_stratifiedTree40 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati40_1,treeList[[i]],check.labels = FALSE)
# P_List_stratifiedTree60 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati60_1,treeList[[i]],check.labels = FALSE)
# P_List_stratifiedTree80 <- foreach(i=1:length(treeList)) %do% path.dist(fullyStrati80_1,treeList[[i]],check.labels = FALSE)
# 
# 
# 
# #unlist the list for downstream analysis
# # Unlist_P_List_stratifiedTree<- unlist(P_List_stratifiedTree)
# 
# Unlist_P_List_stratifiedTree20<- unlist(P_List_stratifiedTree20)
# Unlist_P_List_stratifiedTree40<- unlist(P_List_stratifiedTree40)
# Unlist_P_List_stratifiedTree60<- unlist(P_List_stratifiedTree60)
# Unlist_P_List_stratifiedTree80<- unlist(P_List_stratifiedTree80)


# Read trees generated from RAxML into R 
# read all of the trees at once ensuring we have the correct path to our small_dataset file 
# So for me the path is /home/thanu/Desktop/FishData/missing_pattern/NewBiasPhylipfiles_Oct/ and the pattern is a regular expression to extract out the files we need
# Then I'm using foreach package to just iterate through each file
#install.packages("foreach")
library(foreach)
library(phangorn)
library(ggplot2)
files_bias <- list.files(path="/home/thanu/Desktop/FishData/missing_pattern/NewBiasPhylipfiles_Oct", pattern="*bestTree.bs[0-9]{1,2}", full.names=TRUE, recursive=FALSE)
treeList_bias <- foreach(i=1:length(files_bias)) %do% read.tree(files_bias[i])

# Now I have a list with each element being its own separate tree, so treeList[[1]] is my first tree 20_1_new, treeList[[2]] is 20_2_new etc.

# Also we can name each list element by substituting out our path from the names of each file
names(treeList_bias) <- gsub("/home/thanu/Desktop/FishData/missing_pattern/NewBiasPhylipfiles_Oct/", "", files_bias)

#I compare these 10 bias trees with the complete bias tree ("RAxML_bestTree.biasTree").  
# I call it fullyRandomparseTreeComp, a separate variable since I'm using it for all the analyses below
# biasTree <- read.tree("RAxML_bestTree.biasTree")

# For Robinson-Foulds values, we can use foreach loop to iterate through each tree 
RF_List_biasTree <- foreach(i=1:length(treeList_bias)) %do% RF.dist(parseTreeComp,treeList_bias[[i]],normalize = TRUE,check.labels = FALSE)
RF_List_biasTree <- as.data.frame(t(data.frame(RF_List_biasTree)))


# Same thing for path distances 
P_List_biasTree <- foreach(i=1:length(treeList_bias)) %do% path.dist(parseTreeComp,treeList_bias[[i]],check.labels = FALSE)
P_List_biasTree <- as.data.frame(t(data.frame(P_List_biasTree)))

# For Robinson-Foulds values, we can use foreach loop to iterate through each tree 
RF_List_biasTree_S <- foreach(i=1:length(treeList_bias)) %do% RF.dist(parseTreeComp_strat,treeList_bias[[i]],normalize = TRUE,check.labels = FALSE)
RF_List_biasTree_S <- as.data.frame(t(data.frame(RF_List_biasTree_S)))

# Same thing for path distances 
P_List_biasTree_S <- foreach(i=1:length(treeList_bias)) %do% path.dist(parseTreeComp_strat,treeList_bias[[i]],check.labels = FALSE)
P_List_biasTree_S <- as.data.frame(t(data.frame(P_List_biasTree_S)))


RF_merge <-cbind(RF_metric, RF_metric_S,RF_List_biasTree$V1,RF_List_biasTree_S$V1)

colnames(RF_merge) <- c("RF_20%","RF_40%","RF_60%","RF_80%","S_20%","S_40%","S_60%","S_80%","bias_data_random","bias_data_stratified")

path_merge <- cbind(P_metric,P_metric_S,P_List_biasTree$V1,P_List_biasTree_S$V1)
colnames(path_merge) <- c("RF_20%","RF_40%","RF_60%","RF_80%","S_20%","S_40%","S_60%","S_80%","bias_data_random","bias_data_stratified")

RF_stack <- stack(RF_merge)
path_stack <- stack(path_merge)

colors_10 <- c("#6BAED6" ,"#4292C6", "#2171B5", "#084594","#B8E186", "#7FBC41", "#4D9221", "#276419","red","yellow")

#testing
boxplot(values~ind, data = RF_stack, outline=FALSE,
        xlab = "Sampling Method",
        ylab= "RF distance",
        ylim = rev(range(RF_stack$values)),
        main= "RF distance difference among Random sampling, Stratified sampling and biased sampling",
        col= colors_10,
        names= c("20% Random","40% Random", "60% Random", "80% Random","20% Stratified","40% Stratified", "60% Stratified", "80% Stratified","bias_data_random","bias_data_stratified"))

#testing
boxplot(values~ind, data = path_stack, outline=FALSE,
        xlab = "Sampling Method",
        ylab= "path distance",
        ylim = rev(range(path_stack$values)),
        main= "path distance difference among Random sampling, Stratified sampling and biased sampling",
        col= colors_10,
        names= c("20% Random","40% Random", "60% Random", "80% Random","20% Stratified","40% Stratified", "60% Stratified", "80% Stratified","bias_data_random","bias_data_stratified"))




