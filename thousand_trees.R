###############################################################################
### Testing the Association Between Sexual Dichromatism/Sexual Selection 
###                   and speciation using 1,000 trees       ##################
###############################################################################

##############################
#####       Setup      #######
##############################
library(foreach)
library(doParallel)
library(nlme)
library(phytools)
library(ape)
library(dplyr)
source("functions/diversification_rate_calculation.R")

trees <- list()

do_one_tree <- function(tree){
  xfjvbzdfsg
  write.rds(zxfv)
}

slurm.apply(trees, do_one_tree)


###Generate 1000 trees for the analysis

#thousand_trees <- passerine.trees
#t5.trees <- read.nexus('data/trees/passerine_2500.nex')
# '%ni%' <- Negate('%in%')
# hundred.tree.names <- names(passerine.trees)
# 
# for(x in names(t5.trees)){
#       ifelse(length(thousand_trees) < 1000, 
#              thousand_trees[[x]] <- t5.trees[[x]],
#              break)
#     }

##Save DFs

# full.data.for.model <- restricted.data %>% select(binomial, TipLabel, SDi, range.size.m2, bioclim4, residuals.PC1, PC1.LIG, NPP)
# 
# subset.data.for.model <- SS.subset %>% select(binomial = binomial.x, TipLabel, Sexual_selection_ppca, SDi, range.size.m2, bioclim4, residuals.PC1, PC1.LIG, NPP)

#saveRDS(thousand_trees, file = "data/trees/thousand_trees.rds")
# saveRDS(full.data.for.model, file = "data/full.data.for.model.rds")
# saveRDS(subset.data.for.model, file = "data/subset.data.for.model.rds")

###############################################################
###Running the loop: 

##Read in the trees and data for model

thousand_trees <- readRDS(file = "data/trees/thousand_trees.rds")

full.data.for.model <- readRDS(file = "data/full.data.for.model.rds")
subset.data.for.model <- readRDS(file = "data/subset.data.for.model.rds")

#Set up cluster and register it with R

########################################
######### SPECIFY CORES HERE!!!! #######
########################################


#setup parallel backend to use many processors
cl <- makeCluster(4) #not to overload your computer
registerDoParallel(cl)

##Run loop, where we: 
    #1: Extract DR,ND
    #2: Add the DR and ND to the dfs
    #Run model: DR Sexual dichromatism, ND Sexual dichromatism, DR/ND for SS
    #Save model object

#Get names of trees
thousand_trees <- thousand_trees[1:4]
tree.names <- names(thousand_trees)

###############################
##### DR + SD MODELS #########
##############################
Start.time <- Sys.time() #Time
i_rep <- 0 #iteration rep set to zero

##Create empty logs
resp_log_DR_SD <- ""

### Start loop for the first set of models
DR.SD.models <- foreach(x=thousand_trees,
                        .final = function(x) setNames(x, names(thousand_trees)),
                        .packages = c("nlme", "dplyr", "foreach", "phytools",  "ape"),
                        .inorder = TRUE) %dopar% {

  #start time and rep counter
  i_rep <- i_rep+1
  t1 <- Sys.time()
  
### CALCULATE DR: These functions call the function loaded in at the start of this script
dr.values <- calculate.log.es(x)

#create frame and turn tiplabel to character (avoids warning during left_join)
rate.frame <- data.frame(TipLabel = names(dr.values[[1]]), 
                         DR = unlist(dr.values, use.names = F)) %>% mutate(TipLabel = as.character(TipLabel))

#Join the rates to the existing dfs
model.data <- left_join(full.data.for.model, rate.frame, by = "TipLabel")

#Prune the tree and set rownames to the TipLabel
pruned.tree <- drop.tip(x,x$tip.label[-match(model.data$TipLabel, x$tip.label)])
rownames(model.data) <- model.data$TipLabel

#Initiate the model through tryCatch
  out <- list() #Set up list
  out <- tryCatch({
    out[["model"]]  <- gls(DR ~ SDi 
                         + log(range.size.m2)
                         + bioclim4 #Seasonal variation
                         + residuals.PC1 #Spatial variation
                         + PC1.LIG #Long-term climate variation
                         + NPP,
                         correlation = corPagel(1, phy = pruned.tree, fixed = FALSE), #Fixed = FALSE allows for lambda to be calculated independently
                         data = model.data, 
                         method = "REML")
    out[["log"]] <- paste0(tree.names[i_rep],": Completed", "\n")
    
    return(out)
  },
  
  error = function(e){
    out[["model"]] <- gls(DR ~ SDi 
                        + log(range.size.m2)
                        + bioclim4 #Seasonal variation
                        + residuals.PC1 #Spatial variation
                        + PC1.LIG #Long-term climate variation
                        + NPP,
                        correlation = corPagel(1, phy = pruned.tree, fixed = TRUE), #Fixed = TRUE keeps lambda at 1
                        data = model.data, 
                        method = "REML")
    out[["log"]] <- paste0(tree.names[i_rep], ": ", e, "\n")
    return(out)})
  
  #Add the most recent log line to the list
  resp_log_DR_SD <- paste0(resp_log_DR_SD, "\n", out['log'])
  
  #second time counter for log
  t2 <- Sys.time()
  
  #Print the  log file and save it to the wd
  sink("log_DR_SD.txt")
  cat("Total run time:", round(t2-Start.time,1), attributes(t2-Start.time)$units,"\n")
  cat("Time for last model:", round(t2-t1,1), attributes(t2-t1)$units,"\n")
  cat("-----------------------------------------------------", "\n")
  cat(resp_log_DR_SD)
  sink()
  return(out[["model"]])
  
}

#Save and RDS of the models into a directory called output
if(!dir.exists("output/")){
  dir.create("output/")
}
saveRDS(DR.SD.models, "output/DR.SD.models.rds")

############################################################################################################

##############################
##### ND + SD MODELS ########
#############################

Start.time <- Sys.time() #Time
i_rep <- 0 #iteration rep set to zero

##Create empty logs
resp_log_ND_SD <- ""

### Start loop for the first set of models
ND.SD.models <- foreach(x=thousand_trees,
                        .final = function(x) setNames(x, names(thousand_trees)),
                        .packages = c("nlme", "dplyr", "foreach", "phytools",  "ape"),
                        .inorder = TRUE) %dopar% {
                          
                          #start time and rep counter
                          i_rep <- i_rep+1
                          t1 <- Sys.time()
                          
                          ### CALCULATE ND: These functions call the function loaded in at the start of this script
                          nd.values <- calculate.nd(x)
                          
                          #create frame and turn tiplabel to character (avoids warning during left_join)
                          rate.frame <- data.frame(TipLabel = names(nd.values[[1]]), 
                                                   ND = unlist(nd.values, use.names = F)) %>% mutate(TipLabel = as.character(TipLabel))
                          
                          #Join the rates to the existing dfs
                          model.data <- left_join(full.data.for.model, rate.frame, by = "TipLabel")
                          
                          #Prune the tree and set rownames to the TipLabel
                          pruned.tree <- drop.tip(x,x$tip.label[-match(model.data$TipLabel, x$tip.label)])
                          rownames(model.data) <- model.data$TipLabel
                          
                          #Initiate the model through tryCatch
                          out <- list() #Set up list
                          out <- tryCatch({
                            out[["model"]]  <- gls(ND ~ SDi 
                                                   + log(range.size.m2)
                                                   + bioclim4 #Seasonal variation
                                                   + residuals.PC1 #Spatial variation
                                                   + PC1.LIG #Long-term climate variation
                                                   + NPP,
                                                   correlation = corPagel(1, phy = pruned.tree, fixed = FALSE), #Fixed = FALSE allows for lambda to be calculated independently
                                                   data = model.data, 
                                                   method = "REML")
                            out[["log"]] <- paste0(tree.names[i_rep],": Completed", "\n")
                            return(out)
                          },
                          
                          error = function(e){
                            out[["model"]] <- gls(ND ~ SDi 
                                                  + log(range.size.m2)
                                                  + bioclim4 #Seasonal variation
                                                  + residuals.PC1 #Spatial variation
                                                  + PC1.LIG #Long-term climate variation
                                                  + NPP,
                                                  correlation = corPagel(1, phy = pruned.tree, fixed = TRUE), #Fixed = TRUE keeps lambda at 1
                                                  data = model.data, 
                                                  method = "REML")
                            out[["log"]] <- paste0(tree.names[i_rep], ": ", e, "\n")
                            return(out)})
                          
                          #Add the most recent log line to the list
                          resp_log_ND_SD <- paste0(resp_log_ND_SD, "\n", out['log'])
                          
                          #second time counter for log
                          t2 <- Sys.time()
                          
                          #Print the  log file and save it to the wd
                          sink("log_ND_SD.txt")
                          cat("Total run time:", round(t2-Start.time,1), attributes(t2-Start.time)$units,"\n")
                          cat("Time for last model:", round(t2-t1,1), attributes(t2-t1)$units,"\n")
                          cat("-----------------------------------------------------", "\n")
                          cat(resp_log_ND_SD)
                          sink()
                          return(out[["model"]])
                          
                        }

#Save and RDS of the models into a directory called output
if(!dir.exists("output/")){
  dir.create("output/")
}
saveRDS(ND.SD.models, "output/ND.SD.models.rds")

############################################################################################################

###############################
##### DR + SS MODELS #########
##############################
Start.time <- Sys.time() #Time
i_rep <- 0 #iteration rep set to zero

##Create empty logs
resp_log_DR_SS <- ""

### Start loop for the first set of models
DR.SS.models <- foreach(x=thousand_trees,
                        .final = function(x) setNames(x, names(thousand_trees)),
                        .packages = c("nlme", "dplyr", "foreach", "phytools",  "ape"),
                        .inorder = TRUE) %dopar% {
                          
                          #start time and rep counter
                          i_rep <- i_rep+1
                          t1 <- Sys.time()
                          
                          ### CALCULATE DR: These functions call the function loaded in at the start of this script
                          dr.values <- calculate.log.es(x)
                          
                          #create frame and turn tiplabel to character (avoids warning during left_join)
                          rate.frame <- data.frame(TipLabel = names(dr.values[[1]]), 
                                                   DR = unlist(dr.values, use.names = F)) %>% mutate(TipLabel = as.character(TipLabel))
                          
                          #Join the rates to the existing dfs
                          model.data <- left_join(subset.data.for.model, rate.frame, by = "TipLabel")
                          
                          #Prune the tree and set rownames to the TipLabel
                          pruned.tree <- drop.tip(x,x$tip.label[-match(model.data$TipLabel, x$tip.label)])
                          rownames(model.data) <- model.data$TipLabel
                          
                          #Initiate the model through tryCatch
                          out <- list() #Set up list
                          out <- tryCatch({
                            out[["model"]]  <- gls(DR ~ Sexual_selection_ppca 
                                                   + log(range.size.m2)
                                                   + bioclim4 #Seasonal variation
                                                   + residuals.PC1 #Spatial variation
                                                   + PC1.LIG #Long-term climate variation
                                                   + NPP,
                                                   correlation = corPagel(1, phy = pruned.tree, fixed = FALSE), #Fixed = FALSE allows for lambda to be calculated independently
                                                   data = model.data, 
                                                   method = "REML")
                            out[["log"]] <- paste0(tree.names[i_rep],": Completed", "\n")
                            return(out)
                          },
                          
                          error = function(e){
                            out[["model"]] <- gls(DR ~ Sexual_selection_ppca 
                                                  + log(range.size.m2)
                                                  + bioclim4 #Seasonal variation
                                                  + residuals.PC1 #Spatial variation
                                                  + PC1.LIG #Long-term climate variation
                                                  + NPP,
                                                  correlation = corPagel(1, phy = pruned.tree, fixed = TRUE), #Fixed = TRUE keeps lambda at 1
                                                  data = model.data, 
                                                  method = "REML")
                            out[["log"]] <- paste0(tree.names[i_rep], ": ", e, "\n")
                            return(out)})
                          
                          #Add the most recent log line to the list
                          resp_log_DR_SS <- paste0(resp_log_DR_SS, "\n", out['log'])
                          
                          #second time counter for log
                          t2 <- Sys.time()
                          
                          #Print the  log file and save it to the wd
                          sink("log_DR_SS.txt")
                          cat("Total run time:", round(t2-Start.time,1), attributes(t2-Start.time)$units,"\n")
                          cat("Time for last model:", round(t2-t1,1), attributes(t2-t1)$units,"\n")
                          cat("-----------------------------------------------------", "\n")
                          cat(resp_log_DR_SS)
                          sink()
                          return(out[["model"]])
                          
                        }

#Save and RDS of the models into a directory called output
if(!dir.exists("output/")){
  dir.create("output/")
}
saveRDS(DR.SS.models, "output/DR.SS.models.rds")

############################################################################################################

###############################
##### ND + SS MODELS #########
##############################
Start.time <- Sys.time() #Time
i_rep <- 0 #iteration rep set to zero

##Create empty logs
resp_log_ND_SS <- ""

### Start loop for the first set of models
ND.SS.models <- foreach(x=thousand_trees,
                        .final = function(x) setNames(x, names(thousand_trees)),
                        .packages = c("nlme", "dplyr", "foreach", "phytools",  "ape"),
                        .inorder = TRUE) %dopar% {
                          
                          #start time and rep counter
                          i_rep <- i_rep+1
                          t1 <- Sys.time()
                          
                          ### CALCULATE DR: These functions call the function loaded in at the start of this script
                          dr.values <- calculate.nd(x)
                          
                          #create frame and turn tiplabel to character (avoids warning during left_join)
                          rate.frame <- data.frame(TipLabel = names(dr.values[[1]]), 
                                                   ND = unlist(dr.values, use.names = F)) %>% mutate(TipLabel = as.character(TipLabel))
                          
                          #Join the rates to the existing dfs
                          model.data <- left_join(subset.data.for.model, rate.frame, by = "TipLabel")
                          
                          #Prune the tree and set rownames to the TipLabel
                          pruned.tree <- drop.tip(x,x$tip.label[-match(model.data$TipLabel, x$tip.label)])
                          rownames(model.data) <- model.data$TipLabel
                          
                          #Initiate the model through tryCatch
                          out <- list() #Set up list
                          out <- tryCatch({
                            out[["model"]]  <- gls(ND ~ Sexual_selection_ppca 
                                                   + log(range.size.m2)
                                                   + bioclim4 #Seasonal variation
                                                   + residuals.PC1 #Spatial variation
                                                   + PC1.LIG #Long-term climate variation
                                                   + NPP,
                                                   correlation = corPagel(1, phy = pruned.tree, fixed = FALSE), #Fixed = FALSE allows for lambda to be calculated independently
                                                   data = model.data, 
                                                   method = "REML")
                            out[["log"]] <- paste0(tree.names[i_rep],": Completed", "\n")
                            return(out)
                          },
                          
                          error = function(e){
                            out[["model"]] <- gls(ND ~ Sexual_selection_ppca 
                                                  + log(range.size.m2)
                                                  + bioclim4 #Seasonal variation
                                                  + residuals.PC1 #Spatial variation
                                                  + PC1.LIG #Long-term climate variation
                                                  + NPP,
                                                  correlation = corPagel(1, phy = pruned.tree, fixed = TRUE), #Fixed = TRUE keeps lambda at 1
                                                  data = model.data, 
                                                  method = "REML")
                            out[["log"]] <- paste0(tree.names[i_rep], ": ", e, "\n")
                            return(out)})
                          
                          #Add the most recent log line to the list
                          resp_log_ND_SS <- paste0(resp_log_ND_SS, "\n", out['log'])
                          
                          #second time counter for log
                          t2 <- Sys.time()
                          
                          #Print the  log file and save it to the wd
                          sink("log_ND_SS.txt")
                          cat("Total run time:", round(t2-Start.time,1), attributes(t2-Start.time)$units,"\n")
                          cat("Time for last model:", round(t2-t1,1), attributes(t2-t1)$units,"\n")
                          cat("-----------------------------------------------------", "\n")
                          cat(resp_log_ND_SS)
                          sink()
                          return(out[["model"]])
                          
                        }

#Save and RDS of the models into a directory called output
if(!dir.exists("output/")){
  dir.create("output/")
}
saveRDS(ND.SS.models, "output/ND.SS.models.rds")
