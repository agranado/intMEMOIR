# Author:
# Alejandro Granados Castro, Ph.D.
# Postdoctoral Scholar at Elowitz Lab
# California Institute of Technology
# agranado@caltech.edu
# Compiled: Oct-15-2019

# Please install the following libraries before running the script:

library(phangorn)
library(cluster)
library(stringr)
library(phytools)
library(dplyr)
library(dendextend)

# #USAGE Example:
# #EXECUTE main script to load function and calculate mutation parameters
#source("memoirReconstruction.R")
#READ one of the newick files, which contains the barcodes as leafs' ID
#a_tree = read.newick("training_nwk/sub1_train_1.nwk")
#a_reconstructed_tree = reconstructLineage(a_tree, mu, alpha, return_tree = T)
#For more details, visit: https://agranado.github.io/DataScienceWorkshopMX/DREAM.html

#Notes:
# 1. Some functions translate the readout from the 012 notation to uxr, this does not affect any of the methods

# 2. The main file train_setDREAM2019.txt can be read with read.table and it does contain all ground truth trees (Recommended method)
# 2.1 Newick files can also be read directly using read.newick from phytools package


# 3. We assume that all your train files ( sub1_train_1.txt , sub1_train_2.txt, ...) are in the folder ~/training/10mer_2019/

integrase_folder = "training/"
file.path=paste(integrase_folder,"10mer_2019/",sep="")

# 4. for estimation of parameters ALL training files are merged into a single file which is then read to compute mutation rates (this file is also provided as allBarcodes.txt)
# This file has NO cellID, since it is only used for estimation of edit rates across the whole population
# GLOBAL PARAMETERS
# \mu : total edit rate  \alpha: bias towards state 2
# such that for a given generation the mutation probabilities are:
# Pr(edit) = \mu
# Pr(state = 0| edit) = 1-\alpha
# Pr(state = 2| edit) = \alpha
# indicate the path to the file allBarcodes.txt

estim.params.global <-function(estimG = 3,file = "allBarcodes.txt"){
  all.barcodes<-read.table(file)
  rates<-apply(all.barcodes,2,table)
  rates.pct = rates/apply(rates,2,sum)
  #estim parameters per site:
  estimMu = 1-rates.pct[2,]^(1/estimG)
  #rates.pct[3,] means "X"
  estimAlpha = 1-rates.pct[3,]/(1-rates.pct[2,])

  return(list(estimMu,estimAlpha))
}

# Execute the function so the parameters are available for other functions
# in our case, allBarcodes is in ./parameters/ so: 
params_global = estim.params.global(estimG = 4, file = paste("parameters/allBarcodes.txt",sep=""))
mu = params_global[[1]]
alpha = params_global[[2]]

# Though the edit rate is supposed to be homogenoeus across recording units, this is not observed experimentall and therefore we estimate
# edit rates per site. \mu and \alpha are 10-unit arrays
#returns a list of two vectors corresponding to edits rates and alphas
estim.params.global <-function(estimG = 3,file = "allBarcodes.txt"){
  all.barcodes<-read.table(file)
  rates<-apply(all.barcodes,2,table)
  rates.pct = rates/apply(rates,2,sum)
  #estim parameters per site:
  estimMu = 1-rates.pct[2,]^(1/estimG)
  #rates.pct[3,] means "X"
  estimAlpha = 1-rates.pct[3,]/(1-rates.pct[2,])

  return(list(estimMu,estimAlpha))
}


# this function reconstructs a single lineage
# It takes the barcode from the ground truth tree
reconstructLineage<-function(ground_phylo,mu,alpha,return_tree = F,clust.method = "diana",nGen = 4){

  #get the barcode and cell id (cell id is not neccessarily continuous numbers)
  barcodes = str_split(ground_phylo$tip.label,"_",simplify=T)[,2]
  cell_ids = str_split(ground_phylo$tip.label,"_",simplify=T)[,1]

  # translate the barcode data to the old notation uxr
  barcodes_urx = str_replace_all(barcodes, c("2" = "r", "1" = "u","0"="x"))
  barcodes_urx = str_replace_all(barcodes_urx, c("3" = "u"))

  #recontruct the tree with likelihood method
  matdist_=manualDistML_2(as.character(barcodes_urx),mu = mu ,alpha = alpha,nGen = nGen)



  row.names(matdist_)<- paste(cell_ids,barcodes,sep="_")
  colnames(matdist_)<- paste(cell_ids,barcodes,sep="_")


    hclust.tree = clusterDistMatrix(matdist_,clust.method)


  #here ground_phylo is "the ground truth" which is a simulation but this is what we want
  d = 1- RF.dist(hclust.tree,ground_phylo,normalize =T)
  if(return_tree){
    return(hclust.tree)
  }else{
    return(d)
  }

}

# In principle, we can use any clustering method to recover a dendrogram from the
# distance matrix. Divisive clustering (diana) show higher performance than other clustering method
# D2.ward & complete_linkage methods also perform well
clusterDistMatrix<-function(matdist_,clust.method = "diana"){

  if(clust.method =="diana"){
    hclust.tree = as.phylo(as.hclust( diana(as.dist(t(matdist_)))))
  }else {

    hclust.tree=as.phylo(hclust(as.dist(t(matdist_)),method = clust.method))
  }

    return(hclust.tree)
}



#Functions to calculate the probabilities of different mutations:
# This is an unpublished method for reconstructing lineages
# It does calculate the probability that a particular unit in the array is observed to be in a given state after n generations
# It then calcualtes the probablity that two given cells are siters i.e., share a parent ancestor in the previous generation


Pr_edit <- function (nGen,mu,alpha){
  #probability that at nGen generations there is a mutation
  cumSum =0
  if(nGen>0)
    for (nG in 1:nGen)
      cumSum = cumSum + mu * (1-mu)^(nG-1) * alpha
  else
    cumSum =0

  return(cumSum)
}

Pr_noedit<-function(nGen,mu){
  (1-mu)^nGen
}

#this function gives the pr that a site is in each of the thre states for a given number of generations.

Pr_s0<-function(a,mu,alpha,nGen){
  if(a=="u")
    pr = Pr_noedit(nGen,mu)
  else if(a=="r")
    pr =Pr_edit(nGen,mu,alpha)
  else
    pr =Pr_edit(nGen,mu,1-alpha)

  return(pr)
}

all.probs<-function(mu,alpha,nGen,alphabet,  cousin = F){
  #this is of the form i,j   Pr(i->j | i)
  Tran.pr = t(matrix(c(1-mu, mu*alpha,mu*(1-alpha),0,1,0,0,0,1),nrow=3,ncol=3))
  #may29 NOTE: here we need a list of Tran.pr matrices, one per site

  # if we calculate cousin distance, the parent is g-2 instead of g-1
  if(cousin) parent_time = 2 else parent_time = 1;

  #NULL MODEL: probability of observing sites as independent events:
  Pr = array()
  Pr[1] = (1-mu)^nGen

  Pr[2] = Pr_edit(nGen,mu,alpha)

  Pr[3] = Pr_edit(nGen,mu,1-alpha)

  PrMatrix  = array(0,dim=c(length(Pr),length(Pr)))
  #calcualte probabilistic model: P(x,y)
  #this just means the pr that those two sites have those characters by random. This is the expected pr for each site
  #it assummes independence but does not tell you how likely they are to come from a common ancestor
  for (p1 in 1:length(alphabet)){
    for (p2 in 1:length(alphabet)){
      PrMatrix[p1,p2] = Pr[p1] * Pr[p2]
    }
  }

  #Probability for a cel in the previous generation to have a given state
  Pr_s0_array = array()
  for (i in 1:length(alphabet)){
    Pr_s0_array[i]= Pr_s0(alphabet[i], mu, alpha, nGen-parent_time)
  }

  return(list(PrMatrix,Pr_s0_array,Tran.pr))
}

#MAIN method. Computes the likelihood of being sisters for all pairs of cells in the dataset:
manualDistML_2 <- function(barcodeLeaves,mu,alpha,nGen){

  removeLargeDels = F

  alphabet = array()
  alphabet[1]= "u"
  alphabet[2]= "r"
  alphabet[3]= "x"
  #mu & alpha indicate an explicit probabilistic model
  #Probability of no mutation for nGen cell divisions (mu is rate of edit per cell division)

  #transition probabilities

  #weights for number of sustitutions
  equalU =0
  oneSust = 1
  twoSust = 2

  nBarcodes = length(barcodeLeaves)
  barcodeLength=nchar(barcodeLeaves[1])
  distMat= array(0,dim =c(nBarcodes,nBarcodes))
  ratioMat = array(0,dim=c(nBarcodes,nBarcodes))
  productMat = array(0,dim=c(nBarcodes,nBarcodes))


  #go through all the elements in the barcode array
  for (i in 1:(nBarcodes-1)){

    for (j in (i+1):nBarcodes){

      barcodeLeaves1 = barcodeLeaves[i]
      barcodeLeaves2 = barcodeLeaves[j]

      if(removeLargeDels){

            #New function that will look for large deletions (more than 4 consecutive x)
            #The function will replace the xxxx to match whatever state the other cell has such that
            #those sites have minimum effect on the distance calculation
            #IF large_del in barcodeArray1 THEN copy what barcodeArray2 has
            maskedBarcodes = revomeDels(barcodeLeaves1,barcodeLeaves2) #barcodeArray2 is not modified
            barcodeLeaves1 = maskedBarcodes[1]
            barcodeLeaves2 = maskedBarcodes[2]

            #END section: large scale deletion
      }
      #after deletion check we can split now the barcode into a character array
      #and continue with the normal function
      barcodeArray1 =strsplit(barcodeLeaves1,"")[[1]]
      barcodeArray2 =strsplit(barcodeLeaves2,"")[[1]]


      distSum = 1
      ratio.sum =0
      ratio.product=1
      #for each pairwise comparison
      Pr.s1_s2=array()
      Pr.ind=array()
      for (s in 1:barcodeLength){ #may29 NOTE here is where we go site by site
        #new edit:
        #mu and alpha are vectors, so for each element we estimate the transition matrix and P(x,y)
        res.list =all.probs(as.numeric(mu[s]),as.numeric(alpha[s]),nGen,alphabet )
        PrMatrix = res.list[[1]]
        Pr_s0_array = res.list[[2]] #Pr of different letters in previous generation
        Tran.pr = res.list[[3]]
        #

        # PrMatrix = PrMatrix_list[[s]] #NOTE
        Pr.ind[s] = PrMatrix[which(alphabet ==barcodeArray1[s]),which(alphabet ==barcodeArray2[s])]
        #Pr of observing these characters independtly arising Pr1 * Pr2 : assuming independence at nGen
        # Sum{i=1}{3} Pr(s1|S0=i)Pr(s2|s0=i)Ps(S0=i)
        s1=which(alphabet==barcodeArray1[s])
        s2=which(alphabet==barcodeArray2[s])

        #sum over all possible S0 states
        sumPr=0
        #may29 NOTE here get a Tran.pr
        #Tran.pr = Tran.pr_list[[s]]
        for(a in 1:length(alphabet)){
           sumPr = sumPr+ Tran.pr[a,s1] * Tran.pr[a,s2] * Pr_s0_array[a]   # Pr_s0(alphabet[a],mu,alpha,nGen-1)
        }
        Pr.s1_s2[s] =sumPr
      } #end barcode loop
      #end may 29 NOTE



        ratio.product=sum(log(Pr.s1_s2))/sum(log(Pr.ind))
        # NOTE: remove the following line and everything will be OK
        # ratio.product=sum(log(Pr.s1_s2)) #/sum(log(Pr.ind))


      distMat[i,j]= distSum
      ratioMat[i,j]=1/ratio.sum *distSum
      productMat[i,j] = ratio.product

    }
  }
  return(productMat)
}
