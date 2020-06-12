    #function
    #distance between ground truth and reconstructed trees

    #this function takes the reconstructed tree and the ground truth to calculated;
#     #the RF distance. Normalization comes from the distribution of RF distnaces
#     #for random trees (with the same numner of leaves. )
#
# rand.RF.dist<-function(new.tree){
#
#   n.reps = 100
#   rand.dist =c()
#   for(i in 1:n.reps){
#     random.tree = rcoal(length(new.tree$tip.label),tip.label = new.tree$tip.label)
#     rand.dist[i] = RF.dist(new.tree, random.tree)
#   }
#   return(rand.dist)
# }
#
# normalized.RF.experiment<-function(true.tree, manualTree, barcodes,posInfo){
#     library(ape)
#     barcodeLength = 10
#
#     #Calculate the distance of all leaves to the root of the trees
#     leaves.length= diag(vcv.phylo(true.tree))
#     ground.labels = true.tree$tip.label
#
#     #specific for pos40
#
#     #barcodes = barcodes[-length(barcodes)]
#     method = 2
#     # keep only the cells that stayed until the end of the Movie
#     #In theory, the max value of length should correspond to  the last movie frame
#     #Some cells, though, have max-1 (or something like that)
#     if(method==1){
#         length.threshold = max(leaves.length)-1
#         drop.leaves = leaves.length[ leaves.length < length.threshold]
#         keep.leaves = leaves.length[ leaves.length >= length.threshold]
#         #remove dead cells:
#       #method 1 : more systematic
#       #identify the lenght of the branches then remove dead cells
#         new.tree = true.tree
#         for(i in 1:length(drop.leaves)){
#           new.tree = drop.tip(new.tree,names(drop.leaves)[i])
#       }
#     }else if(method==2){
#       #method 2: more biologist-style
#       #we know the lenght of the branches so just match the name (Grace put that!)
#       last.frame = 216 #we know this from the experiment
#       lf.pattern = paste(toString(last.frame),"_",sep="")
#       #keep.index contains the cells that are alive  (i.e. have FISH readout)
#       #however, still contains cells that have xxxx.. in their data. So we need to remove them
#       keep.index = grep(lf.pattern,ground.labels)
#
#       #remove XXXX
#       good.cells = posInfo$cell[!is.na(posInfo$cell) & !posInfo$state==paste(rep("x",barcodeLength),collapse="")]
#       bad.cells = posInfo$cell[!(!is.na(posInfo$cell) & !posInfo$state==paste(rep("x",barcodeLength),collapse=""))]
#
#       good.cells = as.character(as.integer(good.cells))
#
#       keep.leaves = paste(lf.pattern,good.cells,sep="")
#
#       keep.them = c();for(h in 1:length(keep.leaves)){keep.them[h] = grep(paste(keep.leaves[h],"$",sep=""),ground.labels)}
#       drop.leaves = ground.labels[-keep.them]
#
#       # if(length(bad.cells>0)){
#       #   drop.leaves.x = paste(lf.pattern,bad.cells,sep="")
#       #   drop.leaves = c(ground.labels[-grep(lf.pattern,ground.labels)],drop.leaves.x)
#       # }else{
#       #
#       #     drop.leaves = ground.labels[-grep(lf.pattern,ground.labels)]
#       # }
#
#
#       #keep.leaves = ground.labels[grep(lf.pattern,ground.labels)]
#
#       #EDIT: we can take the cells from the excel file directly since there is only FISH for cells that
#       #make it to the end fo the movie, so the filtering for dead cells was already done by GRACE
#
#       new.tree = true.tree
#       for(i in 1:length(drop.leaves)){
#           new.tree = drop.tip(new.tree,drop.leaves[i])#drop.leaves are the names already
#       }
#     }
#
#     #THE NEW labels for the alive cells
#     alive.cells = new.tree$tip.label
#     new.alive.cells = array(0,length(alive.cells))
#
#     id.n=names(barcodes)
#     for(i in 1:length(barcodes)){
#         this.cell= grep(paste(lf.pattern,id.n[i],"$",sep=""),alive.cells)
#         new.alive.cells[this.cell] = paste(names(barcodes)[i],"_",  toString(barcodes[i]),sep="")
#     }
#
#
#     tip.names = new.alive.cells
#     for(tt in 1:length(tip.names)){
#       this.cell.n = str_extract(tip.names[tt],"\\d+")
#       if(nchar(this.cell.n)==1){
#         tip.names[tt]=paste("0",tip.names[tt],sep="")
#       }
#     }
#
#     tip.names=gsub("r","2",tip.names)
#     tip.names=gsub("x","0",tip.names)
#     tip.names=gsub("u","1",tip.names)
#     #  OLD notation
#     # tip.names=gsub("r","2",tip.names)
#     # tip.names=gsub("x","1",tip.names)
#     # tip.names=gsub("u","0",tip.names)
#
#
#     new.tree$tip.label<-tip.names
#
#     #calculate the distribution of RF distances using random clustering based on the number of leaves:
#     n.reps = 100
#     rand.dist =c()
#     for(i in 1:n.reps){
#       random.tree = rcoal(length(new.tree$tip.label),tip.label = new.tree$tip.label)
#       rand.dist[i] = RF.dist(new.tree, random.tree)
#     }
#     if(mean(rand.dist)==0){
#       #The case for unrooted trees of size 3
#       this.score = 1 #kind of cheating, but is true
#     }else{
#       this.score = 1-RF.dist(manualTree,new.tree)/mean(rand.dist)
#       #this.score = 1-RF.dist(manualTree,new.tree, normalize = T)
#     }
#
#     if(this.score<0){this.score=0}
#     return(list(this.score,new.tree,mean(rand.dist)))
#
#   }
# file.path.membow = file.path
# #file.path.memoir = paste("/Users/alejandrog/MEGA/Caltech/trees/",integrase_folder,"10mer/",sep="")
#
# normalized.RF.membow__old<-function(true.tree,manualTree,barcodes,posInfo){
#
#
#   # 1. Find clusters based on identical genotypes
#       unique.genotypes = unique(substr(manualTree$tip.label,4,14))
#       if(length(unique.genotypes)>2){
#         # 2. Reconstruct based on the new clusters of identical cells
#         #      {What is the right order for 1 & 2?}
#           RF.list = normalized.RF.experiment(true.tree,manualTree,barcodes,posInfo)
#           alive.tree =RF.list[[2]]
#           rand.dist = RF.list[[3]]
#           #simplify the reconstructed tree, by removing duplicated genotypes
#           new.tree = manualTree
#           new.alive.tree = alive.tree
#           for(g in 1:length(unique.genotypes)){
#               drop.leaves = manualTree$tip.label[grep(unique.genotypes[g],manualTree$tip.label)]
#               #delete all except one
#
#               if(length(drop.leaves)>1){
#                   for(i in 1:(length(drop.leaves)-1)){
#                       new.tree = drop.tip(new.tree,drop.leaves[i])
#                       new.alive.tree = drop.tip(new.alive.tree,drop.leaves[i])
#                   }
#               }
#
#           }
#
#           #edit 2019 Apr 1st
#           #after filtering for unique genotypes, we can reconstruct the tree again
#           barcodes_ = new.tree$tip.label
#           barcodes.raw  =  substr(barcodes_,4,14)
#           barcodes.raw = str_replace_all(barcodes.raw, c("2" = "r", "1" = "u","0"="x"))
#           matdist_2 = manualDistML(barcodes.raw,estimMu,estimAlpha,estimG)
#           colnames(matdist_2)<-barcodes_; row.names(matdist_2)<-barcodes_
#           new.tree=as.phylo(hclust(as.dist(t(matdist_2))))
#
#           pdf(paste(file.path,"membow/", file.name,".pdf",sep=""))
#             par(mfrow=c(1,2));plot(new.alive.tree,main = "ground truth"); plot(new.tree,main="reconstruction")
#           dev.off()
#
#
#
#   # 3. Using the ground truth, remove cells that were missclassified, such that _groundTruth
#   # and reconstruction are comparable.
#
#   # 4. Somehow map the ground truth to the dendrogram by excluding cells that
#   # were missclassified from both trees (make the clades identical)
#           this.rand =mean(rand.RF.dist(new.tree))
#           if(this.rand<2){ this.rand=2}
#
#           this.score = 1-RF.dist(new.alive.tree,new.tree)/this.rand
#           this.score = 1-RF.dist(new.alive.tree,new.tree,normalize = T)
#     }else{
#           this.score = 1
#
#     }
#   # 5.
#   return(list(this.score,length(unique.genotypes)))
# }
#
#
# normalized.RF.membow<-function(true.tree,manualTree,barcodes,posInfo,make.plots=F,file.name = "",
#                       clust.method = "complete",global = T,globalG =4,file.path = "/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer_2019/" ){
#
#
#       mem.method = 1
#   # 1. Find clusters based on identical genotypes
#       unique.genotypes = unique(substr(manualTree$tip.label,4,14))
#       if(length(unique.genotypes)>2){
#         # 2. Reconstruct based on the new clusters of identical cells
#         #      {What is the right order for 1 & 2?}
#           RF.list = normalized.RF.experiment(true.tree,manualTree,barcodes,posInfo)
#           #alive.tree is the ground.truth after removing dead cells:
#           alive.tree =RF.list[[2]]
#           rand.dist = RF.list[[3]]
#           #simplify the reconstructed tree, by removing duplicated genotypes
#           new.tree = manualTree
#           new.alive.tree = alive.tree
#           for(g in 1:length(unique.genotypes)){
#               drop.leaves = manualTree$tip.label[grep(unique.genotypes[g],manualTree$tip.label)]
#               #delete all except one
#
#               if(length(drop.leaves)>1){
#                   for(i in 1:(length(drop.leaves)-1)){
#                       new.tree = drop.tip(new.tree,drop.leaves[i])
#                       new.alive.tree = drop.tip(new.alive.tree,drop.leaves[i])
#                   }
#               }
#
#           }
#           if(mem.method==1){
#               #edit 2019 Apr 1st
#               #after filtering for unique genotypes, we can reconstruct the tree again
#               barcodes_ = new.tree$tip.label
#               barcodes.raw  =  substr(barcodes_,4,14)
#               barcodes.raw = str_replace_all(barcodes.raw, c("2" = "r", "1" = "u","0"="x"))
#               param.list = estim.params(barcodes.raw) #estimate all parameters from the unique list
#
#               if(global){
#                 param.list = estim.params.global(estimG = globalG,file = paste(file.path,"editRate/allBarcodes.txt",sep=""))
#                 estimMu = param.list[[1]] #these are vectors!
#                 estimAlpha = param.list[[2]]
#
#                 matdist_2 = manualDistML_2(barcodes.raw,estimMu,estimAlpha,3)
#
#
#               }else{
#                 matdist_2 = manualDistML__(barcodes.raw,param.list[[1]],param.list[[2]],param.list[[3]])
#               }
#
#               #until here, barcodes have names urx, we need to convert back (by renaming the matrix)
#
#               #then compare to the ground truth tree
#
#               #let's rename the reconstructed tree without ID number for cells, just the (unique) barcode with 102 notation
#               colnames(matdist_2)<- substr(barcodes_,4,14); row.names(matdist_2)<- substr(barcodes_,4,14)
#
#
#               if(clust.method=="diana"){
#
#                 new.tree = as.phylo(as.hclust(diana(as.dist(t(matdist_2)))))
#               }else{
#
#                 new.tree=as.phylo(hclust(as.dist(t(matdist_2)), method =clust.method))
#               }
#
#
#
#               pdf(paste(file.path,"membow/", file.name,".pdf",sep=""))
#                 par(mfrow=c(1,2));plot(new.alive.tree,main = "ground truth"); plot(new.tree,main="reconstruction")
#               dev.off()
#
#
#
#               # 3. Using the ground truth, remove cells that were missclassified, such that _groundTruth
#               # and reconstruction are comparable.
#
#
#               new.alive.tree_ = new.alive.tree;
#               new.alive.tree_$tip.label = substr(new.alive.tree$tip.label,4,14)
#               # 4. Somehow map the ground truth to the dendrogram by excluding cells that
#               # were missclassified from both trees (make the clades identical)
#               this.rand =mean(rand.RF.dist(new.tree))
#               if(this.rand<2){ this.rand=2}
#           }else{ #do nothin additional after removing repeated barcodes and just calculate the distance without number ID
#
#               new.alive.tree_ = new.alive.tree
#               new.alive.tree_$tip.label = substr(new.alive.tree$tip.label,4,14)
#
#               new.tree$tip.label = substr(new.tree$tip.label,4,14)
#               this.rand =mean(rand.RF.dist(new.tree))
#
#               if(make.plots){
#
#                   plot(new.alive.tree_, main = "ground truth reduced")
#                   plot(new.tree, main =  "reconstruction reduced")
#
#                   pdf(paste(file.path,"membow/", file.name,".pdf",sep=""))
#                     par(mfrow=c(1,2));plot(new.alive.tree_,main = "ground truth"); plot(new.tree,main="reconstruction")
#                   dev.off()
#               }
#
#           }
#
#
#
#           this.score = 1-RF.dist(new.alive.tree_,new.tree)/this.rand
#           this.score = 1-RF.dist(new.alive.tree_,new.tree,normalize = T) #RF distance without number ID may20
#
#           if(is.na(this.score)){
#             #add a root to the tree, based on the initial state of the barcodes: 1111...1
#             #then calcualte the distance using RF(x, root = T)
#             this.score = 1-RF.dist(root(add.tips(new.tree, "1111111111", 4),"1111111111"),
#                 root(add.tips(new.alive.tree_, "1111111111", 4),"1111111111"),normalize = T,rooted = T)
#           }
#
#     }else{
#           this.score = 1
#
#     }
#   # 5.
#   return(list(this.score,length(unique.genotypes) , new.alive.tree_, new.tree))
# }
#
#
# estim.params <-function(barcodes,lengthBarcode = 10,estimG = 0){
#   nCells=length(as.character(barcodes))
#
#   if(estimG==0){
#     estimG = log(nCells)/log(2)
#   }#by default estimG is estimated by the number of cells, but user can specify it
#
#   avgEditRate = array(); avgEditRate_r=array(); avgEditRate_x=array()
#   for(i in 1:nCells){
#      avgEditRate[i] = sum(strsplit(toString(barcodes[i]),"")[[1]]=="u")/lengthBarcode #fraction of "u" in the array
#      avgEditRate_r[i] = sum(strsplit(toString(barcodes[i]),"")[[1]]=="r")/lengthBarcode
#      avgEditRate_x[i] = sum(strsplit(toString(barcodes[i]),"")[[1]]=="x")/lengthBarcode
#    }
#    pos.vals=!(avgEditRate_x==0 & avgEditRate_r==0) #which cells to take into accoutn
#
#    estimMu = 1-mean(avgEditRate[pos.vals])^(1/estimG)
#    estimAlpha = 1-mean(avgEditRate_x[pos.vals] / (avgEditRate_r[pos.vals] + avgEditRate_x[pos.vals]))
#
#    return(list(estimMu,estimAlpha,estimG))
# }
#
# #estimation of edit rate and alpha per site using all data
#returns a list of two vectors corresponding to edits rates and alphas
#May 28
estim.params.global <-function(estimG = 3,file = "../GraceData/10mer_2019/editRate/allBarcodes.txt"){
  all.barcodes<-read.table(file)
  rates<-apply(all.barcodes,2,table)
  rates.pct = rates/apply(rates,2,sum)
  #estim parameters per site:
  estimMu = 1-rates.pct[2,]^(1/estimG)
  #rates.pct[3,] means "X"
  estimAlpha = 1-rates.pct[3,]/(1-rates.pct[2,])

  return(list(estimMu,estimAlpha))
}
