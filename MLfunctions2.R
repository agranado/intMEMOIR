

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
#works fine as May 8th

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
