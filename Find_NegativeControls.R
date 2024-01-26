#This is a function for identifying validated negative controls for a treatment and outcome from a collection of candidate variables
#Current version uses bnlearn's ci.test function! Be sure to load bnlearn before using this function!
#data is the input data. tname is the name of the treatment variable. oname is the name of the outcome variable.
#alpha is the threshold for significance in conditional independence tests, default is 0.05.
#output is a collection of different validated sets of negative controls for the treatment and outcome
#
#written by Erich Kummerfeld (erichk@umn.edu). last updated June 21 2020
FindNegativeControls <- function(data,tname,oname,alpha = 0.05,vtalpha = NULL,
                                 ConnectednessFilter = FALSE,MeasuredScreenFilter = FALSE,
                                 PsetCIFilter = FALSE, FindBest = FALSE) {
  ###########################################private functions first##################################
  #test of independence that will be used##
  ci <- function(A,Y,Z=c(),alpha,data){
    if (is.null(Z) || is.numeric(Z)){
      temp <- ci.test(A,Y,data = data)
    } else{
      temp <- ci.test(A,Y,z=Z,data = data)
    }
    return(temp$p.value > alpha)
  }
  
  #connectedness function checks if a variable is connected to treatment and outcome###
  connectedness <- function(A,alpha,tname,oname,data){
    return(  !ci(A,tname,alpha=alpha,data=data) && 
               !ci(A,oname,alpha=alpha,data=data) && 
               !ci(A,tname,c(oname),alpha=alpha,data=data) && 
               !ci(A,oname,c(tname),alpha=alpha,data=data) )
  }
  
  #no_measured_screen function checks if a variable is not screened off from the treatment or outcome by any measured variable
  no_measured_screen <- function(A,alpha,tname,oname,data){
    noTOA <- setdiff(names(data),c(tname,oname,A))
    for(Y in noTOA){
      if(ci(A,tname,c(Y),alpha=alpha,data=data) && ci(A,oname,c(Y),alpha=alpha,data=data)){
        #print(paste(A,Y)) #used for bugchecking
        return(FALSE)
      }
    }
    return(TRUE)
  }
  
  #CFilter identifies variables that are connected to the treatment and outcome, and not measurably screened off###
  CFilter <- function(data,tname,oname,alpha){
    noTO <- setdiff(names(data),c(tname,oname))
    connected <- noTO[sapply(noTO,connectedness,USE.NAMES = FALSE,alpha=alpha,data=data,
                             tname=tname,oname=oname)]   #this sapply would be a good place to parallelize
    #if connected is too small, then return NA
    if(length(connected)<3){
      warning("Less than 3 connected variables")
      return(matrix(nrow = 3,ncol = 0))
    }
    return(connected)
  }
  
  #vtet function determines whether or not a tetrad vanishes####
  #  ! ! ! currently this is only implementing the wishart test, so assumes linearity and gaussianity ! ! !
  # x y z w are variables names. Order partly matters, as this tests the tetrad {{x,y},{z,w}}
  # N is the sample size of the data, Co is the correlation matrix, vtalpha is the significance threshold
  vtet <- function(x,y,z,w,N,Co,vtalpha){
    #N <-nrow(data)
    #calculate the determinant of {x,y}{z,w}
    D <- det(Co[c(x,y),c(z,w)])
    
    #calculate the standard deviation of our statistic
    Variance <- ( ((N+1)/(N-1))*det(Co[c(x,y),c(x,y)])*det(Co[c(z,w),c(z,w)])
                  - det(Co[c(x,y,z,w),c(x,y,z,w)]) ) / (N-2)
    
    
    #check if estimated variance is negative. Send a warning if so, but continue anyway after abs()
    if(Variance<0){
      warning(paste("Estimate of Variance was negative for tetrad:",x,y,z,w))
    }
    SD <- sqrt(abs(Variance))
    
    #calcualte the ratio of the determinate to the standard deviation
    ratio <- D / SD
    
    #calculate a 2-sided pvalue for the value of ratio on a normal distribution
    pval <- 2 * pnorm(- abs(ratio))
    
    #return boolean which is TRUE when pvalue is not significant. Uses vtalpha threshold, not alpha.
    res = (pval > vtalpha)
    
    return(list(res = res, pval = pval)) 
  }
  
  #vcheck3 checks if 3 variables satisfy the local negative control validity conditions
  #currently takes both data and Co as inputs. Cheaper to compute Co outside this function, but
  #ci currently runs off the data, so need to use that input as well.
  #input abc is a vector of 3 variable name strings
  vcheck3.all <- function(abc,data,Co,alpha,vtalpha,tname,oname){
    N <-nrow(data)
    tempset <- append(abc,c(tname,oname))
    flag <- TRUE
    #for each pair check ci's on powerset of setdif ####
    if (PsetCIFilter){
      choose2 <- combn(tempset,2)
      for(i in 1:10) {
        A <- choose2[1,i]
        B <- choose2[2,i]
        pset <- powerset(setdiff(tempset,c(A,B)))
        for(j in pset){
          if(ci(A,B,j,alpha,data)){
            flag <- FALSE
            break
          }
        }
        if(! flag){
          break
        }
      }
    }
    #check tempset's vanishing tetrads ####
    if (flag){
      a <- abc[1]
      b <- abc[2]
      c <- abc[3]
      if ( ! (
        vtet(a,b,c,tname,N,Co,vtalpha)$res &&
        vtet(a,c,b,tname,N,Co,vtalpha)$res &&
        vtet(c,b,a,tname,N,Co,vtalpha)$res && #only two of these are independent
        
        vtet(a,b,c,oname,N,Co,vtalpha)$res &&
        vtet(a,c,b,oname,N,Co,vtalpha)$res &&
        vtet(c,b,a,oname,N,Co,vtalpha)$res 
      )){
        flag <- FALSE
      }
    }
    #return flag, so it can be used to filter the set of possible Seeds ####
    return(flag)
  }
  
  vcheck3.best <- function(abc,data,Co,alpha,vtalpha,tname,oname){
    N <-nrow(data)
    tempset <- append(abc,c(tname,oname))
    flag <- TRUE
    pval = 1
    #for each pair check ci's on powerset of setdif ####
    if (PsetCIFilter){
      choose2 <- combn(tempset,2)
      for(i in 1:10) {
        A <- choose2[1,i]
        B <- choose2[2,i]
        pset <- powerset(setdiff(tempset,c(A,B)))
        for(j in pset){
          if(ci(A,B,j,alpha,data)){
            flag <- FALSE
            break
          }
        }
        if(! flag){
          break
        }
      }
    }
    #check tempset's vanishing tetrads ####
    if (flag){
      a <- abc[1]
      b <- abc[2]
      c <- abc[3]
      test1 = vtet(a,b,c,tname,N,Co,vtalpha)
      test2 = vtet(a,c,b,tname,N,Co,vtalpha)
      test3 = vtet(c,b,a,tname,N,Co,vtalpha)
      test4 = vtet(a,b,c,oname,N,Co,vtalpha)
      test5 = vtet(a,c,b,oname,N,Co,vtalpha)
      test6 = vtet(c,b,a,oname,N,Co,vtalpha)
      pval = min(test1$pval, test2$pval, test3$pval, test4$pval, test5$pval, test6$pval)
      if ( ! (
        test1$res && test2$res && test3$res #only two of these are independent
        &&
        test4$res && test5$res && test6$res
      )){
        flag <- FALSE
        
      }
    }
    #return flag, so it can be used to filter the set of possible Seeds ####
    return(list(flag = flag, pval = pval))
  }
  
  #FindSeeds function finds variable triplets that are valid negative controls for the treatment and outcome ####
  FindSeeds <- function(data,ConnectedVariables,alpha,vtalpha,tname,oname,FindBest) {
    Co <- cor(data)
    Seeds <- combn(ConnectedVariables,3)
    numSeeds = ncol(Seeds)
    if(FindBest){
      Seedcheck = unlist(apply(Seeds,2,vcheck3.best,data=data,Co=Co,alpha=alpha,vtalpha=vtalpha,tname=tname,oname=oname))
      SeedPassesCheck <- (Seedcheck[2*(1:numSeeds)-1]==1)  #this apply would be a good place to parallelize
      pvals = Seedcheck[2*(1:numSeeds)]
      maxpval_ind = (pvals == max(pvals)) & (SeedPassesCheck)
      return(Seeds[,maxpval_ind])
    }else{
    SeedPassesCheck <- apply(Seeds,2,vcheck3.all,data=data,Co=Co,alpha=alpha,vtalpha=vtalpha,tname=tname,oname=oname)  #this apply would be a good place to parallelize
    return(Seeds[,SeedPassesCheck])
    }
  }
  
  #powerset function, borrowed from:
  #https://stackoverflow.com/questions/18715580/algorithm-to-calculate-power-set-all-possible-subsets-of-a-set-in-r
  #probably better to write a special function that just does this for my special case of size(s)=3
  #but using this for now for convenience
  powerset = function(s){
    len = length(s)
    l = vector(mode="list",length=2^len) ; l[[1]]=numeric()
    counter = 1L
    for(x in 1L:length(s)){
      for(subset in 1L:counter){
        counter=counter+1L
        l[[counter]] = c(l[[subset]],s[x])
      }
    }
    return(l)
  }
  
  ######Now we can run the public function#############################
  #error handling: make sure tname and oname are in data and not the same. #####
  if(! tname %in% names(data)){
    stop("treatment name not found in data")
  }
  if(! oname %in% names(data)){
    stop("outcome name not found in data")
  }
  if(oname == tname){
    stop("outcome and treatment are the same")
  }
  #If vtalpha is left NULL, then assign it value of inverse sample size, as per experiences with FOFC ####
  if(is.null(vtalpha)) {
    vtalpha <- 1/ nrow(data)
  }
  #run the find  ####
  Variables <- setdiff(names(data),c(tname,oname))
  #check if Connectedness Filter should be used.
  if(ConnectednessFilter){
    Variables <- CFilter(data,tname,oname,alpha)
    #check whether this reduces to less than 3 variables, if so stop and return empty
    if(length(Variables)<3){
      warning("Less than 3 candidate negative controls remain after ConnectednessFilter")
      return(matrix(nrow = 3,ncol = 0))
    }
  }
  #check if measured screening filter should be used.
  if(MeasuredScreenFilter){
    Variables <- Variables[sapply(Variables,no_measured_screen,USE.NAMES = FALSE,alpha=alpha,data=data,tname=tname,oname=oname)]
  }
  #print(ConnectedVariables) #used for bugchecking
  #if ConnectedVariables is too small, then return NA
  if(length(Variables)<3){
    warning("Less than 3 candidate negative controls available for tetrad testing")
    return(matrix(nrow = 3,ncol = 0))
  } else {
    Seeds <- FindSeeds(data,Variables,alpha,vtalpha,tname,oname,FindBest)
    #check if Seeds is a vector. This only happens if there is only one Seed triplet returned. If so, return coerced into a matrix.
    if(is.vector(Seeds)){
      return(matrix(data = Seeds, nrow = 3,ncol = 1))
    } else {
      return(Seeds)
    }
  }
}
