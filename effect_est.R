# Input: df = data frame of data including A, Y, NC candidates
# Output: ATE

effect.est <- function(df){
  rslt.all = NULL
  A = df[,"A"] # A is a binary treatment
  Y = df[,"Y"] # Y is an outcome
  NCf = FindNegativeControls(df, "A", "Y") #NC found
  NC = subset(df, select = -c(A, Y))
  numNC = ncol(NC)
  numP = numNC*(numNC-1)
  weight = rep(0, numP)
  indw = 1
  ### Estimate the ATE for each pair of NC
  par.all = NULL
  if(ncol(NCf) > 0){
    for(i in 1:numNC){
      for(j in setdiff(1:numNC,i)){
        W = NC[,i]
        Z = NC[,j]
        
        for(l in 1:ncol(NCf)){
          if((names(NC)[i] %in% NCf[,l]) & (names(NC)[j] %in% NCf[,l])) weight[indw] = weight[indw]+1
        }
        indw = indw+1
        ###for one pair
        U.continuous = function(par){
          ### each pair gives three parameters, ATE is the third one
          beta0=par[1]
          beta=par[2]
          ate=par[3]## always make the last parameter the ATE
          mu = beta0+beta*W+ate*A
          U = cbind(1,A,Z)*(Y-mu)
          return(U)
        }
        U.sum = function(par){
          return(norm(apply(U.continuous(par),2,mean),type="2"))
        }
        est = nlm(f=U.sum,p=rep(1,3))
        par.all = c(par.all,est$estimate)
      }
    }
    
    weight = weight / sum(weight)
  }
  
  ind.ate = seq(3,numP*3,by=3) ## ATE is always the third one
  ate.adj = as.vector(t(weight)%*%par.all[ind.ate]) ### mean of adj
  return (ate.adj = ate.adj)
}

effect.est(df)
