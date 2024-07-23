# Simulation codes for ROC of negative controls and ATE estimation
# test for fp, fn and ability to find NCs
sim.roc <- function(cont, size, graph, alpha = 0.05, type, c, d, ntry){
  multiple = 2^((1:12)-6)
  na = rep(0, length(multiple)) # probabiliy of no triplet ever found
  answer = answer.key(graph)
  count = rep(list(matrix(0, 3, dim(answer)[2])), length(multiple)) # count for each sample size
  try = rep(list(), length(multiple))
  tp = fp = fn = tn = rep(0, length(multiple))
  tot = base::ifelse(graph == 1, 4, 35)
  noError = 0
  for(j in 1:ntry){
    set.seed(j)
    df <- sim.data(cont=cont, graph = graph, n = size, c = c, d = d)
    varError <<- FALSE
    tryCatch(try <- discover.nc(df, "A", "Y", alpha, multiple), warning = function(w) { varError <<- TRUE }, error = function(e) { varError <<- TRUE })
    
    if(!varError){
      noError = noError + 1
      for(i in 1:length(multiple)){
        if(dim(try[[i]])[2] == 0){
          fn[i] = fn[i] + dim(answer)[2]
          tn[i] = tn[i] + tot - dim(answer)[2]
          na[i] = na[i] + 1
        }else{
          test = test.check(try[[i]], answer)
          tp[i] = tp[i] + test$tp
          fp[i] = fp[i] + test$fp
          fn[i] = fn[i] + dim(answer)[2] - test$tp
          tn[i] = tn[i] + tot - dim(answer)[2] - test$fp
          count[[i]] = count[[i]] + test$count
        }
      }
    }
  }
  
  tpr = tp/(tp+fn)
  fpr = fp/(fp+tn)
  na = round(na / noError, 2)
  
  summary = round(rbind(tpr, fpr, na), 2)
  colnames(summary) <- 2^((1:12)-6)
  
  return (list(tpr = tpr, fpr = fpr, na = na, summary = summary))
}

# ATE estimation for random (all pairs) and DANCE pairs
sim.effect <- function(size, graph, cont = TRUE, alpha = 0.05, c, d, ntry){
  answer = answer.key(graph)
  true <- effect.true(cont, c) # true effect can be calculated prior to for loop
  rslt.all = NULL;eff.all = NULL;
  
  for(k in 1:ntry){
    set.seed(k)
    df <- sim.data(cont=cont, graph = graph, n = size, c = c, d = d)
    A = df[,1]
    Y = df[,2]
    NC = df[,3:dim(df)[2]] #random (all)
    numNC = ncol(NC)
    numP = numNC*(numNC-1)
    ran_ind = sample(3:dim(df)[2], size = 2)-2
    NCf = matrix(names(NC[,ran_ind]))
    weight = rep(0, numP)
    indw = 1
    ### first estimate the ATE for each pair
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
          numDeriv::jacobian(func=U.sum,x=est$estimate)
        }
      }
      
      weight = weight / sum(weight)
      
      # estimating function for all pairs
      U.cont = function(par){
        ate.all=U.all=NULL
        for(i in 1:numNC){
          for(j in setdiff(1:numNC,i)){
            W = NC[,i]
            Z = NC[,j]
            count = 3*length(ate.all)
            beta0.ij=par[count+1]
            beta.ij=par[count+2]
            ate.ij=par[count+3]
            mu.ij = beta0.ij+beta.ij*W+ate.ij*A
            U.ij = cbind(1,A,Z)*(Y-mu.ij)
            U.all = cbind(U.all,U.ij)
            ate.all = c(ate.all,ate.ij)
          }
        }
        return(U.all)
      }
      G.continuous = function(par){return(apply(U.cont(par),2,sum))}
      bread=numDeriv::jacobian(func=G.continuous,x=par.all)
      meat.half=U.cont(par=par.all)
      IF = meat.half%*%t(solve(-bread)) #ginv
      vcov = t(IF)%*%IF ## variance covariance matrix for all parameters
      ind.ate = seq(3,numP*3,by=3) ## ATE is always the third one
      vcov.ate = vcov[ind.ate,ind.ate]
      
      ate.ran = t(rep(1/numP,numP))%*%par.all[ind.ate] ### mean of all ate
      se.ate.ran = sqrt(rep(1/numP,numP)%*%vcov.ate%*%rep(1/numP,numP))  ### se of mean of all ate
      
      ate.adj = t(weight)%*%par.all[ind.ate] ### mean of adj
      se.ate.adj = sqrt(t(weight)%*%vcov.ate%*%weight)  ### se of mean of adj
      
      rslt.all = rbind(rslt.all,c(ate.ran=ate.ran,se.ran=se.ate.ran,ate.adj=ate.adj,se.adj=se.ate.adj))
    }
  }
  
  ate.ran = rslt.all[,1]; se.ran=rslt.all[,2]; ate.adj = rslt.all[,3]; se.adj=rslt.all[,4]
  bias.ran = mean(ate.ran - true); bias.adj = mean(ate.adj - true) # bias
  sd.ran = sd(ate.ran); sd.adj = sd(ate.adj) # MC se
  prop.ran = bias.ran / true; prop.adj = bias.adj / true # proportional bias
  mse.ran = mean((ate.ran - true)^2); mse.adj = mean((ate.adj - true)^2) # MSE
  se.ran = mean(se.ran); se.adj = mean(se.adj) # estimated se
  
  cover.ran = mean(ate.ran-qnorm(0.975)*se.ran<true & true<ate.ran+qnorm(0.975)*se.ran) ## coverage probability of ran
  cover.adj = mean(ate.adj-qnorm(0.975)*se.adj<true & true<ate.adj+qnorm(0.975)*se.adj) ## coverage probability of adj
  
  ran.out = c(bias.ran, prop.ran, sd.ran, se.ran, mse.ran, cover.ran)
  adj.out = c(bias.adj, prop.adj, sd.adj, se.adj, mse.adj, cover.adj)
  output = rbind(ran.out, adj.out)
  colnames(output) = c("Bias", "Prop", "MC SE", "Estimated SE", "MSE", "Coverage")
  
  return (output)
}

# Raw effect estimation for graphs
sim.raw.effect <- function(cont, size, graph, alpha = 0.05, c, d, ntry){
  try = rep(list(), 12)
  answer = answer.key(graph)
  
  true <- effect.true(cont, c) # true effect can be calculated prior to for loop
  
  cover = ran_exp_sum = ran_var_sum = ran_cover_sum = adj_exp_sum = adj_var_sum = adj_cover_sum = noError = 0
  raw_exp_v = ran_exp_v = ran_var_v = ran_cover_v = NULL
  adj_exp_v = adj_var_v = adj_cover_v = NULL
  
  for(k in 1:ntry){
    set.seed(k)
    df <- sim.data(cont=cont, graph = graph, n = size, c = c, d = d)
    cand <- seq(from = 3, to = length(df))
    raw_exp=ran_exp=ran_var=cover=ran_cover=count1=0
    
    if(cont){
      raw_exp = lm(df$Y ~ df$A)$coefficients[2]
      raw_exp_v = c(raw_exp_v, raw_exp)
    }
    else{
      coeff.A = glm(df$Y ~ df$A, family = "binomial")$coefficients
      raw_exp = plogis(coeff.A[1]+coeff.A[2]*1) - plogis(coeff.A[1])
      #raw_exp = plogis(coeff.A[2]) - plogis(0)
      raw_exp_v = c(raw_exp_v, raw_exp)
    }
  }
  
  sd.raw = sd(raw_exp_v)
  
  raw_lb = quantile((raw_exp_v-true)/true, probs = 0.025, na.rm = TRUE)
  raw_ub = quantile((raw_exp_v-true)/true, probs = 0.975, na.rm = TRUE)
  bias_raw = (mean(raw_exp_v) - true)/true
  return (c(raw_mean = bias_raw, raw_lb = raw_lb, raw_ub = raw_ub))
}


answer.key <- function(graph){
  if(graph == 1){
    answer = matrix(c("Z1","Z3","Z4","Z2","Z3","Z4"), 3, 2)
  }
  else{
    answer = matrix(c("Z1", "Z3", "Z6", "Z1", "Z3", "Z7", "Z1", "Z4", "Z6", "Z1", "Z4", "Z7", "Z1", "Z5", "Z6", "Z1", "Z5", "Z7", "Z2", "Z3", "Z6", "Z2", "Z3", "Z7", "Z2", "Z4", "Z6", "Z2", "Z4", "Z7", "Z2", "Z5", "Z6", "Z2", "Z5", "Z7"), 3, 12)
  }
  return (answer)
}

test.check <- function(result, answer){ #compare our result to the answer
  dim1 = dim(result)[2]
  dim2 = dim(answer)[2]
  
  tp = 0
  fp = 0
  count = rep(0, dim2)
  
  if(dim1 == 0){
    return (list(tp = tp, fp = fp, count = count))
  }
  for(i in 1:dim1){
    for(j in 1:dim2){
      if (length(setdiff(answer[,j], result[,i]))==0){
        tp = tp + 1
        count[j] = count[j] + 1
      }
    }
  }
  fp = dim1 - tp
  
  return (list(tp = tp, fp = fp, count = count))
}

# discover NCs
discover.nc <- function(df, tname, oname, alpha, multiple){
  res <- rep(list(), length(multiple))
  n = nrow(df)
  for(i in (1:length(multiple))){
    vtalpha = multiple[i] / n
    res[[i]] <- FindNegativeControls(data = df, tname = tname, oname = oname, vtalpha = vtalpha, alpha = alpha)
  }
  return (res)
}

# simulate data
sim.data <- function(cont, graph, n=n, c, d){ #cont = 0(binomial),1(normal), graph = 1~2, d = edge btw cand, c = other edge
  if(graph == 1 & cont == FALSE){
    U = rbinom(n, 1, 0.5)
    A = rbinom(n, 1, plogis(-1+c[2]*U))
    Y = rbinom(n, 1, plogis(-1+c[3]*U+c[1]*A+c[4]*U*A))
    Z1 = rbinom(n, 1, plogis(-1+c[5]*U))
    Z2 = rbinom(n, 1, plogis(-1+c[6]*U+d[1]*Z1+d[2]*U*Z1))
    Z3 = rbinom(n, 1, plogis(-1+c[7]*U))
    Z4 = rbinom(n, 1, plogis(-1+c[8]*U))
    df = data.frame(A,Y,Z1,Z2,Z3,Z4)
  }
  else if(graph == 1 & cont == TRUE){ 
    U = rnorm(n,mean=0,sd=2)
    A = c[2]*U + rnorm(n,0,1)
    Y = c[1]*A + c[3]*U + rnorm(n,0,1)
    Z1 = c[4]*U + rnorm(n,0,1)
    Z2 = c[5]*U + d[1]*Z1 + rnorm(n,0,1)
    Z3 = c[6]*U + rnorm(n,0,1)
    Z4 = c[7]*U + rnorm(n,0,1)
    df = data.frame(A,Y,Z1,Z2,Z3,Z4)
  }
  else if(graph == 2 & cont == FALSE){
    U = rbinom(n, 1, 0.5)
    A = rbinom(n, 1, plogis(-1+c[2]*U))
    Y = rbinom(n, 1, plogis(-1+c[3]*U+c[1]*A+c[4]*U*A))
    Z1 = rbinom(n, 1, plogis(-1+c[5]*U))
    Z2 = rbinom(n, 1, plogis(-1+c[6]*U+d[1]*Z1+d[2]*U*Z1))
    Z3 = rbinom(n, 1, plogis(-1+c[7]*U))
    Z4 = rbinom(n, 1, plogis(-1+c[8]*U+d[3]*Z3+d[4]*U*Z3))
    Z5 = rbinom(n, 1, plogis(-1+c[9]*U+d[5]*Z3+d[6]*Z4+d[7]*U*Z3+d[8]*U*Z4+d[9]*Z3*Z4+d[10]*U*Z3*Z4))
    Z6 = rbinom(n, 1, plogis(-1+c[10]*U))
    Z7 = rbinom(n, 1, plogis(-1+c[11]*U+d[11]*Z6+d[12]*U*Z6))
    
    df = data.frame(A,Y,Z1,Z2,Z3,Z4,Z5,Z6,Z7)
  }
  else{
    U = rnorm(n,mean=0,sd=2)
    A = c[2]*U + rnorm(n,0,1)
    Y = c[1]*A + c[3]*U + rnorm(n,0,1)
    Z1 = c[4]*U + rnorm(n,0,1)
    Z2 = c[5]*U + d[1]*Z1 + rnorm(n,0,1)
    Z3 = c[6]*U + rnorm(n,0,1)
    Z4 = c[7]*U + d[2]*Z3 + rnorm(n,0,1)
    Z5 = c[8]*U + d[3]*Z3 + d[4]*Z4 + rnorm(n,0,1)
    Z6 = c[9]*U + rnorm(n,0,1)
    Z7 = c[10]*U + d[5]*Z6 + rnorm(n,0,1)
    df = data.frame(A,Y,Z1,Z2,Z3,Z4,Z5,Z6,Z7)
  }
}

expit = function(x){exp(x)/(1+exp(x))} 

# Edge strength
c.wk = c(0.6262472, 0.3817754, 0.6795356, 0.3837477, 0.3565017, 0.4012941, 0.3450560, 0.5648041, 0.3724304, 0.5551002, 0.4074619, 0.3746760, 0.4442358, 0.4431560, 0.3439058, 0.5940274, 0.6612201, 0.3517061, 0.3082211, 0.5345576)
d.wk = c(1.413934, 1.245254, 1.270487, 1.479389, 1.652445, 1.836359, 1.466640, 1.542527, 1.048164, 1.589371, 1.960653, 1.841541, 1.371016, 1.488179, 1.361607, 1.816975, 1.833493, 1.337110, 1.171592, 1.558444)
c.str = c(0.8800869, 0.8755867, 0.8124299, 0.6633984, 0.6223731, 0.7401199, 0.9992366, 0.7725701, 0.6719096, 0.8698893, 0.7598463, 0.6975442, 0.7107355, 0.7997793, 0.9001751, 0.6447044, 0.6683586, 0.8351065, 0.6506573, 0.6790947)
d.str = c(3.749785, 3.370727, 2.306108, 2.881296, 2.347439, 3.622491, 2.573287, 2.680017, 2.166247, 3.917405, 2.155028, 3.238999, 3.962036, 3.788603, 3.248937, 2.461236, 2.460915, 2.980951, 2.005269, 2.548388)
c.bin = c(1.380936, 1.488693, 1.650510, 1.315629, 1.465796, 1.928048, 1.597286, 1.364279, 1.049533, 1.485732, 1.064882, 1.005135, 1.401728, 1.382683, 1.389307, 1.499075, 1.243281, 1.732134, 1.546704, 1.167175, 1.975921, 1.160478, 1.270171, 1.082588, 1.908118, 1.869197, 1.626949, 1.442009, 1.158431, 1.539326)
d.bin = c(1.226158, 1.969087, 1.901274, 1.847480, 1.946074, 1.583777, 1.314048, 1.967713, 1.434478, 1.527845, 1.155267, 1.267736, 1.044821, 1.084514, 1.404189, 1.067947, 1.868043, 1.500284, 1.093246, 1.784776)

#--------------------------------------------------------------------------------------------------
# ROC curve simulation example (weak, simple)
s = c(10, 30, 100, 300, 1000)
auc <- rep(0, 5)
roc = list(list(), list(), list(), list(), list())
for(i in 1:5){
  roc[[i]] = sim.roc(cont=TRUE, size=s[i], graph = 1, ntry = 200, c = c.wk, d = d.wk)
  auc[i] = round(MESS::auc(roc[[i]]$fpr, roc[[i]]$tpr, rule = 2, from = 0, to = 1), 3)
}


col = brewer.pal(6, "Set1")

#tiff("Plot1.tiff", width = 4, height = 3, units = 'in', res = 200)
par(mar = c(3, 3, 1, 1))
plot(0, 0, xlim=c(0, 1), ylim=c(0, 1), xlab = "", ylab = "", type = "n", font.lab = 1, cex.lab = 0.7, cex.axis = 0.7)
title(ylab="True positive rate", cex.lab=0.7, mgp = c(2, 2, 0))
title(xlab="False positive rate", cex.lab=0.7, mgp = c(2, 2, 0))
lines(c(0,1), c(0,1), lwd = 2)

for(i in 1:5){
  lines(c(1, roc[[i]]$fpr, 0), c(1, roc[[i]]$tpr, 0), lty = i, col = col[i], lwd = 2)
  points(x = roc[[i]]$fpr[6], y = roc[[i]]$tpr[6], col = col[i], lwd = 1, pch = 19)
  points(x = roc[[i]]$fpr[7], y = roc[[i]]$tpr[7], col = col[i], lwd = 1, pch = 1)
}
legend("bottomright", inset = 0.05, legend = c(paste("AUC: ", auc[1]), paste("AUC: ", auc[2]), paste("AUC: ", auc[3]), paste("AUC: ", "1.000"), paste("AUC: ", "1.000")), lty = 1:5, lwd = 3, col = col, bty = 'n', cex = 0.6, pt.cex = 0.5)

#dev.off()
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
# ATE simulation example (weak, simple)
# Continuous effect estimation (weak)
eff.raw.wk.100 <- sim.raw.effect(cont=TRUE, 100, 1, 0.05, c=c.wk, d=d.wk, ntry=200)
eff.raw.wk.300 <- sim.raw.effect(cont=TRUE, 300, 1, 0.05, c=c.wk, d=d.wk, ntry=200)
eff.raw.wk.1000 <- sim.raw.effect(cont=TRUE, 1000, 1, 0.05, c=c.wk, d=d.wk, ntry=200)
eff.raw.wk.3000 <- sim.raw.effect(cont=TRUE, 3000, 1, 0.05, c=c.wk, d=d.wk, ntry=200)
eff.raw.wk.10000 <- sim.raw.effect(cont=TRUE, 10000, 1, 0.05, c=c.wk, d=d.wk, ntry=200)
eff.raw.wk.30000 <- sim.raw.effect(cont=TRUE, 30000, 1, 0.05, c=c.wk, d=d.wk, ntry=200)
eff.raw.wk = rbind(eff.raw.wk.300, eff.raw.wk.1000, eff.raw.wk.3000, eff.raw.wk.10000, eff.raw.wk.30000)

#Simple graph (raw bias, raw MC error for random pair, random all, dance total and best)/true
true.eff = effect.true(TRUE, c.wk)

#output
eff.100 = sim.effect(100, 1, TRUE, 0.05, c=c.wk, d=d.wk, ntry = 200)
eff.300 = sim.effect(300, 1, TRUE, 0.05, c=c.wk, d=d.wk, ntry = 200)
eff.1000 = sim.effect(100, 1, TRUE, 0.05, c=c.wk, d=d.wk, ntry = 200)
eff.3000 = sim.effect(100, 1, TRUE, 0.05, c=c.wk, d=d.wk, ntry = 200)
eff.10000 = sim.effect(100, 1, TRUE, 0.05, c=c.wk, d=d.wk, ntry = 200)
eff.30000 = sim.effect(100, 1, TRUE, 0.05, c=c.wk, d=d.wk, ntry = 200)
#eff.100 = c(0.12288034, 1.3158623, 0.07725033, 0.7507134, 0.01703966, 0.7438612, -0.05774123, 0.3565645)
#eff.300 = c(0.06681372, 0.21960917, 0.04575083, 0.09387759, -0.02507188, 0.10692813, -0.02414502, 0.11105707)
#eff.1000 = c(0.06750306, 0.1710122, 0.062625805, 0.04895231, -0.006678859, 0.05738008, -0.006960372, 0.05876792)
#eff.3000 = c(0.08190447, 0.1604826, 6.795682e-02, 0.03050280, 7.733872e-05, 0.03428526, 0.0006818565,  0.03482281)
#eff.10000 = c(0.07329124, 0.1628806, 0.065832926, 0.01611102, -0.001472503, 0.01862019, -0.001535787, 0.01929533)
#eff.30000 = c(0.07888576, 0.1613969, 6.747192e-02, 0.008537057, 1.252259e-05, 0.009835970, 4.716967e-05, 0.010196606)
eff <- rbind(eff.300, eff.1000, eff.3000, eff.10000, eff.30000)/true.eff
eff.ub1 <- eff[,1] + 1.96*eff[,2]
eff.lb1 <- eff[,1] - 1.96*eff[,2]
eff.ub2 <- eff[,3] + 1.96*eff[,4]
eff.lb2 <- eff[,3] - 1.96*eff[,4]
eff.ub3 <- eff[,5] + 1.96*eff[,6]
eff.lb3 <- eff[,5] - 1.96*eff[,6]
eff.ub4 <- eff[,7] + 1.96*eff[,8]
eff.lb4 <- eff[,7] - 1.96*eff[,8]
effect <- cbind(eff.raw.wk, eff[,1], eff.lb1, eff.ub1, eff[,3], eff.lb2, eff.ub2, eff[,5], eff.lb3, eff.ub3, eff[,7], eff.lb4, eff.ub4)
effect

x0 = c(1, 3.5, 6, 8.5, 11)
x1 = c(1.5, 4, 6.5, 9, 11.5)
x2 = c(2, 4.5, 7, 9.5, 12)
x3 = c(2.5, 5, 7.5, 10, 12.5)
x4 = c(3, 5.5, 8, 10.5, 13)

pdf("Plot1.pdf",h=3,w=4)
par(mar = c(3, 3, 1, 1), cex.lab = 0.7, cex.axis = 0.7)
plot(x0, y = effect[,1], xlim = c(0, 14), ylim = c(-0.6, 1.5), type = "p", pch = 19, col = col[1], xaxt = "n", xlab = "", ylab = "", font.lab = 1, cex.lab = 0.7, cex.axis = 0.7, cex = 0.7)
#plot(x0, y = effect[,1], xlim = c(0, 14), ylim = c(-0.6, 1.0), type = "p", pch = 19, col = col[1], xaxt = "n", xlab = "", ylab = "", font.lab = 1, cex.lab = 0.7, cex.axis = 0.7, cex = 0.7)
#plot(x0, y = effect[,1], xlim = c(0, 14), ylim = c(-0.4, 0.6), type = "p", pch = 19, col = col[1], xaxt = "n", xlab = "", ylab = "", font.lab = 1, cex.lab = 0.7, cex.axis = 0.7, cex = 0.7)

#gap.plot(x0, y = effect[,1]+0.5, xlim = c(0, 14), ylim = c(-1.3, 1.3), type = "p", pch = 19, col = col[1], font.lab = 1, cex.lab = 0.7, cex.axis = 0.7, cex = 0.7, gap = c(-1.2, -0.7, 0.65, 1.1), xtics = c(2, 4.5, 7, 9.5, 12), ytics = c(-1.3, -0.5, 0, 0.5, 1.3), xaxt = "n", xlab = "", ylab = "", xticlab = FALSE)

axis(1, c(2, 4.5, 7, 9.5, 12), c("300", "1000", "3000", "10000", "30000"), cex.axis = 0.7) #cont
#axis(1, c(2, 4.5, 7, 9.5, 12), c("1000", "3000", "10000", "30000", "100000"), cex.axis = 0.7) #bin
title(ylab="Proportion bias", cex.lab=0.7, mgp = c(2, 2, 0))
title(xlab="Size", cex.lab=0.7, mgp = c(2, 2, 0))

points(x1, y = effect[,4], col = col[2], pch = 19, cex = 0.7)
points(x2, y = effect[,7], col = col[3], pch = 19, cex = 0.7)
points(x3, y = effect[,10], col = col[4], pch = 19, cex = 0.7)
points(x4, y = effect[,13], col = col[5], pch = 19, cex = 0.7)

abline(h=0)
#abline(h=-0.5)
arrows(x0, effect[,2], x0, effect[,3], length=0.05, angle=90, code=3, col = col[1], lty = 1, lwd = 2)
arrows(x1, effect[,5], x1, effect[,6], length=0.05, angle=90, code=3, col = col[2], lty = 2, lwd = 2)
arrows(x2, effect[,8], x2, effect[,9], length=0.05, angle=90, code=3, col = col[3], lty = 3, lwd = 2)
arrows(x3, effect[,11], x3, effect[,12], length=0.05, angle=90, code=3, col = col[4], lty = 4, lwd = 2)
arrows(x4, effect[,14], x4, effect[,15], length=0.05, angle=90, code=3, col = col[5], lty = 5, lwd = 2)
#dev.off()

#tiff("Plot2.jpg", width = 4, height = 3, units = 'in', res = 200)
par(mar = c(3, 3, 1, 1), cex.lab = 0.7, cex.axis = 0.7)
plot(0, xlim = c(0, 14), ylim = c(-0.6, 1.5), type = "n", pch = 19, col = col[1], xaxt = "n", xlab = "", ylab = "", font.lab = 1, cex.lab = 0.7, cex.axis = 0.7, cex = 0.7)
legend(x = 7, y = 0.3, legend = c("Naive", "No validation (pair)", "No validation (avg)", "DANCE (all)", "DANCE (best)"), bty = 'n', lty = 1:5, lwd = 2, col = col[1:5], cex = 0.7)
#dev.off()

#--------------------------------------------------------------------------------------------------