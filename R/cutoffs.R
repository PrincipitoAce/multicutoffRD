#' Creates NOT DONE, purpose of function
#'
#' @param X_in X value
#' @param Y_in Y value
#' @param C_in C value
#' @param c.vec A vector of cutoff values
#' @param kk Multiplicative factor
#' @param cost cost
#' @param K K-folds, default is 20
#' @param Lip_0temp Smoothness
#' @param Lip_1temp Smoothness
#' @param B.0m B0
#' @param B.1m B1
#' @return NOT DONE: Smoothness parameters
#'
#' @examples
#' smoother(NOT DONE)
#'
#' @export

cutoffs = function(X_in, Y_in, C_in, c.vec, kk, cost, K=20, Lip_0temp, Lip_1temp, B.0m, B.1m){

  G_in = match(C_in,c.vec)  # Group index
  D_in = as.numeric(X_in>=C_in) # Treatment
  n=length(Y_in) # sample size
  q=length(c.vec) # number of groups

  data_all = data.frame(X=X_in, Y=Y_in, C=C_in, D=D_in, G=G_in)

  ##############################################################################

  lip_extra = function(x.train,group,g,g.prim){ # extrapolation function

    if(group=="B1"){ #B1 G=1
      d=1;Lip=Lip_1temp[g,g.prim]
      B.m=B.1m[g,g.prim];
      eval.main = unique(C[G == max(g,g.prim)])
    }
    if(group=="B0"){ #B1 G=1
      d=0;Lip=Lip_0temp[g,g.prim]
      B.m=B.0m[g,g.prim];
      eval.main = unique(C[G == min(g,g.prim)])
    }

    B.up=  B.m
    B.low= B.m

    upper = sapply(x.train,function(x_prime) min(1, min( B.m   +  Lip * abs(x_prime - eval.main)) ))
    lower= sapply(x.train,function(x_prime) max(-1, max( B.m - Lip * abs(x_prime - eval.main) )) )
    return(list(upper = upper, lower = lower ))
  }

  ############################################
  ########  Learning optimal cutoffs
  ############################################

    Lip_1 = kk*Lip_1temp ; Lip_0 = kk*Lip_0temp

    c.all= rep(0,length(c.vec))

    for(g in seq(1,q,1)){

      eval.dat1 = c(data_all %>% filter(G==g, X>=c.vec[1], X<c.vec[q],X<c.vec[g]) %>% select(X))$X #d(1)
      IND.1 = sapply(eval.dat1, function(x) sum(c.vec<x))
      eval.dat0 = c(data_all %>% filter(G==g,  X>=c.vec[1], X<c.vec[q],X>=c.vec[g]) %>% select(X))$X #d(0)
      IND.0 = sapply(eval.dat0, function(x) sum(c.vec<x))

      tryCatch(
        {  data_all[data_all$G==g &  data_all$X>=c.vec[1] & data_all$X<c.vec[q] & data_all$X<c.vec[g],paste0("d",1)]=
          apply( cbind(eval.dat1,IND.1),1, function(x) sum(unlist(sapply(x[2]:x[2],function(g.temp) lip_extra(x.train=x[1],group="B1",g=g,g.prim = g.temp))[2,])))
        },error=function(e) return(0))
      tryCatch(
        {  data_all[data_all$G==g &  data_all$X>=c.vec[1] & data_all$X<c.vec[q] & data_all$X>=c.vec[g],paste0("d",0)]=
          apply( cbind(eval.dat0,IND.0),1, function(x) sum(unlist(sapply((x[2]+1):(x[2]+1),function(g.temp) lip_extra(x.train=x[1],group="B0",g=g,g.prim = g.temp))[2,])) )
        },error=function(e) return(0))

    }

    data_mid = data_all %>% filter(X>=min(c.vec),X<max(c.vec))

    regret_sum=NULL

    for(g in seq(1,q,1)){
      regret=NULL
      for( c.alt in unique(X[X>=c.vec[1]&X<c.vec[q]]) ){
        if(c.alt>=c.vec[g]){
          temp1= tryCatch(-sum(data_mid[data_mid $X>=c.vec[g] & data_mid $X<c.alt & data_mid $G==g,"Y"])/n, error=function(e) return(0))
          ###########
          dat.temp = data_mid %>% filter(G==g, X<c.alt, X>=c.vec[g])

          tempDB1=
            tryCatch( sum( dat.temp[,"mu.m"] )  /n, error=function(e) return(0))

          tempd =   tryCatch( sum( dat.temp[,paste0("d",0)])/n, error=function(e) return(0))
          ###########
          dat.temp = data_mid %>% filter( X<c.alt, X>=c.vec[g], X>=c.vec[ifelse(G==1,1,G-1)],X<c.vec[G] )  #& X>=c.vec[G-1] & X<c.vec[G]


          tempDB2=
            tryCatch( sum( with(dat.temp, eval(parse(text =paste0("pseudo.ps",g)))/eval(parse(text =paste0("pseudo.ps",G)))*(Y-eval(parse(text ="mu.aug"))) )
            )/n, error=function(e) return(0))

          tempcost= tryCatch(cost*dim(data_mid[data_mid $X>=c.vec[g] & data_mid $X<c.alt & data_mid $G==g,"Y"])[1]/n, error=function(e) return(0))

          temp.reg=temp1+tempDB1+tempd+tempDB2 +tempcost
        }
        if(c.alt<c.vec[g]){
          temp1= tryCatch(-sum(data_mid[data_mid $X<c.vec[g] & data_mid $X>=c.alt & data_mid $G==g,"Y"])/n, error=function(e) return(0))
          ###########
          dat.temp = data_mid %>% filter(G==g, X>=c.alt, X<c.vec[g])


          tempDB1=
            tryCatch( sum( dat.temp[,"mu.m"] )  /n, error=function(e) return(0))

          tempd =   tryCatch( sum( dat.temp[,paste0("d",1)])/n, error=function(e) return(0))

          dat.temp = data_mid %>% filter( X>=c.alt, X<c.vec[g], X>=c.vec[G],X<c.vec[ifelse(G==q,q,G+1)] )

          tempDB2=
            tryCatch( sum( with(dat.temp, eval(parse(text =paste0("pseudo.ps",g)))/eval(parse(text =paste0("pseudo.ps",G)))*(Y-eval(parse(text ="mu.aug"))) ) )/n, error=function(e) return(0))

          tempcost= tryCatch(cost*dim(data_mid[data_mid $X<c.vec[g] & data_mid $X>=c.alt & data_mid $G==g,"Y"])[1]/n, error=function(e) return(0))

          temp.reg=temp1+tempDB1+tempd+tempDB2 -tempcost

        }
        regret=c(regret,temp.reg)

      }

      if(max(regret)==0){
        c.all[g]=c.vec[g]
      }else{
        c.all[g]= unique(X[X>=c.vec[1]&X<c.vec[q]])[which(regret==max(regret))[1]]
      }
      regret_sum=c(regret_sum,max(regret))
    }

  ##############################################################################
  # Returns data frame of original cutoff values, new cutoff values, & difference
  diffs = c.vec - c.all
  result = data.frame(unlist(c.vec), unlist(c.all), unlist(diffs))
  colnames(result) <- c("Original Cutoffs", "New Cutoffs", "Difference")
  return(result)
}
