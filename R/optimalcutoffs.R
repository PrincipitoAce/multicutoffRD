#' Creates NOT DONE, purpose of function
#'
#' @param X Uni-variate running variable
#' @param Y Outcome variable
#' @param C A vector of cutoff values
#' @param c.vec A vector of BASELINE cutoff values
#' @param K Folds
#' @param kk Multiplicative factor on the smoothness parameter, set cost=0 when varying kk
#' @param cost Cost value
#' @return Data frame of new, optimal cutoff values alongside original cutoff values
#'
#' @examples
#' optimalcutoffs(X=X, Y=Y, C=C, c.vec=original_cutoffs, K=20, kk=1, cost=0.2)
#'
#' @export

optimalcutoffs = function(X, Y, C, c.vec, K=20, kk=1, cost=0){
  # Class checks
  if (!is.numeric(X)){
    stop("Error: X must be numeric")
  }
  if (!is.numeric(Y)){
    stop("Error: Y must be numeric")
  }
  if (!is.numeric(c.vec)){
    stop("Error: c.vec must be numeric")
  }
  if (!is.numeric(K)){
    stop("Error: K must be numeric")
  }
  if (!is.numeric(kk)){
    stop("Error: kk must be numeric")
  }
  if (!is.numeric(cost)){
    stop("Error: cost must be numeric")
  }

  # If c.vec is unsorted, then sort
  if (is.unsorted(c.vec)){
    # Sort c.vec
    c.vec <- sort(c.vec, decreasing = TRUE)
  }

  G = match(C, c.vec)  # Group index
  D = as.numeric(X>=C) # Treatment
  n = length(Y) # Sample size
  q = length(c.vec) # Number of groups
  mu.fit = NULL

  datall =  data.frame(Y=Y,X=X,C=C,D=D,G=G)
  data_split = datall %>% dplyr::mutate(fold_id=sample(1:K,size=dim(datall)[1],replace=T)) %>% group_by(fold_id) %>%
    tidyr::nest() %>% arrange(fold_id)
  data_all = data_split  %>% tidyr::unnest(data) %>% ungroup()

  mu.fit=NULL
  for(k in 1:K){
    data_train = data_split %>% filter(fold_id!=k) %>% tidyr::unnest(data) %>% ungroup() %>% dplyr::select(-fold_id)
    data_test = data_split %>% filter(fold_id==k) %>% tidyr::unnest(data) %>% ungroup() %>% dplyr::select(-fold_id)

    # conditional prob
    gamfit= nnet::multinom(formula =G ~ X, data =  data_train)
    for(g in seq(1,q,1)){
      eval.dat1.m = c( data_test %>% filter(X>=c.vec[g],X<c.vec[min(g+1,q)],D==0) %>% dplyr::select(X))$X
      eval.dat1.aug = c( data_test %>% filter(X>=c.vec[g],X<c.vec[min(g+1,q)],G==g) %>% dplyr::select(X))$X
      eval.dat1.pseudo = c( data_test %>% filter(D==1, X>=c.vec[g]) %>% dplyr::select(X))$X
      eval.dat1.all = c(eval.dat1.m,eval.dat1.aug, eval.dat1.pseudo)

      mu.fit1= nprobust::lprobust(data_train$Y[data_train$D==1 & data_train$G==g],data_train$X[data_train$D==1 & data_train$G==g],eval = eval.dat1.all,bwselect="imse-dpi")$Estimate[,5]


      tryCatch( { data_all[  data_all$fold_id==k & data_all$X>=c.vec[g] & data_all$X<c.vec[min(g+1,q)] & data_all$D==0, paste0("mu",".m")]= mu.fit1[1:length(eval.dat1.m)] },error=function(e) return(0))
      tryCatch( {  data_all[  data_all$fold_id==k & data_all$X>=c.vec[g] & data_all$X<c.vec[min(g+1,q)] & data_all$G==g, paste0("mu",".aug")]= mu.fit1[(length(eval.dat1.m)+1):(length(eval.dat1.m)+length(eval.dat1.aug))] },error=function(e) return(0))
      tryCatch( { data_all[  data_all$fold_id==k & data_all$D==1 & data_all$X>=c.vec[g], paste0("pseudo.",g)]= mu.fit1[c((length(eval.dat1.m)+length(eval.dat1.aug)+1):length(eval.dat1.all))]} ,error=function(e) return(0))

      eval.dat0.m = c( data_test %>% filter(X>=c.vec[max(g-1,1)],X<c.vec[g],D==1) %>% dplyr::select(X))$X
      eval.dat0.aug = c( data_test %>% filter(X>=c.vec[max(g-1,1)],X<c.vec[g],G==g) %>% dplyr::select(X))$X
      eval.dat0.pseudo = c( data_test %>% filter(D==0, X<c.vec[g]) %>% dplyr::select(X))$X
      eval.dat0.all = c(eval.dat0.m,eval.dat0.aug, eval.dat0.pseudo)

      mu.fit0= nprobust::lprobust(data_train$Y[data_train$D==0 & data_train$G==g],data_train$X[data_train$D==0 & data_train$G==g],eval = eval.dat0.all,bwselect="imse-dpi")$Estimate[,5]

      tryCatch( { data_all[  data_all$fold_id==k & data_all$X>=c.vec[max(g-1,1)] & data_all$X<c.vec[g] & data_all$D==1, paste0("mu",".m")]= mu.fit0[1:length(eval.dat0.m)] },error=function(e) return(0))
      tryCatch( {  data_all[  data_all$fold_id==k & data_all$X>=c.vec[max(g-1,1)] & data_all$X<c.vec[g] & data_all$G==g, paste0("mu",".aug")]= mu.fit0[(length(eval.dat0.m)+1):(length(eval.dat0.m)+length(eval.dat0.aug))] },error=function(e) return(0))
      tryCatch( { data_all[  data_all$fold_id==k & data_all$D==0 & data_all$X<c.vec[g], paste0("pseudo.",g)]= mu.fit0[c((length(eval.dat0.m)+length(eval.dat0.aug)+1):length(eval.dat0.all))]} ,error=function(e) return(0))
    }
    ###### pseudo estimate ########

    eval.dat1.p = data_test %>% filter(D==1)
    tryCatch(
      { pred =predict(gamfit, newdata=   eval.dat1.p, "probs")

      if(dim( pred)[1]==1){
        data_all[data_all$fold_id==k & data_all$D==1,paste0("pseudo.ps",seq(1,q,1))]= t(as.matrix(pred ,byrow=F))
      }
      if(dim( pred)[1]!=1){
        data_all[data_all$fold_id==k & data_all$D==1,paste0("pseudo.ps",seq(1,q,1))]=pred
      }},error=function(e) return(0))
    eval.dat0.p = data_test %>% filter(D==0)
    tryCatch(
      { pred =predict(gamfit, newdata=   eval.dat0.p, "probs")

      if(dim( pred)[1]==1){
        data_all[data_all$fold_id==k & data_all$D==0,paste0("pseudo.ps",seq(1,q,1))]= t(as.matrix(pred ,byrow=F))
      }
      if(dim( pred)[1]!=1){
        data_all[data_all$fold_id==k & data_all$D==0,paste0("pseudo.ps",seq(1,q,1))]=pred
      }},error=function(e) return(0))

  }

  ## Second, (1) estimate cross-group differences and (2) choose the value of smoothness parameter


  psd_dat1 = psd_dat0 =NULL;
  Lip_1=Lip_0=matrix(0,q,q) # storing the value of smoothness parameter;  1/0: treatment/control
  B.1m=B.0m=matrix(0,nrow=q,ncol=q) # storing the value of estimated cross-group differences at cutoff point

  for(g in seq(1,q-1,1))
    for(g.pr in seq(g+1,q,1)){


      temp.dat = data_all %>% filter(D==1 & X>=max(c.vec[g.pr],c.vec[g]))
      temp.vc =
        data.frame( temp.dat[,paste0("pseudo.",g)]-temp.dat[,paste0("pseudo.",g.pr)]+
                      with(temp.dat,I(G==g)*(Y-eval(parse(text =paste0("pseudo.",g))))/eval(parse(text =paste0("pseudo.ps",g))) )-
                      with(temp.dat,I(G==g.pr)*(Y-eval(parse(text =paste0("pseudo.",g.pr))))/eval(parse(text =paste0("pseudo.ps",g.pr))) )
                    ,  temp.dat $X,g,g.pr)

      names(temp.vc)[1:2]=c("psout","X")
      psd_dat1=rbind(psd_dat1,   temp.vc )

      # Section 4.3
      Lip_1[g,g.pr]=abs(nprobust::lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = max(c.vec[g.pr],c.vec[g]),deriv = 1,p=2,bwselect="mse-dpi")$Estimate[,5])
      # Algorithm 1
      B.1m[g,g.pr]=nprobust::lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = max(c.vec[g.pr],c.vec[g]),bwselect="mse-dpi")$Estimate[,5]

      temp.dat = data_all %>% filter(D==0 & X<min(c.vec[g.pr],c.vec[g]))

      temp.vc = data.frame("psout"=temp.dat[,paste0("pseudo.",g)]-temp.dat[,paste0("pseudo.",g.pr)]+
                             with(temp.dat,I(G==g)*(Y-eval(parse(text =paste0("pseudo.",g))))/eval(parse(text =paste0("pseudo.ps",g))) )-
                             with(temp.dat,I(G==g.pr)*(Y-eval(parse(text =paste0("pseudo.",g.pr))))/eval(parse(text =paste0("pseudo.ps",g.pr))) ),   temp.dat $X,g,g.pr)
      names(temp.vc)[1:2]=c("psout","X")
      psd_dat0=rbind(psd_dat0,   temp.vc )

      Lip_0[g,g.pr]=abs(nprobust::lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = min(c.vec[g.pr],c.vec[g]),deriv = 1,p=2,bwselect="mse-dpi")$Estimate[,5])
      B.0m[g,g.pr]=nprobust::lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = min(c.vec[g.pr],c.vec[g]),bwselect="mse-dpi")$Estimate[,5]

    }
  Lip_1=Lip_1+t(Lip_1);Lip_0=Lip_0+t(Lip_0)
  B.1m= B.1m + t(-B.1m) ; B.0m= B.0m + t(-B.0m)

  # initial values for smoothness parameter
  Lip_1temp=Lip_1;Lip_0temp=Lip_0


  lip_extra = function(x.train,group,g,g.prim){ # extrapolation function

    if(group=="B1"){ #B1 G=1
      d=1;Lip=Lip_1[g,g.prim]
      B.m=B.1m[g,g.prim];
      eval.main = unique(C[G == max(g,g.prim)])
    }
    if(group=="B0"){ #B1 G=1
      d=0;Lip=Lip_0[g,g.prim]
      B.m=B.0m[g,g.prim];
      eval.main = unique(C[G == min(g,g.prim)])
    }

    B.up=  B.m
    B.low= B.m

    upper = sapply(x.train,function(x_prime) min(1, min( B.m   +  Lip * abs(x_prime - eval.main)) ))
    lower= sapply(x.train,function(x_prime) max(-1, max( B.m - Lip * abs(x_prime - eval.main) )) )
    return(list(upper = upper, lower = lower ))
  }

  Lip_1 = kk*Lip_1temp ; Lip_0 = kk*Lip_0temp

  c.all= rep(0,length(c.vec))

  for(g in seq(1,q,1)){

    eval.dat1 = c(data_all %>% filter(G==g, X>=c.vec[1], X<c.vec[q],X<c.vec[g]) %>% dplyr::select(X))$X #d(1)
    IND.1 = sapply(eval.dat1, function(x) sum(c.vec<x))
    eval.dat0 = c(data_all %>% filter(G==g,  X>=c.vec[1], X<c.vec[q],X>=c.vec[g]) %>% dplyr::select(X))$X #d(0)
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

  # Returns data frame of original cutoff values, new cutoff values, & difference
  diffs = c.vec - c.all
  result = data.frame(unlist(c.vec), unlist(c.all), unlist(diffs))
  colnames(result) <- c("Original Cutoffs", "New Cutoffs", "Difference")
  lip_b <- list(Lip_0temp, Lip_1temp, B.0m, B.1m)
  return(lip_b)

}
