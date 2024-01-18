#' Creates NOT DONE, purpose of function
#'
#' @param X X value
#' @param Y Y value
#' @param C C value
#' @param c.vec A vector of cutoff values
#' @param K K folds
#' @return NOT DONE: Smoothness parameters
#'
#' @examples
#' smooth(NOT DONE)
#'
#' @export

smooth = function(X, Y, C, c.vec, K){

  G = match(C,c.vec)  # Group index
  D = as.numeric(X>=C) # Treatment

  n=length(Y) # sample size
  q=length(c.vec) # number of groups

  datall =  data.frame(Y=Y,X=X,C=C,D=D,G=G)
  data_split = datall %>% mutate(fold_id=sample(1:K,size=dim(datall)[1],replace=T)) %>% group_by(fold_id) %>%
    nest() %>% arrange(fold_id)
  data_all = data_split  %>% unnest(data) %>% ungroup()
  ##############################################################################

  for(k in 1:K){
    data_train = data_split %>% filter(fold_id!=k) %>% unnest(data) %>% ungroup() %>% select(-fold_id)
    data_test = data_split %>% filter(fold_id==k) %>% unnest(data) %>% ungroup() %>% select(-fold_id)

    # conditional prob
    gamfit= multinom(formula =G ~ X, data =  data_train)
    for(g in seq(1,q,1)){
      eval.dat1.m = c( data_test %>% filter(X>=c.vec[g],X<c.vec[min(g+1,q)],D==0) %>% select(X))$X
      eval.dat1.aug = c( data_test %>% filter(X>=c.vec[g],X<c.vec[min(g+1,q)],G==g) %>% select(X))$X
      eval.dat1.pseudo = c( data_test %>% filter(D==1, X>=c.vec[g]) %>% select(X))$X
      eval.dat1.all = c(eval.dat1.m,eval.dat1.aug, eval.dat1.pseudo)

      mu.fit1= lprobust(data_train$Y[data_train$D==1 & data_train$G==g],data_train$X[data_train$D==1 & data_train$G==g],eval = eval.dat1.all,bwselect="imse-dpi")$Estimate[,5]


      tryCatch( { data_all[  data_all$fold_id==k & data_all$X>=c.vec[g] & data_all$X<c.vec[min(g+1,q)] & data_all$D==0, paste0("mu",".m")]= mu.fit1[1:length(eval.dat1.m)] },error=function(e) return(0))
      tryCatch( {  data_all[  data_all$fold_id==k & data_all$X>=c.vec[g] & data_all$X<c.vec[min(g+1,q)] & data_all$G==g, paste0("mu",".aug")]= mu.fit1[(length(eval.dat1.m)+1):(length(eval.dat1.m)+length(eval.dat1.aug))] },error=function(e) return(0))
      tryCatch( { data_all[  data_all$fold_id==k & data_all$D==1 & data_all$X>=c.vec[g], paste0("pseudo.",g)]= mu.fit1[c((length(eval.dat1.m)+length(eval.dat1.aug)+1):length(eval.dat1.all))]} ,error=function(e) return(0))

      eval.dat0.m = c( data_test %>% filter(X>=c.vec[max(g-1,1)],X<c.vec[g],D==1) %>% select(X))$X
      eval.dat0.aug = c( data_test %>% filter(X>=c.vec[max(g-1,1)],X<c.vec[g],G==g) %>% select(X))$X
      eval.dat0.pseudo = c( data_test %>% filter(D==0, X<c.vec[g]) %>% select(X))$X
      eval.dat0.all = c(eval.dat0.m,eval.dat0.aug, eval.dat0.pseudo)

      mu.fit0= lprobust(data_train$Y[data_train$D==0 & data_train$G==g],data_train$X[data_train$D==0 & data_train$G==g],eval = eval.dat0.all,bwselect="imse-dpi")$Estimate[,5]

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
      Lip_1[g,g.pr]=abs(lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = max(c.vec[g.pr],c.vec[g]),deriv = 1,p=2,bwselect="mse-dpi")$Estimate[,5])
      # Algorithm 1
      B.1m[g,g.pr]=lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = max(c.vec[g.pr],c.vec[g]),bwselect="mse-dpi")$Estimate[,5]

      temp.dat = data_all %>% filter(D==0 & X<min(c.vec[g.pr],c.vec[g]))

      temp.vc = data.frame("psout"=temp.dat[,paste0("pseudo.",g)]-temp.dat[,paste0("pseudo.",g.pr)]+
                             with(temp.dat,I(G==g)*(Y-eval(parse(text =paste0("pseudo.",g))))/eval(parse(text =paste0("pseudo.ps",g))) )-
                             with(temp.dat,I(G==g.pr)*(Y-eval(parse(text =paste0("pseudo.",g.pr))))/eval(parse(text =paste0("pseudo.ps",g.pr))) ),   temp.dat $X,g,g.pr)
      names(temp.vc)[1:2]=c("psout","X")
      psd_dat0=rbind(psd_dat0,   temp.vc )

      Lip_0[g,g.pr]=abs(lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = min(c.vec[g.pr],c.vec[g]),deriv = 1,p=2,bwselect="mse-dpi")$Estimate[,5])
      B.0m[g,g.pr]=lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = min(c.vec[g.pr],c.vec[g]),bwselect="mse-dpi")$Estimate[,5]

    }
  Lip_1=Lip_1+t(Lip_1);Lip_0=Lip_0+t(Lip_0)
  B1= B.1m + t(-B.1m) ; B0= B.0m + t(-B.0m)

  results <- c(Lip_0, Lip_1, B0, B1)
  return(results)
}
