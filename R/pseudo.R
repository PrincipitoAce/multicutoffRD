#' Creates NOT DONE, purpose of function
#'
#' @param data_set A data set with longitude and latitude variables
#'
#' @return NOT DONE
#'
#' @examples
#' function1(NOT DONE)
#'
#' @export

pseudo = function(data_train, data_test, c.vec, G){
  for(k in 1:K){
    data_train = data_split %>%
      filter(fold_id!=k) %>%
      unnest(data) %>%
      ungroup() %>%
      select(-fold_id)
    data_test = data_split %>%
      filter(fold_id==k) %>%
      unnest(data) %>%
      ungroup() %>%
      select(-fold_id)

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

    ################################################################################

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
}
