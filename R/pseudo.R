#' Creates NOT DONE, purpose of function
#'
#' @param datall A data set with variables X, Y, C, W, G
#' @param K The number of folds desired, the default is K=20
#' @param c.vec A vector of cutoff values
#' @return NOT DONE
#'
#' @examples
#' pseudo(NOT DONE)
#'
#' @export

pseudo = function(datall, K=20, c.vec){
  # Numeric class checks
  if (!is.numeric(K)){
    stop("Error: K must be numeric")
  }

  data_split = datall |>
    dplyr::mutate(fold_id=sample(1:K,size=dim(datall)[1],replace=T)) |>
    dplyr::group_by(fold_id) |>
    tidyr::nest() |>
    dplyr::arrange(fold_id)

  data_all = data_split |>
    tidyr::unnest(data) |>
    dplyr::ungroup()

  mu.fit = NULL
  q = length(c.vec)

  for(k in 1:K){
    data_train = data_split |>
      dplyr::filter(fold_id!=k) |>
      tidyr::unnest(data) |>
      dplyr::ungroup() |>
      dplyr::select(-fold_id)
    data_test = data_split |>
      dplyr::filter(fold_id==k) |>
      tidyr::unnest(data) |>
      dplyr::ungroup() |>
      dplyr::select(-fold_id)

      # conditional prob
      gamfit = nnet::multinom(formula = G ~ X, data =  data_train)
      for(g in seq(1,q,1)){
        eval.dat1.m = c( data_test |> dplyr::filter(X>=c.vec[g],X<c.vec[min(g+1,q)],W==0) |> dplyr::select(X))$X
        eval.dat1.aug = c( data_test |> dplyr::filter(X>=c.vec[g],X<c.vec[min(g+1,q)],G==g) |> dplyr::select(X))$X
        eval.dat1.pseudo = c( data_test |> dplyr::filter(W==1, X>=c.vec[g]) |> dplyr::select(X))$X
        eval.dat1.all = c(eval.dat1.m,eval.dat1.aug, eval.dat1.pseudo)

        # Calonico, Cattaneo and Farrell (2019): https://nppackages.github.io/
        mu.fit1= nprobust::lprobust(data_train$Y[data_train$W==1 & data_train$G==g],data_train$X[data_train$W==1 & data_train$G==g],eval = eval.dat1.all,bwselect="imse-dpi")$Estimate[,5]

        tryCatch( { data_all[  data_all$fold_id==k & data_all$X>=c.vec[g] & data_all$X<c.vec[min(g+1,q)] & data_all$W==0, paste0("mu",".m")]= mu.fit1[1:length(eval.dat1.m)] },error=function(e) return(0))
        tryCatch( {  data_all[  data_all$fold_id==k & data_all$X>=c.vec[g] & data_all$X<c.vec[min(g+1,q)] & data_all$G==g, paste0("mu",".aug")]= mu.fit1[(length(eval.dat1.m)+1):(length(eval.dat1.m)+length(eval.dat1.aug))] },error=function(e) return(0))
        tryCatch( { data_all[  data_all$fold_id==k & data_all$W==1 & data_all$X>=c.vec[g], paste0("pseudo.",g)]= mu.fit1[c((length(eval.dat1.m)+length(eval.dat1.aug)+1):length(eval.dat1.all))]} ,error=function(e) return(0))

        eval.dat0.m = c( data_test |> filter(X>=c.vec[max(g-1,1)],X<c.vec[g],W==1) |> select(X))$X
        eval.dat0.aug = c( data_test |> filter(X>=c.vec[max(g-1,1)],X<c.vec[g],G==g) |> select(X))$X
        eval.dat0.pseudo = c( data_test |> filter(W==0, X<c.vec[g]) |> select(X))$X
        eval.dat0.all = c(eval.dat0.m,eval.dat0.aug, eval.dat0.pseudo)

        mu.fit0= nprobust::lprobust(data_train$Y[data_train$W==0 & data_train$G==g],data_train$X[data_train$W==0 & data_train$G==g],eval = eval.dat0.all,bwselect="imse-dpi")$Estimate[,5]

        tryCatch( { data_all[  data_all$fold_id==k & data_all$X>=c.vec[max(g-1,1)] & data_all$X<c.vec[g] & data_all$W==1, paste0("mu",".m")]= mu.fit0[1:length(eval.dat0.m)] },error=function(e) return(0))
        tryCatch( {  data_all[  data_all$fold_id==k & data_all$X>=c.vec[max(g-1,1)] & data_all$X<c.vec[g] & data_all$G==g, paste0("mu",".aug")]= mu.fit0[(length(eval.dat0.m)+1):(length(eval.dat0.m)+length(eval.dat0.aug))] },error=function(e) return(0))
        tryCatch( { data_all[  data_all$fold_id==k & data_all$W==0 & data_all$X<c.vec[g], paste0("pseudo.",g)]= mu.fit0[c((length(eval.dat0.m)+length(eval.dat0.aug)+1):length(eval.dat0.all))]} ,error=function(e) return(0))
      }

    ###### pseudo estimate ########

    eval.dat1.p = data_test |> dplyr::filter(W==1)
    tryCatch(
      { pred =predict(gamfit, newdata=   eval.dat1.p, "probs")

      if(dim( pred)[1]==1){
        data_all[data_all$fold_id==k & data_all$W==1,paste0("pseudo.ps",seq(1,q,1))]= t(as.matrix(pred ,byrow=F))
      }
      if(dim( pred)[1]!=1){
        data_all[data_all$fold_id==k & data_all$W==1,paste0("pseudo.ps",seq(1,q,1))]=pred
      }},error=function(e) return(0))
    eval.dat0.p = data_test |> dplyr::filter(W==0)
    tryCatch(
      { pred =predict(gamfit, newdata=   eval.dat0.p, "probs")

      if(dim( pred)[1]==1){
        data_all[data_all$fold_id==k & data_all$W==0,paste0("pseudo.ps",seq(1,q,1))]= t(as.matrix(pred ,byrow=F))
      }
      if(dim( pred)[1]!=1){
        data_all[data_all$fold_id==k & data_all$W==0,paste0("pseudo.ps",seq(1,q,1))]=pred
      }},error=function(e) return(0))

  }
}
