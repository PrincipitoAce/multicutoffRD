#' Creates NOT DONE, purpose of function
#'
#' @param X Uni-variate running variable
#' @param Y Outcome variable
#' @param C A vector of cutoff values
#' @param c.vec A vector of BASELINE cutoff values
#' @param K Folds
#' @param kk Multiplicative factor on the smoothness parameter, set cost=0 when varying kk
#' @param cost Cost value
#' @return NOT DONE
#'
#' @examples
#' optimalcutoffs(NOT DONE)
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
  W = as.numeric(X>=C) # Treatment

  n = length(Y) # sample size
  q = length(c.vec) # number of groups
  mu.fit = NULL

  datall =  data.frame(Y=Y,X=X,C=C,W=W,G=G)
  data_split = datall |>
    dplyr::mutate(fold_id=sample(1:K,size=dim(datall)[1],replace=T)) |>
    dplyr::group_by(fold_id) |>
    tidyr::nest() |>
    dplyr::arrange(fold_id)
  data_all = data_split |>
    tidyr::unnest(data) |>
    dplyr::ungroup()

  for(k in 1:K){

    # Train/test split
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
    gamfit = nnet::multinom(formula =G ~ X, data =  data_train)
    for(g in seq(1,q,1)){
      eval.dat1.m = c( data_test |>
                         dplyr::filter(X>=c.vec[g],X<c.vec[min(g+1,q)],W==0) |>
                         dplyr::select(X))$X
      eval.dat1.aug = c( data_test |>
                           dplyr::filter(X>=c.vec[g],X<c.vec[min(g+1,q)],G==g) |>
                           dplyr::select(X))$X
      eval.dat1.pseudo = c( data_test |>
                              dplyr::filter(W==1, X>=c.vec[g]) |>
                              dplyr::select(X))$X
      eval.dat1.all = c(eval.dat1.m,eval.dat1.aug, eval.dat1.pseudo)

      mu.fit1= nprobust::lprobust(data_train$Y[data_train$W==1 & data_train$G==g],data_train$X[data_train$W==1 & data_train$G==g],eval = eval.dat1.all,bwselect="imse-dpi")$Estimate[,5]


      tryCatch( { data_all[  data_all$fold_id==k & data_all$X>=c.vec[g] & data_all$X<c.vec[min(g+1,q)] & data_all$W==0, paste0("mu",".m")]= mu.fit1[1:length(eval.dat1.m)] },error=function(e) return(0))
      tryCatch( {  data_all[  data_all$fold_id==k & data_all$X>=c.vec[g] & data_all$X<c.vec[min(g+1,q)] & data_all$G==g, paste0("mu",".aug")]= mu.fit1[(length(eval.dat1.m)+1):(length(eval.dat1.m)+length(eval.dat1.aug))] },error=function(e) return(0))
      tryCatch( { data_all[  data_all$fold_id==k & data_all$W==1 & data_all$X>=c.vec[g], paste0("pseudo.",g)]= mu.fit1[c((length(eval.dat1.m)+length(eval.dat1.aug)+1):length(eval.dat1.all))]} ,error=function(e) return(0))

      eval.dat0.m = c( data_test |>
                         dplyr::filter(X>=c.vec[max(g-1,1)],X<c.vec[g],W==1) |>
                         dplyr::select(X))$X
      eval.dat0.aug = c( data_test |>
                           dplyr::filter(X>=c.vec[max(g-1,1)],X<c.vec[g],G==g) |>
                           dplyr::select(X))$X
      eval.dat0.pseudo = c( data_test |>
                              dplyr::filter(W==0, X<c.vec[g]) |>
                              dplyr::select(X))$X
      eval.dat0.all = c(eval.dat0.m,eval.dat0.aug, eval.dat0.pseudo)

      mu.fit0= nprobust::lprobust(data_train$Y[data_train$W==0 & data_train$G==g],data_train$X[data_train$W==0 & data_train$G==g],eval = eval.dat0.all,bwselect="imse-dpi")$Estimate[,5]

      tryCatch( { data_all[  data_all$fold_id==k & data_all$X>=c.vec[max(g-1,1)] & data_all$X<c.vec[g] & data_all$W==1, paste0("mu",".m")]= mu.fit0[1:length(eval.dat0.m)] },error=function(e) return(0))
      tryCatch( {  data_all[  data_all$fold_id==k & data_all$X>=c.vec[max(g-1,1)] & data_all$X<c.vec[g] & data_all$G==g, paste0("mu",".aug")]= mu.fit0[(length(eval.dat0.m)+1):(length(eval.dat0.m)+length(eval.dat0.aug))] },error=function(e) return(0))
      tryCatch( { data_all[  data_all$fold_id==k & data_all$W==0 & data_all$X<c.vec[g], paste0("pseudo.",g)]= mu.fit0[c((length(eval.dat0.m)+length(eval.dat0.aug)+1):length(eval.dat0.all))]} ,error=function(e) return(0))
    }

    ###### pseudo estimate ########

    eval.dat1.p = data_test |>
      dplyr::filter(W==1)

    tryCatch(
      { pred = predict(gamfit, newdata=   eval.dat1.p, "probs")

      if(dim( pred)[1]==1){
        data_all[data_all$fold_id==k & data_all$W==1,paste0("pseudo.ps",seq(1,q,1))]= t(as.matrix(pred ,byrow=F))
      }
      if(dim( pred)[1]!=1){
        data_all[data_all$fold_id==k & data_all$W==1,paste0("pseudo.ps",seq(1,q,1))]=pred
      }},error=function(e) return(0))

    eval.dat0.p = data_test |>
      dplyr::filter(W==0)

    tryCatch(
      { pred =predict(gamfit, newdata=   eval.dat0.p, "probs")

      if(dim( pred)[1]==1){
        data_all[data_all$fold_id==k & data_all$W==0,paste0("pseudo.ps",seq(1,q,1))]= t(as.matrix(pred ,byrow=F))
      }
      if(dim( pred)[1]!=1){
        data_all[data_all$fold_id==k & data_all$W==0,paste0("pseudo.ps",seq(1,q,1))]=pred
      }},error=function(e) return(0))

  }

  # Obtaining smoothness parameters: Lip_1temp, Lip_0temp, B.1m, & B.0m
  results <- smooth(data = data_all, c.vec_initial = c.vec)
  Lip_1temp <- results[1]
  Lip_0temp <- results[2]
  B.1m <- results[3]
  B.0m <- results[4]

  # ### NOT DONE: Option 1 Code for if Lip_1 & Lip_0 are arguments:
  # # Bounds B.1m & B.0m
  # bounds <- ulbounds(data = data_all, c.vec_initial = c.vec)
  # B.1m <- bounds[1]
  # B.0m <- bounds[2]
  #
  # # Smoothness parameter(s)
  # if (is.null(Lip_1_input)){
  #   results <- smooth(data = data_all, c.vec_initial = c.vec)
  #   Lip_1temp <- results[1]
  #   Lip_0temp <- results[2]
  # }
  # else{
  #   Lip_1temp <- Lip_1_input
  #   Lip_0temp <- Lip_0_input
  # }

  ############################################
  ########  Learning optimal cutoffs
  ############################################

    Lip_1 = kk*Lip_1temp ; Lip_0 = kk*Lip_0temp

    c.all= rep(0, length(c.vec))

    for(g in seq(1,q,1)){

      eval.dat1 = c(data_all |>
                      dplyr::filter(G==g, X>=c.vec[1], X<c.vec[q],X<c.vec[g]) |>
                      dplyr::select(X))$X #d(1)
      IND.1 = sapply(eval.dat1, function(x) sum(c.vec<x))
      eval.dat0 = c(data_all |>
                      dplyr::filter(G==g,  X>=c.vec[1], X<c.vec[q],X>=c.vec[g]) |>
                      dplyr::select(X))$X #d(0)
      IND.0 = sapply(eval.dat0, function(x) sum(c.vec<x))

      tryCatch(
        {  data_all[data_all$G==g &  data_all$X>=c.vec[1] & data_all$X<c.vec[q] & data_all$X<c.vec[g],paste0("d",1)]=
          apply( cbind(eval.dat1,IND.1),1, function(x) sum(unlist(sapply(x[2]:x[2],function(g.temp) lip_extra(x.train=x[1], group="B1", g=g, g.prim=g.temp, Lip_1_initial = Lip_1, Lip_0_initial = Lip_0, B.1m_initial = B.1m, B.0m_initial = B.0m))[2,])))
        },error=function(e) return(0))
      tryCatch(
        {  data_all[data_all$G==g &  data_all$X>=c.vec[1] & data_all$X<c.vec[q] & data_all$X>=c.vec[g],paste0("d",0)]=
          apply( cbind(eval.dat0,IND.0),1, function(x) sum(unlist(sapply((x[2]+1):(x[2]+1),function(g.temp) lip_extra(x.train=x[1],group="B0",g=g,g.prim = g.temp, Lip_1_initial = Lip_1, Lip_0_initial = Lip_0, B.1m_initial = B.1m, B.0m_initial = B.0m))[2,])) )
        },error=function(e) return(0))

    }

    data_mid = data_all |>
      dplyr::filter(X>=min(c.vec),X<max(c.vec))

    regret_sum=NULL
    for(g in seq(1,q,1)){
      regret=NULL
      for( c.alt in unique(X[X>=c.vec[1]&X<c.vec[q]]) ){
        if(c.alt>=c.vec[g]){
          temp1= tryCatch(-sum(data_mid[data_mid $X>=c.vec[g] & data_mid $X<c.alt & data_mid $G==g,"Y"])/n, error=function(e) return(0))
          ###########
          dat.temp = data_mid |>
            dplyr::filter(G==g, X<c.alt, X>=c.vec[g])

          tempDB1=
            tryCatch( sum( dat.temp[,"mu.m"] )  /n, error=function(e) return(0))

          tempd =   tryCatch( sum( dat.temp[,paste0("d",0)])/n, error=function(e) return(0))
          ###########
          dat.temp = data_mid |>
            dplyr::filter( X<c.alt, X>=c.vec[g], X>=c.vec[ifelse(G==1,1,G-1)],X<c.vec[G] )  #& X>=c.vec[G-1] & X<c.vec[G]


          tempDB2=
            tryCatch( sum( with(dat.temp, eval(parse(text =paste0("pseudo.ps",g)))/eval(parse(text =paste0("pseudo.ps",G)))*(Y-eval(parse(text ="mu.aug"))) )
            )/n, error=function(e) return(0))

          tempcost= tryCatch(cost*dim(data_mid[data_mid $X>=c.vec[g] & data_mid $X<c.alt & data_mid $G==g,"Y"])[1]/n, error=function(e) return(0))

          temp.reg=temp1+tempDB1+tempd+tempDB2 +tempcost
        }
        if(c.alt<c.vec[g]){
          temp1= tryCatch(-sum(data_mid[data_mid $X<c.vec[g] & data_mid $X>=c.alt & data_mid $G==g,"Y"])/n, error=function(e) return(0))
          ###########
          dat.temp = data_mid |>
            dplyr::filter(G==g, X>=c.alt, X<c.vec[g])


          tempDB1=
            tryCatch( sum( dat.temp[,"mu.m"] )  /n, error=function(e) return(0))

          tempd =   tryCatch( sum( dat.temp[,paste0("d",1)])/n, error=function(e) return(0))

          dat.temp = data_mid |>
            dplyr::filter( X>=c.alt, X<c.vec[g], X>=c.vec[G],X<c.vec[ifelse(G==q,q,G+1)] )

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

    # NOT DONE: Check if c.all aligned with correct Groups
    # Returns data frame of original and new cutoff values
    result = data.frame(unlist(c.vec), unlist(c.all))
    colnames(result) <- c("Original Cutoff", "New Cutoff")
    return(result)

}
