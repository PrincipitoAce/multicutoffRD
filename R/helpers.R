#' Creates NOT DONE, purpose of function
#'
#' @param data A data set with variables X, Y, C, W, G, and pseudo outcomes
#' @param c.vec_initial A vector of cutoff values
#' @return NOT DONE: Smoothness parameters
#'
#' @examples
#' smooth(NOT DONE)
#'
#' @export

smooth = function(data, c.vec_initial){
  q = length(c.vec_initial)
  Lip_1 = Lip_0 = matrix(0,q,q) # storing the value of smoothness parameter;  1/0: treatment/control
  B.1m = B.0m = matrix(0,nrow=q,ncol=q) # storing the value of estimated cross-group differences at cutoff point

  for(g in seq(1,q-1,1))
    for(g.pr in seq(g+1,q,1)){

      temp.dat = data |> dplyr::filter(W==1 & X>=max(c.vec_initial[g.pr], c.vec_initial[g]))
      temp.vc =
        data.frame( temp.dat[,paste0("pseudo.",g)]-temp.dat[,paste0("pseudo.",g.pr)]+
                      with(temp.dat,I(G==g)*(Y-eval(parse(text =paste0("pseudo.",g))))/eval(parse(text =paste0("pseudo.ps",g))) )-
                      with(temp.dat,I(G==g.pr)*(Y-eval(parse(text =paste0("pseudo.",g.pr))))/eval(parse(text =paste0("pseudo.ps",g.pr))) )
                    ,  temp.dat $X,g,g.pr)

      names(temp.vc)[1:2]=c("psout","X")

      # Section 4.3
      Lip_1[g,g.pr] = abs(nprobust::lprobust(temp.vc[,"psout"], temp.vc[,"X"], eval = max(c.vec_initial[g.pr], c.vec_initial[g]), deriv = 1, p=2, bwselect="mse-dpi")$Estimate[,5])

      # Algorithm 1
      B.1m[g,g.pr] = nprobust::lprobust(temp.vc[,"psout"], temp.vc[,"X"], eval = max(c.vec_initial[g.pr], c.vec_initial[g]), bwselect="mse-dpi")$Estimate[,5]

      temp.dat = data |> dplyr::filter(W==0 & X<min(c.vec_initial[g.pr],c.vec_initial[g]))

      temp.vc = data.frame("psout"=temp.dat[,paste0("pseudo.",g)]-temp.dat[,paste0("pseudo.",g.pr)]+
                             with(temp.dat,I(G==g)*(Y-eval(parse(text =paste0("pseudo.",g))))/eval(parse(text =paste0("pseudo.ps",g))) )-
                             with(temp.dat,I(G==g.pr)*(Y-eval(parse(text =paste0("pseudo.",g.pr))))/eval(parse(text =paste0("pseudo.ps",g.pr))) ),   temp.dat $X,g,g.pr)
      names(temp.vc)[1:2]=c("psout","X")

      Lip_0[g,g.pr] = abs(nprobust::lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = min(c.vec_initial[g.pr],c.vec_initial[g]),deriv = 1,p=2,bwselect="mse-dpi")$Estimate[,5])
      B.0m[g,g.pr] = nprobust::lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = min(c.vec_initial[g.pr],c.vec_initial[g]),bwselect="mse-dpi")$Estimate[,5]

    }
  Lip_1=Lip_1+t(Lip_1);Lip_0=Lip_0+t(Lip_0)
  B.1m= B.1m + t(-B.1m) ; B.0m= B.0m + t(-B.0m)
  result <- c(Lip_1, Lip_0, B.1m, B.0m)
  return(result)
}

#' Creates NOT DONE, purpose of function e.g. extrapolation function
#'
#' @param x.train ...
#' @param group ...
#' @param g ...
#' @param g.prim ...
#' @param Lip_1_initial ...
#' @param Lip_0_initial ...
#' @param B.1m_initial ...
#' @param B.0m_initial ...
#' @return NOT DONE
#'
#' @examples
#' lip_extra(NOT DONE)
#'
#' @export

lip_extra = function(x.train, group, g, g.prim, Lip_1_initial, Lip_0_initial, B.1m_initial, B.0m_initial){

  if(group=="B1"){ #B1 G=1
    d = 1; Lip = Lip_1_initial[g,g.prim]
    B.m = B.1m_initial[g,g.prim];
    eval.main = unique(C[G == max(g,g.prim)])
  }
  if(group=="B0"){ #B0 G=1
    d = 0; Lip = Lip_0_initial[g,g.prim]
    B.m = B.0m_initial[g,g.prim];
    eval.main = unique(C[G == min(g,g.prim)])
  }

  B.up =  B.m
  B.low = B.m

  upper = sapply(x.train,function(x_prime) min(1, min( B.m   +  Lip * abs(x_prime - eval.main)) ))
  lower= sapply(x.train,function(x_prime) max(-1, max( B.m - Lip * abs(x_prime - eval.main) )) )
  return(list(upper = upper, lower = lower ))
}

#' Creates NOT DONE, purpose of function
#'
#' @param data A data set with variables X, Y, C, W, G, and pseudo outcomes
#' @param c.vec_initial A vector of cutoff values
#' @return NOT DONE: Smoothness parameters
#'
#' @examples
#' ulbounds(NOT DONE)
#'
#' @export

ulbounds <- function(data, c.vec_initial){
  q = length(c.vec_initial)
  B.1m = B.0m = matrix(0,nrow=q,ncol=q) # storing the value of estimated cross-group differences at cutoff point

  for(g in seq(1,q-1,1))
    for(g.pr in seq(g+1,q,1)){

      temp.dat = data |> dplyr::filter(W==1 & X>=max(c.vec_initial[g.pr], c.vec_initial[g]))
      temp.vc =
        data.frame( temp.dat[,paste0("pseudo.",g)]-temp.dat[,paste0("pseudo.",g.pr)]+
                      with(temp.dat,I(G==g)*(Y-eval(parse(text =paste0("pseudo.",g))))/eval(parse(text =paste0("pseudo.ps",g))) )-
                      with(temp.dat,I(G==g.pr)*(Y-eval(parse(text =paste0("pseudo.",g.pr))))/eval(parse(text =paste0("pseudo.ps",g.pr))) )
                    ,  temp.dat $X,g,g.pr)

      names(temp.vc)[1:2]=c("psout","X")

      # Algorithm 1
      B.1m[g,g.pr] = nprobust::lprobust(temp.vc[,"psout"], temp.vc[,"X"], eval = max(c.vec_initial[g.pr], c.vec_initial[g]), bwselect="mse-dpi")$Estimate[,5]

      temp.dat = data |> dplyr::filter(W==0 & X<min(c.vec_initial[g.pr],c.vec_initial[g]))

      temp.vc = data.frame("psout"=temp.dat[,paste0("pseudo.",g)]-temp.dat[,paste0("pseudo.",g.pr)]+
                             with(temp.dat,I(G==g)*(Y-eval(parse(text =paste0("pseudo.",g))))/eval(parse(text =paste0("pseudo.ps",g))) )-
                             with(temp.dat,I(G==g.pr)*(Y-eval(parse(text =paste0("pseudo.",g.pr))))/eval(parse(text =paste0("pseudo.ps",g.pr))) ),   temp.dat $X,g,g.pr)
      names(temp.vc)[1:2]=c("psout","X")

      B.0m[g,g.pr] = nprobust::lprobust(temp.vc[,"psout"],temp.vc[,"X"],eval = min(c.vec_initial[g.pr],c.vec_initial[g]),bwselect="mse-dpi")$Estimate[,5]

    }
  B.1m= B.1m + t(-B.1m) ; B.0m= B.0m + t(-B.0m)
  result <- c(B.1m, B.0m)
  return(result)
}
