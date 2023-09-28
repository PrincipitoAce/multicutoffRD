#' Creates NOT DONE, purpose of function
#'
#' @param data_all A data set with variables X, Y, C, W, G, and pseudo outcomes
#' @param c.vec A vector of cutoff values
#' @return NOT DONE
#'
#' @examples
#' smooth(NOT DONE)
#'
#' @export

smooth = function(data_all, c.vec){
  q=length(c.vec)
  psd_dat1 = psd_dat0 = NULL; # empty vector for pseudo
  Lip_1 = Lip_0=matrix(0,q,q) # storing the value of smoothness parameter;  1/0: treatment/control
  B.1m=B.0m=matrix(0,nrow=q,ncol=q) # storing the value of estimated cross-group differences at cutoff point

  for(g in seq(1,q-1,1))
    for(g.pr in seq(g+1,q,1)){

      temp.dat = data_all |> dplyr::filter(W==1 & X>=max(c.vec[g.pr],c.vec[g]))
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

      temp.dat = data_all |> dplyr::filter(W==0 & X<min(c.vec[g.pr],c.vec[g]))

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
  return(cat("Lip_1:", Lip_1, "\n", "Lip_0:", Lip_0, "\n", "B.1m:", B.1m, "\n", "B.0m:", B.0m))
  # NOT DONE: Ask Yi if other labels would be more user-friendly, confusion on notation
}
